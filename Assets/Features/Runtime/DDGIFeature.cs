using System;
using System.Collections.Generic;
using System.Linq;
using UnityEditor.SceneManagement;
using UnityEngine;
using UnityEngine.Experimental.Rendering;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;
using UnityEngine.SceneManagement;
using ProfilingScope = UnityEngine.Rendering.ProfilingScope;
using Random = System.Random;

public sealed class DDGIFeature : ScriptableRendererFeature
{
    public sealed class DDGIPass : ScriptableRenderPass
    {
        // ----------------------------------------------------
        //              Members and Defines
        // ----------------------------------------------------
        private DDGI mddgiOverride;
        private DDGICustomBounds mCustomGIVolume;
        
        private bool mIsInitialized = false;

        private bool mNeedToResetProbeHistory = true;
        private bool mNeedToResetProbeRelocation = true;
        private bool mNeedToResetProbeClassification = true;

        private static readonly int PROBE_NUM_IRRADIANCE_INTERIOR_TEXELS = 6;
        private static readonly int PROBE_NUM_IRRADIANCE_TEXELS = PROBE_NUM_IRRADIANCE_INTERIOR_TEXELS + 2; // Including 1 texel for up border and 1 texel for down border;
        private static readonly int PROBE_NUM_DISTANCE_INTERIOR_TEXELS = 14;
        private static readonly int PROBE_NUM_DISTANCE_TEXELS = PROBE_NUM_DISTANCE_INTERIOR_TEXELS + 2;
        
        struct DDGIVolumeCpu
        {
            public Vector3 Origin;
            public Vector3 Extents;
            public Vector3Int NumProbes;
            public int MaxNumRays;
            public int NumRays;
        }
        private DDGIVolumeCpu mDDGIVolumeCpu;

        struct DDGIVolumeGpu
        {
            public Vector4 _ProbeRotation;
            public Vector3 _StartPosition;
            public int _RaysPerProbe;
            public Vector3 _ProbeSize;
            public int _MaxRaysPerProbe;
            public Vector3Int _ProbeCount;
            public float _NormalBias;
            public Vector3 _RandomVector;
            public float _EnergyPreservation;
            
            public float _RandomAngle;
            public float _HistoryBlendWeight;
            public float _IndirectIntensity;
            public float _NormalBiasMultiplier;
            
            public float _ViewBiasMultiplier;
            public int DDGI_PROBE_CLASSIFICATION;
            public int DDGI_PROBE_RELOCATION;
            public float _ProbeFixedRayBackfaceThreshold;
            
            public float _ProbeMinFrontfaceDistance;
            public int _DirectionalLightCount;
            public int _PunctualLightCount;
            public int DDGI_SKYLIGHT_MODE;

            public Vector4 _SkyboxTintColor;
            public Vector4 _SkyColor;
            public Vector4 _EquatorColor;
            public Vector4 _GroundColor;
            public Vector4 _AmbientColor;
            
            public int DDGI_PROBE_REDUCTION;
            public float _SkyboxIntensityMultiplier;
            public float _SkyboxExposure;
            public float _Pad0;
        }
        private DDGIVolumeGpu mDDGIVolumeGpu;
        private ConstantBuffer<DDGIVolumeGpu> mDDGIVolumeGpuCB;

        #region [Shader Resources]

        private RayTracingShader mDDGIRayTraceShader;
        private RayTracingAccelerationStructure mAccelerationStructure;
        
        private ComputeBuffer mRayBuffer;

        private readonly ComputeShader mUpdateIrradianceCS;
        private readonly int mUpdateIrradianceKernel;
        private readonly ComputeShader mUpdateDistanceCS;
        private readonly int mUpdateDistanceKernel;
        private readonly ComputeShader mProbeClassificationCS;
        private readonly int mResetClassificationKernel;
        private readonly int mProbeClassificationKernel;
        private readonly ComputeShader mRelocateProbeCS;
        private readonly int mResetRelocationKernel;
        private readonly int mRelocateProbeKernel;
        private readonly ComputeShader mProbeReductionCS;
        private readonly int mReductionKernel;
        private readonly int mExtraReductionKernel;

        private readonly Shader mCubemapSkyPS;
        
        #endregion

        #region [Probe Volume Textures]
        
        private enum DDGIVolumeTextureType
        {
            RayData = 0,
            Irradiance = 1,
            Distance = 2,
            ProbeData = 3,
            Variability = 4,
            VariabilityAverage = 5,
            Count
        }
        
        private RenderTexture mProbeIrradiance;
        private RenderTargetIdentifier mProbeIrradianceId;
        private RenderTexture mProbeDistance;
        private RenderTargetIdentifier mProbeDistanceId;
        private RenderTexture mProbeIrradianceHistory;
        private RenderTargetIdentifier mProbeIrradianceHistoryId;
        private RenderTexture mProbeDistanceHistory;
        private RenderTargetIdentifier mProbeDistanceHistoryId;
        private RenderTexture mProbeData;                           // For Probe Relocation
        private RenderTargetIdentifier mProbeDataId;
        private RenderTexture mProbeVariability;                    // For Probe Variability
        private RenderTargetIdentifier mProbeVariabilityId;
        private RenderTexture mProbeVariabilityAverage;             // For Probe Variability
        private RenderTargetIdentifier mProbeVariabilityAverageId;
        
        #endregion

        #region [Probe Variability]
        
        private bool mIsConverged;
        private readonly uint mMinimumVariabilitySamples = 16u;
        private bool mClearProbeVariability;
        private uint mNumVolumeVariabilitySamples = 0u;
        
        #endregion

        #region [Light Update and Change Dectect]
        
        private ComputeBuffer mDirectionalLightBuffer;
        private ComputeBuffer mPunctualLightBuffer;

        // Used to collect directional light data during the Build Light Structured Buffer process
        private struct DirectionalLight
        {
            public Vector4 direction;
            public Vector4 color;
        }

        // Used to collect punctual light data (point lights, spotlights, area lights) during the Build Light Structured Buffer process
        // Reference: RealtimeLights.hlsl 153 | Note: We consider point lights, spotlights, and area lights as punctual lights, while directional lights are not included here
        private struct PunctualLight
        {
            public Vector4 position;
            public Vector4 color;
            public Vector4 distanceAndSpotAttenuation;
            public Vector4 spotDirection;
        }

        // Used to determine the sky light mode during the Build Light Structured Buffer process
        // Raytrace shader does not support multi_compile, so we use an int define to determine the sky light mode
        private enum SkyLightMode
        {
            DDGI_SKYLIGHT_MODE_SKYBOX_CUBEMAP = 0,
            DDGI_SKYLIGHT_MODE_GRADIENT = 1,
            DDGI_SKYLIGHT_MODE_COLOR = 2,
            DDGI_SKYLIGHT_MODE_UNSUPPORTED = 3
        }

        // For URP's default Skybox parameter IDs, used during the Build Light Structured Buffer process and Probe Variability light comparison
        private static class SkyboxParam
        {
            public static readonly int _Tint = Shader.PropertyToID("_Tint");
            public static readonly int _Exposure = Shader.PropertyToID("_Exposure");
            public static readonly int _Rotation = Shader.PropertyToID("_Rotation");
            public static readonly int _Tex = Shader.PropertyToID("_Tex");
        }

        // Used only to determine if the Sky Light settings have changed in the latest frame, unrelated to the Build Light Structured Buffer process, only used for Probe Variability
        private class SkyLight
        {
            public SkyLight(Material skybox, AmbientMode ambientMode, float ambientIntensity,
                Color skyColor, Color equatorColor, Color groundColor)
            {
                if (skybox != null)
                {
                    skyboxTint = skybox.GetColor(SkyboxParam._Tint);
                    skyboxExposure = skybox.GetFloat(SkyboxParam._Exposure);
                    skyboxRotation = skybox.GetFloat(SkyboxParam._Rotation);
                    skyboxTex = skybox.GetTexture(SkyboxParam._Tex);
                }
                this.ambientMode = ambientMode;
                this.ambientIntensity = ambientIntensity;
                this.skyColor = skyColor;
                this.equatorColor = equatorColor;
                this.groundColor = groundColor;
            }

            public bool Equals(SkyLight skyLight)
            {
                bool result = true;
                if (ambientMode == skyLight.ambientMode && ambientMode == AmbientMode.Skybox)
                {
                    result &= skyboxTint == skyLight.skyboxTint;
                    result &= FloatEqual(skyboxExposure, skyLight.skyboxExposure);
                    result &= FloatEqual(skyboxRotation, skyLight.skyboxRotation);
                    result &= skyboxTex == skyLight.skyboxTex;
                }
                result &= ambientMode == skyLight.ambientMode;
                result &= FloatEqual(ambientIntensity, skyLight.ambientIntensity);
                result &= skyColor == skyLight.skyColor;
                result &= equatorColor == skyLight.equatorColor;
                result &= groundColor == skyLight.groundColor;
                return result;
            }

            private static bool FloatEqual(float a, float b) => Mathf.Abs(a - b) < 0.0001f;
            
            private Color skyboxTint = Color.black;
            private float skyboxExposure = 0.0f;
            private float skyboxRotation = 0.0f;
            private Texture skyboxTex = null;
            private AmbientMode ambientMode;
            private float ambientIntensity;
            private Color skyColor;
            private Color equatorColor;
            private Color groundColor;
        }

        // Lighting data cache, used for the Probe Variability phase
        private List<DirectionalLight> mCachedDirectionalLights = new List<DirectionalLight>();
        private List<PunctualLight> mCachedPunctualLights = new List<PunctualLight>();
        private SkyLight mCachedSkyLight = new SkyLight(null, AmbientMode.Flat, 0.0f, Color.black, Color.black, Color.black);
        private bool mAnyLightChanged;
        private bool mSkyChanged;
        
        #endregion
        
        private static class GpuParams
        {
            // For Probe Tracing and Updating
            public static readonly string RayGenShaderName = "DDGI_RayGen";
            
            public static readonly int RayBuffer = Shader.PropertyToID("RayBuffer");
            public static readonly int DirectionalLightBuffer = Shader.PropertyToID("DirectionalLightBuffer");
            public static readonly int PunctualLightBuffer = Shader.PropertyToID("PunctualLightBuffer");
            
            public static readonly int DDGIVolumeGpu = Shader.PropertyToID("DDGIVolumeGpu");

            public static readonly int _ProbeIrradiance = Shader.PropertyToID("_ProbeIrradiance");
            public static readonly int _ProbeIrradianceHistory = Shader.PropertyToID("_ProbeIrradianceHistory");
            public static readonly int _ProbeDistance = Shader.PropertyToID("_ProbeDistance");
            public static readonly int _ProbeDistanceHistory = Shader.PropertyToID("_ProbeDistanceHistory");

            public static readonly int _AccelerationStructure = Shader.PropertyToID("_AccelerationStructure");

            // For Sky Light Sampling
            public static readonly string DDGI_SKYLIGHT_MODE = "DDGI_SKYLIGHT_MODE";
            public static readonly int _SkyboxCubemap = Shader.PropertyToID("_SkyboxCubemap");

            // For Probe Relocation
            public static readonly string DDGI_PROBE_RELOCATION = "DDGI_PROBE_RELOCATION";
            public static readonly int _ProbeData = Shader.PropertyToID("_ProbeData");

            // For Probe Debug
            public static readonly string DDGI_SHOW_INDIRECT_ONLY = "DDGI_SHOW_INDIRECT_ONLY";
            public static readonly string DDGI_SHOW_PURE_INDIRECT_RADIANCE = "DDGI_SHOW_PURE_INDIRECT_RADIANCE";
            
            // For Probe Reduction (Variability)
            public static readonly int _ReductionInputSize = Shader.PropertyToID("_ReductionInputSize");
            public static readonly int _ProbeVariability = Shader.PropertyToID("_ProbeVariability");
            public static readonly int _ProbeVariabilityAverage = Shader.PropertyToID("_ProbeVariabilityAverage");
        }

        
        // ----------------------------------------------------
        //           Core Functions and Render Loops
        // ----------------------------------------------------
        public DDGIPass()
        {
            renderPassEvent = RenderPassEvent.BeforeRenderingOpaques;

            mDDGIRayTraceShader = Resources.Load<RayTracingShader>("Shaders/DDGIRayTracing");
            mUpdateIrradianceCS = Resources.Load<ComputeShader>("Shaders/DDGIUpdateIrradiance");
            mUpdateIrradianceKernel = mUpdateIrradianceCS.FindKernel("DDGIUpdateIrradiance");
            mUpdateDistanceCS = Resources.Load<ComputeShader>("Shaders/DDGIUpdateDistance");
            mUpdateDistanceKernel = mUpdateDistanceCS.FindKernel("DDGIUpdateDistance");
            mProbeClassificationCS = Resources.Load<ComputeShader>("Shaders/DDGIProbeClassification");
            mResetClassificationKernel = mProbeClassificationCS.FindKernel("DDGIProbeClassificationResetCS");
            mProbeClassificationKernel = mProbeClassificationCS.FindKernel("DDGIProbeClassificationCS");
            mRelocateProbeCS = Resources.Load<ComputeShader>("Shaders/DDGIRelocateProbe");
            mResetRelocationKernel = mRelocateProbeCS.FindKernel("DDGIResetRelocation");
            mRelocateProbeKernel = mRelocateProbeCS.FindKernel("DDGIRelocateProbe");
            mProbeReductionCS = Resources.Load<ComputeShader>("Shaders/DDGIReduction");
            mReductionKernel = mProbeReductionCS.FindKernel("DDGIReductionCS");
            mExtraReductionKernel = mProbeReductionCS.FindKernel("DDGIExtraReductionCS");

            RayTracingAccelerationStructure.RASSettings setting = new RayTracingAccelerationStructure.RASSettings
                (RayTracingAccelerationStructure.ManagementMode.Automatic, RayTracingAccelerationStructure.RayTracingModeMask.Everything,  255);
            mAccelerationStructure = new RayTracingAccelerationStructure(setting);

            mDDGIVolumeGpuCB = new ConstantBuffer<DDGIVolumeGpu>();

            // Shader.Find is not reliable, as shaders may be missing after packaging, making Find ineffective
            // For demonstration purposes, using a quick workaround here
            mCubemapSkyPS = Shader.Find("Skybox/Cubemap");
        }

        private void Initialize()
        {
            if (mIsInitialized || mddgiOverride == null) return;

            // ---------------------------------------
            // Initialize cpu-side volume parameters
            // ---------------------------------------
            var sceneBoundingBox = GenerateSceneMeshBounds();
            if (sceneBoundingBox.extents == Vector3.zero) return;   // The zero bounding box means the scene has no geometry, no GI significance
            mDDGIVolumeCpu.Origin = sceneBoundingBox.center;
            mDDGIVolumeCpu.Extents = 1.1f * sceneBoundingBox.extents;
            mDDGIVolumeCpu.NumProbes = new Vector3Int(mddgiOverride.probeCountX.value, mddgiOverride.probeCountY.value, mddgiOverride.probeCountZ.value);
            mDDGIVolumeCpu.NumRays = mddgiOverride.raysPerProbe.value;
            mDDGIVolumeCpu.MaxNumRays = 512;

            // ---------------------------------------
            // Initialize Ray Data Buffer
            // ---------------------------------------
            if (mRayBuffer != null) { mRayBuffer.Release(); mRayBuffer = null; }
            int numProbesFlat = mDDGIVolumeCpu.NumProbes.x * mDDGIVolumeCpu.NumProbes.y * mDDGIVolumeCpu.NumProbes.z;
            mRayBuffer = new ComputeBuffer(numProbesFlat * mDDGIVolumeCpu.MaxNumRays, 16 /* float4 */, ComputeBufferType.Default);

            // Note: Try to use GraphicsFormat to provide explicit identification of floating-point / fixed-point numbers
            // For example, Distance Texture, previously used RenderTextureFormat.RG32, which uses 16-bit unsigned fixed-point numbers, but distances need to be floating-point
            // Allocating RG32 as Distance Texture would ignore the decimal places of the distance, causing Chebyshev visibility tests to have Edge Clamp Artifacts.
            // ---------------------------------------
            // Radiance and Distance Texture2DArray
            // ---------------------------------------
            if (mProbeIrradiance != null) { mProbeIrradiance.Release(); mProbeIrradiance = null; }
            var probeIrradianceDimensions = GetDDGIVolumeTextureDimensions(mDDGIVolumeCpu, DDGIVolumeTextureType.Irradiance);
            mProbeIrradiance = new RenderTexture(probeIrradianceDimensions.x, probeIrradianceDimensions.y, 0, GraphicsFormat.R16G16B16A16_SFloat);
            mProbeIrradiance.filterMode = FilterMode.Bilinear;
            mProbeIrradiance.useMipMap = false;
            mProbeIrradiance.autoGenerateMips = false;
            mProbeIrradiance.enableRandomWrite = true;
            mProbeIrradiance.name = "DDGI Probe Irradiance";
            mProbeIrradiance.dimension = TextureDimension.Tex2DArray;
            mProbeIrradiance.volumeDepth = probeIrradianceDimensions.z;
            mProbeIrradiance.Create();
            mProbeIrradianceId = new RenderTargetIdentifier(mProbeIrradiance);
            
            if(mProbeDistance != null) { mProbeDistance.Release(); mProbeDistance = null; }
            var probeDistanceDimensions = GetDDGIVolumeTextureDimensions(mDDGIVolumeCpu, DDGIVolumeTextureType.Distance);
            mProbeDistance = new RenderTexture(probeDistanceDimensions.x, probeDistanceDimensions.y, 0, GraphicsFormat.R16G16_SFloat);
            mProbeDistance.filterMode = FilterMode.Bilinear;
            mProbeDistance.useMipMap = false;
            mProbeDistance.autoGenerateMips = false;
            mProbeDistance.enableRandomWrite = true;
            mProbeDistance.name = "DDGI Probe Distance";
            mProbeDistance.dimension = TextureDimension.Tex2DArray;
            mProbeDistance.volumeDepth = probeDistanceDimensions.z;
            mProbeDistance.Create();
            mProbeDistanceId = new RenderTargetIdentifier(mProbeDistance);
            
            if(mProbeIrradianceHistory != null) { mProbeIrradianceHistory.Release(); mProbeIrradianceHistory = null; }
            mProbeIrradianceHistory = new RenderTexture(mProbeIrradiance.descriptor);
            mProbeIrradianceHistory.name = "DDGI Probe Irradiance History";
            mProbeIrradianceHistory.Create();
            mProbeIrradianceHistoryId = new RenderTargetIdentifier(mProbeIrradianceHistory);
            
            if(mProbeDistanceHistory != null) { mProbeDistanceHistory.Release(); mProbeDistanceHistory = null; }
            mProbeDistanceHistory = new RenderTexture(mProbeDistance.descriptor);
            mProbeDistanceHistory.name = "DDGI Probe Distance History";
            mProbeDistanceHistory.Create();
            mProbeDistanceHistoryId = new RenderTargetIdentifier(mProbeDistanceHistory);

            // ---------------------------------------
            // Create Probe Data
            // ---------------------------------------
            if(mProbeData != null) { mProbeData.Release(); mProbeData = null; }
            var probeDataDimensions = GetDDGIVolumeTextureDimensions(mDDGIVolumeCpu, DDGIVolumeTextureType.ProbeData);
            mProbeData = new RenderTexture(probeDataDimensions.x, probeDataDimensions.y, 0, GraphicsFormat.R16G16B16A16_SFloat);
            mProbeData.filterMode = FilterMode.Bilinear;
            mProbeData.useMipMap = false;
            mProbeData.autoGenerateMips = false;
            mProbeData.enableRandomWrite = true;
            mProbeData.name = "DDGI Probe Data";
            mProbeData.dimension = TextureDimension.Tex2DArray;
            mProbeData.volumeDepth = probeDataDimensions.z;
            mProbeData.Create();
            mProbeDataId = new RenderTargetIdentifier(mProbeData);
            
            // ---------------------------------------
            // Create Probe Variability
            // ---------------------------------------
            if (mProbeVariability != null) { mProbeVariability.Release(); mProbeVariability = null; }
            var probeVariabilityDimensions = GetDDGIVolumeTextureDimensions(mDDGIVolumeCpu, DDGIVolumeTextureType.Variability);
            mProbeVariability = new RenderTexture(probeVariabilityDimensions.x, probeVariabilityDimensions.y, 0, GraphicsFormat.R32_SFloat);
            mProbeVariability.filterMode = FilterMode.Bilinear;
            mProbeVariability.useMipMap = false;
            mProbeVariability.autoGenerateMips = false;
            mProbeVariability.enableRandomWrite = true;
            mProbeVariability.name = "DDGI Probe Variability";
            mProbeVariability.dimension = TextureDimension.Tex2DArray;
            mProbeVariability.volumeDepth = probeVariabilityDimensions.z;
            mProbeVariability.Create();
            mProbeVariabilityId = new RenderTargetIdentifier(mProbeVariability);

            if (mProbeVariabilityAverage != null) { mProbeVariabilityAverage.Release(); mProbeVariabilityAverage = null; }
            var probeVariabilityAverageDimensions = GetDDGIVolumeTextureDimensions(mDDGIVolumeCpu, DDGIVolumeTextureType.VariabilityAverage);
            mProbeVariabilityAverage = new RenderTexture(mProbeVariability.descriptor);
            mProbeVariabilityAverage.graphicsFormat = GraphicsFormat.R32G32_SFloat;
            mProbeVariabilityAverage.width = probeVariabilityAverageDimensions.x;
            mProbeVariabilityAverage.height = probeVariabilityAverageDimensions.y;
            mProbeVariabilityAverage.volumeDepth = probeVariabilityAverageDimensions.z;
            mProbeVariabilityAverage.name = "DDGI Probe Variability Average";
            mProbeVariabilityAverage.Create();
            mProbeVariabilityAverageId = new RenderTargetIdentifier(mProbeVariabilityAverage);

            mIsInitialized = true;
        }

        public void Reinitialize()
        {
            mIsInitialized = false;
            mNeedToResetProbeHistory = true;
            mNeedToResetProbeRelocation = true;
            mNeedToResetProbeClassification = true;
            mClearProbeVariability = true;
            mIsConverged = false;
        }
        
        public override void OnCameraSetup(CommandBuffer cmd, ref RenderingData renderingData)
        {
            base.OnCameraSetup(cmd, ref renderingData);
            
            mddgiOverride = VolumeManager.instance.stack.GetComponent<DDGI>();

            Initialize();
        }

        public override void Execute(ScriptableRenderContext context, ref RenderingData renderingData)
        {
            if (mddgiOverride == null || !mIsInitialized) return;
            if (!mddgiOverride.IsActive()) return;

            var cmd = CommandBufferPool.Get("DDGI Pass");
            var camera = renderingData.cameraData.camera;

            ResetHistoryInfoIfNeeded(cmd);

            // Ensure this function runs before PushGpuConstants, otherwise light-related constants will not be pushed
            UpdateSceneLights(cmd);

            // Note: Each time this function is called, the random numbers are updated, which causes _RandomVector and _RandomAngle to change
            // If the random number update logic is not separated, this function can only be called once per frame!
            PushGpuConstants(cmd);

            int numProbesFlat = mDDGIVolumeCpu.NumProbes.x * mDDGIVolumeCpu.NumProbes.y * mDDGIVolumeCpu.NumProbes.z;
            
            using (new ProfilingScope(cmd, new ProfilingSampler("DDGI Ray Trace Pass")))
            {
                if (!mIsConverged)
                {
                    cmd.BuildRayTracingAccelerationStructure(mAccelerationStructure);
                    cmd.SetRayTracingAccelerationStructure(mDDGIRayTraceShader, GpuParams._AccelerationStructure, mAccelerationStructure);
                
                    cmd.SetRayTracingShaderPass(mDDGIRayTraceShader, "DDGIRayTracing");
                    cmd.SetGlobalTexture(GpuParams._ProbeIrradianceHistory, mProbeIrradianceHistoryId);
                    cmd.SetGlobalTexture(GpuParams._ProbeDistanceHistory, mProbeDistanceHistoryId);
                    cmd.SetGlobalTexture(GpuParams._ProbeData, mProbeDataId);

                    cmd.SetRayTracingBufferParam(mDDGIRayTraceShader, GpuParams.RayBuffer, mRayBuffer);
                    cmd.SetGlobalBuffer(GpuParams.DirectionalLightBuffer, mDirectionalLightBuffer);     // We will use it in closest hit shader, not in actual .raytrace shader
                    cmd.SetGlobalBuffer(GpuParams.PunctualLightBuffer, mPunctualLightBuffer);           // We will use it in closest hit shader, not in actual .raytrace shader
                    
                    cmd.DispatchRays(mDDGIRayTraceShader, GpuParams.RayGenShaderName, (uint)mDDGIVolumeCpu.NumRays, (uint)numProbesFlat, 1, camera);
                }
            }

            using (new ProfilingScope(cmd, new ProfilingSampler("DDGI Update Irradiance Pass")))
            {
                if (!mIsConverged)
                {
                    cmd.SetComputeBufferParam(mUpdateIrradianceCS, mUpdateIrradianceKernel, GpuParams.RayBuffer, mRayBuffer);
                    cmd.SetComputeTextureParam(mUpdateIrradianceCS, mUpdateIrradianceKernel, GpuParams._ProbeIrradiance, mProbeIrradianceId);
                    cmd.SetComputeTextureParam(mUpdateIrradianceCS, mUpdateIrradianceKernel, GpuParams._ProbeIrradianceHistory, mProbeIrradianceHistoryId);
                    cmd.SetComputeTextureParam(mUpdateIrradianceCS, mUpdateIrradianceKernel, GpuParams._ProbeVariability, mProbeVariabilityId);

                    // Note that we use Y-UP, so the Dispatch needs to be reversed
                    cmd.DispatchCompute(mUpdateIrradianceCS, mUpdateIrradianceKernel, mDDGIVolumeCpu.NumProbes.x, mDDGIVolumeCpu.NumProbes.z, mDDGIVolumeCpu.NumProbes.y);
                }
            }

            using (new ProfilingScope(cmd, new ProfilingSampler("DDGI Update Distance Pass")))
            {
                if (!mIsConverged)
                {
                    cmd.SetComputeBufferParam(mUpdateDistanceCS, mUpdateDistanceKernel, GpuParams.RayBuffer, mRayBuffer);
                    cmd.SetComputeTextureParam(mUpdateDistanceCS, mUpdateDistanceKernel, GpuParams._ProbeDistance, mProbeDistanceId);
                    cmd.SetComputeTextureParam(mUpdateDistanceCS, mUpdateDistanceKernel, GpuParams._ProbeDistanceHistory, mProbeDistanceHistoryId);

                    // Note that we use Y-UP, so the Dispatch needs to be reversed
                    cmd.DispatchCompute(mUpdateDistanceCS, mUpdateDistanceKernel, mDDGIVolumeCpu.NumProbes.x, mDDGIVolumeCpu.NumProbes.z, mDDGIVolumeCpu.NumProbes.y);
                }
            }


            if (mddgiOverride.enableProbeRelocation.value)
            {
                using (new ProfilingScope(cmd, new ProfilingSampler("DDGI Relocate Probe Pass")))
                {
                    var numGroupsX = Mathf.CeilToInt(numProbesFlat / 32.0f /*relocationGroupSizeX*/);
                    
                    if (mNeedToResetProbeRelocation)
                    {
                        cmd.SetComputeTextureParam(mRelocateProbeCS, mResetRelocationKernel, GpuParams._ProbeData, mProbeDataId);
                        cmd.DispatchCompute(mRelocateProbeCS, mResetRelocationKernel, numGroupsX, 1, 1);
                        mNeedToResetProbeRelocation = false;
                    }
                    
                    cmd.SetComputeTextureParam(mRelocateProbeCS, mRelocateProbeKernel, GpuParams._ProbeData, mProbeDataId);
                    cmd.SetComputeBufferParam(mRelocateProbeCS, mRelocateProbeKernel, GpuParams.RayBuffer, mRayBuffer);
                    cmd.DispatchCompute(mRelocateProbeCS, mRelocateProbeKernel, numGroupsX, 1, 1);
                }
            }
            else
            {
                if (!mNeedToResetProbeRelocation)
                {
                    var numGroupsX = Mathf.CeilToInt(numProbesFlat / 32.0f /*relocationGroupSizeX*/);
                    
                    cmd.SetComputeTextureParam(mRelocateProbeCS, mResetRelocationKernel, GpuParams._ProbeData, mProbeDataId);
                    cmd.DispatchCompute(mRelocateProbeCS, mResetRelocationKernel, numGroupsX, 1, 1);
                    mNeedToResetProbeRelocation = true;
                }
            }

            if (mddgiOverride.enableProbeClassification.value)
            {
                using (new ProfilingScope(cmd, new ProfilingSampler("DDGI Classify Probe Pass")))
                {
                    var numGroupsX = Mathf.CeilToInt(numProbesFlat / 32.0f /*relocationGroupSizeX*/);

                    if (mNeedToResetProbeClassification)
                    {
                        cmd.SetComputeTextureParam(mProbeClassificationCS, mResetClassificationKernel, GpuParams._ProbeData, mProbeDataId);
                        cmd.DispatchCompute(mProbeClassificationCS, mResetClassificationKernel, numGroupsX, 1, 1);
                        mNeedToResetProbeClassification = false;
                    }
                    
                    cmd.SetComputeTextureParam(mProbeClassificationCS, mProbeClassificationKernel, GpuParams._ProbeData, mProbeDataId);
                    cmd.SetComputeBufferParam(mProbeClassificationCS, mProbeClassificationKernel, GpuParams.RayBuffer, mRayBuffer);
                    cmd.DispatchCompute(mProbeClassificationCS, mProbeClassificationKernel, numGroupsX, 1, 1);
                }
            }
            else
            {
                if (!mNeedToResetProbeClassification)
                {
                    var numGroupsX = Mathf.CeilToInt(numProbesFlat / 32.0f /*relocationGroupSizeX*/);
                    
                    cmd.SetComputeTextureParam(mProbeClassificationCS, mResetClassificationKernel, GpuParams._ProbeData, mProbeDataId);
                    cmd.DispatchCompute(mProbeClassificationCS, mResetClassificationKernel, numGroupsX, 1, 1);
                    mNeedToResetProbeClassification = true;
                }
            }
            
            if (mddgiOverride.enableProbeVariability.value)
            {
                using (new ProfilingScope(cmd, new ProfilingSampler("DDGI Variability Pass")))
                {
                    // TODO: Y-UP Probe Volume hardcoding, if you need to modify the Volume axis, branches will be required
                    var inputTexels = new Vector3Int(mDDGIVolumeCpu.NumProbes.x * PROBE_NUM_IRRADIANCE_INTERIOR_TEXELS,
                        mDDGIVolumeCpu.NumProbes.z * PROBE_NUM_IRRADIANCE_INTERIOR_TEXELS,
                        mDDGIVolumeCpu.NumProbes.y);
                    var NumThreadsInGroup = new Vector3Int(4, 8, 4);
                    var ThreadSampleFootprint = new Vector2Int(4, 2);
                    
                    // -------------------------
                    // First Reduction Pass
                    // -------------------------
                    {
                        cmd.SetComputeTextureParam(mProbeReductionCS, mReductionKernel, GpuParams._ProbeVariability, mProbeVariabilityId);
                        cmd.SetComputeTextureParam(mProbeReductionCS, mReductionKernel, GpuParams._ProbeVariabilityAverage, mProbeVariabilityAverageId);
                        cmd.SetComputeVectorParam(mProbeReductionCS, GpuParams._ReductionInputSize, new Vector4(inputTexels.x, inputTexels.y, inputTexels.z, 0.0f));

                        var outputTexelsX = Mathf.CeilToInt((float)inputTexels.x / (float)(NumThreadsInGroup.x * ThreadSampleFootprint.x));
                        var outputTexelsY = Mathf.CeilToInt((float)inputTexels.y / (float)(NumThreadsInGroup.y * ThreadSampleFootprint.y));
                        var outputTexelsZ = Mathf.CeilToInt((float)inputTexels.z / (float)(NumThreadsInGroup.z));
                    
                        cmd.DispatchCompute(mProbeReductionCS, mReductionKernel, outputTexelsX, outputTexelsY, outputTexelsZ);

                        inputTexels = new Vector3Int(outputTexelsX, outputTexelsY, outputTexelsZ);
                    }
                    
                    // -------------------------
                    // Extra Reduction Pass
                    // -------------------------
                    {
                        while (inputTexels.x > 1 || inputTexels.y > 1 || inputTexels.z > 1)
                        {
                            var outputTexelsX = Mathf.CeilToInt((float)inputTexels.x / (float)(NumThreadsInGroup.x * ThreadSampleFootprint.x));
                            var outputTexelsY = Mathf.CeilToInt((float)inputTexels.y / (float)(NumThreadsInGroup.y * ThreadSampleFootprint.y));
                            var outputTexelsZ = Mathf.CeilToInt((float)inputTexels.z / (float)(NumThreadsInGroup.z));
                            
                            cmd.SetComputeTextureParam(mProbeReductionCS, mExtraReductionKernel, GpuParams._ProbeVariabilityAverage, mProbeVariabilityAverageId);
                            cmd.SetComputeVectorParam(mProbeReductionCS, GpuParams._ReductionInputSize, new Vector4(inputTexels.x, inputTexels.y, inputTexels.z, 0.0f));
                            
                            cmd.DispatchCompute(mProbeReductionCS, mExtraReductionKernel, outputTexelsX, outputTexelsY, outputTexelsZ);
                            
                            inputTexels = new Vector3Int(outputTexelsX, outputTexelsY, outputTexelsZ);
                        }
                    }
                    
                    // ---------------------------------
                    // Readback From Variability Average
                    // ---------------------------------
                    // Grab First Pixel of Variability Average
                    AsyncGPUReadback.Request(mProbeVariabilityAverage, 0, 0, 1, 0, 1, 0, 1, VariabilityEstimate);
                }
            }
            else
            {
                // If the variability feature is not enabled, we assume that the integration process will never converge
                mIsConverged = false;
                mClearProbeVariability = true;
                mNumVolumeVariabilitySamples = 0u;
            }

            cmd.CopyTexture(mProbeIrradianceId, mProbeIrradianceHistoryId);
            cmd.CopyTexture(mProbeDistanceId, mProbeDistanceHistoryId);

            context.ExecuteCommandBuffer(cmd);
            CommandBufferPool.Release(cmd);
        }

        public void Release()
        {
            if (mRayBuffer != null) { mRayBuffer.Release(); mRayBuffer = null; }
            if (mDirectionalLightBuffer != null) { mDirectionalLightBuffer.Release(); mDirectionalLightBuffer = null; }
            if (mPunctualLightBuffer != null) { mPunctualLightBuffer.Release(); mPunctualLightBuffer = null; }
            
            if (mAccelerationStructure != null) { mAccelerationStructure.Release(); mAccelerationStructure = null; }
            
            if (mProbeIrradiance != null) { mProbeIrradiance.Release(); mProbeIrradiance = null; }
            if (mProbeDistance != null) { mProbeDistance.Release(); mProbeDistance = null; }
            if (mProbeIrradianceHistory != null) { mProbeIrradianceHistory.Release(); mProbeIrradianceHistory = null; }
            if (mProbeDistanceHistory != null) { mProbeDistanceHistory.Release(); mProbeDistanceHistory = null; }
            if (mProbeData != null) { mProbeData.Release(); mProbeData = null; }
            if (mProbeVariability != null) { mProbeVariability.Release(); mProbeVariability = null; }
            if (mProbeVariabilityAverage != null) { mProbeVariabilityAverage.Release(); mProbeVariabilityAverage = null; }

            if (mDDGIVolumeGpuCB != null) { mDDGIVolumeGpuCB.Release(); mDDGIVolumeGpuCB = null; }
        }
        
        
        // ----------------------------------------------------
        //           Wrapper and Utility Functions
        // ----------------------------------------------------
        public Vector3Int GetNumProbes() => mDDGIVolumeCpu.NumProbes; // For Visualization Pass
        public RenderTargetIdentifier GetProbeData() => mProbeDataId; // For Visualization Pass

        private void ResetHistoryInfoIfNeeded(CommandBuffer cmd)
        {
            if (mNeedToResetProbeHistory)
            {
                CoreUtils.SetRenderTarget(cmd, mProbeIrradianceHistoryId, ClearFlag.Color, new Color(0,0,0,0));
                CoreUtils.SetRenderTarget(cmd, mProbeDistanceHistoryId, ClearFlag.Color, new Color(0,0,0,0));
                mNeedToResetProbeHistory = false;
            }
        }
        
        private void PushGpuConstants(CommandBuffer cmd)
        {
            var random = (float)NextDouble(new Random(), 0.0f, 1.0f, 5); // Generate a random number between 0 and 1, with 5 decimal places
            var randomVec = Vector3.Normalize(new Vector3(2.0f * random - 1.0f, 2.0f * random - 1.0f, 2.0f * random - 1.0f));
            var randomAngle = random * Mathf.PI * 2.0f;

            // -------------------------------------------------
            // Fill GPU constants (lighting-related constants are updated in UpdateSceneLights)
            // -------------------------------------------------
            Quaternion rotation;
            if (mddgiOverride.useCustomBounds.value && mCustomGIVolume != null) { rotation = mCustomGIVolume.transform.rotation; }
            else { rotation = Quaternion.Euler(mddgiOverride.probeRotationDegrees.value); }
            mDDGIVolumeGpu._ProbeRotation = new Vector4(rotation.x, rotation.y, rotation.z, rotation.w);
            mDDGIVolumeGpu._StartPosition = mDDGIVolumeCpu.Origin - mDDGIVolumeCpu.Extents;
            mDDGIVolumeGpu._RaysPerProbe = mDDGIVolumeCpu.NumRays;
            var a = 2.0f * mDDGIVolumeCpu.Extents;
            var b = new Vector3(mDDGIVolumeCpu.NumProbes.x, mDDGIVolumeCpu.NumProbes.y, mDDGIVolumeCpu.NumProbes.z) - Vector3.one;
            mDDGIVolumeGpu._ProbeSize = new Vector3(a.x / b.x, a.y / b.y, a.z / b.z);
            mDDGIVolumeGpu._MaxRaysPerProbe = mDDGIVolumeCpu.MaxNumRays;
            mDDGIVolumeGpu._ProbeCount = new Vector3Int(mDDGIVolumeCpu.NumProbes.x, mDDGIVolumeCpu.NumProbes.y, mDDGIVolumeCpu.NumProbes.z);
            mDDGIVolumeGpu._NormalBias = 0.25f;
            mDDGIVolumeGpu._RandomVector = randomVec;
            mDDGIVolumeGpu._EnergyPreservation = 0.85f;
            mDDGIVolumeGpu._RandomAngle = randomAngle;
            mDDGIVolumeGpu._HistoryBlendWeight = 0.98f;
            mDDGIVolumeGpu._IndirectIntensity = mddgiOverride.indirectIntensity.value;
            mDDGIVolumeGpu._NormalBiasMultiplier = mddgiOverride.normalBiasMultiplier.value;
            mDDGIVolumeGpu._ViewBiasMultiplier = mddgiOverride.viewBiasMultiplier.value;
            mDDGIVolumeGpu.DDGI_PROBE_CLASSIFICATION = mddgiOverride.enableProbeClassification.value ? 1 : 0;
            mDDGIVolumeGpu.DDGI_PROBE_RELOCATION = mddgiOverride.enableProbeRelocation.value ? 1 : 0;
            mDDGIVolumeGpu._ProbeFixedRayBackfaceThreshold = mddgiOverride.probeFixedRayBackfaceThreshold.value;
            mDDGIVolumeGpu._ProbeMinFrontfaceDistance = mddgiOverride.probeMinFrontfaceDistance.value;
            mDDGIVolumeGpu.DDGI_PROBE_REDUCTION = mddgiOverride.enableProbeVariability.value ? 1 : 0;
            mDDGIVolumeGpu._Pad0 = 0.0f;
            mDDGIVolumeGpuCB?.PushGlobal(cmd, mDDGIVolumeGpu, GpuParams.DDGIVolumeGpu); // TODO: How can this be null? Debug!

            // -------------------------------------------------
            // Shader Keywords.
            // -------------------------------------------------
            cmd.DisableShaderKeyword(GpuParams.DDGI_SHOW_INDIRECT_ONLY);
            cmd.DisableShaderKeyword(GpuParams.DDGI_SHOW_PURE_INDIRECT_RADIANCE);
            if (mddgiOverride.debugIndirect.value)
            {
                switch (mddgiOverride.indirectDebugMode.value)
                {
                    case IndirectDebugMode.FullIndirectRadiance:
                        cmd.EnableShaderKeyword(GpuParams.DDGI_SHOW_INDIRECT_ONLY);
                        break;
                    case IndirectDebugMode.PureIndirectRadiance:
                        cmd.EnableShaderKeyword(GpuParams.DDGI_SHOW_PURE_INDIRECT_RADIANCE);
                        break;
                    default:
                        break;
                }
            }
        }

        #region [Light Update]

        // Update all lights and monitor changes in light data
        private void UpdateSceneLights(CommandBuffer cmd)
        {
            BuildLightStructuredBuffer(cmd);
            UpdateSkyLight(cmd);
            mClearProbeVariability = mAnyLightChanged || mSkyChanged;
        }

        // Unity by default culls extra lights in the scene, which affects our ability to collect global lighting information for the scene.
        // Therefore, we have to manually collect this on the CPU side.
        private void BuildLightStructuredBuffer(CommandBuffer cmd)
        {
            var cpuLights = FindObjectsOfType<Light>();

            var gpuDirectionalLights = new List<DirectionalLight>();
            var gpuPunctualLights = new List<PunctualLight>();
            foreach (var cpuLight in cpuLights)
            {
                if (cpuLight.lightmapBakeType == LightmapBakeType.Baked) continue;

                // Dynamic global lighting for area lights is not supported yet...
                if (cpuLight.type == LightType.Point || cpuLight.type == LightType.Spot)

                {
                    var position = cpuLight.transform.position;
                    var color = cpuLight.color * cpuLight.intensity;
                    var lightAttenuation = new Vector4(0.0f, 1.0f, 0.0f, 1.0f);
                    var lightSpotDir = new Vector4(0.0f, 0.0f, 1.0f, 0.0f);
                    
                    GetPunctualLightDistanceAttenuation(cpuLight.range, ref lightAttenuation);

                    if (cpuLight.type == LightType.Spot)
                    {
                        GetSpotDirection(cpuLight.transform.forward, out lightSpotDir);
                        GetSpotAngleAttenuation(cpuLight.spotAngle, cpuLight.innerSpotAngle, ref lightAttenuation);
                    }
                    
                    PunctualLight punctualLight;
                    punctualLight.position = new Vector4(position.x, position.y, position.z, 1.0f);
                    punctualLight.color = color;
                    punctualLight.distanceAndSpotAttenuation = lightAttenuation;
                    punctualLight.spotDirection = lightSpotDir;
                    
                    gpuPunctualLights.Add(punctualLight);
                }
                else if (cpuLight.type == LightType.Directional)
                {
                    var lightForward = cpuLight.transform.forward;
                    
                    DirectionalLight directionalLight;
                    directionalLight.direction = new Vector4(-lightForward.x, -lightForward.y, -lightForward.z, 0.0f);
                    directionalLight.color = cpuLight.color;
                    
                    gpuDirectionalLights.Add(directionalLight);
                }
            }

            // If the light array size is 0, only allocate a buffer with 1 element. Creating a ComputeBuffer with size 0 will cause an error.
            if (mDirectionalLightBuffer != null) { mDirectionalLightBuffer.Release(); mDirectionalLightBuffer = null; }
            mDirectionalLightBuffer = new ComputeBuffer(Mathf.Max(gpuDirectionalLights.Count, 1), 2 * 16, ComputeBufferType.Default);

            if (mPunctualLightBuffer != null) { mPunctualLightBuffer.Release(); mPunctualLightBuffer = null; }
            mPunctualLightBuffer = new ComputeBuffer(Mathf.Max(gpuPunctualLights.Count, 1), 4 * 16, ComputeBufferType.Default);

            mDirectionalLightBuffer.SetData(gpuDirectionalLights.ToArray());
            mPunctualLightBuffer.SetData(gpuPunctualLights.ToArray());

            mDDGIVolumeGpu._DirectionalLightCount = gpuDirectionalLights.Count;
            mDDGIVolumeGpu._PunctualLightCount = gpuPunctualLights.Count;
            /*cmd.SetGlobalInt(GpuParams._DirectionalLightCount, gpuDirectionalLights.Count);
            cmd.SetGlobalInt(GpuParams._PunctualLightCount, gpuPunctualLights.Count);*/



            // -----------------------------
            // Any Light Changed Determine
            // -----------------------------
            mAnyLightChanged = (!mCachedDirectionalLights.SequenceEqual(gpuDirectionalLights)) || (!mCachedPunctualLights.SequenceEqual(gpuPunctualLights));
            if (mAnyLightChanged)
            {
                mCachedDirectionalLights = new List<DirectionalLight>(gpuDirectionalLights);
                mCachedPunctualLights = new List<PunctualLight>(gpuPunctualLights);
            }
        }
        
        // Reference: UniversalRenderPipelineCore.cs 1634
        private static void GetPunctualLightDistanceAttenuation(float lightRange, ref Vector4 lightAttenuation)
        {
            // Light attenuation in universal matches the unity vanilla one (HINT_NICE_QUALITY).
            // attenuation = 1.0 / distanceToLightSqr
            // The smoothing factor makes sure that the light intensity is zero at the light range limit.
            // (We used to offer two different smoothing factors.)

            // The current smoothing factor matches the one used in the Unity lightmapper.
            // smoothFactor = (1.0 - saturate((distanceSqr * 1.0 / lightRangeSqr)^2))^2
            float lightRangeSqr = lightRange * lightRange;
            float fadeStartDistanceSqr = 0.8f * 0.8f * lightRangeSqr;
            float fadeRangeSqr = (fadeStartDistanceSqr - lightRangeSqr);
            float lightRangeSqrOverFadeRangeSqr = -lightRangeSqr / fadeRangeSqr;
            float oneOverLightRangeSqr = 1.0f / Mathf.Max(0.0001f, lightRangeSqr);

            // On all devices: Use the smoothing factor that matches the GI.
            lightAttenuation.x = oneOverLightRangeSqr;
            lightAttenuation.y = lightRangeSqrOverFadeRangeSqr;
        }

        // Reference: UniversalRenderPipelineCore.cs 1654
        private static void GetSpotAngleAttenuation(float spotAngle, float? innerSpotAngle, ref Vector4 lightAttenuation)
        {
            // Spot Attenuation with a linear falloff can be defined as
            // (SdotL - cosOuterAngle) / (cosInnerAngle - cosOuterAngle)
            // This can be rewritten as
            // invAngleRange = 1.0 / (cosInnerAngle - cosOuterAngle)
            // SdotL * invAngleRange + (-cosOuterAngle * invAngleRange)
            // If we precompute the terms in a MAD instruction
            float cosOuterAngle = Mathf.Cos(Mathf.Deg2Rad * spotAngle * 0.5f);
            // We need to do a null check for particle lights
            // This should be changed in the future
            // Particle lights will use an inline function
            float cosInnerAngle;
            if (innerSpotAngle.HasValue)
                cosInnerAngle = Mathf.Cos(innerSpotAngle.Value * Mathf.Deg2Rad * 0.5f);
            else
                cosInnerAngle = Mathf.Cos((2.0f * Mathf.Atan(Mathf.Tan(spotAngle * 0.5f * Mathf.Deg2Rad) * (64.0f - 18.0f) / 64.0f)) * 0.5f);
            float smoothAngleRange = Mathf.Max(0.001f, cosInnerAngle - cosOuterAngle);
            float invAngleRange = 1.0f / smoothAngleRange;
            float add = -cosOuterAngle * invAngleRange;

            lightAttenuation.z = invAngleRange;
            lightAttenuation.w = add;
        }
        
        // Reference: UniversalRenderPipelineCore.cs 1681
        private static void GetSpotDirection(Vector3 forward, out Vector4 lightSpotDir)
        {
            lightSpotDir = new Vector4(-forward.x, -forward.y, -forward.z, 0.0f);
        }

        // Update the skylight for sampling in the Miss Shader (Window->Rendering->Lighting)
        private void UpdateSkyLight(CommandBuffer cmd)
        {
            // -----------------------------
            // Sky Light Changed Determine
            // -----------------------------
            var currSkyLight = new SkyLight(RenderSettings.skybox, RenderSettings.ambientMode, RenderSettings.ambientIntensity,
                RenderSettings.ambientSkyColor, RenderSettings.ambientEquatorColor, RenderSettings.ambientGroundColor);

            mSkyChanged = !mCachedSkyLight.Equals(currSkyLight);
            if (mSkyChanged) { mCachedSkyLight = currSkyLight; }
            
            switch (RenderSettings.ambientMode)
            {
                case AmbientMode.Skybox:
                    UpdateSkyLightAsSkybox(cmd);
                    break;
                case AmbientMode.Trilight:
                    UpdateSkyLightAsGradient(cmd);
                    break;
                case AmbientMode.Flat:
                    UpdateSkyLightAsColor(cmd);
                    break;
            }
        }

        private void UpdateSkyLightAsSkybox(CommandBuffer cmd)
        {
            var skybox = RenderSettings.skybox;
            if (skybox == null)
            {
                // If the skybox material is not set correctly, fallback to a solid color (Ambient Color), consistent with Unity's behavior
                UpdateSkyLightAsColor(cmd);
                return;
            }

            if (mCubemapSkyPS == null)
            {
                Debug.LogWarning("DDGIFeature failed to find the built-in skybox shader in URP, please check.");
                UpdateSkyLightAsBlack(cmd);
                return;
            }

            if (skybox.shader == mCubemapSkyPS)
            {
                mDDGIVolumeGpu.DDGI_SKYLIGHT_MODE = (int)SkyLightMode.DDGI_SKYLIGHT_MODE_SKYBOX_CUBEMAP;
                mDDGIVolumeGpu._SkyboxIntensityMultiplier = RenderSettings.ambientIntensity;
                mDDGIVolumeGpu._SkyboxTintColor = skybox.GetColor(SkyboxParam._Tint);
                mDDGIVolumeGpu._SkyboxExposure = skybox.GetFloat(SkyboxParam._Exposure);
                /*cmd.SetRayTracingIntParam(mDDGIRayTraceShader, GpuParams.DDGI_SKYLIGHT_MODE, (int)SkyLightMode.DDGI_SKYLIGHT_MODE_SKYBOX_CUBEMAP);
                cmd.SetRayTracingFloatParam(mDDGIRayTraceShader, GpuParams._SkyboxIntensityMultiplier, RenderSettings.ambientIntensity);
                cmd.SetRayTracingVectorParam(mDDGIRayTraceShader, GpuParams._SkyboxTintColor, skybox.GetColor(SkyboxParam._Tint));
                cmd.SetRayTracingFloatParam(mDDGIRayTraceShader, GpuParams._SkyboxExposure, skybox.GetFloat(SkyboxParam._Exposure));*/
                cmd.SetRayTracingTextureParam(mDDGIRayTraceShader, GpuParams._SkyboxCubemap, skybox.GetTexture(SkyboxParam._Tex));
            }
            else
            {
                // Currently, we only support the most commonly used Cubemap-style skyboxes. Other types of skyboxes are not supported and will fallback to pure black.
                UpdateSkyLightAsBlack(cmd);
            }
        }

        private void UpdateSkyLightAsGradient(CommandBuffer cmd)
        {
            mDDGIVolumeGpu.DDGI_SKYLIGHT_MODE = (int)SkyLightMode.DDGI_SKYLIGHT_MODE_GRADIENT;
            mDDGIVolumeGpu._SkyColor = RenderSettings.ambientSkyColor;
            mDDGIVolumeGpu._EquatorColor = RenderSettings.ambientEquatorColor;
            mDDGIVolumeGpu._GroundColor = RenderSettings.ambientGroundColor;
            /*cmd.SetRayTracingIntParam(mDDGIRayTraceShader, GpuParams.DDGI_SKYLIGHT_MODE, (int)SkyLightMode.DDGI_SKYLIGHT_MODE_GRADIENT);
            cmd.SetRayTracingVectorParam(mDDGIRayTraceShader, GpuParams._SkyColor, RenderSettings.ambientSkyColor);
            cmd.SetRayTracingVectorParam(mDDGIRayTraceShader, GpuParams._EquatorColor, RenderSettings.ambientEquatorColor);
            cmd.SetRayTracingVectorParam(mDDGIRayTraceShader, GpuParams._GroundColor, RenderSettings.ambientGroundColor);*/
        }

        private void UpdateSkyLightAsColor(CommandBuffer cmd)
        {
            mDDGIVolumeGpu.DDGI_SKYLIGHT_MODE = (int)SkyLightMode.DDGI_SKYLIGHT_MODE_COLOR;
            mDDGIVolumeGpu._AmbientColor = RenderSettings.ambientSkyColor;
            /*cmd.SetRayTracingIntParam(mDDGIRayTraceShader, GpuParams.DDGI_SKYLIGHT_MODE, (int)SkyLightMode.DDGI_SKYLIGHT_MODE_COLOR);
            cmd.SetRayTracingVectorParam(mDDGIRayTraceShader, GpuParams._AmbientColor, RenderSettings.ambientSkyColor);*/
        }

        private void UpdateSkyLightAsBlack(CommandBuffer cmd)
        {
            mDDGIVolumeGpu.DDGI_SKYLIGHT_MODE = (int)SkyLightMode.DDGI_SKYLIGHT_MODE_UNSUPPORTED;
            //cmd.SetRayTracingIntParam(mDDGIRayTraceShader, GpuParams.DDGI_SKYLIGHT_MODE, (int)SkyLightMode.DDGI_SKYLIGHT_MODE_UNSUPPORTED);
        }

        #endregion
        
        private void VariabilityEstimate(AsyncGPUReadbackRequest request)
        {
            if (request.hasError)
            {
                Debug.LogError("DDGI: An error occurred while reading back the Variability Average!");
            }
            else if (request.done)
            {
                // Our Variability Average uses the R32G32_SFLOAT format, so reading a float on the CPU corresponds exactly to 32 bits
                // At this point, the size of readbackPixels should be 2, corresponding to the R and G channels, and we only take the R channel
                var readbackPixels = request.GetData<float>().ToArray();
                if (readbackPixels.Length > 0)
                {
                    var volumeAverageVariability = readbackPixels[0];

                    if (mClearProbeVariability) mNumVolumeVariabilitySamples = 0;
                    
                    mIsConverged = (mNumVolumeVariabilitySamples++ > mMinimumVariabilitySamples) &&
                                   (volumeAverageVariability < mddgiOverride.probeVariabilityThreshold.value);
                }
                else
                {
                    Debug.LogError("DDGI: Variability Average readback completed, but unexpectedly returned empty data, please check.");
                }
            }
        }

        private Bounds GenerateSceneMeshBounds()
        {
            Bounds bounds = new Bounds(Vector3.zero, Vector3.zero);

            if (mddgiOverride != null && mddgiOverride.useCustomBounds.value)
            {
                // Currently, only a single custom bounding box is supported
                mCustomGIVolume = FindFirstObjectByType<DDGICustomBounds>();
                var boxCollider = mCustomGIVolume.GetComponent<BoxCollider>();
                if (boxCollider != null) bounds = boxCollider.bounds;
            }
            else
            {
                // Automatically generate the bounding box based on the scene's meshes
                foreach (var meshRenderer in FindObjectsOfType<MeshRenderer>())
                {
                    bounds.Encapsulate(meshRenderer.bounds);
                }

                // In theory, we won't update the bounding box frame by frame, so it's not necessary to force include skinned meshes. Removing the following part is also fine.
                foreach (var skinnedMeshRenderer in FindObjectsOfType<SkinnedMeshRenderer>())
                {
                    bounds.Encapsulate(skinnedMeshRenderer.bounds);
                }
            }

            return bounds;
        }

        private static Vector3Int GetDDGIVolumeTextureDimensions(DDGIVolumeCpu volumeDescCpu, DDGIVolumeTextureType type)
        {
            // TODO: Y-UP Probe Volume hardcoded, if you need to modify the Volume axis direction, branching is required
            // In Unity, we use Y-UP for the DDGI Volume
            // Our texture encoding principle is: the axis that points up corresponds to the array size
            var width = volumeDescCpu.NumProbes.x;
            var height = volumeDescCpu.NumProbes.z;
            var arraySize = volumeDescCpu.NumProbes.y;

            // Since ProbeData is one-by-one, there's no need for additional branching, just return the code above directly
            switch (type)
            {
                case DDGIVolumeTextureType.RayData:
                    {
                        // Each row represents all ray info for a single probe; the number of rows is the total number of probes in a plane
                        height = width * height;
                        width = volumeDescCpu.NumRays;
                        break;
                    }
                case DDGIVolumeTextureType.Irradiance:
                    {
                        width *= PROBE_NUM_IRRADIANCE_TEXELS;
                    height *= PROBE_NUM_IRRADIANCE_TEXELS;
                    break;    
                }
                case DDGIVolumeTextureType.Distance:
                {
                    width *= PROBE_NUM_DISTANCE_TEXELS;
                    height *= PROBE_NUM_DISTANCE_TEXELS;
                    break;
                }
                case DDGIVolumeTextureType.Variability:
                {
                    width *= PROBE_NUM_IRRADIANCE_INTERIOR_TEXELS;
                    height *= PROBE_NUM_IRRADIANCE_INTERIOR_TEXELS;
                    break;
                }
                case DDGIVolumeTextureType.VariabilityAverage:
                {
                    width *= PROBE_NUM_IRRADIANCE_INTERIOR_TEXELS;
                    height *= PROBE_NUM_IRRADIANCE_INTERIOR_TEXELS;
                    
                    var NumThreadsInGroup = new Vector3Int(4, 8, 4);
                    var DimensionScale = new Vector3Int(NumThreadsInGroup.x * 4, NumThreadsInGroup.y * 2, NumThreadsInGroup.z);
                    
                    width = (width + DimensionScale.x - 1) / DimensionScale.x;
                    height = (height + DimensionScale.y - 1) / DimensionScale.y;
                    arraySize = (arraySize + DimensionScale.z - 1) / DimensionScale.z;
                    
                    break;
                }
            }

            return new Vector3Int(width, height, arraySize);
        }

        // Random Generator
        private static double NextDouble(Random ran, double minValue, double maxValue, int decimalPlace)
        {
            double randNum = ran.NextDouble() * (maxValue - minValue) + minValue;
            return Convert.ToDouble(randNum.ToString("f" + decimalPlace));
        }
    }

    public sealed class DDGIVisualizePass : ScriptableRenderPass
    {
        private DDGI mddgiOverride;
        
        private Shader mVisualizeShader;
        private Material mVisualizeMaterial;
        private Mesh mVisualizeMesh;

        private DDGIPass mDDGIPass;

        private ComputeBuffer mArgsBuffer;

        private static class GpuParams
        {
            public static readonly string DDGI_DEBUG_IRRADIANCE = "DDGI_DEBUG_IRRADIANCE";
            public static readonly string DDGI_DEBUG_DISTANCE = "DDGI_DEBUG_DISTANCE";
            public static readonly string DDGI_DEBUG_OFFSET = "DDGI_DEBUG_OFFSET";

            public static readonly int _ProbeData = Shader.PropertyToID("_ProbeData");
            public static readonly int _ddgiSphere_ObjectToWorld = Shader.PropertyToID("_ddgiSphere_ObjectToWorld");
        }
        
        public DDGIVisualizePass()
        {
            renderPassEvent = RenderPassEvent.AfterRenderingOpaques;

            mVisualizeShader = Resources.Load<Shader>("Shaders/DDGIVisualize");
            mVisualizeMaterial = CoreUtils.CreateEngineMaterial(mVisualizeShader);
            mVisualizeMaterial.enableInstancing = true;
        }

        public void Setup(Mesh debugMesh, DDGIPass ddgiPass)
        {
            mVisualizeMesh = debugMesh;
            mDDGIPass = ddgiPass;
        }
        
        public override void Configure(CommandBuffer cmd, RenderTextureDescriptor cameraTextureDescriptor)
        {
            base.Configure(cmd, cameraTextureDescriptor);

            mddgiOverride = VolumeManager.instance.stack.GetComponent<DDGI>();
        }

        public override void Execute(ScriptableRenderContext context, ref RenderingData renderingData)
        {
            if (mVisualizeMesh == null || mDDGIPass == null) return;
            if (mddgiOverride == null || !mddgiOverride.IsActive()) return;

            var ddgiOverride = VolumeManager.instance.stack.GetComponent<DDGI>();
            if (ddgiOverride == null) return;

            if (!ddgiOverride.debugProbe.value) return;
            
            var cmd = CommandBufferPool.Get("DDGI Visualize");
            var camera = renderingData.cameraData.camera;
            var renderer = renderingData.cameraData.renderer;

            // Configure Debug Mode.
            {
                cmd.DisableShaderKeyword(GpuParams.DDGI_DEBUG_IRRADIANCE);
                cmd.DisableShaderKeyword(GpuParams.DDGI_DEBUG_DISTANCE);
                cmd.DisableShaderKeyword(GpuParams.DDGI_DEBUG_OFFSET);

                if (ddgiOverride.debugProbe.value)
                {
                    switch (ddgiOverride.probeDebugMode.value)
                    {
                        case ProbeDebugMode.Irradiance:
                            cmd.EnableShaderKeyword(GpuParams.DDGI_DEBUG_IRRADIANCE);
                            break;
                        case ProbeDebugMode.Distance:
                            cmd.EnableShaderKeyword(GpuParams.DDGI_DEBUG_DISTANCE);
                            break;
                        case ProbeDebugMode.RelocationOffset:
                            cmd.EnableShaderKeyword(GpuParams.DDGI_DEBUG_OFFSET);
                            break;
                    }
                }
            }

            // Prepare Rendering Parameters
            {
                var matrix = Matrix4x4.TRS(Vector3.zero, Quaternion.identity, Vector3.one * mddgiOverride.probeRadius.value);
                cmd.SetGlobalMatrix(GpuParams._ddgiSphere_ObjectToWorld, matrix);
                
                cmd.SetGlobalTexture(GpuParams._ProbeData, mDDGIPass.GetProbeData());
            }

            // Construct Indirect Draw Arguments.
            {
                var numProbes = mDDGIPass.GetNumProbes();
                var numProbesFlat = numProbes.x * numProbes.y * numProbes.z; 
                
                uint[] args = new uint[5] { 0, 0, 0, 0, 0 };
                args[0] = (uint)mVisualizeMesh.GetIndexCount(0);
                args[1] = (uint)numProbesFlat;
                args[2] = (uint)mVisualizeMesh.GetIndexStart(0);
                args[3] = (uint)mVisualizeMesh.GetBaseVertex(0);
                
                if(mArgsBuffer != null) { mArgsBuffer.Release(); mArgsBuffer = null; }
                mArgsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
                mArgsBuffer.SetData(args);
            }

            // Draw Spheres.
            {
                cmd.SetRenderTarget(renderer.cameraColorTargetHandle, renderer.cameraDepthTargetHandle);
                // cmd.DrawMeshInstanced has a limit of 1024 instances per pass, so indirect drawing is used instead
                cmd.DrawMeshInstancedIndirect(mVisualizeMesh, 0, mVisualizeMaterial, 0, mArgsBuffer);
            }

            context.ExecuteCommandBuffer(cmd);
            CommandBufferPool.Release(cmd);
        }

        public void Release()
        {
            CoreUtils.Destroy(mVisualizeMaterial);

            if (mArgsBuffer != null) { mArgsBuffer.Release(); mArgsBuffer = null; }
        }
    }

    private DDGIPass mDDGIPass;
    private DDGIVisualizePass mDDGIVisualizePass;

    private Mesh mDDGIVisualizeSphere;
    private bool mIsRayTracingSupported;

    public override void Create()
    {
        if (!isActive)
        {
            mDDGIPass?.Release();
            mDDGIVisualizePass?.Release();
            return;
        }
        
        mIsRayTracingSupported = SystemInfo.supportsRayTracing;
        if (!mIsRayTracingSupported) return;

        mDDGIPass ??= new DDGIPass();
        mDDGIVisualizePass ??= new DDGIVisualizePass();

        mDDGIVisualizeSphere = Resources.Load<Mesh>("Meshes/DDGIVisualizationSphere");

    #if UNITY_EDITOR
        EditorSceneManager.sceneOpened += OnSceneOpened;
    #endif
    }

    public override void AddRenderPasses(ScriptableRenderer renderer, ref RenderingData renderingData)
    {
        if (renderingData.cameraData.isPreviewCamera) return;
        if (!mIsRayTracingSupported) return;
        
        renderer.EnqueuePass(mDDGIPass);

        mDDGIVisualizePass.Setup(mDDGIVisualizeSphere, mDDGIPass);
        renderer.EnqueuePass(mDDGIVisualizePass);
    }

    protected override void Dispose(bool disposing)
    {
        base.Dispose(disposing);

        mDDGIPass?.Release(); 
        mDDGIVisualizePass?.Release();

    #if UNITY_EDITOR
        EditorSceneManager.sceneOpened -= OnSceneOpened;
    #endif
    }

    public void Reinitialize()
    {
        mDDGIPass.Reinitialize();
    }
    
    private void OnSceneOpened(Scene scene, OpenSceneMode mode)
    {
        Reinitialize();
    }
}

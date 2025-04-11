using System;
using System.Collections.Generic;
using System.Linq;
using UnityEditor.SceneManagement;
using UnityEngine;
using UnityEngine.Experimental.Rendering;
using UnityEngine.Rendering;
using UnityEngine.Rendering.RenderGraphModule;
using UnityEngine.Rendering.Universal;
using UnityEngine.SceneManagement;
using ProfilingScope = UnityEngine.Rendering.ProfilingScope;
using Random = System.Random;

namespace DDGI
{
    public sealed class DDGIFeature : ScriptableRendererFeature
    {
        public sealed class DDGIPass : ScriptableRenderPass
        {
            // ----------------------------------------------------
            //              Members and Defines
            // ----------------------------------------------------
            private DDGI _ddgiOverride;
            private DDGICustomBounds _customGIVolume;

            private bool _isInitialized = false;

            private bool _needToResetProbeHistory = true;
            private bool _needToResetProbeRelocation = true;
            private bool _needToResetProbeClassification = true;

            private static readonly int ProbeNumIrradianceInteriorTexels = 6;
            private static readonly int ProbeNumIrradianceTexels = ProbeNumIrradianceInteriorTexels + 2; // Including 1 texel for up border and 1 texel for down border;
            private static readonly int ProbeNumDistanceInteriorTexels = 14;
            private static readonly int ProbeNumDistanceTexels = ProbeNumDistanceInteriorTexels + 2;

            struct DDGIVolumeCpu
            {
                public Vector3 origin;
                public Vector3 extents;
                public Vector3Int numProbes;
                public int maxNumRays;
                public int numRays;
            }
            private DDGIVolumeCpu _ddgiVolumeCpu;

            struct DDGIVolumeGpu
            {
                public Vector4 probeRotation;
                public Vector3 startPosition;
                public int raysPerProbe;
                public Vector3 probeSize;
                public int maxRaysPerProbe;
                public Vector3Int probeCount;
                public float normalBias;
                public Vector3 randomVector;
                public float energyPreservation;

                public float randomAngle;
                public float historyBlendWeight;
                public float indirectIntensity;
                public float normalBiasMultiplier;

                public float viewBiasMultiplier;
                public int ddgiProbeClassification;
                public int ddgiProbeRelocation;
                public float probeFixedRayBackfaceThreshold;

                public float probeMinFrontfaceDistance;
                public int directionalLightCount;
                public int punctualLightCount;
                public int ddgiSkylightMode;

                public Vector4 skyboxTintColor;
                public Vector4 skyColor;
                public Vector4 equatorColor;
                public Vector4 groundColor;
                public Vector4 ambientColor;

                public int ddgiProbeReduction;
                public float skyboxIntensityMultiplier;
                public float skyboxExposure;
                public float padding0;
            }
            private DDGIVolumeGpu _ddgiVolumeGpu;
            private ConstantBuffer<DDGIVolumeGpu> _ddgiVolumeGpuCB;

            #region [Shader Resources]

            private RayTracingShader _ddgiRayTraceShader;
            private RayTracingAccelerationStructure _accelerationStructure;

            //private GraphicsBuffer _rayBuffer;

            private readonly ComputeShader _updateIrradianceCS;
            private readonly int _updateIrradianceKernel;
            private readonly ComputeShader _updateDistanceCS;
            private readonly int _updateDistanceKernel;
            private readonly ComputeShader _probeClassificationCS;
            private readonly int _resetClassificationKernel;
            private readonly int _probeClassificationKernel;
            private readonly ComputeShader _relocateProbeCS;
            private readonly int _resetRelocationKernel;
            private readonly int _relocateProbeKernel;
            private readonly ComputeShader _probeReductionCS;
            private readonly int _reductionKernel;
            private readonly int _extraReductionKernel;

            private readonly Shader _cubemapSkyPS;

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

            #endregion

            #region [Probe Variability]

            private bool _isConverged;
            private readonly uint _minimumVariabilitySamples = 16u;
            private bool _clearProbeVariability;
            private uint _numVolumeVariabilitySamples = 0u;

            #endregion

            #region [Light Update and Change Dectect]

            private GraphicsBuffer _directionalLightBuffer;
            private GraphicsBuffer _punctualLightBuffer;

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
                Cubemap = 0,
                Gradient = 1,
                Color = 2,
                Unsupported = 3
            }

            // For URP's default Skybox parameter IDs, used during the Build Light Structured Buffer process and Probe Variability light comparison
            private static class SkyboxParam
            {
                public static readonly int tint = Shader.PropertyToID("_Tint");
                public static readonly int exposure = Shader.PropertyToID("_Exposure");
                public static readonly int rotation = Shader.PropertyToID("_Rotation");
                public static readonly int tex = Shader.PropertyToID("_Tex");
            }

            // Used only to determine if the Sky Light settings have changed in the latest frame, unrelated to the Build Light Structured Buffer process, only used for Probe Variability
            private class SkyLight
            {
                public SkyLight(Material skybox, AmbientMode ambientMode, float ambientIntensity,
                    Color skyColor, Color equatorColor, Color groundColor)
                {
                    if (skybox != null)
                    {
                        _skyboxTint = skybox.GetColor(SkyboxParam.tint);
                        _skyboxExposure = skybox.GetFloat(SkyboxParam.exposure);
                        _skyboxRotation = skybox.GetFloat(SkyboxParam.rotation);
                        _skyboxTex = skybox.GetTexture(SkyboxParam.tex);
                    }
                    _ambientMode = ambientMode;
                    _ambientIntensity = ambientIntensity;
                    _skyColor = skyColor;
                    _equatorColor = equatorColor;
                    _groundColor = groundColor;
                }

                public bool Equals(SkyLight skyLight)
                {
                    bool result = true;
                    if (_ambientMode == skyLight._ambientMode && _ambientMode == AmbientMode.Skybox)
                    {
                        result &= _skyboxTint == skyLight._skyboxTint;
                        result &= Mathf.Approximately(_skyboxExposure, skyLight._skyboxExposure);
                        result &= Mathf.Approximately(_skyboxRotation, skyLight._skyboxRotation);
                        result &= _skyboxTex == skyLight._skyboxTex;
                    }
                    result &= _ambientMode == skyLight._ambientMode;
                    result &= Mathf.Approximately(_ambientIntensity, skyLight._ambientIntensity);
                    result &= _skyColor == skyLight._skyColor;
                    result &= _equatorColor == skyLight._equatorColor;
                    result &= _groundColor == skyLight._groundColor;
                    return result;
                }

                private Color _skyboxTint = Color.black;
                private float _skyboxExposure = 0.0f;
                private float _skyboxRotation = 0.0f;
                private Texture _skyboxTex = null;
                private AmbientMode _ambientMode;
                private float _ambientIntensity;
                private Color _skyColor;
                private Color _equatorColor;
                private Color _groundColor;
            }

            // Lighting data cache, used for the Probe Variability phase
            private List<DirectionalLight> _cachedDirectionalLights = new List<DirectionalLight>();
            private List<PunctualLight> _cachedPunctualLights = new List<PunctualLight>();
            private SkyLight _cachedSkyLight = new SkyLight(null, AmbientMode.Flat, 0.0f, Color.black, Color.black, Color.black);
            private bool _anyLightChanged;
            private bool _skyChanged;

            #endregion

            public Texture _skyboxTex;
            public RTHandle _skyboxTexHandle;

            public ComputeBuffer _rayBuffer;

            public RenderTexture _probeData;                           // For Probe Relocation
            public RenderTexture _probeDistance;
            public RenderTexture _probeIrradiance;
            public RenderTexture _probeVariability;                    // For Probe Variability
            public RenderTexture _probeDistanceHistory;
            public RenderTexture _probeIrradianceHistory;
            public RenderTexture _probeVariabilityAverage;             // For Probe Variability
            
            public RenderTargetIdentifier _probeDataId;
            public RenderTargetIdentifier _probeDistanceId;
            public RenderTargetIdentifier _probeIrradianceId;
            public RenderTargetIdentifier _probeVariabilityId;
            public RenderTargetIdentifier _probeDistanceHistoryId;
            public RenderTargetIdentifier _probeIrradianceHistoryId;
            public RenderTargetIdentifier _probeVariabilityAverageId;

            public RTHandle _probeDataHandle;
            public RTHandle _probeDistanceHandle;
            public RTHandle _probeIrradianceHandle;
            public RTHandle _probeVariabilityHandle;
            public RTHandle _probeDistanceHistoryHandle;
            public RTHandle _probeIrradianceHistoryHandle;
            public RTHandle _probeVariabilityAverageHandle;

            private static class GpuParams
            {
                // For Probe Tracing and Updating
                public static readonly string RayGenShaderName = "DDGI_RayGen";

                public static readonly int RayBuffer = Shader.PropertyToID("RayBuffer");
                public static readonly int DirectionalLightBuffer = Shader.PropertyToID("DirectionalLightBuffer");
                public static readonly int PunctualLightBuffer = Shader.PropertyToID("PunctualLightBuffer");

                public static readonly int DDGIVolumeGpu = Shader.PropertyToID("DDGIVolumeGpu");

                public static readonly int ProbeIrradiance = Shader.PropertyToID("_ProbeIrradiance");
                public static readonly int ProbeIrradianceHistory = Shader.PropertyToID("_ProbeIrradianceHistory");
                public static readonly int ProbeDistance = Shader.PropertyToID("_ProbeDistance");
                public static readonly int ProbeDistanceHistory = Shader.PropertyToID("_ProbeDistanceHistory");

                public static readonly int AccelerationStructure = Shader.PropertyToID("_AccelerationStructure");

                // For Sky Light Sampling
                public static readonly string DDGISkylightMode = "DDGI_SKYLIGHT_MODE";
                public static readonly int SkyboxCubemap = Shader.PropertyToID("_SkyboxCubemap");

                // For Probe Relocation
                public static readonly string DDGIProbeRelocation = "DDGI_PROBE_RELOCATION";
                public static readonly int ProbeData = Shader.PropertyToID("_ProbeData");

                // For Probe Debug
                public static readonly string DDGIShowIndirectOnly = "DDGI_SHOW_INDIRECT_ONLY";
                public static readonly string DDGIShowPureIndirectRadiance = "DDGI_SHOW_PURE_INDIRECT_RADIANCE";

                // For Probe Reduction (Variability)
                public static readonly int ReductionInputSize = Shader.PropertyToID("_ReductionInputSize");
                public static readonly int ProbeVariability = Shader.PropertyToID("_ProbeVariability");
                public static readonly int ProbeVariabilityAverage = Shader.PropertyToID("_ProbeVariabilityAverage");
            }


            // ----------------------------------------------------
            //           Core Functions and Render Loops
            // ----------------------------------------------------
            public DDGIPass()
            {
                renderPassEvent = RenderPassEvent.BeforeRenderingOpaques;

                _ddgiRayTraceShader = Resources.Load<RayTracingShader>("Shaders/DDGIRayTracing");
                _updateIrradianceCS = Resources.Load<ComputeShader>("Shaders/DDGIUpdateIrradiance");
                _updateIrradianceKernel = _updateIrradianceCS.FindKernel("DDGIUpdateIrradiance");
                _updateDistanceCS = Resources.Load<ComputeShader>("Shaders/DDGIUpdateDistance");
                _updateDistanceKernel = _updateDistanceCS.FindKernel("DDGIUpdateDistance");
                _probeClassificationCS = Resources.Load<ComputeShader>("Shaders/DDGIProbeClassification");
                _resetClassificationKernel = _probeClassificationCS.FindKernel("DDGIProbeClassificationResetCS");
                _probeClassificationKernel = _probeClassificationCS.FindKernel("DDGIProbeClassificationCS");
                _relocateProbeCS = Resources.Load<ComputeShader>("Shaders/DDGIRelocateProbe");
                _resetRelocationKernel = _relocateProbeCS.FindKernel("DDGIResetRelocation");
                _relocateProbeKernel = _relocateProbeCS.FindKernel("DDGIRelocateProbe");
                _probeReductionCS = Resources.Load<ComputeShader>("Shaders/DDGIReduction");
                _reductionKernel = _probeReductionCS.FindKernel("DDGIReductionCS");
                _extraReductionKernel = _probeReductionCS.FindKernel("DDGIExtraReductionCS");

                RayTracingAccelerationStructure.Settings setting = new RayTracingAccelerationStructure.Settings
                    (RayTracingAccelerationStructure.ManagementMode.Automatic, RayTracingAccelerationStructure.RayTracingModeMask.Everything, 255);
                _accelerationStructure = new RayTracingAccelerationStructure(setting);

                _ddgiVolumeGpuCB = new ConstantBuffer<DDGIVolumeGpu>();

                // Shader.Find is not reliable, as shaders may be missing after packaging, making Find ineffective
                // For demonstration purposes, using a quick workaround here
                _cubemapSkyPS = Shader.Find("Skybox/Cubemap");
            }

            private void Initialize(RenderGraph renderGraph)
            {
                void ReleaseIfNotNull(ref RenderTexture texture)
                {
                    if (texture != null)
                    {
                        texture.Release();
                        texture = null;
                    }
                }

                if (_isInitialized || _ddgiOverride == null) return;

                // ---------------------------------------
                // Initialize cpu-side volume parameters
                // ---------------------------------------
                var sceneBoundingBox = GenerateSceneMeshBounds();
                if (sceneBoundingBox.extents == Vector3.zero) return;
                if (sceneBoundingBox.extents == Vector3.zero) return;
                _ddgiVolumeCpu.origin = sceneBoundingBox.center;
                _ddgiVolumeCpu.extents = 1.1f * sceneBoundingBox.extents;
                _ddgiVolumeCpu.numProbes = new Vector3Int(_ddgiOverride.probeCountX.value, _ddgiOverride.probeCountY.value, _ddgiOverride.probeCountZ.value);
                _ddgiVolumeCpu.numRays = _ddgiOverride.raysPerProbe.value;
                _ddgiVolumeCpu.maxNumRays = 512;

                // ---------------------------------------
                // Initialize Ray Data Buffer
                // ---------------------------------------
                if (_rayBuffer != null) { _rayBuffer.Release(); _rayBuffer = null; }
                int numProbesFlat = _ddgiVolumeCpu.numProbes.x * _ddgiVolumeCpu.numProbes.y * _ddgiVolumeCpu.numProbes.z;
                _rayBuffer = new ComputeBuffer(numProbesFlat * _ddgiVolumeCpu.maxNumRays, 16 /* float4 */, ComputeBufferType.Default);

                // Note: Try to use GraphicsFormat to provide explicit identification of floating-point / fixed-point numbers
                // For example, Distance Texture, previously used RenderTextureFormat.RG32, which uses 16-bit unsigned fixed-point numbers, but distances need to be floating-point
                // Allocating RG32 as Distance Texture would ignore the decimal places of the distance, causing Chebyshev visibility tests to have Edge Clamp Artifacts.
                // ---------------------------------------
                // Radiance and Distance Texture2DArray
                // ---------------------------------------
                ReleaseIfNotNull(ref _probeIrradiance);
                var probeIrradianceDimensions = GetDDGIVolumeTextureDimensions(_ddgiVolumeCpu, DDGIVolumeTextureType.Irradiance);
                _probeIrradiance = new RenderTexture(probeIrradianceDimensions.x, probeIrradianceDimensions.y, 0, GraphicsFormat.R16G16B16A16_SFloat)
                {
                    filterMode = FilterMode.Bilinear,
                    useMipMap = false,
                    autoGenerateMips = false,
                    enableRandomWrite = true,
                    name = "DDGI Probe Irradiance",
                    dimension = TextureDimension.Tex2DArray,
                    volumeDepth = probeIrradianceDimensions.z
                };
                if (_probeIrradiance.Create())
                {
                    _probeIrradianceId = new RenderTargetIdentifier(_probeIrradiance);
                    _probeIrradianceHandle = RTHandles.Alloc(_probeIrradiance);
                }

                ReleaseIfNotNull(ref _probeDistance);
                var probeDistanceDimensions = GetDDGIVolumeTextureDimensions(_ddgiVolumeCpu, DDGIVolumeTextureType.Distance);
                _probeDistance = new RenderTexture(probeDistanceDimensions.x, probeDistanceDimensions.y, 0, GraphicsFormat.R16G16_SFloat)
                {
                    filterMode = FilterMode.Bilinear,
                    useMipMap = false,
                    autoGenerateMips = false,
                    enableRandomWrite = true,
                    name = "DDGI Probe Distance",
                    dimension = TextureDimension.Tex2DArray,
                    volumeDepth = probeDistanceDimensions.z
                };
                if (_probeDistance.Create())
                {
                    _probeDistanceId = new RenderTargetIdentifier(_probeDistance);
                    _probeDistanceHandle = RTHandles.Alloc(_probeDistance);
                }

                ReleaseIfNotNull(ref _probeIrradianceHistory);
                _probeIrradianceHistory = new RenderTexture(_probeIrradiance.descriptor)
                {
                    name = "DDGI Probe Irradiance History"
                };
                if (_probeIrradianceHistory.Create())
                {
                    _probeIrradianceHistoryId = new RenderTargetIdentifier(_probeIrradianceHistory);
                    _probeIrradianceHistoryHandle = RTHandles.Alloc(_probeIrradianceHistory);
                }

                ReleaseIfNotNull(ref _probeDistanceHistory);
                _probeDistanceHistory = new RenderTexture(_probeDistance.descriptor)
                {
                    name = "DDGI Probe Distance History"
                };
                if (_probeDistanceHistory.Create())
                {
                    _probeDistanceHistoryId = new RenderTargetIdentifier(_probeDistanceHistory);
                    _probeDistanceHistoryHandle = RTHandles.Alloc(_probeDistanceHistory);
                }

                // ---------------------------------------
                // Create Probe Data
                // ---------------------------------------
                ReleaseIfNotNull(ref _probeData);
                var probeDataDimensions = GetDDGIVolumeTextureDimensions(_ddgiVolumeCpu, DDGIVolumeTextureType.ProbeData);
                _probeData = new RenderTexture(probeDataDimensions.x, probeDataDimensions.y, 0, GraphicsFormat.R16G16B16A16_SFloat)
                {
                    filterMode = FilterMode.Bilinear,
                    useMipMap = false,
                    autoGenerateMips = false,
                    enableRandomWrite = true,
                    name = "DDGI Probe Data",
                    dimension = TextureDimension.Tex2DArray,
                    volumeDepth = probeDataDimensions.z
                };
                if (_probeData.Create())
                {
                    _probeDataId = new RenderTargetIdentifier(_probeData);
                    _probeDataHandle = RTHandles.Alloc(_probeData);
                }

                // ---------------------------------------
                // Create Probe Variability
                // ---------------------------------------
                ReleaseIfNotNull(ref _probeVariability);
                var probeVariabilityDimensions = GetDDGIVolumeTextureDimensions(_ddgiVolumeCpu, DDGIVolumeTextureType.Variability);
                _probeVariability = new RenderTexture(probeVariabilityDimensions.x, probeVariabilityDimensions.y, 0, GraphicsFormat.R32_SFloat)
                {
                    filterMode = FilterMode.Bilinear,
                    useMipMap = false,
                    autoGenerateMips = false,
                    enableRandomWrite = true,
                    name = "DDGI Probe Variability",
                    dimension = TextureDimension.Tex2DArray,
                    volumeDepth = probeVariabilityDimensions.z
                };
                if (_probeVariability.Create())
                {
                    _probeVariabilityId = new RenderTargetIdentifier(_probeVariability);
                    _probeVariabilityHandle = RTHandles.Alloc(_probeVariability);
                }

                ReleaseIfNotNull(ref _probeVariabilityAverage);
                var probeVariabilityAverageDimensions = GetDDGIVolumeTextureDimensions(_ddgiVolumeCpu, DDGIVolumeTextureType.VariabilityAverage);
                _probeVariabilityAverage = new RenderTexture(_probeVariability.descriptor)
                {
                    graphicsFormat = GraphicsFormat.R32G32_SFloat,
                    width = probeVariabilityAverageDimensions.x,
                    height = probeVariabilityAverageDimensions.y,
                    volumeDepth = probeVariabilityAverageDimensions.z,
                    name = "DDGI Probe Variability Average"
                };
                if (_probeVariabilityAverage.Create())
                {
                    _probeVariabilityAverageId = new RenderTargetIdentifier(_probeVariabilityAverage);
                    _probeVariabilityAverageHandle = RTHandles.Alloc(_probeVariabilityAverage);
                }

                _isInitialized = true;
            }

            public void Reinitialize()
            {
                _isInitialized = false;
                _needToResetProbeHistory = true;
                _needToResetProbeRelocation = true;
                _needToResetProbeClassification = true;
                _clearProbeVariability = true;
                _isConverged = false;
            }

            public class DDGIPassData
            {
                internal Camera camera;

                internal TextureHandle skyboxTex;

                internal TextureHandle probeDataHandle;
                internal TextureHandle probeDistanceHandle;
                internal TextureHandle probeIrradianceHandle;
                internal TextureHandle probeVariabilityHandle;
                internal TextureHandle probeDistanceHistoryHandle;
                internal TextureHandle probeIrradianceHistoryHandle;
                internal TextureHandle probeVariabilityAverageHandle;

                internal BufferHandle directionalBufferHandle;
                internal BufferHandle punctualBufferHandle;
            }

            public static RenderTargetInfo GetRenderTargetInfo(Texture texture)
            {
                return new RenderTargetInfo
                {
                    width = texture.width,
                    height = texture.height,
                    format = texture.graphicsFormat,
                    volumeDepth = 1,
                };
            }

            public static RenderTargetInfo GetRenderTargetInfo(RenderTexture texture)
            {
                return new RenderTargetInfo
                {
                    width = texture.width,
                    height = texture.height,
                    format = texture.graphicsFormat,
                    volumeDepth = texture.volumeDepth,
                };
            }

            public override void RecordRenderGraph(RenderGraph renderGraph, ContextContainer frameData)
            {
                _ddgiOverride = VolumeManager.instance.stack.GetComponent<DDGI>();
                if (_ddgiOverride == null || !_ddgiOverride.IsActive()) return;

                Initialize(renderGraph);
                if (!_isInitialized || _accelerationStructure == null) return;

                UpdateSceneLights();

                var skyboxTexHandle                 = renderGraph.ImportTexture(_skyboxTexHandle, GetRenderTargetInfo(_skyboxTex));
                var probeDataHandle                 = renderGraph.ImportTexture(_probeDataHandle, GetRenderTargetInfo(_probeData));
                var probeDistanceHandle             = renderGraph.ImportTexture(_probeDistanceHandle, GetRenderTargetInfo(_probeDistance));
                var probeIrradianceHandle           = renderGraph.ImportTexture(_probeIrradianceHandle, GetRenderTargetInfo(_probeIrradiance));
                var probeVariabilityHandle          = renderGraph.ImportTexture(_probeVariabilityHandle, GetRenderTargetInfo(_probeVariability));
                var probeDistanceHistoryHandle      = renderGraph.ImportTexture(_probeDistanceHistoryHandle, GetRenderTargetInfo(_probeDistanceHistory));
                var probeIrradianceHistoryHandle    = renderGraph.ImportTexture(_probeIrradianceHistoryHandle, GetRenderTargetInfo(_probeIrradianceHistory));
                var probeVariabilityAverageHandle   = renderGraph.ImportTexture(_probeVariabilityAverageHandle, GetRenderTargetInfo(_probeVariabilityAverage));

                var punctualBufferHandle            = renderGraph.ImportBuffer(_punctualLightBuffer);
                var directionalBufferHandle         = renderGraph.ImportBuffer(_directionalLightBuffer);

                using (var builder = renderGraph.AddUnsafePass<DDGIPassData>("DDGI Pass", out var passData))
                {
                    UniversalCameraData cameraData = frameData.Get<UniversalCameraData>();
                    UniversalResourceData resourceData = frameData.Get<UniversalResourceData>();

                    passData.camera = cameraData.camera;

                    builder.EnableAsyncCompute(false);

                    builder.UseTexture(passData.skyboxTex = skyboxTexHandle);
                    builder.UseTexture(passData.probeDataHandle = probeDataHandle, AccessFlags.ReadWrite);
                    builder.UseTexture(passData.probeDistanceHandle = probeDistanceHandle, AccessFlags.ReadWrite);
                    builder.UseTexture(passData.probeIrradianceHandle = probeIrradianceHandle, AccessFlags.ReadWrite);
                    builder.UseTexture(passData.probeVariabilityHandle = probeVariabilityHandle, AccessFlags.ReadWrite);
                    builder.UseTexture(passData.probeDistanceHistoryHandle = probeDistanceHistoryHandle, AccessFlags.ReadWrite);
                    builder.UseTexture(passData.probeIrradianceHistoryHandle = probeIrradianceHistoryHandle, AccessFlags.ReadWrite);
                    builder.UseTexture(passData.probeVariabilityAverageHandle = probeVariabilityAverageHandle, AccessFlags.ReadWrite);

                    builder.UseBuffer(passData.punctualBufferHandle = punctualBufferHandle, AccessFlags.ReadWrite);
                    builder.UseBuffer(passData.directionalBufferHandle = directionalBufferHandle, AccessFlags.ReadWrite);

                    builder.SetRenderFunc((DDGIPassData data, UnsafeGraphContext context) => ExecuteUnsafePass(data, context));
                }
            }

            public void ExecuteUnsafePass(DDGIPassData data, UnsafeGraphContext frameContext)
            {
                var camera = data.camera;

                ResetHistoryInfoIfNeeded(data, frameContext.cmd);

                if (_ddgiVolumeGpu.ddgiSkylightMode == (int)SkyLightMode.Cubemap)
                {
                    frameContext.cmd.SetRayTracingTextureParam(_ddgiRayTraceShader, GpuParams.SkyboxCubemap, data.skyboxTex);
                }

                // Note: Each time this function is called, the random numbers are updated, which causes _RandomVector and _RandomAngle to change
                // If the random number update logic is not separated, this function can only be called once per frame!
                PushGpuConstants(frameContext.cmd);

                int numProbesFlat = _ddgiVolumeCpu.numProbes.x * _ddgiVolumeCpu.numProbes.y * _ddgiVolumeCpu.numProbes.z;
            
                using (new ProfilingScope(frameContext.cmd, new ProfilingSampler("DDGI Ray Trace Pass")))
                {
                    if (!_isConverged)
                    {
                        frameContext.cmd.BuildRayTracingAccelerationStructure(_accelerationStructure);
                        frameContext.cmd.SetRayTracingAccelerationStructure(_ddgiRayTraceShader, GpuParams.AccelerationStructure, _accelerationStructure);

                        CommandBufferHelpers.GetNativeCommandBuffer(frameContext.cmd).SetRayTracingShaderPass(_ddgiRayTraceShader, "DDGIRayTracing");
                        
                        frameContext.cmd.SetGlobalTexture(GpuParams.ProbeIrradianceHistory, data.probeIrradianceHistoryHandle);
                        frameContext.cmd.SetGlobalTexture(GpuParams.ProbeDistanceHistory, data.probeDistanceHistoryHandle);
                        frameContext.cmd.SetGlobalTexture(GpuParams.ProbeData, data.probeDataHandle);

                        frameContext.cmd.SetRayTracingBufferParam(_ddgiRayTraceShader, GpuParams.RayBuffer, _rayBuffer);
                        frameContext.cmd.SetGlobalBuffer(GpuParams.DirectionalLightBuffer, data.directionalBufferHandle);     // We will use it in closest hit shader, not in actual .raytrace shader
                        frameContext.cmd.SetGlobalBuffer(GpuParams.PunctualLightBuffer, data.punctualBufferHandle);           // We will use it in closest hit shader, not in actual .raytrace shader
                    
                        frameContext.cmd.DispatchRays(_ddgiRayTraceShader, GpuParams.RayGenShaderName, (uint)_ddgiVolumeCpu.numRays, (uint)numProbesFlat, 1, camera);
                    }
                }

                using (new ProfilingScope(frameContext.cmd, new ProfilingSampler("DDGI Update Irradiance Pass")))
                {
                    if (!_isConverged)
                    {
                        frameContext.cmd.SetComputeBufferParam(_updateIrradianceCS, _updateIrradianceKernel, GpuParams.RayBuffer, _rayBuffer);
                        frameContext.cmd.SetComputeTextureParam(_updateIrradianceCS, _updateIrradianceKernel, GpuParams.ProbeIrradiance, data.probeIrradianceHandle);
                        frameContext.cmd.SetComputeTextureParam(_updateIrradianceCS, _updateIrradianceKernel, GpuParams.ProbeIrradianceHistory, data.probeIrradianceHistoryHandle);
                        frameContext.cmd.SetComputeTextureParam(_updateIrradianceCS, _updateIrradianceKernel, GpuParams.ProbeVariability, data.probeVariabilityHandle);

                        // Note that we use Y-UP, so the Dispatch needs to be reversed
                        frameContext.cmd.DispatchCompute(_updateIrradianceCS, _updateIrradianceKernel, _ddgiVolumeCpu.numProbes.x, _ddgiVolumeCpu.numProbes.z, _ddgiVolumeCpu.numProbes.y);
                    }
                }

                using (new ProfilingScope(frameContext.cmd, new ProfilingSampler("DDGI Update Distance Pass")))
                {
                    if (!_isConverged)
                    {
                        frameContext.cmd.SetComputeBufferParam(_updateDistanceCS, _updateDistanceKernel, GpuParams.RayBuffer, _rayBuffer);
                        frameContext.cmd.SetComputeTextureParam(_updateDistanceCS, _updateDistanceKernel, GpuParams.ProbeDistance, data.probeDistanceHandle);
                        frameContext.cmd.SetComputeTextureParam(_updateDistanceCS, _updateDistanceKernel, GpuParams.ProbeDistanceHistory, data.probeDistanceHistoryHandle);

                        // Note that we use Y-UP, so the Dispatch needs to be reversed
                        frameContext.cmd.DispatchCompute(_updateDistanceCS, _updateDistanceKernel, _ddgiVolumeCpu.numProbes.x, _ddgiVolumeCpu.numProbes.z, _ddgiVolumeCpu.numProbes.y);
                    }
                }


                if (_ddgiOverride.enableProbeRelocation.value)
                {
                    using (new ProfilingScope(frameContext.cmd, new ProfilingSampler("DDGI Relocate Probe Pass")))
                    {
                        var numGroupsX = Mathf.CeilToInt(numProbesFlat / 32.0f /*relocationGroupSizeX*/);
                    
                        if (_needToResetProbeRelocation)
                        {
                            frameContext.cmd.SetComputeTextureParam(_relocateProbeCS, _resetRelocationKernel, GpuParams.ProbeData, data.probeDataHandle);
                            frameContext.cmd.DispatchCompute(_relocateProbeCS, _resetRelocationKernel, numGroupsX, 1, 1);
                            _needToResetProbeRelocation = false;
                        }
                    
                        frameContext.cmd.SetComputeTextureParam(_relocateProbeCS, _relocateProbeKernel, GpuParams.ProbeData, data.probeDataHandle);
                        frameContext.cmd.SetComputeBufferParam(_relocateProbeCS, _relocateProbeKernel, GpuParams.RayBuffer, _rayBuffer);
                        frameContext.cmd.DispatchCompute(_relocateProbeCS, _relocateProbeKernel, numGroupsX, 1, 1);
                    }
                }
                else
                {
                    if (!_needToResetProbeRelocation)
                    {
                        var numGroupsX = Mathf.CeilToInt(numProbesFlat / 32.0f /*relocationGroupSizeX*/);
                    
                        frameContext.cmd.SetComputeTextureParam(_relocateProbeCS, _resetRelocationKernel, GpuParams.ProbeData, data.probeDataHandle);
                        frameContext.cmd.DispatchCompute(_relocateProbeCS, _resetRelocationKernel, numGroupsX, 1, 1);
                        _needToResetProbeRelocation = true;
                    }
                }

                if (_ddgiOverride.enableProbeClassification.value)
                {
                    using (new ProfilingScope(frameContext.cmd, new ProfilingSampler("DDGI Classify Probe Pass")))
                    {
                        var numGroupsX = Mathf.CeilToInt(numProbesFlat / 32.0f /*relocationGroupSizeX*/);

                        if (_needToResetProbeClassification)
                        {
                            frameContext.cmd.SetComputeTextureParam(_probeClassificationCS, _resetClassificationKernel, GpuParams.ProbeData, data.probeDataHandle);
                            frameContext.cmd.DispatchCompute(_probeClassificationCS, _resetClassificationKernel, numGroupsX, 1, 1);
                            _needToResetProbeClassification = false;
                        }
                    
                        frameContext.cmd.SetComputeTextureParam(_probeClassificationCS, _probeClassificationKernel, GpuParams.ProbeData, data.probeDataHandle);
                        frameContext.cmd.SetComputeBufferParam(_probeClassificationCS, _probeClassificationKernel, GpuParams.RayBuffer, _rayBuffer);
                        frameContext.cmd.DispatchCompute(_probeClassificationCS, _probeClassificationKernel, numGroupsX, 1, 1);
                    }
                }
                else
                {
                    if (!_needToResetProbeClassification)
                    {
                        var numGroupsX = Mathf.CeilToInt(numProbesFlat / 32.0f /*relocationGroupSizeX*/);
                    
                        frameContext.cmd.SetComputeTextureParam(_probeClassificationCS, _resetClassificationKernel, GpuParams.ProbeData, data.probeDataHandle);
                        frameContext.cmd.DispatchCompute(_probeClassificationCS, _resetClassificationKernel, numGroupsX, 1, 1);
                        _needToResetProbeClassification = true;
                    }
                }
            
                if (_ddgiOverride.enableProbeVariability.value)
                {
                    using (new ProfilingScope(frameContext.cmd, new ProfilingSampler("DDGI Variability Pass")))
                    {
                        // TODO: Y-UP Probe Volume hardcoding, if you need to modify the Volume axis, branches will be required
                        var inputTexels = new Vector3Int(_ddgiVolumeCpu.numProbes.x * ProbeNumIrradianceInteriorTexels,
                            _ddgiVolumeCpu.numProbes.z * ProbeNumIrradianceInteriorTexels,
                            _ddgiVolumeCpu.numProbes.y);
                        var NumThreadsInGroup = new Vector3Int(4, 8, 4);
                        var ThreadSampleFootprint = new Vector2Int(4, 2);
                    
                        // -------------------------
                        // First Reduction Pass
                        // -------------------------
                        {
                            frameContext.cmd.SetComputeTextureParam(_probeReductionCS, _reductionKernel, GpuParams.ProbeVariability, data.probeVariabilityHandle);
                            frameContext.cmd.SetComputeTextureParam(_probeReductionCS, _reductionKernel, GpuParams.ProbeVariabilityAverage, data.probeVariabilityAverageHandle);
                            frameContext.cmd.SetComputeVectorParam(_probeReductionCS, GpuParams.ReductionInputSize, new Vector4(inputTexels.x, inputTexels.y, inputTexels.z, 0.0f));

                            var outputTexelsX = Mathf.CeilToInt((float)inputTexels.x / (float)(NumThreadsInGroup.x * ThreadSampleFootprint.x));
                            var outputTexelsY = Mathf.CeilToInt((float)inputTexels.y / (float)(NumThreadsInGroup.y * ThreadSampleFootprint.y));
                            var outputTexelsZ = Mathf.CeilToInt((float)inputTexels.z / (float)(NumThreadsInGroup.z));
                    
                            frameContext.cmd.DispatchCompute(_probeReductionCS, _reductionKernel, outputTexelsX, outputTexelsY, outputTexelsZ);

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
                            
                                frameContext.cmd.SetComputeTextureParam(_probeReductionCS, _extraReductionKernel, GpuParams.ProbeVariabilityAverage, data.probeVariabilityAverageHandle);
                                frameContext.cmd.SetComputeVectorParam(_probeReductionCS, GpuParams.ReductionInputSize, new Vector4(inputTexels.x, inputTexels.y, inputTexels.z, 0.0f));
                            
                                frameContext.cmd.DispatchCompute(_probeReductionCS, _extraReductionKernel, outputTexelsX, outputTexelsY, outputTexelsZ);
                            
                                inputTexels = new Vector3Int(outputTexelsX, outputTexelsY, outputTexelsZ);
                            }
                        }
                    
                        // ---------------------------------
                        // Readback From Variability Average
                        // ---------------------------------
                        // Grab First Pixel of Variability Average
                        AsyncGPUReadback.Request(data.probeVariabilityAverageHandle, 0, 0, 1, 0, 1, 0, 1, VariabilityEstimate);
                    }
                }
                else
                {
                    // If the variability feature is not enabled, we assume that the integration process will never converge
                    _isConverged = false;
                    _clearProbeVariability = true;
                    _numVolumeVariabilitySamples = 0u;
                }

                using (new ProfilingScope(frameContext.cmd, new ProfilingSampler("DDGI Copy Irradiance/Distance Pass")))
                {
                    CommandBufferHelpers.GetNativeCommandBuffer(frameContext.cmd).CopyTexture(_probeIrradianceId, _probeIrradianceHistoryId);
                    CommandBufferHelpers.GetNativeCommandBuffer(frameContext.cmd).CopyTexture(_probeDistanceId, _probeDistanceHistoryId);
                }
            }

            public void Release()
            {
                if (_accelerationStructure != null) { _accelerationStructure.Release(); _accelerationStructure = null; }
                if (_ddgiVolumeGpuCB != null) { _ddgiVolumeGpuCB.Release(); _ddgiVolumeGpuCB = null; }

                ReleaseBufferSafe(_punctualLightBuffer);
                ReleaseBufferSafe(_directionalLightBuffer);

                _rayBuffer?.Release();

                _probeData?.Release();
                _probeDistance?.Release();
                _probeIrradiance?.Release();
                _probeVariability?.Release();
                _probeDistanceHistory?.Release();
                _probeIrradianceHistory?.Release();
                _probeVariabilityAverage?.Release();

                RTHandles.Release(_skyboxTexHandle);
                RTHandles.Release(_probeDataHandle);
                RTHandles.Release(_probeDistanceHandle);
                RTHandles.Release(_probeIrradianceHandle);
                RTHandles.Release(_probeVariabilityHandle);
                RTHandles.Release(_probeDistanceHistoryHandle);
                RTHandles.Release(_probeIrradianceHistoryHandle);
                RTHandles.Release(_probeVariabilityAverageHandle);

                _isInitialized = false;
            }
        
        
            // ----------------------------------------------------
            //           Wrapper and Utility Functions
            // ----------------------------------------------------
            public Vector3Int GetNumProbes() => _ddgiVolumeCpu.numProbes; // For Visualization Pass

            private void ResetHistoryInfoIfNeeded(DDGIPassData data, UnsafeCommandBuffer cmd)
            {
                if (_needToResetProbeHistory)
                {
                    cmd.SetRenderTarget(data.probeIrradianceHistoryHandle);
                    cmd.ClearRenderTarget(false, true, Color.clear);

                    cmd.SetRenderTarget(data.probeDistanceHistoryHandle);
                    cmd.ClearRenderTarget(false, true, Color.clear);

                    _needToResetProbeHistory = false;
                }
            }
        
            private void PushGpuConstants(UnsafeCommandBuffer cmd)
            {
                var random = (float)NextDouble(new Random(), 0.0f, 1.0f, 5); // Generate a random number between 0 and 1, with 5 decimal places
                var randomVec = Vector3.Normalize(new Vector3(2.0f * random - 1.0f, 2.0f * random - 1.0f, 2.0f * random - 1.0f));
                var randomAngle = random * Mathf.PI * 2.0f;

                // -------------------------------------------------
                // Fill GPU constants (lighting-related constants are updated in UpdateSceneLights)
                // -------------------------------------------------
                Quaternion rotation;
                if (_ddgiOverride.useCustomBounds.value && _customGIVolume != null) { rotation = _customGIVolume.transform.rotation; }
                else { rotation = Quaternion.Euler(_ddgiOverride.probeRotationDegrees.value); }
                _ddgiVolumeGpu.probeRotation = new Vector4(rotation.x, rotation.y, rotation.z, rotation.w);
                _ddgiVolumeGpu.startPosition = _ddgiVolumeCpu.origin - _ddgiVolumeCpu.extents;
                _ddgiVolumeGpu.raysPerProbe = _ddgiVolumeCpu.numRays;
                var a = 2.0f * _ddgiVolumeCpu.extents;
                var b = new Vector3(_ddgiVolumeCpu.numProbes.x, _ddgiVolumeCpu.numProbes.y, _ddgiVolumeCpu.numProbes.z) - Vector3.one;
                _ddgiVolumeGpu.probeSize = new Vector3(a.x / b.x, a.y / b.y, a.z / b.z);
                _ddgiVolumeGpu.maxRaysPerProbe = _ddgiVolumeCpu.maxNumRays;
                _ddgiVolumeGpu.probeCount = new Vector3Int(_ddgiVolumeCpu.numProbes.x, _ddgiVolumeCpu.numProbes.y, _ddgiVolumeCpu.numProbes.z);
                _ddgiVolumeGpu.normalBias = 0.25f;
                _ddgiVolumeGpu.randomVector = randomVec;
                _ddgiVolumeGpu.energyPreservation = 0.85f;
                _ddgiVolumeGpu.randomAngle = randomAngle;
                _ddgiVolumeGpu.historyBlendWeight = 0.98f;
                _ddgiVolumeGpu.indirectIntensity = _ddgiOverride.indirectIntensity.value;
                _ddgiVolumeGpu.normalBiasMultiplier = _ddgiOverride.normalBiasMultiplier.value;
                _ddgiVolumeGpu.viewBiasMultiplier = _ddgiOverride.viewBiasMultiplier.value;
                _ddgiVolumeGpu.ddgiProbeClassification = _ddgiOverride.enableProbeClassification.value ? 1 : 0;
                _ddgiVolumeGpu.ddgiProbeRelocation = _ddgiOverride.enableProbeRelocation.value ? 1 : 0;
                _ddgiVolumeGpu.probeFixedRayBackfaceThreshold = _ddgiOverride.probeFixedRayBackfaceThreshold.value;
                _ddgiVolumeGpu.probeMinFrontfaceDistance = _ddgiOverride.probeMinFrontfaceDistance.value;
                _ddgiVolumeGpu.ddgiProbeReduction = _ddgiOverride.enableProbeVariability.value ? 1 : 0;
                _ddgiVolumeGpu.padding0 = 0.0f;

                _ddgiVolumeGpuCB.PushGlobal(CommandBufferHelpers.GetNativeCommandBuffer(cmd), _ddgiVolumeGpu, GpuParams.DDGIVolumeGpu); 

                // -------------------------------------------------
                // Shader Keywords.
                // -------------------------------------------------
                cmd.DisableShaderKeyword(GpuParams.DDGIShowIndirectOnly);
                cmd.DisableShaderKeyword(GpuParams.DDGIShowPureIndirectRadiance);
                if (_ddgiOverride.debugIndirect.value)
                {
                    switch (_ddgiOverride.indirectDebugMode.value)
                    {
                        case IndirectDebugMode.FullIndirectRadiance:
                            cmd.EnableShaderKeyword(GpuParams.DDGIShowIndirectOnly);
                            break;
                        case IndirectDebugMode.PureIndirectRadiance:
                            cmd.EnableShaderKeyword(GpuParams.DDGIShowPureIndirectRadiance);
                            break;
                        default:
                            break;
                    }
                }
            }

            #region [Light Update]

            // Update all lights and monitor changes in light data
            private void UpdateSceneLights()
            {
                BuildLightStructuredBuffer();
                UpdateSkyLight();
                _clearProbeVariability = _anyLightChanged || _skyChanged;
            }

            // Unity by default culls extra lights in the scene, which affects our ability to collect global lighting information for the scene.
            // Therefore, we have to manually collect this on the CPU side.
            private void BuildLightStructuredBuffer()
            {
                var cpuLights = FindObjectsByType<Light>(FindObjectsInactive.Exclude, FindObjectsSortMode.None);

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
                ReleaseBufferSafe(_directionalLightBuffer);
                _directionalLightBuffer = new GraphicsBuffer(GraphicsBuffer.Target.Structured, Mathf.Max(gpuDirectionalLights.Count, 1), 2 * 16);

                ReleaseBufferSafe(_punctualLightBuffer);
                _punctualLightBuffer = new GraphicsBuffer(GraphicsBuffer.Target.Structured, Mathf.Max(gpuPunctualLights.Count, 1), 4 * 16);

                _directionalLightBuffer.SetData(gpuDirectionalLights.ToArray());
                _punctualLightBuffer.SetData(gpuPunctualLights.ToArray());

                _ddgiVolumeGpu.directionalLightCount = gpuDirectionalLights.Count;
                _ddgiVolumeGpu.punctualLightCount = gpuPunctualLights.Count;

                // -----------------------------
                // Any Light Changed Determination
                // -----------------------------
                _anyLightChanged = (!_cachedDirectionalLights.SequenceEqual(gpuDirectionalLights)) || (!_cachedPunctualLights.SequenceEqual(gpuPunctualLights));
                if (_anyLightChanged)
                {
                    _cachedDirectionalLights = new List<DirectionalLight>(gpuDirectionalLights);
                    _cachedPunctualLights = new List<PunctualLight>(gpuPunctualLights);
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
            private void UpdateSkyLight()
            {
                // -----------------------------
                // Sky Light Changed Determination
                // -----------------------------
                _skyboxTex = RenderSettings.skybox.GetTexture(SkyboxParam.tex);
                if (_skyboxTexHandle == null) _skyboxTexHandle = RTHandles.Alloc(_skyboxTex);

                var currSkyLight = new SkyLight(RenderSettings.skybox, RenderSettings.ambientMode, RenderSettings.ambientIntensity,
                    RenderSettings.ambientSkyColor, RenderSettings.ambientEquatorColor, RenderSettings.ambientGroundColor);

                _skyChanged = !_cachedSkyLight.Equals(currSkyLight);
                if (_skyChanged) { _cachedSkyLight = currSkyLight; }
            
                switch (RenderSettings.ambientMode)
                {
                    case AmbientMode.Skybox:
                    {
                        var skybox = RenderSettings.skybox;
                        if (skybox == null)
                        {
                            _ddgiVolumeGpu.ddgiSkylightMode = (int)SkyLightMode.Color;
                            _ddgiVolumeGpu.ambientColor = RenderSettings.ambientSkyColor;
                            return;
                        }

                        if (_cubemapSkyPS == null)
                        {
                            Debug.LogWarning("DDGIFeature failed to find the built-in skybox shader in URP, please check.");
                            _ddgiVolumeGpu.ddgiSkylightMode = (int)SkyLightMode.Unsupported;
                            return;
                        }

                        if (skybox.shader == _cubemapSkyPS)
                        {
                            _ddgiVolumeGpu.ddgiSkylightMode = (int)SkyLightMode.Cubemap;
                            _ddgiVolumeGpu.skyboxIntensityMultiplier = RenderSettings.ambientIntensity;
                            _ddgiVolumeGpu.skyboxTintColor = skybox.GetColor(SkyboxParam.tint);
                            _ddgiVolumeGpu.skyboxExposure = skybox.GetFloat(SkyboxParam.exposure);
                        }
                        else
                        {
                            // Currently, we only support the most commonly used Cubemap-style skyboxes. Other types of skyboxes are not supported and will fallback to pure black.
                            _ddgiVolumeGpu.ddgiSkylightMode = (int)SkyLightMode.Unsupported;
                        }
                    }
                        break;
                    case AmbientMode.Trilight:
                    {
                        _ddgiVolumeGpu.ddgiSkylightMode = (int)SkyLightMode.Gradient;
                        _ddgiVolumeGpu.skyColor = RenderSettings.ambientSkyColor;
                        _ddgiVolumeGpu.equatorColor = RenderSettings.ambientEquatorColor;
                        _ddgiVolumeGpu.groundColor = RenderSettings.ambientGroundColor;
                    }
                        break;
                    case AmbientMode.Flat:
                    {
                        _ddgiVolumeGpu.ddgiSkylightMode = (int)SkyLightMode.Color;
                        _ddgiVolumeGpu.ambientColor = RenderSettings.ambientSkyColor;
                    }
                        break;
                }
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

                        if (_clearProbeVariability) _numVolumeVariabilitySamples = 0;
                    
                        _isConverged = (_numVolumeVariabilitySamples++ > _minimumVariabilitySamples) &&
                                       (volumeAverageVariability < _ddgiOverride.probeVariabilityThreshold.value);
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

                if (_ddgiOverride != null && _ddgiOverride.useCustomBounds.value)
                {
                    // Currently, only a single custom bounding box is supported
                    _customGIVolume = FindFirstObjectByType<DDGICustomBounds>();
                    var boxCollider = _customGIVolume.GetComponent<BoxCollider>();
                    if (boxCollider != null) bounds = boxCollider.bounds;
                }
                else
                {
                    // Automatically generate the bounding box based on the scene's meshes
                    foreach (var meshRenderer in FindObjectsByType<MeshRenderer>(FindObjectsSortMode.None))
                    {
                        bounds.Encapsulate(meshRenderer.bounds);
                    }

                    // In theory, we won't update the bounding box frame by frame, so it's not necessary to force include skinned meshes. Removing the following part is also fine.
                    foreach (var skinnedMeshRenderer in FindObjectsByType<SkinnedMeshRenderer>(FindObjectsSortMode.None))
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
                var width = volumeDescCpu.numProbes.x;
                var height = volumeDescCpu.numProbes.z;
                var arraySize = volumeDescCpu.numProbes.y;

                // Since ProbeData is one-by-one, there's no need for additional branching, just return the code above directly
                switch (type)
                {
                    case DDGIVolumeTextureType.RayData:
                    {
                        // Each row represents all ray info for a single probe; the number of rows is the total number of probes in a plane
                        height = width * height;
                        width = volumeDescCpu.numRays;
                        break;
                    }
                    case DDGIVolumeTextureType.Irradiance:
                    {
                        width *= ProbeNumIrradianceTexels;
                        height *= ProbeNumIrradianceTexels;
                        break;    
                    }
                    case DDGIVolumeTextureType.Distance:
                    {
                        width *= ProbeNumDistanceTexels;
                        height *= ProbeNumDistanceTexels;
                        break;
                    }
                    case DDGIVolumeTextureType.Variability:
                    {
                        width *= ProbeNumIrradianceInteriorTexels;
                        height *= ProbeNumIrradianceInteriorTexels;
                        break;
                    }
                    case DDGIVolumeTextureType.VariabilityAverage:
                    {
                        width *= ProbeNumIrradianceInteriorTexels;
                        height *= ProbeNumIrradianceInteriorTexels;
                    
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

                //double randNum = ran.NextDouble() * (maxValue - minValue) + minValue;
                //double multiplier = Math.Pow(10, decimalPlace);
                //return Math.Truncate(randNum * multiplier) / multiplier;
            }
        }

        public sealed class DDGIVisualizePass : ScriptableRenderPass
        {
            private DDGI _ddgiOverride;
        
            private Shader _visualizeShader;
            private Material _visualizeMaterial;
            private Mesh _visualizeMesh;

            private DDGIPass _ddgiPass;

            private ComputeBuffer _argsBuffer;

            private static class GpuParams
            {
                public static readonly string DDGIDebugIrradiance = "DDGI_DEBUG_IRRADIANCE";
                public static readonly string DDGIDebugDistance = "DDGI_DEBUG_DISTANCE";
                public static readonly string DDGIDebugOffset = "DDGI_DEBUG_OFFSET";

                public static readonly int ProbeData = Shader.PropertyToID("_ProbeData");
                public static readonly int DDGISphereObjectToWorld = Shader.PropertyToID("_ddgiSphere_ObjectToWorld");
            }

            public class DDGIVisualizePassData
            {
                internal TextureHandle cameraActiveColorTexture; // UnsafeRenderPass -> Write Only
                internal TextureHandle cameraActiveDepthTexture; // UnsafeRenderPass -> Write Only

                internal TextureHandle probeData;   // UnsafeRenderPass -> Read Only

                internal bool debugProbe;
                internal ProbeDebugMode probeDebugMode;
                internal float probeRadius;
            }

            public DDGIVisualizePass()
            {
                renderPassEvent = RenderPassEvent.AfterRenderingOpaques;

                _visualizeShader = Resources.Load<Shader>("Shaders/DDGIVisualize");
                _visualizeMaterial = CoreUtils.CreateEngineMaterial(_visualizeShader);
                _visualizeMaterial.enableInstancing = true;
            }

            public void Setup(Mesh debugMesh, DDGIPass ddgiPass)
            {
                _visualizeMesh = debugMesh;
                _ddgiPass = ddgiPass;
            }

            public override void RecordRenderGraph(RenderGraph renderGraph, ContextContainer frameData)
            {
                _ddgiOverride = VolumeManager.instance.stack.GetComponent<DDGI>();

                if (_visualizeMesh == null || _ddgiPass == null) return;
                if (_ddgiOverride == null || !_ddgiOverride.IsActive()) return;
                if (!_ddgiOverride.debugProbe.value) return;

                using (var builder = renderGraph.AddUnsafePass<DDGIVisualizePassData>("DDGI Visualize", out var passData))
                {
                    UniversalResourceData resourceData = frameData.Get<UniversalResourceData>();

                    passData.debugProbe = _ddgiOverride.debugProbe.value;
                    passData.probeRadius = _ddgiOverride.probeRadius.value;
                    passData.probeDebugMode = _ddgiOverride.probeDebugMode.value;

                    passData.probeData = renderGraph.ImportTexture(_ddgiPass._probeDataHandle, DDGIPass.GetRenderTargetInfo(_ddgiPass._probeData));

                    builder.UseTexture(passData.probeData, AccessFlags.Read);
                    builder.UseTexture(passData.cameraActiveColorTexture = resourceData.activeColorTexture, AccessFlags.Write);
                    builder.UseTexture(passData.cameraActiveDepthTexture = resourceData.activeDepthTexture, AccessFlags.Write);

                    builder.SetRenderFunc((DDGIVisualizePassData data, UnsafeGraphContext context) => ExecuteUnsafePass(data, context));
                }
            }

            public void ExecuteUnsafePass(DDGIVisualizePassData data, UnsafeGraphContext frameContext)
            {
                // Configure Debug Mode.
                {
                    frameContext.cmd.DisableShaderKeyword(GpuParams.DDGIDebugIrradiance);
                    frameContext.cmd.DisableShaderKeyword(GpuParams.DDGIDebugDistance);
                    frameContext.cmd.DisableShaderKeyword(GpuParams.DDGIDebugOffset);

                    if (data.debugProbe)
                    {
                        switch (data.probeDebugMode)
                        {
                            case ProbeDebugMode.Irradiance:
                                frameContext.cmd.EnableShaderKeyword(GpuParams.DDGIDebugIrradiance);
                                break;
                            case ProbeDebugMode.Distance:
                                frameContext.cmd.EnableShaderKeyword(GpuParams.DDGIDebugDistance);
                                break;
                            case ProbeDebugMode.RelocationOffset:
                                frameContext.cmd.EnableShaderKeyword(GpuParams.DDGIDebugOffset);
                                break;
                        }
                    }
                }

                // Prepare Rendering Parameters
                {
                    var matrix = Matrix4x4.TRS(Vector3.zero, Quaternion.identity, Vector3.one * data.probeRadius);
                    frameContext.cmd.SetGlobalMatrix(GpuParams.DDGISphereObjectToWorld, matrix);
                    frameContext.cmd.SetGlobalTexture(GpuParams.ProbeData, data.probeData);
                }

                // Construct Indirect Draw Arguments.
                {
                    var numProbes = _ddgiPass.GetNumProbes();
                    var numProbesFlat = numProbes.x * numProbes.y * numProbes.z; 
                
                    uint[] args = new uint[5] { 0, 0, 0, 0, 0 };
                    args[0] = (uint)_visualizeMesh.GetIndexCount(0);
                    args[1] = (uint)numProbesFlat;
                    args[2] = (uint)_visualizeMesh.GetIndexStart(0);
                    args[3] = (uint)_visualizeMesh.GetBaseVertex(0);
                
                    if(_argsBuffer != null) { _argsBuffer.Release(); _argsBuffer = null; }
                    _argsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
                    _argsBuffer.SetData(args);
                }

                // Draw Spheres.
                {
                    frameContext.cmd.SetRenderTarget(data.cameraActiveColorTexture, data.cameraActiveDepthTexture);
                    frameContext.cmd.DrawMeshInstancedIndirect(_visualizeMesh, 0, _visualizeMaterial, 0, _argsBuffer);
                }
            }

            public void Release()
            {
                CoreUtils.Destroy(_visualizeMaterial);

                if (_argsBuffer != null) { _argsBuffer.Release(); _argsBuffer = null; }
            }
        }

        public static void ReleaseBufferSafe(GraphicsBuffer buffer)
        {
            if (buffer != null) { if (buffer.IsValid()) { buffer.Release(); } buffer = null; }
        }

        public static RTHandle ConvertRenderTargetIdentifierToRTHandle(RenderTargetIdentifier rtId)
        {
            RTHandleStaticHelpers.SetRTHandleStaticWrapper(rtId);
            return RTHandleStaticHelpers.s_RTHandleWrapper;
        }

        private DDGIPass _ddgiPass;
        private DDGIVisualizePass _ddgiVisualizePass;

        private Mesh _ddgiVisualizeSphere;
        private bool _isRayTracingSupported;

        public override void Create()
        {
            if (!isActive)
            {
                _ddgiPass?.Release();
                _ddgiVisualizePass?.Release();
                return;
            }
        
            _isRayTracingSupported = SystemInfo.supportsRayTracing;
            if (!_isRayTracingSupported) return;

            _ddgiPass ??= new DDGIPass();
            _ddgiVisualizePass ??= new DDGIVisualizePass();

            _ddgiVisualizeSphere = Resources.Load<Mesh>("Meshes/DDGIVisualizationSphere");

        #if UNITY_EDITOR
            EditorSceneManager.sceneOpened += OnSceneOpened;
        #endif
        }

        public override void AddRenderPasses(ScriptableRenderer renderer, ref RenderingData renderingData)
        {
            if (renderingData.cameraData.isPreviewCamera) return;
            if (!_isRayTracingSupported) return;
        
            renderer.EnqueuePass(_ddgiPass);

            _ddgiVisualizePass.Setup(_ddgiVisualizeSphere, _ddgiPass);
            renderer.EnqueuePass(_ddgiVisualizePass);
        }

        protected override void Dispose(bool disposing)
        {
            base.Dispose(disposing);

            _ddgiPass?.Release(); 
            _ddgiVisualizePass?.Release();

        #if UNITY_EDITOR
            EditorSceneManager.sceneOpened -= OnSceneOpened;
        #endif
        }

        public void Reinitialize()
        {
            _ddgiPass?.Reinitialize();
        }
    
        private void OnSceneOpened(Scene scene, OpenSceneMode mode)
        {
            Reinitialize();
        }
    }
}

using System;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;
using UnityEngine.Serialization;

public enum ProbeDebugMode
{
    Irradiance,
    Distance,
    RelocationOffset
}

public enum IndirectDebugMode
{
    FullIndirectRadiance = 0,
    PureIndirectRadiance = 1
}

[Serializable, VolumeComponentMenuForRenderPipeline("Lighting/Dynamic Diffuse Global Illumination", typeof(UniversalRenderPipeline))]
public class DDGI : VolumeComponent, IPostProcessComponent
{
    // This parameter is to resolve an issue where the scene volume does not have DDGI attached, but VolumeManager.instance.stack still accesses other scene DDGI components.
    // The exact cause of this problem is not yet clear.
    public BoolParameter enableDDGI = new BoolParameter(false);

    #region Dynamic Lighting Settings

    [Header("Dynamic Lighting Settings")]
    public ClampedFloatParameter indirectIntensity = new ClampedFloatParameter(1.0f, 0.0f, 3.0f);

    public ClampedFloatParameter normalBiasMultiplier = new ClampedFloatParameter(0.2f, 0.0f, 1.0f);

    public ClampedFloatParameter viewBiasMultiplier = new ClampedFloatParameter(0.8f, 0.0f, 1.0f);

    public Vector3Parameter probeRotationDegrees = new Vector3Parameter(Vector3.zero);

    //public ClampedFloatParameter biasMultiplier = new ClampedFloatParameter(0.25f, 0.0f, 1.0f);
    //public ClampedFloatParameter axialDistanceMultiplier = new ClampedFloatParameter(0.75f, 0.0f, 1.0f);

    #endregion


    #region Probe Feature Settings

    [Header("Probe Feature Settings")]
    [Tooltip("Reallocate the position of the Probe in real-time to avoid Probes being stuck inside geometry.")]
    public BoolParameter enableProbeRelocation = new BoolParameter(true);

    public ClampedFloatParameter probeMinFrontfaceDistance = new ClampedFloatParameter(0.3f, 0.0f, 2.0f);

    [Tooltip("Disable Probes in voxel blocks with no valid geometry, allowing them to emit only a small number of fixed rays to save performance.")]
    public BoolParameter enableProbeClassification = new BoolParameter(true);

    [Tooltip("This parameter is used for both Relocation and Classification. For Relocation, when the backface hit ray ratio exceeds this value, the Probe will be offset; for Classification, the Probe will be forced to update when the ratio exceeds this value.")]
    public ClampedFloatParameter probeFixedRayBackfaceThreshold = new ClampedFloatParameter(0.25f, 0.0f, 1.0f);

    // ! Experimental Feature !
    [Tooltip("When enabled, the system will analyze the variability of GI results. If the variability is low, indicating convergence of the integral results, Probe Texture updates will stop to avoid indirect light flickering and save costs. [Note: This feature does not support emissive objects!]")]
    public BoolParameter enableProbeVariability = new BoolParameter(false);

    [Tooltip("When variability falls below this value, GI will stop updating. A higher value causes the integral to converge earlier. [Note: Too high a value may result in insufficient GI convergence, while too low a value won't prevent indirect light flickering.]")]
    public ClampedFloatParameter probeVariabilityThreshold = new ClampedFloatParameter(0.025f, 0.0f, 1.0f);

    #endregion


    #region Debug Options

    [Header("Debug Options")]
    public BoolParameter debugProbe = new BoolParameter(false);

    public ProbeDebugParameter probeDebugMode = new ProbeDebugParameter(ProbeDebugMode.Irradiance);

    public ClampedFloatParameter probeRadius = new ClampedFloatParameter(11.0f, 0.01f, 20.0f);

    [Tooltip("Only display the indirect lighting results (considering the geometry surface albedo).")]
    public BoolParameter debugIndirect = new BoolParameter(false);

    public IndirectDebugParameter indirectDebugMode = new IndirectDebugParameter(IndirectDebugMode.FullIndirectRadiance);

    #endregion


    #region Reinitialize Settings

    [Header("Reinitialize Settings (Need Refresh)")]
    [Tooltip("Use a custom DDGI bounding box. This feature requires creating a DDGI Custom Bounds in the scene.")]
    public BoolParameter useCustomBounds = new BoolParameter(false);

    public ClampedIntParameter probeCountX = new ClampedIntParameter(22, 1, 25);
    public ClampedIntParameter probeCountY = new ClampedIntParameter(22, 1, 25);
    public ClampedIntParameter probeCountZ = new ClampedIntParameter(22, 1, 25);

    [Tooltip("The number of rays each Probe emits.")]
    public ClampedIntParameter raysPerProbe = new ClampedIntParameter(144, 32, 256);

    #endregion


    public bool IsActive() => indirectIntensity.value > 0.0f && enableDDGI.value;

    public bool IsTileCompatible() => false;
}

[Serializable]
public sealed class ProbeDebugParameter : VolumeParameter<ProbeDebugMode>
{
    public ProbeDebugParameter(ProbeDebugMode value, bool overrideState = false) : base(value, overrideState) { }
}

[Serializable]
public sealed class IndirectDebugParameter : VolumeParameter<IndirectDebugMode>
{
    public IndirectDebugParameter(IndirectDebugMode value, bool overrideState = false) : base(value, overrideState) { }
}

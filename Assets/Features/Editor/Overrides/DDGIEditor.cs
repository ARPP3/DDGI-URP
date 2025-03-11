using System;
using System.Reflection;
using UnityEditor;
using UnityEditor.Rendering;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;

[CustomEditor(typeof(DDGI))]
public sealed class DDGIEditor : VolumeComponentEditor
{
    private SerializedDataParameter mEnableDDGI;
    private SerializedDataParameter mIndirectIntensity;
    private SerializedDataParameter mNormalBiasMultiplier;
    private SerializedDataParameter mViewBiasMultiplier;
    private SerializedDataParameter mProbeRotationDegrees;
    private SerializedDataParameter mDebugProbe;
    private SerializedDataParameter mProbeDebugMode;
    private SerializedDataParameter mProbeRadius;
    private SerializedDataParameter mDebugIndirect;
    private SerializedDataParameter mIndirectDebugMode;
    private SerializedDataParameter mEnableProbeRelocation;
    private SerializedDataParameter mProbeMinFrontfaceDistance;
    private SerializedDataParameter mEnableProbeClassification;
    private SerializedDataParameter mProbeFixedRayBackfaceThreshold;
    private SerializedDataParameter mEnableProbeVariability;
    private SerializedDataParameter mProbeVariabilityThreshold;
    private SerializedDataParameter mUseCustomBounds;
    private SerializedDataParameter mProbeCountX;
    private SerializedDataParameter mProbeCountY;
    private SerializedDataParameter mProbeCountZ;
    private SerializedDataParameter mRaysPerProbe;

    public override void OnEnable()
    {
        var o = new PropertyFetcher<DDGI>(serializedObject);

        mEnableDDGI = Unpack(o.Find(x => x.enableDDGI));
        mIndirectIntensity = Unpack(o.Find(x => x.indirectIntensity));
        mNormalBiasMultiplier = Unpack(o.Find(x => x.normalBiasMultiplier));
        mViewBiasMultiplier = Unpack(o.Find(x => x.viewBiasMultiplier));
        mProbeRotationDegrees = Unpack(o.Find(x => x.probeRotationDegrees));
        mDebugProbe = Unpack(o.Find(x => x.debugProbe));
        mProbeDebugMode = Unpack(o.Find(x => x.probeDebugMode));
        mProbeRadius = Unpack(o.Find(x => x.probeRadius));
        mDebugIndirect = Unpack(o.Find(x => x.debugIndirect));
        mIndirectDebugMode = Unpack(o.Find(x => x.indirectDebugMode));
        mEnableProbeRelocation = Unpack(o.Find(x => x.enableProbeRelocation));
        mProbeMinFrontfaceDistance = Unpack(o.Find(x => x.probeMinFrontfaceDistance));
        mEnableProbeClassification = Unpack(o.Find(x => x.enableProbeClassification));
        mProbeFixedRayBackfaceThreshold = Unpack(o.Find(x => x.probeFixedRayBackfaceThreshold));
        mEnableProbeVariability = Unpack(o.Find(x => x.enableProbeVariability));
        mProbeVariabilityThreshold = Unpack(o.Find(x => x.probeVariabilityThreshold));
        mUseCustomBounds = Unpack(o.Find(x => x.useCustomBounds));
        mProbeCountX = Unpack(o.Find(x => x.probeCountX));
        mProbeCountY = Unpack(o.Find(x => x.probeCountY));
        mProbeCountZ = Unpack(o.Find(x => x.probeCountZ));
        mRaysPerProbe = Unpack(o.Find(x => x.raysPerProbe));
    }

    public override void OnInspectorGUI()
    {
        if (!SystemInfo.supportsRayTracing)
        {
            EditorGUILayout.HelpBox("DDGI relies on hardware ray tracing and is only supported on DX12, Playstation 5, and Xbox Series X.", MessageType.Warning);
            return;
        }
        
        PropertyField(mEnableDDGI);

    #region Dynamic Lighting Settings

        PropertyField(mIndirectIntensity);
        PropertyField(mNormalBiasMultiplier);
        PropertyField(mViewBiasMultiplier);
        PropertyField(mProbeRotationDegrees);
        EditorGUILayout.Space(5.0f);

    #endregion

    
    #region Probe Feature Settings

        PropertyField(mEnableProbeRelocation);
        if (mEnableProbeRelocation.value.boolValue)
        {
            PropertyField(mProbeMinFrontfaceDistance);
        }
        EditorGUILayout.Space(3.0f);
        
        PropertyField(mEnableProbeClassification);
        EditorGUILayout.Space(3.0f);
        
        if (mEnableProbeRelocation.value.boolValue || mEnableProbeClassification.value.boolValue)
        {
            PropertyField(mProbeFixedRayBackfaceThreshold);
            EditorGUILayout.Space(3.0f);
        }
        
        PropertyField(mEnableProbeVariability);
        if (mEnableProbeVariability.value.boolValue)
        {
            PropertyField(mProbeVariabilityThreshold);
            EditorGUILayout.HelpBox("Probe Variability is currently an experimental feature and does not support emissive objects. Please consider using it with caution.", MessageType.Info);
        }
        EditorGUILayout.Space(5.0f);

    #endregion

    
    #region Debug Options

        PropertyField(mDebugProbe);
        if (mDebugProbe.value.boolValue)
        {
            EditorGUI.indentLevel++;
            PropertyField(mProbeDebugMode);
            PropertyField(mProbeRadius);
            EditorGUI.indentLevel--;
        }
        PropertyField(mDebugIndirect);
        if (mDebugIndirect.value.boolValue)
        {
            EditorGUI.indentLevel++;
            PropertyField(mIndirectDebugMode);
            EditorGUI.indentLevel--;
        }
        EditorGUILayout.Space(5.0f);

    #endregion


    #region Reinitialize Settings
    
        PropertyField(mUseCustomBounds);
        if (mUseCustomBounds.value.boolValue)
        {
            var customBounds = FindFirstObjectByType<DDGICustomBounds>();
            if (customBounds == null)
            {
                EditorGUILayout.HelpBox("No valid DDGI Custom Bounds detected in the current scene. You may have never created it or have it disabled; " +
                                        "To create it, you can right-click in the Hierarchy -> Light -> DDGI Custom Bounds",
                    MessageType.Warning);
            }
        }
            
        PropertyField(mProbeCountX);
        PropertyField(mProbeCountY);
        PropertyField(mProbeCountZ);
        PropertyField(mRaysPerProbe);
    
    #endregion
        
    
        if (GUILayout.Button("Refresh DDGI Settings"))
        {
            var urpAsset = GraphicsSettings.renderPipelineAsset;
            
            if (urpAsset != null && urpAsset is UniversalRenderPipelineAsset)
            {
                Type urpAssetType = urpAsset.GetType();
                FieldInfo scriptableRendererDataListField = urpAssetType.GetField("m_RendererDataList", BindingFlags.Instance | BindingFlags.NonPublic);
                    
                if (scriptableRendererDataListField != null)
                {
                    ScriptableRendererData[] rendererDataList = scriptableRendererDataListField.GetValue(urpAsset) as ScriptableRendererData[];

                    if (rendererDataList == null) return;
                        
                    foreach (var rendererData in rendererDataList)
                    {
                        var ddgiFeature = (DDGIFeature)rendererData.rendererFeatures.Find(x => x.GetType() == typeof(DDGIFeature));
                        if(ddgiFeature != null) ddgiFeature.Reinitialize();
                    }
                }
            }
        }
    }
}

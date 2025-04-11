using System;
using System.Reflection;
using UnityEditor;
using UnityEditor.Rendering;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;

namespace DDGI.Editor
{ 
    [CustomEditor(typeof(DDGI))]
    public sealed class DDGIEditor : VolumeComponentEditor
    {
        private SerializedDataParameter _enableDDGI;
        private SerializedDataParameter _indirectIntensity;
        private SerializedDataParameter _normalBiasMultiplier;
        private SerializedDataParameter _viewBiasMultiplier;
        private SerializedDataParameter _probeRotationDegrees;
        private SerializedDataParameter _debugProbe;
        private SerializedDataParameter _probeDebugMode;
        private SerializedDataParameter _probeRadius;
        private SerializedDataParameter _debugIndirect;
        private SerializedDataParameter _indirectDebugMode;
        private SerializedDataParameter _enableProbeRelocation;
        private SerializedDataParameter _probeMinFrontfaceDistance;
        private SerializedDataParameter _enableProbeClassification;
        private SerializedDataParameter _probeFixedRayBackfaceThreshold;
        private SerializedDataParameter _enableProbeVariability;
        private SerializedDataParameter _probeVariabilityThreshold;
        private SerializedDataParameter _useCustomBounds;
        private SerializedDataParameter _probeCountX;
        private SerializedDataParameter _probeCountY;
        private SerializedDataParameter _probeCountZ;
        private SerializedDataParameter _raysPerProbe;

        public override void OnEnable()
        {
            var o = new PropertyFetcher<DDGI>(serializedObject);

            _enableDDGI = Unpack(o.Find(x => x.enableDDGI));
            _indirectIntensity = Unpack(o.Find(x => x.indirectIntensity));
            _normalBiasMultiplier = Unpack(o.Find(x => x.normalBiasMultiplier));
            _viewBiasMultiplier = Unpack(o.Find(x => x.viewBiasMultiplier));
            _probeRotationDegrees = Unpack(o.Find(x => x.probeRotationDegrees));
            _debugProbe = Unpack(o.Find(x => x.debugProbe));
            _probeDebugMode = Unpack(o.Find(x => x.probeDebugMode));
            _probeRadius = Unpack(o.Find(x => x.probeRadius));
            _debugIndirect = Unpack(o.Find(x => x.debugIndirect));
            _indirectDebugMode = Unpack(o.Find(x => x.indirectDebugMode));
            _enableProbeRelocation = Unpack(o.Find(x => x.enableProbeRelocation));
            _probeMinFrontfaceDistance = Unpack(o.Find(x => x.probeMinFrontfaceDistance));
            _enableProbeClassification = Unpack(o.Find(x => x.enableProbeClassification));
            _probeFixedRayBackfaceThreshold = Unpack(o.Find(x => x.probeFixedRayBackfaceThreshold));
            _enableProbeVariability = Unpack(o.Find(x => x.enableProbeVariability));
            _probeVariabilityThreshold = Unpack(o.Find(x => x.probeVariabilityThreshold));
            _useCustomBounds = Unpack(o.Find(x => x.useCustomBounds));
            _probeCountX = Unpack(o.Find(x => x.probeCountX));
            _probeCountY = Unpack(o.Find(x => x.probeCountY));
            _probeCountZ = Unpack(o.Find(x => x.probeCountZ));
            _raysPerProbe = Unpack(o.Find(x => x.raysPerProbe));
        }

        public override void OnInspectorGUI()
        {
            if (!SystemInfo.supportsRayTracing)
            {
                EditorGUILayout.HelpBox("DDGI relies on hardware ray tracing and is only supported on DX12, Playstation 5, and Xbox Series X.", MessageType.Warning);
                return;
            }
        
            PropertyField(_enableDDGI);

        #region Dynamic Lighting Settings

            PropertyField(_indirectIntensity);
            PropertyField(_normalBiasMultiplier);
            PropertyField(_viewBiasMultiplier);
            PropertyField(_probeRotationDegrees);
            EditorGUILayout.Space(5.0f);

        #endregion

    
        #region Probe Feature Settings

            PropertyField(_enableProbeRelocation);
            if (_enableProbeRelocation.value.boolValue)
            {
                PropertyField(_probeMinFrontfaceDistance);
            }
            EditorGUILayout.Space(3.0f);
        
            PropertyField(_enableProbeClassification);
            EditorGUILayout.Space(3.0f);
        
            if (_enableProbeRelocation.value.boolValue || _enableProbeClassification.value.boolValue)
            {
                PropertyField(_probeFixedRayBackfaceThreshold);
                EditorGUILayout.Space(3.0f);
            }
        
            PropertyField(_enableProbeVariability);
            if (_enableProbeVariability.value.boolValue)
            {
                PropertyField(_probeVariabilityThreshold);
                EditorGUILayout.HelpBox("Probe Variability is currently an experimental feature and does not support emissive objects. Please consider using it with caution.", MessageType.Info);
            }
            EditorGUILayout.Space(5.0f);

        #endregion

    
        #region Debug Options

            PropertyField(_debugProbe);
            if (_debugProbe.value.boolValue)
            {
                EditorGUI.indentLevel++;
                PropertyField(_probeDebugMode);
                PropertyField(_probeRadius);
                EditorGUI.indentLevel--;
            }
            PropertyField(_debugIndirect);
            if (_debugIndirect.value.boolValue)
            {
                EditorGUI.indentLevel++;
                PropertyField(_indirectDebugMode);
                EditorGUI.indentLevel--;
            }
            EditorGUILayout.Space(5.0f);

        #endregion


        #region Reinitialize Settings
    
            PropertyField(_useCustomBounds);
            if (_useCustomBounds.value.boolValue)
            {
                var customBounds = FindFirstObjectByType<DDGICustomBounds>();
                if (customBounds == null)
                {
                    EditorGUILayout.HelpBox("No valid DDGI Custom Bounds detected in the current scene. You may have never created it or have it disabled; " +
                                            "To create it, you can right-click in the Hierarchy -> Light -> DDGI Custom Bounds",
                        MessageType.Warning);
                }
            }
            
            PropertyField(_probeCountX);
            PropertyField(_probeCountY);
            PropertyField(_probeCountZ);
            PropertyField(_raysPerProbe);
    
        #endregion
        
    
            if (GUILayout.Button("Refresh DDGI Settings"))
            {
                DDGI.RefreshDDGISettings();
            }
        }
    }
}

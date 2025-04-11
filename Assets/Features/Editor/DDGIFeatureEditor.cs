using UnityEditor;
using UnityEngine;
using UnityEngine.Rendering;

namespace DDGI.Editor
{ 
    [CustomEditor(typeof(DDGIFeature))]
    public class DDGIFeatureEditor : UnityEditor.Editor
    {
        public override void OnInspectorGUI()
        {
            if (!SystemInfo.supportsRayTracing)
            {
                EditorGUILayout.HelpBox("DDGI relies on hardware ray tracing and is only supported on DX12, Playstation 5, and Xbox Series X", MessageType.Warning);
                return;
            }
        }

        [InitializeOnLoadMethod]
        private static void OnAssemblyReloaded()
        {
            DDGI.RefreshDDGISettings();
        }
    }
}

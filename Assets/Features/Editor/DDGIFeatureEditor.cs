using UnityEditor;
using UnityEngine;

[CustomEditor(typeof(DDGIFeature))]
public class DDGIFeatureEditor : Editor
{
    private void OnEnable()
    {

    }

    public override void OnInspectorGUI()
    {
        if (!SystemInfo.supportsRayTracing)
        {
            EditorGUILayout.HelpBox("DDGI relies on hardware ray tracing and is only supported on DX12, Playstation 5, and Xbox Series X", MessageType.Warning);
            return;
        }
    }
}

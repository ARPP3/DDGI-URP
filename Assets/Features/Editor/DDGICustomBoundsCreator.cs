using UnityEditor;
using UnityEngine;

namespace DDGI.Editor
{
    public class DDGICustomBoundsCreator : UnityEditor.Editor
    {
        [MenuItem("GameObject/Light/DDGI Custom Bounds")]
        public static void CreateDDGICustomBounds()
        {
            if (FindObjectsByType<DDGICustomBounds>(FindObjectsInactive.Exclude, FindObjectsSortMode.None).Length > 0)
            {
                EditorUtility.DisplayDialog("Duplicate DDGI Custom Bounds Not Allowed", "DDGI Custom Bounds already exists in the scene", "OK");
                return;
            }
        
            var ddgiBounds = new GameObject("DDGI Custom Bounds");
            ddgiBounds.AddComponent<BoxCollider>();
            ddgiBounds.AddComponent<DDGICustomBounds>();
        }
    }
}

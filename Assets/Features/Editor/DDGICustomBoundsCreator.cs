using System.Collections;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;

public class DDGICustomBoundsCreator : Editor
{
    [MenuItem("GameObject/Light/DDGI Custom Bounds")]
    public static void CreateDDGICustomBounds()
    {
        if (FindObjectsOfType<DDGICustomBounds>().Length > 0)
        {
            EditorUtility.DisplayDialog("Duplicate DDGI Custom Bounds Not Allowed", "DDGI Custom Bounds already exists in the scene", "OK");
            return;
        }
        
        var ddgiBounds = new GameObject("DDGI Custom Bounds");
        ddgiBounds.AddComponent<BoxCollider>();
        ddgiBounds.AddComponent<DDGICustomBounds>();
    }
}

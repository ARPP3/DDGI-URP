using UnityEngine;

namespace DDGI
{
    [RequireComponent(typeof(BoxCollider)), ExecuteAlways]
    public class DDGICustomBounds : MonoBehaviour
    {
        [SerializeField]
        private bool _drawOutline = true;

        private BoxCollider _boxCollider;

        // This class is empty and is only used to identify a GameObject for calling the FindObjectByType method.
        private void OnEnable()
        {
            if (FindObjectsByType<DDGICustomBounds>(FindObjectsInactive.Include, FindObjectsSortMode.None).Length > 1)
            {
                Debug.LogError("Duplicate DDGI Custom Bounds Not Allowed");
#if UNITY_EDITOR
                if (Application.isPlaying) Destroy(this);
                else                       DestroyImmediate(this);
#else
                Destroy(this);
#endif
                return;
            }

            _boxCollider = GetComponent<BoxCollider>();
            _boxCollider.isTrigger = false;
        }


        private void OnDrawGizmos()
        {
            if (!_drawOutline) return;

            Gizmos.matrix = transform.localToWorldMatrix;
            Gizmos.color = Color.green;

            Gizmos.DrawWireCube(_boxCollider.center, _boxCollider.size);

            Gizmos.color = Color.white;
            Gizmos.matrix = Matrix4x4.identity;
        }
    }
}

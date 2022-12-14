using UnityEngine;

namespace Mirror.Cloud.Examples
{
    /// <summary>
    /// Instantiate a new NetworkManager if one does not already exist
    /// </summary>
    public class InstantiateNetworkManager : MonoBehaviour
    {
        public GameObject prefab;

        public void myinstantiate()
        {
            Instantiate(prefab);
        }
    }
}

using UnityEngine;
using UnityEngine.UI;
using UnityEngine.XR;

using System.Collections.Generic;
using System.Collections;
using System.IO;
using System;
using System.Linq;
using System.Reflection;

using UnityEngine.UI.Extensions;
using UnityEngine.UI.Extensions.ColorPicker;
using UnityEngine.EventSystems;

using UMol.API;
using Mirror;
public class StartManager : MonoBehaviour
{
    public void UIStartHost()
    {
        GameObject.Find("NetworkManager").GetComponent<NetworkManagerHUD>().UIStartHost();
    }
    public void UIStartClient(GameObject t)
    {
        GameObject.Find("NetworkManager").GetComponent<NetworkManagerHUD>().UIStartClient(t);
    }
    public void UIStartServer()
    {
        GameObject.Find("NetworkManager").GetComponent<NetworkManagerHUD>().UIStartServer();
    }
}

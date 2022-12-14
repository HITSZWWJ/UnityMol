using System;
using System.Collections.Generic;
using Mirror;
using UnityEngine;
using UnityEngine.UI;
using Random = UnityEngine.Random;
using UMol;
using UMol.API;


public class PlayerScript : NetworkBehaviour
{
    public GameObject floatingInfo;
    public TextMesh nameText;
    private Material playerMaterialClone;

    private SceneScript sceneScript;

    [SyncVar(hook = nameof(OnPlayerNameChanged))]
    private string playerName;

    [SyncVar(hook = nameof(OnPlayerColorChanged))]
    private Color playerColor;

    private GameObject VRTK_SDKManager;
    private GameObject VRTK_SDKSetups;
    private GameObject SteamVR;
    private GameObject VRCameraRig;
    private GameObject VRCamera;
    public bool isSender = false;
    public bool isRecipter = false;
    public GameObject CanvasMainUI;
    public float rotateSpeed = 0.05f;

    private void OnPlayerNameChanged(string oldStr, string newStr)
    {
        nameText.text = newStr;
    }

    private void OnPlayerColorChanged(Color oldCol, Color newCol)
    {
        nameText.color = newCol;

        playerMaterialClone = new Material(GetComponent<Renderer>().material);
        playerMaterialClone.SetColor("_EmissionColor", newCol);

        GetComponent<Renderer>().material = playerMaterialClone;
    }

    [Command]
    private void CmdSetupPlayer(string nameValue, Color colorValue)
    {
        playerName = nameValue;
        playerColor = colorValue;
        sceneScript.statusText = $"{playerName} joined.";
        RpcSetupPlayer(nameValue, colorValue);
    }

    [ClientRpc]
    private void RpcSetupPlayer(string nameValue, Color colorValue)
    {
        playerName = nameValue;
        playerColor = colorValue;
    }

    [Command]
    public void CmdSendPlayerMessage()
    {
        if (sceneScript)
        {
            sceneScript.statusText = $"{playerName} says hello {Random.Range(1, 99)}";
        }
    }

    //执行命令
    [Command(ignoreAuthority = true)]
    public void CmdExecuteCommand(string command, NetworkConnectionToClient sender = null)
    {
        Debug.Log("Cmd");
        RpcExecuteCommand(command);
    }
    [ClientRpc]
    public void RpcExecuteCommand(string command)
    {
        Debug.Log("Rpc: " + command);
        if (isSender)
        {
            Debug.Log("is sender");
            isSender = false;
            return;
        }
        GameObject.Find("Player").GetComponent<PlayerScript>().isRecipter = true;
        GameObject.Find("ConsolePython_Autocomplete").GetComponent<PythonConsole2>().ExecuteCommand(command);
    }


    public override void OnStartLocalPlayer()
    {
        base.OnStartLocalPlayer();

        sceneScript = FindObjectOfType<SceneScript>();
        //依次获取VR设备
        VRTK_SDKManager = GameObject.Find("[VRTK_SDKManager]").gameObject;
        VRTK_SDKSetups = VRTK_SDKManager.transform.Find("[VRTK_SDKSetups]").gameObject;
        SteamVR = VRTK_SDKSetups.transform.Find("SteamVR").gameObject;
        VRCameraRig = SteamVR.transform.Find("[CameraRig]").gameObject;
        VRCamera = VRCameraRig.transform.Find("UICamera").gameObject;

        this.name = "Player";//本地生成的是player，其他是player(clone)
    }


    private void Update()
    {
        if (!isLocalPlayer) return;
        var VRCameraPosition = VRCamera.transform.localPosition;
        var VRCameraRotation = VRCamera.transform.localRotation;

        transform.localPosition = VRCameraPosition;

        RotateFollow();

        //if (Input.GetKeyDown(KeyCode.C))
        //{
        // ChangePlayerNameAndColor();
        //}
    }

    private void RotateFollow()
    {//让一个物体始终与相机的y轴方向上旋转保持一致
        float angle;
        angle = transform.eulerAngles.y - VRCamera.transform.eulerAngles.y;
        if (angle < 0)
        {
            angle += 360;
        }
        if (angle > rotateSpeed && angle < 180)
        {
            transform.Rotate(0f, -rotateSpeed, 0f);
        }
        else if (angle > 180 && angle < 360 - rotateSpeed)
        {
            transform.Rotate(0f, rotateSpeed, 0f);
        }
    }


    private void ChangePlayerNameAndColor()
    {
        var tempName = $"Player {Random.Range(1, 999)}";
        var tempColor = new Color
        (
            Random.Range(0f, 1f),
            Random.Range(0f, 1f),
            Random.Range(0f, 1f),
            1
        );

        CmdSetupPlayer(tempName, tempColor);
    }
}

using System;
using System.Collections.Generic;
using Mirror;
using UnityEngine;
using UnityEngine.UI;
using Random = UnityEngine.Random;
using UMol;
using UMol.API;

public class CmdManager : NetworkBehaviour
{
    public bool isSender = false;
    public bool isRecipter = false;
    //执行命令
    [Command]
    public void CmdExecuteCommand(string command, NetworkConnectionToClient sender = null)
    {
        //if (!isLocalPlayer) return;
        Debug.Log("Cmd");
        isSender = true;
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
        GameObject.Find("ConsolePython_Autocomplete").GetComponent<PythonConsole2>().ExecuteCommand(command);
        isRecipter = true;
    }
}

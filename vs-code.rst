How to use VS Code connected (or not) to HiperGator
***********

Setup
================

Downloading VS Code
----------------------------

Go to https://code.visualstudio.com/download, and click on the download button related to your computer (Mac, Windows and Linux options).

Installing VS Code
----------------------------

Go to https://code.visualstudio.com/docs/setup/setup-overview. On the left, below 'SETUP', click on your operational system (Mac, Windows or Linux) for instructions on how to install. The installing process is pretty straight forward though, similar to other applications.


Get started
================

Here are some cool short videos on how to get started on a bunch of thing in VS Code: https://code.visualstudio.com/docs/getstarted/introvideos.

The things that I find most important are accesed from the left hand side options:

<img src="images/VSCode_docs.jpg" alt="drawing" width="30"/> Explorer: where you find all your folders and files. You can open an specific folder or clone from github.

<img src="images/VSCode_github.jpg" alt="drawing" width="30"/> Source Control: where you can see all github (or any other source control platform) commands; you can pull, push, commit, etc.

<img src="images/VSCode_extensions.jpg" alt="drawing" width="30"/> Extensions: where you find other softwares and interfaces that you need to be installed in your VS Code, such as: Python, Jupyter, and Remote - SSH.

<img src="images/VSCode_settings.jpg" alt="drawing" width="31
"/> Manage: click on this and then on "Command Palette" to access and search commands. Type "terminal" and click on the first option to open up a terminal window inside VS Code, for example.


Connecting VS Code to HiperGator
================

Make VS Code SSH to HPG
----------------------------

Detailed instructions are here:
https://code.visualstudio.com/docs/remote/ssh.

In summary:

1) Go to "Extensions" on the left-side icon list, search for "Remote - SSH" and download it;
2) Go to <img src="images/VSCode_compal.jpg" alt="drawing" width="20"/> on the bottom left side of the window and a command palette will show up.
3) Click on "Connect to Host" -> "Configure SSH Hosts" -> "/Users/username/.ssh/config" and the config file will open up
4)  Write this in the file:
    Host hpg
    HostName hpg2.rc.ufl.edu
    User yourusername
    ControlMaster auto
    ControlPath ~/.ssh/%r@%h:%p
Save and close it.
5) Try to go to <img src="images/VSCode_compal.jpg" alt="drawing" width="20"/> -> "Connect to Host" -> "hpg" -> login with passord and authetication process -> open a folder (if you got until here it worked)

Now what you're going to do everytime you'll ssh through VSCode:
1) Open a terminal and ssh into hpg
2) Open VSCode and go to <img src="images/VSCode_compal.jpg" alt="drawing" width="20"/> -> "Connect to Host" -> "hpg"
3) Open the folder you want to on hpg
Doing these steps will prevent you from being asked your passowrd a bunch of times.

For more about this password fixing: https://stackoverflow.com/questions/69277631/2fa-with-vs-code-remote-ssh
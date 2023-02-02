Common Software Tools and Procedures
**********


Anaconda Environments
============

Using new anaconda environments for individual projects or major code
initiatives can be super helpful, especially when testing and
debugging.  In what follows, we outline how to go about this.

After installing Anaconda, you can actiavte a new environment (and activate it) via::

  conda create --name myenv python=3.9
  conda activate myenv



When you do this, you are in your new environment, and any package you
install via conda will remain local to this conda environment, and not
see (or be touched) by any other environment.  Why would you want to
do this?

#. You can then delete an entire environment, and all of it's packages, without messing up other environments via::

     conda remove --name myenv --all


#. Imagine you're running a ton of (e.g.,) powderday simulations, and
   you want to code some tests up in a new branch of powderday and
   debug them without messing up your major simulation campaign you
   have going on.  Then you can install the new powderday test branch
   to a new environment, and run locally there, without worrying about
   your other environments getting messed up.
		
You can find many other options/nuances at
https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#removing-an-environment
, though the above should be enough to get you going.
		


Tmux
============
Tmux stands for the terminal multiplexer. Don't worry about the fancy terminology. The main use of tmux for us at least is to keep the single/multiple terminal sessions running even after you log out from hipergator. To load tmux just do::

    	module load tmux

and to start a tmux session just type:: 

    	tmux
    
This will start a tmux session which will look like a new terminal window. You can run your codes in this window and they will keep running even after you log out from hipergator. Once you have your code running, you can exit the tmux window by hitting::

  	ctrl–b–d
	
To list all the tmux sessions you have runnning use tmux ls and to open an existing tmux session use:: 

	tmux attach -t session-name


Keep in mind that you cannot scroll by default inside a tmmux session. To scroll you must hit::
	
	ctrl–b–[
	
to go into scrolling mode. To get out of it just hit *q*
	
A list of the most useful tmux commands can be found here: https://danielmiessler.com/study/tmux/

One common use for tmux sessions is to ask for a development node inside a tmux session. You can do so via::

	srun --pty --partition=hpg-dev --time=10:00:00 --nodes=1 --ntasks=10  --mem=60gb --qos=narayanan bash -i
	
Once the node is allocated you can close the tmux session and then you can access this node anywhere using::
	
	ssh -XC \$(squeue -u $USER -n bash --states=R --noheader --Format=nodelist | head -n1)
	
It is recommended to alias the above two commands in your .bashrc file. This can be done by adding::

	alias devnode="ssh -XC \$(squeue -u $USER -n bash --states=R --noheader --Format=nodelist | head -n1)"
	alias dev='srun --pty --partition=hpg-dev --time=10:00:00 --nodes=1 --ntasks=10  --mem=60gb --qos=narayanan bash -i'
	
to your .bashrc file which can be found in your home directory. Once you have edited the bashrc file type::

	source .bashrc
	
to reload the .bashrc file and now you should be able access these commands just by using their aliases. This can be done for any commands you frequently use.

Git
============

Microsoft Virtual Studio Code
============

How to use VS Code connected (or not) to HiperGator

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



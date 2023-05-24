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

Git (todo)
============

.. include:: vs-code.rst


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



Git
============

Microsoft Virtual Studio Code
============


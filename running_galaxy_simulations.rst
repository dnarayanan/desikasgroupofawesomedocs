Running Galaxy Simulations
**********


Cosmological Simulations (Large-Box)
============

Cosmological Simulations (Zoom-in)
============

Idealized Galaxy Simulations
============

First, all of the relevant files live on a private GitHub space at
https://github.com/dnarayanan/arepo_ics - be sure to ask if you don't
have permissions here.

Initial Conditions
-----------------

Make New Disk:

The first thing we'll need to do is run MakeNewDisk in order to make
an idealized galaxy disk IC.  You can find some examples on HiPerGator
here::

  /home/paul.torrey/InitialConditions/PhilsSpecificICs/MW
  /home/paul.torrey/InitialConditions/PhilsSpecificICs/SMC
  /home/paul.torrey/InitialConditions/PhilsSpecificICs/Sbc
  /home/paul.torrey/InitialConditions/PhilsSpecificICs/HiZ
  
(though note because of the size only the MW example is on GitHub).  You have to edit the resolution in main.c and then::

  make clean
  make

This is light weight code and can usually run from even a login node.  You can just run from something like::

  ./MakeHubbleType ./output/MW.dat

which produces your basic initial condition file.

**Adding Backgrounds for Arepo Simulations** If we were running a
gizmo simulation, we could stop right here.  Note, however, if we are
running an arepo simulation, this will require a background grid to be
added to the IC.  To do this, we use sbatch_makeIC.sh which is a job
script that calls arepo from the arepo_addbg directory.  Note, you'll
need a working executable for arepo in that directory.  There are
param files that are necessary for adding the background -- these
parameter files will have the params you'll use at simulation runtime,
and there are examples in::

  /blue/narayanan/desika.narayanan/MakeGalaxy/arepo_addbg

as well as the GitHub repo.  An example of this .sh file (for posterity) is::

  (pd4env_gcc) [desika.narayanan@login2 MakeGalaxy]$ more sbatch_makeIC.sh
  #!/bin/bash
  #SBATCH --job-name=makeIC
  #SBATCH --mail-type=ALL
  #SBATCH --mail-user=desika.narayanan@gmail.com
  #SBATCH --time=6:00:00
  #SBATCH --nodes=1
  #SBATCH --tasks-per-node=8
  #SBATCH --ntasks-per-socket=8
  #SBATCH --cpus-per-task=1
  #SBATCH --distribution=cyclic:cyclic
  #SBATCH --mem-per-cpu=8gb
  #SBATCH --partition=hpg2-compute
  #SBATCH --account=narayanan
  #SBATCH --qos=narayanan-b

  
  module purge
  #module load ddt/18.0.2
  module load intel/2018
  module load gsl
  module load openmpi/3.1.2
  module load hdf5
  #module load grackle
  
  DATADIR=$SLURM_SUBMIT_DIR
  
  export OMPI_MCA_pml="ucx"
  export OMPI_MCA_btl="^vader,tcp,openib"
  export OMPI_MCA_oob_tcp_listen_mode="listen_thread"
  
  srun --mpi=pmix_v2     ./arepo_addbg/Arepo   arepo_addbg/param_MW_ultra_lowres.txt 0        1> output_makeIC/OUTPUT  2> output_makeIC/ERROR



When you run this background addition, it
will automatically make a new file that has the appendate
"--with-grid.hdf5".  For example if your input file was called
"MW_lr.dat", your output file (assuming output type 3 is used) will be
called "MW_lr.dat-with-grid.hdf5"

**Further Modifying ICs for Dust**

We have one final step which is to move all type 2 and type 3 particles to type 4.  We do this with a script written by Qi Li called modifyIC.py -- there's an example on GitHub as well as at::

  /blue/narayanan/desika.narayanan/MakeGalaxy

Run this and we should have a new IC that is ready for Arepo as well!

Run Time
-----------------

You can find exmaples of the parameter files in the GitHub subdirectory idealized_repo_example or at::

  /blue/narayanan/desika.narayanan/arepo_runs/idealized/MW_ultra_lowres

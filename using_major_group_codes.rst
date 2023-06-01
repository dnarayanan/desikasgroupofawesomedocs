Major Codes Used by our Group
**********

.. contents:: Section Contents
    :local:

Note --- some of these may be developed by our group, and some may
have authors outside our group.  What follows is a quick instruction
manual for any of these codes to get started quickly.  This is not
intended to be a substitute for the main code manuals/websites, though
those are documented below as well.



Prospector
============

Overview
------------------
Prospector is an SED fitting tool developed by Ben Johnson, Joel Leja, and Charlie Conroy at CfA. The git repo can be found at::

  https://github.com/bd-j/sedpy

Prospector uses fsps to model stellar populations and dynesty to optimize the fsps models to the input observed data. Below are instructions on how to install prospector and its dependencies.

Installation
------------------

sedpy : this manages the details of our observations and allows us to interface with the properties of the photometric filters used. You can find sedpy here: https://github.com/bd-j/sedpy

Clone the sedpy repo by doing::

  >cd $HOME                                                                                                                                    
  >git clone https://github.com/bd-j/sedpy.git

Then install using::

  >pip install .

dynesty : this is the backbone of propsector that handles the actual fitting methods. It is a form of Bayesian modeling similar to MCMC codes like 'emcee' but a bit more sophisticated. You can find dynesty here: https://github.com/joshspeagle/dynesty


Clone the dynesty repo by doing::

  >cd $HOME                                                                                                                                                       
  >git clone https://github.com/joshspeagle/dynesty.git

Then install using::

  >python setup.py install

fsps : this handles the stellar modeling for prospector. by itself, it's a super useful tool for generating stellar spectra. The core of fsps is a Fortran code but these days, the python bindings for fsps now come with its own fsps source code, meaning we no longer have to compile the Fortran code first then install the python wrapper. All we need is to clone the fortran fsps::

    >export SPS_HOME="/path/where/you/want/to/download/fsps"
    
    >git clone https://github.com/cconroy20/fsps.git $SPS_HOME

And then pip install python-fsps::

    >python -m pip install fsps


With the dependencies installed, Prospector can now be cloned and installed::

  >cd $HOME                                                                                                                                     >git clone https://github.com/bd-j/prospector.git
  >cd prospector
  >python -m pip install .

It also has a pretty decent demo on how to use prospector: https://github.com/bd-j/prospector/blob/main/demo/InteractiveDemo.ipynb


For visualization purposes, we'll also want to install 2 packages: corner, which helps us plot corner plots::

  >cd $HOME               
  >git clone https://github.com/dfm/corner.py.git
  >cd corner
  >python -m pip install .


And arviz, which allows us to do 'advanced' things with corner::

  >conda install arviz

The expanded tutorial on how to use Prospector can be found at::

  https://github.com/dnarayanan/desikasgroupofawesome/blob/main/seds_with_prospector.rst

Slick (including Despotic)
============

Overview
------------------

The main public repository and manual is here:

https://github.com/karolinagarcia/slick


Installing DESPOTIC
-------------------

DESPOTIC is written and maintained by Mark Krumholz (ANU), with code
paper, repository, and manual here.  It is a fairly extensive and
powerful code, and we recommend you read through the examples that
ship with DESPOTIC to have at least a surface-level understanding of
the code before embarking on SLICK.

https://bitbucket.org/krumholz/despotic

https://ui.adsabs.harvard.edu/abs/2014MNRAS.437.1662K/abstract



Clone the DESPOTIC repository if you haven't already::

  >cd $HOME
  >git clone git@bitbucket.org:krumholz/despotic.git

There are a few changes which need to be made in order to install despotic with the dependencies we have now.
These changes can be made from a patch file::

  >cd $HOME/despotic
  >curl -L https://gist.github.com/smsutherland/11556695e210ed8ffd18e0ffbed91963/raw/ | git apply

Finally install DESPOTIC using::

  >python setup.py install

Installing SLICK

::

  >cd $HOME
  >git clone https://github.com/karolinagarcia/slick.git
  >cd slick
  >bash ./setup.sh
  >source $HOME/.bashrc

Congratulations! You've now set up SLICK.

The SLICK Conda Environment
----------------------
Create a new conda environment and install the following packages into it::

  >conda create --name slick python=3.10.4
  >conda activate slick
  >conda install pandas=1.4.2 numpy=1.22.3 tqdm=4.64.0 astropy=5.0.4 scipy=1.7.3
  >pip install Cython==0.29.32

Install yt and caesar using the instructions in Code Installation Notes.

Basic Usage
----------------------

To create a template slick project, run `slick_new.sh`. This creates [parameters.ini](#parameters-file) and `run.sh`.
`run.sh` calls slick_frontend.sh. It assumes slick_frontend.sh is on the path, so this will need to be changed if you haven't run `setup.sh` previously.
Running `sbatch run.sh` will enqueue the job for the slick_init step and the slick_run step (unless the skip_run option is set).

Parameters File
---------------

The behavior of slick is configured by the parameters.ini file. The following describes the options currently available.

.. code-block:: ini

  [snap]
  ; This is used when naming the clouds_per_core file.
  boxsize=[int]
  ; The full path to the yt file being operated on.
  ytfilename=[str]
  ; The full path to the caesar file being operated on.
  caesarfilename=[str]
  [sample]
  ; Either 'total' or 'randomized'.
  ; Determines whether slick should operate on all clouds or a sample of clouds.
  mode=[str]
  ; How many galaxies should slick operate on.
  ; Only required if mode = 'randomized'
  n_galaxies_sample=[int]
  ; Minimum galaxy mass bound to operate on.
  ; Only required if mode = 'randomized'
  min_mass=[float]
  ; Maximum galaxy mass bound to operate on.
  ; Only required if mode = 'randomized'
  max_mass=[float]
  [sbatch]
  ; This section takes any key value pairs which can be used in a job script as 
  #SBATCH --key=value
  ; These are used to configure the slurm job for the slick_run step.
  ; The only parameter not configurable is array, which is set by slick to match the number of runs being prepared.
  ; The following is the default configuration
  nodes=1
  tasks-per-node=1
  cpus-per-task=1
  mem-per-cpu=8gb
  time=96:00:00
  output=/dev/null
  [run]
  ; The directory which slick should output its files.
  ; This does not include any logs generated by slurm.
  ; To change the output of those, use the output parameter in [sbatch]
  ; Defaults to Output_Files
  output_dir=[str]
  ; If true, only the slick_init step is run.
  ; The slick_run step can be triggered manually via `sbatch slick_run_jobscript.sh`
  ; Defaults to false
  skip_run=[bool]




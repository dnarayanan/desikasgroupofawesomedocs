Code Installation Notes
**********


Powderday on HiPerGator 
============



Manual Installation With gcc Compilers
-----------------

[1] As of August 2022, this is the preferred way to install (i.e., least headache).  Thanks to Prerak Garg for these.::

First, load up the compilers that we'll use throughout::

  >module load gcc/9.3.0 openmpi/4.1.1 libz/1.2.11 hdf5/1.10.1 git/2.30.1

  
yt::

  >cd $HOME
  >git clone https://github.com/yt-project/yt
  >cd yt
  >pip install -e .



fsps and python-fsps

The development version of python-fsps now includes the Fortran FSPS source code. You can get both via::

fsps::

  >cd $HOME
  >git clone https://github.con/cconroy20/fsps

in the Makefile set F90=$(FC) and this will ensure that the compilers
`fsps <https://code.google.com/p/fsps/source/checkout>`_ uses are what
you have module loaded.::
  
  >make clean
  >make

then in your .bashrc set the analog to::
  
  >export SPS_HOME=/Users/desika/fsps


python fsps::

>CC=gcc F90=gfortran F77=gfortran python setup.py install



hyperion

As of commit 7cae6d0, a bug has been introduced with the __version__ module. Once cloned, checkout stable commit 4170c6c::

  >cd $HOME
  >git clone https://github.com/hyperion-rt/hyperion.git
  >cd hyperion
  >git checkout 4170c6cc3009893e2b591e133baeb9927122aef1
  >python setup.py install
  >git submodule init
  >git submodule update

  >./configure --prefix=$HOME/local

  >make
  >make install

hyperion dust::

  >cd $HOME
  >wget http://pypi.python.org/packages/source/h/hyperion-dust/hyperion-dust-0.1.0.tar.gz
  >tar -xzvf hyperion-dust-0.1.0.tar.gz
  >cd hyperion-dust-0.1.0
  >python setup.py build_dust

  
powderday::

  >git clone https://github.com/dnarayanan/powderday.git
  >conda install numpy scipy cython h5py matplotlib psutil joblib six astropy scikit-learn ipython
  >cd powderday
  >python setup.py install

  


Manual Installation With Intel Compilers
-----------------

[2] The first set of instructions for the University of Florida
HiPerGator3.0 facility is to employ intel compilers, and to compile
everything manually.  This allows the greatest flexibility, as well as
the ability to use private forks of individual codes.

First, load up the compilers that we'll use throughout::

  >module load intel/2018.1.163
  >module load openmpi/4.0.3
  >module load hdf5/1.10.1
  >module load git

yt::

  >cd $HOME
  >git clone https://github.com/yt-project/yt
  >cd yt
  >pip install -e .


fsps::

  >cd $HOME
  >git clone https://github.con/cconroy20/fsps

in the Makefile set F90=$(FC) and this will ensure that the compilers
`fsps <https://code.google.com/p/fsps/source/checkout>`_ uses are what
you have module loaded.::
  
  >make clean
  >make

then in your .bashrc set the analog to::
  
  >export SPS_HOME=/Users/desika/fsps


python fsps::

  >cd $HOME
  >git clone --recursive https://github.com/dfm/python-fsps.git
  >cd python-fsps
  >CC=icc F90=ifort python setup.py install


hyperion::

  >cd $HOME
  >git clone https://github.com/hyperion-rt/hyperion.git
  >cd hyperion
  >python setup.py install
  >git submodule init
  >git submodule update
  >./configure --prefix=$HOME/local
  >make
  >make install

hyperion dust::

  >cd $HOME
  >wget http://pypi.python.org/packages/source/h/hyperion-dust/hyperion-dust-0.1.0.tar.gz
  >tar -xzvf hyperion-dust-0.1.0.tar.gz
  >cd hyperion-dust-0.1.0
  >python setup.py build_dust
  
powderday::

  >git clone https://github.com/dnarayanan/powderday.git
  >cd powderday
  >python setup.py install







  

Conda Installation
-----------------
  
[3] The final set of instructions use gcc, and the conda installation
of `Hyperion <http://www.hyperion-rt.org>`_.  Thanks to Paul Torrey
for these.::

  >module load openmpi/4.1.1 libz/1.2.11 hdf5/1.10.1 conda/4.12.0 git/2.30.1 gcc
  >conda install -c conda-forge hyperion
  >python -c "import hyperion" (just to ensure no errors thrown)
  >hyperion (just to ensure command is found)
  >python -m pip install fsps
  >[set $SPS_HOME variable in .bashrc)
  >cd $HOME
  >git clone https://github.com/dnarayanan/powderday.git
  >cd powderday
  >python setup.py install

then fix import six line in the equivalent of all of these::

  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/model/model.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/util/validator.py 
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/conf/conf_files.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/filter/filter.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/dust/dust_type.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/model/model_output.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/flared_disk.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/alpha_disk.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/bipolar_cavity.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/ulrich_envelope.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/power_law_envelope.py 
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/ambient_medium.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/model/sed.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/model/image.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/grid/yt3_wrappers.py



Caesar on HiPerGator
============
	    

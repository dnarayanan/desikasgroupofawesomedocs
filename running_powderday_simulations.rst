Running POWDERDAY Simulations
**********

Is POWDERDAY Working?
============

First, if you have not already, you'll need to verify that your POWDERDAY is working properly.
To do so, follow the checks at the POWDERDAY installation docs (https://powderday.readthedocs.io/en/latest/).
My usual test is the SKIRT GIZMO Milky Way zoom simulation. Briefly, you'll want to follow these commands after loading the appropriate modules::

	cd powderday
	python pd_front_end.py tests/SKIRT/gizmo_mw_zoom parameters_master_gizmo parameters_model_gizmo

If all goes will, this should finish without a problem. You should also check the output SED following the POWDERDAY installation docs. 
The 'parameters_master' file you use when running POWDERDAY will dictate exactly what important parameters you are giving to POWDERDAY (e.g. whether to have AGN on, IMF being used, etc.).
The 'parameters_model' file includes information such as where the results will be stored, and where to find the galaxy in the simulation.
If you're unsure about whether everything in these two files is good, post in the Slack.


Workflow Setup
============
Retreive some version of Sidney's scripts for automating (https://github.com/smlower/sl_simulation_tools) the POWDERDAY runs. 
Dhruv's modified versions can be found in ``/home/d.zimmerman/sl_simulation_tools/`` (and hopefully on GitHub soon). These are already set up to help
automate everything and will make your life significantly easier. 

Initial Setup Work
============


Filtering Galaxies
-----------------

At this point, you hopefully have a CAESAR file to reference for galaxies for the appropriate snapshots.
If not, refer to the CAESAR docs (https://caesar.readthedocs.io/en/latest/) for how to get that done.
The first thing you will want to do is filter the galaxy or galaxies you would like to run POWDERDAY on. 
This essentially entails extracting the star and gas particles from full simulation that are associated with the galaxy and copying them over to a new file with only those particles.
This has a few very useful benefits; (1) the file sizes become significantly smaller by virtue of only needing the 
particles necessary for radiative transfer and the ones specifically associated with the galaxy, 
which will make it much faster to load everything and easier to not hog group space with the files, especially if you need to use a galaxy from a large box, and 
(2) you don't need to be as careful with setting the POWDERDAY range because all potential background sources will no longer be included in the simulation file.

The key script for filtering SIMBA simulations is (unsurprisingly) ``filter_simba.py``. This script is built to filter a single galaxy. You'll need to provide it inline with the
path to the complete set of snapshots (i.e. the full simulation snapshot at a specific time), the specific snapshot you want, the specific galaxy number you want, and where to store
the resulting filtered galaxy (this should always be somewhere in your folder on the group's ``/orange/narayanan/`` drive). Notably, you may also need to modify the default CAESAR
file location if you don't want to draw from the old galaxy catalogs. My (Dhruv's) current version of the script, ``/home/d.zimmerman/sl_simulation_tools/filter_simba_all.py`` 
has been modified to filter all the galaxies from a particular snapshot. Without the noise in the script, the essential parts of it look like this::

	import h5py
	import caesar
	import sys
	import glob
	import numpy as np
	import tqdm

	sim = 'm100n1024' # which box you want to pull from
	snapshot_path = '/orange/narayanan/[...]' #where is the snapshot?
	###########
	# Line arguments
	###########
	snap_num = sys.argv[1] # only input left in this version is snapshot number
	##############

	output_path = '/orange/narayanan/d.zimmerman/simba/'+sim+'/snap'+str(snap_num)+'/filtered/' # where should it go?
	ds = snapshot_path
	caesar_file = '/orange/narayanan/d.zimmerman/simba/'+sim+'/caesar_cats/caesar_simba_'+str(snap_num)+'.hdf5' #where's your CAESAR file?
	obj = caesar.load(caesar_file)

	# holdover from naming conventions of snapshots
	if(int(snap_num) < 100):
	        snap_str = "0"+str(snap_num)
	else:
	        snap_str = str(snap_num)

	input_file = h5py.File(ds+str(snap_str)+'.hdf5', 'r')

	galcount = len(obj.galaxies)
	# various modifications can be done here to not run over everything
	for galaxy in range(galcount):
      	  print()
	        print("GALAXY NUM:",str(galaxy))
	        print()
	        glist = obj.galaxies[int(galaxy)].glist
	        slist = obj.galaxies[int(galaxy)].slist

	        with h5py.File(output_path+'galaxy_'+str(galaxy)+'.hdf5', 'w') as output_file:
	            output_file.copy(input_file['Header'], 'Header')
	            print('starting with gas attributes now')
	            output_file.create_group('PartType0')
	            for k in tqdm.tqdm(input_file['PartType0']):
	                output_file['PartType0'][k] = input_file['PartType0'][k][:][glist]
	            print('moving to star attributes now')
	            output_file.create_group('PartType4')
	            for k in tqdm.tqdm(input_file['PartType4']):
	                output_file['PartType4'][k] = input_file['PartType4'][k][:][slist]


	        print('done copying attributes, going to edit header now')
	        outfile_reload = output_path+'galaxy_'+str(galaxy)+'.hdf5'

	        re_out = h5py.File(outfile_reload,'r+')
	        re_out['Header'].attrs.modify('NumPart_ThisFile', np.array([len(glist), 0, 0, 0, len(slist), 0]))
	        re_out['Header'].attrs.modify('NumPart_Total', np.array([len(glist), 0, 0, 0, len(slist), 0]))

	        re_out.close()



Galaxy Positions
-----------------

The next, relatively minor, part of the setup process requires running the ``galaxy_positions.py`` script. The purpose of this script is to use the
newly generated filtered snapshots and simply generate a list of the positions of the center of the galaxies.
Again, Dhruv's current version looks like this::

	import h5py
	import numpy as np
	import sys, os
	import numpy as np
	import glob
	import tqdm

	##############
	# Line arguments
	###############
	snap = int(sys.argv[1])
	snap_dir = '/orange/narayanan/[...]' #where are the filtered galaxies?
	outfile = '/orange/narayanan/[...]'+'_gal_positions.npz' #where do you want the output to go?
	################

	pos = {}
	ngalaxies = {}
	infiles = sorted(glob.glob(snap_dir+'/galaxy_*.hdf5'))
	for i in tqdm.tqdm(range(len(infiles))):
	    try:
	        infile = h5py.File(snap_dir+'/galaxy_'+str(i)+'.hdf5', 'r')
	    except:
	        print(str(i))
	        continue
	    pos['galaxy'+str(i)] = {}


	    gas_masses = infile['PartType0']['Masses']
	    gas_coords = infile['PartType0']['Coordinates']
	    star_masses = infile['PartType4']['Masses']
	    star_coords = infile['PartType4']['Coordinates']
	    total_mass = np.sum(gas_masses) + np.sum(star_masses)

	    x_pos = (np.sum(gas_masses * gas_coords[:,0]) + np.sum(star_masses * star_coords[:,0])) / total_mass
	    y_pos = (np.sum(gas_masses * gas_coords[:,1]) + np.sum(star_masses * star_coords[:,1])) / total_mass
	    z_pos = (np.sum(gas_masses * gas_coords[:,2]) + np.sum(star_masses * star_coords[:,2])) / total_mass
	
	
	    pos['galaxy'+str(i)]['snap'+str(snap)] = np.array([x_pos, y_pos, z_pos])
	    infile.close()
	ngalaxies['snap'+str(snap)] = count


	print("SAVING")
	np.savez(outfile, ngalaxies=ngalaxies, pos=pos)




Getting POWDERDAY to Run on a Single Filtered Galaxy
-----------------



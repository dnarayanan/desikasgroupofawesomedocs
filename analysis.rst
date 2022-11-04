Analysis Procedures
**********


Generating SFHs from Snapshots
============

(so far) There's two ways of calculating the star formation rate for simba galaxies which differ by the timescales the SFR is averaged over. The SFR calculated by caesar for galaxies is based on the current rate of star formation of the gas particles. This is an instantaneous measure of SFR -- literally, how many stars are being produced from gas particles at the moment the snapshot is written. A more general way to quantify the SFR is by averaging the stellar mass produced over some time period. This is more akin to the SFR that observers measure from galaxies, either by relating emission line luminosities to the SFR or by modeling galaxy SEDs with a SFH. Caesar does not currently calculate this SFR for us so I've written a script to do just that. 

The basic idea of the script is to take all the star particles for a particular object (i.e., galaxy) and bin them up according to their formation time and initial mass. While the formation time is an attribute of the star particle that is written upon its creation, the initial mass is not recorded. Thus, we need to use FSPS to estimate the initial mass based on the star's current mass, metallicity, and age using isochrone/evolution tracks. 

The output of this script is an object containing the formation times and initial masses of each star particle in the galaxy/galaxies of interest. To turn this into a SFH, I use a binning function from scipy that returns the SFR in arbitrary spaced bins; this allows us to make SFHs/calculate the recent SFR over whatever timescales we want. 

SFH Script
--------------

The script I use for generating SFHs for simba (or any cosmo sim), using FSPS to calculate the initial star particle masses, can be found at this link::
  
  <https://github.com/smlower/sl_simulation_tools/blob/main/caesar_sfh.py>

Which is just a repo where I keep generally useful scripts for dealing with simulations. The current setup for the script is for it to calculate the SFHs for all galaxies in the simba m25n512 z~5 snapshot. Lines 50-52 are where I make the list of galaxies I want SFHs for. The script can be edited to accept, e.g., command line arguments for snapshot, galaxy number, etc. 

The output is a pickled object containing a list of galaxy IDs and those galaxies' initial star particle masses and formation times. One thing to note is that FSPS takes some time to start up...once the isochrone grids are sampled, the initial mass interpolation speeds up siginificantly. 


Binning and plotting the SFH
--------------

The list output from the above script is not immediately useful, as we need to bin those masses according to their formation times to get the SFH. This can be done tons of ways but the easiest (at least to me) way is using the scipy binned_statistic function. The code is (thankfully) documented a bit so it should be straightfoward::

  <https://github.com/smlower/sl_simulation_tools/blob/main/bin_sfh.py>


The good thing about binning this way is that we don't have to re-run FSPS everytime we want a differently sampled SFH -- just choose a different bin size!



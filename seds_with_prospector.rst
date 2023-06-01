Intro to SED Fitting with Prospector
**************************************
A brief intro to SED modeling and using prospector. I cover two things in this tutorial: (1) basics of SED modeling that apply to any SED fitting code and (2) Bayesian inference (to the best of my ability). 

.. contents:: Section Contents
    :local:

A (decently written) tutorial is available on github that walks through the basics of installing propsector and its dependencies, setting up a prospector run, and doing some analysis with the output. You can find that repo, which contains a tutorial notebook and the needed data files, here::

  https://github.com/smlower/prospector_tutorial

Besides the tutorial on github, the installation instructions can also be found at::

  https://github.com/dnarayanan/desikasgroupofawesome/blob/main/using_major_group_codes.rst#prospector

Step 1: What is SED fitting?
============================

SED modeling basically entails mapping the observed light of a galaxy at different wavelengths to their physical properties like stellar mass, star formation rate, and dust mass. There's a bunch of different codes that do this with varying levels of sophistication: FAST, CIGALE, MAGPHYS, ProSpect, Prospector, BAGPIPES, and many, many more.

To model a galaxy spectrum and infer the galaxy's properties, there's a bunch of assumptions we have to make about how and when the stars formed in a galaxy (the star formation history), how much and what kind of dust is in that galaxy (dust attenuation and emisision), the properties of regions like HII nebulae and AGN. For the purposes of this tutorial, we're going to focus on just the components pertaining to stars and dust. Below is a great summary figure from Charlie Conroy's 2013 review article on the technique of modeling stellar populations. 

.. image:: images/conroy_2013.png
	   :width: 600


Modeling Stellar Populations
------------------------------
The top row concerns modeling stellar evolution. The IMF represents the stellar mass distribution of stars in a galaxy -- literally how many stars of each mass are born in a galaxy. There's a lot of debate over what the IMF is in various regions of our own Galaxy, let alone other galaxies, but for our purposes, we just assume one of the typical distributions from the literature like Kroupa 2002 or Chabrier 2003. 

The isochrones are models detailing how stars move about the HR diagram as a function of time. Each line in that plot represents the distribution of stars for a fixed time. Stars that are born at the same time, but of different masses, will occupy different regions in color-magnitude space, impacting the colors of a galaxy. 

The right-hand plot shows spectra for a few different spectral types. Combine these three ingredients and we get something called a simple stellar population, or SSP, shown in the middle panel of the middle row. A simple stellar population represents the spectrum of a group of stars that are the same age. Think of it as taking a picture of a globular cluster: those stars formed at ~the same time but are of different masses and so the combined spectrum of that cluster will depend on how many high/low mass (hot/cold temp) stars there are and the age of that cluster.

But a galaxy is definitely not just a globular cluster. It's more of a bunch of globular clusters of different ages put together. So how do we describe that spectrum? We combine it with a model for the star formation history of the galaxy -- literally, when did those globular clusters form? The left panel in the middle row shows two examples of a star formation history model (let's ignore the bottom plot showing metallicity). SFHs of real galaxies are diverse from galaxy-to-galaxy and have different short- and long-term variability. Convolving a model for SFH with an SSP gives us a compsite stellar spectrum combining the effects of stellar spectral types and stellar age. We can see an example of this spectrum in the blue curve plotted in the bottom row.


Dust Attenuation and Emission
-------------------------------
Finally, to the shagrin of most UV and optical astronomers, a galaxy is not just stars: it also contains dust. Dust modulates the spectrum of a galaxy by absorbing the UV and optical light of stars and re-emitting that light as a blackbody in the far-infrared. This absorption (+ scattering) is called dust attenuation and is modeled with a dust attenuation curve. In the right panel of the middle row, we see an example of this where tau (optical depth) is plotted as a function of wavelength. Dust preferentially absorbs light (i.e., higher optical depths) in the UV and not as much in the red / near-IR. There's several features to a dust attenuation curve but for now, we'll stick with using one of the most commonly used attenuation curves in the literature: the Calzetti 2001 curve, which is shown in red in the plot. The plot right underneath shows the resulting dust emission spectra. Focusing on the peak of the curve around 100 micron, this spectrum is the result of the dust absorbing light in the UV, heating up, and re-emitting that light as a thermal blackbody. The weird features around 10 micron are the result of a special kind of dust called PAHs (polycyclic aromatic hydrocarbon) -- these are fascinating but we'll ignore them for now. The shape of the blackbody curve is dependent on the amount of dust and the dust temperature.

Putting The Pieces Together
-----------------------------
Combining the two rows, we get a 'composite' galaxy spectra, including the contributions of stars of all ages and dust, shown in red in the bottom plot. This is the general shape of a galaxy SED: the stars dominate in the UV and optical while the dust emission dominates in the mid- to far-infrared. 

To extract information from a galaxy SED, we essentially do the above process backwards: what star and dust spectra are the best fit to the observed SED, from which we can derive the properties of the galaxy? To figure that out, we select models for the stellar evolution, star formation history, and dust. The stellar evolution models are typically fixed (i.e., we choose one model set and stick with it), but the parameters of the star formation history model and the dust attenuation/emission models can vary -- this represents the basis of our MCMC problem: what combination of model parameters give us a best fit to our data? To perform this fit, we'll use fsps+dynesty+prospector.

Step 2: How to fit data?
========================
The basic idea anytime we want to fit a model to data is to literally minimize the difference between the model and the data. In its most basic form, this means generating a model SED and calculating the chi square statistic, with the 'best fit' model having the lowest chi square. In a more sophisticated form, this involves Bayesian inference. I'm never going to do an explanaition of Bayesian statistics justice, so if you're super interested in learning the mechanics of this, I suggest doing outside readins. Regardless, the basis of Bayesian inference is that we have some 'prior' knowledge that we can use to construct the probability distribution of model parameters, which can in turn be used to construct a sort-of best fit model SED. Bayesian inference comes from Bayes theorem (which I'll point to the wiki page for more info: https://en.wikipedia.org/wiki/Bayes%27_theorem) which says the probability distribution of a model parameter (called the posterior distribution) is related to the likelihood of that model parameter * the prior distribution of that model parameter. If you've ever read a paper about model fitting or listened to a colloquium about deriving properties of something from a model, this is where 'prior,' 'posterior,' and 'likelihood' come from.

What this means in practice is that for any set of models we choose for our SED components (star formation history, dust), the model SED is evaluated based on the prior knowledge of the distribution of model parameters and the likelihood of that model parameter representing the true data. For our purposes, the likelihood function is taken care of in the internals of prospector/dynesty. Thus for each variable model parameter, we will choose a prior distribution based on our knowledge of that parameter. Literally, what are the physical or known values this model parameter can take? An example is the age of a galaxy: we know that a galaxy has to have an age greater than zero and less than the age of the universe. Now, priors can have any degree of complexity but most of the time we will use an 'uninformative' prior, i.e., a prior that does not impose a lot of weight on the posterior distribution of the model paramter. An example of an uninformative prior is a uniform distribution, and for the age of the galaxy, the prior would range from 0 to 14 Gyr with every value in between having equal probability within the prior space. An example of an 'informative' prior would be a Gaussian, where galaxy ages around the mean of the Gaussian would have greater weight than ages close to the wings. Neither prior distribution is necessarily 'wrong,' (and believe me, there's tons of discussion on the intricacies of choosing priors) and generally is entirely dependent on the information/data we have and the problem we are trying to solve.

Specifically for prospector, which we'll see below, we don't interface with the actual Bayesian inference at all, besides the initial selection of models and the choices for model priors. After the data has been fit, what we'll have as a result are posterior distributions for each model parameter. In cases where the data is not constraining or is not fit very well, these posterior distributions will resemble the prior distributions, basically a null result. But most of the time, we'll get back posterior distributions that resemble a Gaussian from which we can report the median value +/- the variability -- this is usually what's reported in publications. From this point, we can discuss things like maximum likelihood estimates vs. medians and degeneracies but that's probably outside the scope of just getting started with prospector.

Step 3: Demo
=============
With the above in mind and prospector and its dependencies successfully installed, we're ready to test out our setup with some data! From here, you can follow through the tutorial at https://github.com/smlower/prospector_tutorial/blob/main/tutorial.ipynb, which first fits a prospector model to a toy galaxy model dataset and then fits a more 'realistic' galaxy SED generated by simba and powderday.   There is also example code to model the physical properties from a powderday run here.  This run uses a non-parametric SFH::

  from prospect.models import priors, transforms #helper functions for specifying priors
  import numpy as np

  from prospect.models import sedmodel
  from prospect.io import write_results as writer
  from prospect.fitting import fit_model
  
  import sys


  def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx

  def build_model(**kwargs):

    """
    Function to build model components for SFH and dust.
    The model params are defined by their name, whether they are a free parameter
    their initial value, and their prior distribution if they are variable. The model
    params are then fed to the prospector SedModel class

    All parameters except 'mass' correspond to fsps model parameters, the definitions of which you can find here:
    https://dfm.io/python-fsps/current/stellarpop_api/

    """

    model_params = []
    #redshift of galaxy -- can be a free parameter but we'll fix it for now
    model_params.append({'name': "zred", "N": 1, "isfree": False,"init": 0.01, 'prior': None})
    #IMF model which will be used by the simple stellar population model
    model_params.append({'name': 'imf_type', 'N': 1,'isfree': False,'init': 2, 'prior': None})
    #stellar mass of a galaxy -- what we're interested in! So we'll set it as a free parameter
    model_params.append({'name': 'logmass', 'N': 1,'isfree': True, 'init': 10,'prior': priors.TopHat(mini=8, maxi=12)})
    #stellar metallicity, in units of log(Z/Z_sun)
    model_params.append({'name': 'logzsol', 'N': 1,'isfree': True,'init': -0.5,'prior': priors.TopHat(mini=-1.6, maxi=0.1)})

    #SFH model. here, we are choosing the 'delayed-tau' model and has two free parameters: the age and the e-folding time
    #model_params.append({'name': "sfh", "N": 1, "isfree": False, "init": 4, 'prior': None})
    #age of the galaxy
    #model_params.append({'name': "tage", 'N': 1, 'isfree': True, 'init': 5., 'units': 'Gyr', 'prior': priors.TopHat(mini=0.001, maxi=13.8)})
    #e-folding time
    #model_params.append({'name': "tau", 'N': 1, 'isfree': True,'init': 1., 'units': 'Gyr', 'prior': priors.LogUniform(mini=0.1, maxi=30)})


    #non parametric SFH
    n = [p['name'] for p in model_params]
    tuniv = 14.0
    nbins=6
    tbinmax = (tuniv * 0.85) * 1e9
    lim1, lim2 = 8.0, 8.52 #100 Myr and 330 Myr
    agelims = [0,lim1] + np.linspace(lim2,np.log10(tbinmax),nbins-2).tolist() + [np.log10(tuniv*1e9)]
    agebins = np.array([agelims[:-1], agelims[1:]])

    #model_params.append({'name': "agebins",'agebins':agebins})

    ncomp = nbins
    alpha_sfh = 0.7
    alpha = np.repeat(alpha_sfh,nbins-1)
    tilde_alpha = np.array([alpha[i-1:].sum() for i in range(1,ncomp)])
    zinit = np.array([(i-1)/float(i) for i in range(ncomp, 1, -1)])
    zprior = priors.Beta(alpha=tilde_alpha, beta=alpha, mini=0.0, maxi=1.0)


    model_params.append({'name': "sfh", "N": 1, "isfree": False, "init": 3})
    model_params.append({'name': "mass", 'N': 3, 'isfree': False, 'init': 1., 'depends_on':zfrac_to_masses_log})
    model_params.append({'name': "agebins", 'N': 1, 'isfree': False,'init': []})
    model_params.append({'name': "z_fraction", "N": 2, 'isfree': True, 'init': [0, 0],'prior': priors.Beta(alpha=1.0, beta=1.0, mini=0.0, maxi=1.0)})

    n = [p['name'] for p in model_params]
    model_params[n.index('mass')]['N'] = ncomp
    model_params[n.index('agebins')]['N'] = ncomp
    model_params[n.index('agebins')]['init'] = agebins.T
    model_params[n.index('z_fraction')]['N'] = len(zinit)
    model_params[n.index('z_fraction')]['init'] = zinit
    model_params[n.index('z_fraction')]['prior'] = zprior


    #model_params[n.index('massmet')]['prior'] = MassMet(z_mini=-1.6, z_maxi=0.2, mass_mini=8.0, mass_maxi=12.)

    #dust attenuation model, cardelli
    model_params.append({'name': 'dust_type', 'N': 1,'isfree': False,'init': 1,'prior': None})
    model_params.append({'name': 'mwr', 'N': 1,'isfree': False, 'init': 3.1,'prior': None})
    model_params.append({'name': 'uvb', 'N': 1,'isfree': False,'init': 0,'prior': None})
    #the attenuation (in magnitudes) in the V-band
    model_params.append({'name': 'dust2', 'N': 1,'isfree': True, 'init': 1.0,'prior': priors.Uniform(mini=0.0, maxi=3.0)})

    #for calzetti
    #model_params.append({'name': 'dust2', 'N': 1,'isfree': True, 'init': 0.1,'prior': priors.ClippedNormal(mini=0.0, maxi=2.0, mean=0.0, sigma=0.3)})
    #dust emission model -- only 1 choice, from Draine & Li 2007
    model_params.append({'name': 'add_dust_emission', 'N': 1,'isfree': False,'init': 1,'prior': None})
    #mass fraction of warm dust
    model_params.append({'name': 'duste_gamma', 'N': 1,'isfree': True,'init': 0.01,'prior': priors.TopHat(mini=0.0, maxi=1.0)})
    #minimum radiation field
    model_params.append({'name': 'duste_umin', 'N': 1,'isfree': True,'init': 1.0,'prior': priors.TopHat(mini=0.1, maxi=20.0)})
    #mass fraction of dust in PAHs
    model_params.append({'name': 'duste_qpah', 'N': 1,'isfree': True,'init': 3.0,'prior': priors.TopHat(mini=0.0, maxi=6.0)})


    model = sedmodel.SedModel(model_params)
    return model


  def zfrac_to_masses_log(logmass=None, z_fraction=None, agebins=None, **extras):
    sfr_fraction = np.zeros(len(z_fraction) + 1)
    sfr_fraction[0] = 1.0 - z_fraction[0]
    for i in range(1, len(z_fraction)):
        sfr_fraction[i] = np.prod(z_fraction[:i]) * (1.0 - z_fraction[i])
    sfr_fraction[-1] = 1 - np.sum(sfr_fraction[:-1])
    time_per_bin = np.diff(10**agebins, axis=-1)[:, 0]
    mass_fraction = sfr_fraction * np.array(time_per_bin)
    mass_fraction /= mass_fraction.sum()

    masses = 10**logmass * mass_fraction
    return masses


  from prospect.sources import CSPSpecBasis
  def build_sps(zcontinuous=1,compute_vega_mags=False,**kwargs):
    """
    This is our stellar population model which generates the spectra for stars of a given age and mass.
    Most of the time, you aren't going to need to pay attention to this.
    """

    from prospect.sources import FastStepBasis
    sps = FastStepBasis(zcontinuous=zcontinuous,
                        compute_vega_mags=compute_vega_mags)
    #sps = CSPSpecBasis(zcontinuous=1)
    return sps
  
    #---------------------
    # Setup Observations
   
    #---------------------
    
    galex = ['galex_FUV', 'galex_NUV']
    hst_wfc3_uv  = ['wfc3_uvis_f275w', 'wfc3_uvis_f336w', 'wfc3_uvis_f475w','wfc3_uvis_f555w', 'wfc3_uvis_f606w', 'wfc3_uvis_f814w']
    hst_wfc3_ir = ['wfc3_ir_f105w', 'wfc3_ir_f125w', 'wfc3_ir_f140w', 'wfc3_ir_f160w']
    spitzer_mips = ['spitzer_mips_24']
    herschel_pacs = ['herschel_pacs_70', 'herschel_pacs_100', 'herschel_pacs_160']
    herschel_spire = ['herschel_spire_250', 'herschel_spire_350', 'herschel_spire_500']
        
    jwst_nircam = ['jwst_f115w', 'jwst_f150w', 'jwst_f200w', 'jwst_f277w', 'jwst_f356w', 'jwst_f444w']
    jwst_miri = ['jwst_f560w', 'jwst_f770w', 'jwst_f1000w', 'jwst_f1280w', 'jwst_f1500w', 'jwst_f1800w', 'jwst_f2100w']
    filternames = hst_wfc3_uv+hst_wfc3_ir


  #------------------
  # Build Observations
  #-------------------
  
  def build_obs(pd_dir,**kwargs):

    from sedpy.observate import load_filters

    print('loading obs')
    import sedpy
    from astropy import units as u
    from astropy import constants
    from astropy.cosmology import FlatLambdaCDM
    from hyperion.model import ModelOutput
    cosmo = FlatLambdaCDM(H0=68, Om0=0.3, Tcmb0=2.725)
    m = ModelOutput(pd_dir)
    wav,flux = m.get_sed(inclination=0,aperture=-1)
    wav  = np.asarray(wav)*u.micron #wav is in micron
    

    wav = wav.to(u.AA)
    flux = np.asarray(flux)*u.erg/u.s
    dl = cosmo.luminosity_distance(z_spec).to('cm')
    flux /= (4.*3.14*dl**2.)
    nu = constants.c.cgs/(wav.to(u.cm))
    nu = nu.to(u.Hz)
    flux /= nu
    flux = flux.to(u.Jy)
    maggies = flux / 3631.

    filters_unsorted = load_filters(filternames)
    waves_unsorted = [x.wave_mean for x in filters_unsorted]
    filters = [x for _,x in sorted(zip(waves_unsorted,filters_unsorted))]
    flx = []
    flxe = []
    redshifted_wav = wav*(1.+z_spec)
    for i in range(len(filters)):
        flux_range = []
        wav_range = []
        for j in filters[i].wavelength:
            flux_range.append(maggies[find_nearest(redshifted_wav.value,j)].value)
            wav_range.append(redshifted_wav[find_nearest(redshifted_wav.value,j)].value)
        a = np.trapz(wav_range * filters[i].transmission* flux_range, wav_range, axis=-1)
        b = np.trapz(wav_range * filters[i].transmission, wav_range)
        flx.append(a/b)
        flxe.append(0.03* flx[i])
    flx = np.asarray(flx)
    flxe = np.asarray(flxe)
    flux_mag = flx
    unc_mag = flxe

    obs = {}
    obs['filters'] = filters
    obs['maggies'] = flux_mag
    obs['maggies_unc'] = unc_mag
    obs['phot_mask'] = np.isfinite(flux_mag)
    obs['wavelength'] = None
    obs['spectrum'] = None
    obs['pd_sed'] = maggies
    obs['pd_wav'] = redshifted_wav

    return obs


  def get_sfr100(res, mod):
    agebins = mod.params['agebins']
    thetas = mod.theta_labels()
    agebins_yrs = 10**agebins.T
    dt = agebins_yrs[1, :] - agebins_yrs[0, :]

    zfrac_idx = [i for i, s in enumerate(thetas) if 'z_fraction' in s]
    zfrac_chain = res['chain'][:,zfrac_idx[0]:zfrac_idx[-1]+1]
    try:
        total_mass_chain = res['chain'][:,thetas.index('massmet_1')]
    except:
        total_mass_chain = res['chain'][:,thetas.index('logmass')]
    sfr_chain = []
    for i in range(len(zfrac_chain)):
        #this is the important part of this function -- it takes the proxy model parameters we sampled
        # and transforms them into the masses formed in each timebin
        # from there, we just divide these masses by the time in each bin to get the SFR in each bin
        masses_chain = transforms.zfrac_to_masses(10**total_mass_chain[i], zfrac_chain[i], agebins)
        sfr = masses_chain / dt
        sfr_chain.append(sfr[0])
    return sfr_chain

  if __name__ == '__main__':
    pd_sed = '/orange/narayanan/s.lower/simba/pd_runs/snap305/snap305.galaxy'+str(sys.argv[1])+'.rtout.sed'
    print(pd_sed)
    z_spec = 0.01


    #First, set some runtime parameters. These define important parameters for dynesty like how many live points (akin)
    #to emcee's walkers, and the stopping condition dlogz. We'll mostly leave these alone for now
    run_params = {'output_pickles': False, #our output will be in hdf5 format
                  # dynesty Fitter parameters
                  
                  'nested_bound': 'multi',
		  'nested_sample': 'auto',                                                                                                                    

                  'nested_nlive_init': 400,
                  'nested_nlive_batch': 200,
                  'nested_bootstrap': 0,
                  'nested_dlogz_init': 0.05,
                  'nested_weight_kwargs': {"pfrac": 1.0},
              }

    mod = build_model(z_spec=z_spec)
    print(mod)
    sps = build_sps()
    obs = build_obs(pd_sed)
    output = fit_model(obs, mod, sps, **run_params)

    galaxy_num = sys.argv[1]

    out_file = 'outfiles/dum.uvb0.gal'+str(galaxy_num)+'.h5'
    writer.write_hdf5(out_file, run_params, mod, obs,
                      output["sampling"][0], output["optimization"][0],
                      tsample=output["sampling"][1],
                      toptimize=output["optimization"][1])


Okay but how might we actually use these files to get a physical
property?  This depends a bit on the physical property.  First, let's
do the most complicated one -- SFH.  There are different functions for different priors on the SFH unfortunately (I personally stuff all the below into prospector_sfh_functions.py)::

  import prospect.io.read_results as pread
  from prospector_test import build_model
  from prospect.models import priors, transforms #helper functions for specifying priors
  from prospect.models import sedmodel
  from corner import quantile
  import matplotlib.pyplot as plt
  import numpy as np
  from glob2 import glob
  
  import numpy as np
  from scipy.special import gamma, gammainc

  def get_sfh(res, mod):
    agebins = mod.params['agebins']
    thetas = mod.theta_labels()
    agebins_yrs = 10**agebins.T
    bin_edges = np.unique(agebins_yrs)
    dt = agebins_yrs[1, :] - agebins_yrs[0, :]
    epsilon = 1e-4 #fudge factor used to define the fraction time separation of adjacent points at the bin edges
    t = np.concatenate((bin_edges * (1.-epsilon), bin_edges * (1+epsilon)))
    t.sort()
    t = t[1:-1]
    zfrac_idx = [i for i, s in enumerate(thetas) if 'z_fraction' in s]
    zfrac_chain = res['chain'][:,zfrac_idx[0]:zfrac_idx[-1]+1]
    try:
        total_mass_chain = res['chain'][:,thetas.index('massmet_1')]
    except:
        total_mass_chain = res['chain'][:,thetas.index('logmass')]
    sfr_chain = []
    for i in range(len(zfrac_chain)):
        masses_chain = transforms.zfrac_to_masses(10**total_mass_chain[i], zfrac_chain[i], agebins)
        sfr = masses_chain / dt
        sfrout = np.zeros_like(t)
        sfrout[::2] = sfr
        sfrout[1::2] = sfr
        sfr_chain.append(sfrout)
    return (14. - t/1e9), sfr_chain


  def SL_psb_logsfr_ratios_to_agebins(logsfr_ratios=None, agebins=None,
                                 tlast=None, tflex=None, nflex=None, nfixed=None, **extras):
    """This is a modified version of logsfr_ratios_to_agebins above. This now
    assumes that there are nfixed fixed-edge timebins at the beginning of
    the universe, followed by nflex flexible timebins that each form an equal
    stellar mass. The final bin has variable width and variable SFR; the width
    of the bin is set by the parameter tlast.
    For the flexible bins, we again use the equation:
        delta(t1) = tuniv  / (1 + SUM(n=1 to n=nbins-1) PROD(j=1 to j=n) Sn)
        where Sn = SFR(n) / SFR(n+1) and delta(t1) is width of youngest bin
    """


    # numerical stability
    logsfr_ratios = np.clip(logsfr_ratios, -7, 7)

    # flexible time is t_flex - youngest bin (= tlast, which we fit for)
    # this is also equal to tuniv - upper_time - lower_time
    tf = (tflex - tlast) * 1e9

    # figure out other bin sizes
    n_ratio = logsfr_ratios.shape[0]
    sfr_ratios = 10**logsfr_ratios
    dt1 = tf / (1 + np.sum([np.prod(sfr_ratios[:(i+1)]) for i in range(n_ratio)]))

    # translate into agelims vector (time bin edges)
    agelims = [1, (tlast*1e9), dt1+(tlast*1e9)]
    for i in range(n_ratio):
        agelims += [dt1*np.prod(sfr_ratios[:(i+1)]) + agelims[-1]]
    agelims += list(10**agebins[-nfixed:,1])
    abins = np.log10([agelims[:-1], agelims[1:]]).T

    return abins



  def SL_logsfr_ratios_to_masses_psb(logmass=None, logsfr_ratios=None,
                                 logsfr_ratio_young=None, logsfr_ratio_old=None,
                                 tlast=None, tflex=None, nflex=None, nfixed=None,
                                 agebins=None, **extras):
    """This is a modified version of logsfr_ratios_to_masses_flex above. This now
    assumes that there are nfixed fixed-edge timebins at the beginning of
    the universe, followed by nflex flexible timebins that each form an equal
    stellar mass. The final bin has variable width and variable SFR; the width
    of the bin is set by the parameter tlast.
    The major difference between this and the transform above is that
    logsfr_ratio_old is a vector.
    """

    # clip for numerical stability
    logsfr_ratio_young = np.clip(logsfr_ratio_young[0], -7, 7)
    logsfr_ratio_old = np.clip(logsfr_ratio_old, -7, 7)
    syoung, sold = 10**logsfr_ratio_young, 10**logsfr_ratio_old
    sratios = 10.**np.clip(logsfr_ratios, -7, 7) # numerical issues...

    # get agebins
    abins = SL_psb_logsfr_ratios_to_agebins(logsfr_ratios=logsfr_ratios,
            agebins=agebins, tlast=tlast, tflex=tflex, nflex=nflex, nfixed=nfixed, **extras)

    # get find mass in each bin
    dtyoung, dt1 = (10**abins[:2, 1] - 10**abins[:2, 0])
    dtold = 10**abins[-nfixed-1:, 1] - 10**abins[-nfixed-1:, 0]
    old_factor = np.zeros(nfixed)
    for i in range(nfixed):
        old_factor[i] = (1. / np.prod(sold[:i+1]) * np.prod(dtold[1:i+2]) / np.prod(dtold[:i+1]))
    mbin = 10**logmass / (syoung*dtyoung/dt1 + np.sum(old_factor) + nflex)
    myoung = syoung * mbin * dtyoung / dt1
    mold = mbin * old_factor
    n_masses = np.full(nflex, mbin)

    return np.array([myoung] + n_masses.tolist() + mold.tolist())


  
  def get_sfr_psb(res, mod):
    logmass=res['chain'][:,0]
    logsfr_ratios_chain = res['chain'][:,7:11]
    logsfr_ratios_young_chain = res['chain'][:,3]
    logsfr_ratios_old_chain = res['chain'][:,4:7]
    tlast_chain = res['chain'][:,2]
    tflex = 0.37112653
    nflex=5
    nfixed=3
    agebins = mod.params['agebins']
    sfr_chain = []
    time_chain = []

    weights = res.get('weights',None)
    idx = np.argsort(weights)[-3000:]

    for i in idx:
        #print(logsfr_ratios_young_chain[i])
        masses = SL_logsfr_ratios_to_masses_psb(logmass=logmass[i], logsfr_ratios=logsfr_ratios_chain[i],
                                 logsfr_ratio_young=[logsfr_ratios_young_chain[i]],
                                logsfr_ratio_old=logsfr_ratios_old_chain[i],
                                 tlast=tlast_chain[i], tflex=tflex, nflex=nflex, nfixed=nfixed,
                                 agebins=agebins)


        fit_agebins = SL_psb_logsfr_ratios_to_agebins(logsfr_ratios=logsfr_ratios_chain[i], agebins=agebins,
                                 tlast=tlast_chain[i], tflex=tflex, nflex=nflex, nfixed=nfixed)

        #print(10**fit_agebins)
        #dt = (10**fit_agebins[:, 1] - 10**fit_agebins[:, 0])
        agebins_yrs = 10**fit_agebins.T
        bin_edges = np.unique(agebins_yrs)
        dt = agebins_yrs[1, :] - agebins_yrs[0, :]
        epsilon = 1e-4 #fudge factor used to define the fraction time separation of adjacent points at the bin edges
        t = np.concatenate((bin_edges * (1.-epsilon), bin_edges * (1+epsilon)))
        t.sort()
        t = t[1:-1]

        #print(f'age bin years {agebins_yrs}')
        #print(f'bin edges {bin_edges}')
        #print(f't {t}')
        #print(f'dt {dt}')
        #print(f'masses {masses}')
        sfr = masses/dt
        sfrout = np.zeros_like(t)
        sfrout[::2] = sfr
        sfrout[1::2] = sfr
        sfr_chain.append(sfrout)
        time_chain.append((t[-1] - t)/1e6)
    return t/1.e9,[item for item in sfr_chain]


  def get_sfr_beta(res, mod): #rising SFH from Wang, Labbe et al.
    #agebins = mod.params['agebins']                                                                                              

    agebins = transforms.zred_to_agebins_pbeta(np.array([7.2]))
    thetas = mod.theta_labels()
    agebins_yrs = 10**agebins.T
    bin_edges = np.unique(agebins_yrs)
    dt = agebins_yrs[1, :] - agebins_yrs[0, :]
    epsilon = 1e-4 #fudge factor used to define the fraction time separation of adjacent points at the bin edges 

    t = np.concatenate((bin_edges * (1.-epsilon), bin_edges * (1+epsilon)))
    t.sort()
    t = t[1:-1]
    weights = res.get('weights',None)
    idx = np.argsort(weights)[-3000:]
    sfh_bin_idx = [i for i, s in enumerate(thetas) if 'nzsfh' in s]
    sfr_chain = []
    for i in idx:
        sfh_bin_chain = res['chain'][i,3:sfh_bin_idx[-1]+1]
        total_mass_chain = res['chain'][i,1]
        sfr = transforms.logsfr_ratios_to_sfrs(total_mass_chain, sfh_bin_chain, agebins)
        sfrout = np.zeros_like(t)
        sfrout[::2] = sfr
        sfrout[1::2] = sfr
        #sfrout[::2] = sfr                                                                                                            
        #sfrout[1::2] = sfr                                                                                                           
        #sfr_chain.append(sfr[0])
        sfr_chain.append(sfrout)

    return (t[-1] - t[::-1])/1e9,sfr_chain

    

  def get_sfr_dirichlet(res, mod):
    agebins = mod.params['agebins']
    thetas = mod.theta_labels()
    agebins_yrs = 10**agebins.T
    bin_edges = np.unique(agebins_yrs)
    dt = agebins_yrs[1, :] - agebins_yrs[0, :]
    epsilon = 1e-4 #fudge factor used to define the fraction time separation of adjacent points at the bin edges
    t = np.concatenate((bin_edges * (1.-epsilon), bin_edges * (1+epsilon)))
    t.sort()
    t = t[1:-1]
    zfrac_idx = [i for i, s in enumerate(thetas) if 'z_fraction' in s]
    zfrac_chain = res['chain'][:,zfrac_idx[0]:zfrac_idx[-1]+1]
    try:
        total_mass_chain = res['chain'][:,thetas.index('massmet_1')]
    except:
        total_mass_chain = res['chain'][:,thetas.index('logmass')]
    sfr_chain = []
    for i in range(len(zfrac_chain)):
        masses_chain = transforms.zfrac_to_masses(10**total_mass_chain[i], zfrac_chain[i], agebins)
        sfr = masses_chain / dt
        sfrout = np.zeros_like(t)
        sfrout[::2] = sfr
        sfrout[1::2] = sfr
        sfr_chain.append(sfrout)
    return (t[-1] - t[::-1])/1e9, sfr_chain

  def delayed_tau_sfr(t,tau):
    return (t/tau) * np.exp(-t/tau)


  def get_sfr_tau(res,mod):

    mass = res['chain'][:,model_thetas.index('mass')]
    tage = res['chain'][:,model_thetas.index('tage')]
    tau = res['chain'][:,model_thetas.index('tau')]

    print(np.median(mass),np.median(tage),np.median(tau))

    # for delay tau this function gives the (unnormalized) SFR
    # for any t, tau combo in M_sun/Gyr


    #sfr = lambda t,tau: (t/tau) * np.exp(-t/tau)

    #sfr = lambda t,tau: return (t/tau) * np.exp(-t/tau)
    # now we numerically integrate this SFH from 0 to tage to get the mass formed
    times = np.linspace(0, tage, 1000)
    A = np.trapz(delayed_tau_sfr(times, tau), times)
    # But this could also be done using an incomplete gamma function (integral of xe^{-x})
    A = tau * gamma(2) * gammainc(2, tage/tau)
    # and now we renormalize the formed mass to the actual mass value
    # to get the the SFR in M_sun per Gyr
    psi = mass * delayed_tau_sfr(tage, tau) / A
    # if we want SFR in Msun/year
    psi /= 1e9

    return psi

And then with these functions in hand, we can do (for example)::

  sfr = []
  res,obs,mod = pread.results_from(uvb0file)
  t_sfr_prosp,sfh = get_sfr_dirichlet(res,mod)
  sfh_50, sfh_16, sfh_84 = [], [], []


  for time in range(len(t_sfr_prosp)):
        sfr_quan = quantile([item[time] for item in sfh], [0.16, 0.5, 0.84])
        sfh_50.append(sfr_quan[1])
        sfh_16.append(sfr_quan[0])
        sfh_84.append(sfr_quan[2])

        sfr.append(sfh_50[-1])


We can get other physical quantities such as stellar mass much easier::

  res,obs,mod = pread.results_from(uvb1file)
  z_spec=0.01
  mod = build_model(z_spec=z_spec)

  model_thetas = mod.theta_labels()
  stellar_mass = res['chain'][:,model_thetas.index('logmass')]


Generally, you can print the model object to see what other physical properties we can derive from the model_thetas::

  Free Parameters: (name: prior) 
  -----------
  logmass: <class 'prospect.models.priors.TopHat'>(mini=9.0,maxi=12.0)
  logzsol: <class 'prospect.models.priors.TopHat'>(mini=-1.0,maxi=0.2)
  z_fraction: <class 'prospect.models.priors.Beta'>(mini=0.0,maxi=1.0,alpha=[6.3 5.6 4.9 4.2 3.5 2.8 2.1 1.4 0.7],beta=[1. 1. 1. 1. 1. 1. 1. 1. 1.])
  mwr: <class 'prospect.models.priors.TopHat'>(mini=1.5,maxi=5.0)
  uvb: <class 'prospect.models.priors.TopHat'>(mini=0.0,maxi=4.0)
  dust2: <class 'prospect.models.priors.TopHat'>(mini=0.0,maxi=2.0)
  duste_gamma: <class 'prospect.models.priors.TopHat'>(mini=0.0,maxi=1.0)
  duste_umin: <class 'prospect.models.priors.TopHat'>(mini=0.1,maxi=30.0)
  duste_qpah: <class 'prospect.models.priors.TopHat'>(mini=0.0,maxi=10.0)

  Fixed Parameters: (name: value [, depends_on]) 
  -----------
  lumdist: [1.e-05] 
  pmetals: [-99] 
  imf_type: [2] 
  sfh: [3] 
  mass: [1.] <function zfrac_to_masses_log at 0x2ac0b0e02160>
  agebins: [[ 0.          7.        ]
  [ 7.          8.        ]
  [ 8.          8.29650671]
  [ 8.29650671  8.59301342]
  [ 8.59301342  8.88952013]
  [ 8.88952013  9.18602684]
  [ 9.18602684  9.48253354]
  [ 9.48253354  9.77904025]
  [ 9.77904025 10.07554696]
  [10.07554696 10.14612804]] 
  dust_type: [1] 
  add_dust_emission: [1] 
  add_agb_dust_model: [0] 




Finally, you may wish to compare the true SFH to the model SFH.  This is pretty easy now::

  t_sfr_prosp,sfh = get_sfr_psb(res,mod)
  sfh_50, sfh_16, sfh_84 = [], [], []
    for time in range(len(t_sfr_prosp)):
        sfr_quan = quantile([item[time] for item in sfh], [0.16, 0.5, 0.84])
        sfh_50.append(sfr_quan[1])
        sfh_16.append(sfr_quan[0])
        sfh_84.append(sfr_quan[2])
  ax.plot(t_sfr_prosp[::-1], sfh_50) #this is the median inferred SFH
  ax.fill_between(t_sfr_prosp[::-1], y1=sfh_16, y2=sfh_84, alpha=0.3, zorder=0) #this covers the 1sigma width

  

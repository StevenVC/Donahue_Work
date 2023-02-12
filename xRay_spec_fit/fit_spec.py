import numpy as np
import matplotlib.pyplot as plt
import pandas as pds
import csv
import configparser
import sys
import logging

import astropy.units as au
import sherpa
import sherpa.astro.ui as sh
import ciao_contrib.all

from models import \
    gAbsorber_powerlaw, \
    gAbsorber_pAbsorber_powerlaw, \
    gAbsorber_zAbsorber_emissionSpec_powerlaw, \
    gAbsorber_zAbsorber_powerlaw

from astropy.cosmology import Planck15

def checkQuit(userInputRaw):
    '''
    Parameters
    ----------
    Input : input from user; (str)

    Returns : Quit the script if specified (quit,quit(),exit,exit()). Otherwise return the input.
    -------
    '''
    userInput = userInputRaw.lower().strip( ) #change to lower case, and strip leading and trailing whitespace
    if ((userInput == 'quit') or (userInput == 'quit()') or (userInput == 'exit') or (userInput == 'exit()') or (userInput == 'q')):
        plt.clf()
        plt.close()
        sh.clean()
        sys.exit(0)
    
    return(userInputRaw)


'''
Load Config File Info
'''
conf_file = configparser.ConfigParser()
conf_file.read('config.ini')

data_path = conf_file['PATHS']['data_path']
stat = conf_file['FITTING']['stat']

obslo = float(conf_file['FLUX']['obslo'])
obshi = float(conf_file['FLUX']['obshi'])


'''
Load data info
'''
data_info_path = f'data_info_dir_{stat}.csv'
data_info_csv = np.loadtxt(data_info_path, delimiter=',', dtype=str, skiprows=1)

# unpack data
data_dir_names = data_info_csv[:,0]
data_E_filter_lo = data_info_csv[:,1].astype(float)
data_E_filter_hi = data_info_csv[:,2].astype(float)
data_binning_value = data_info_csv[:,3].astype(int)
data_model = data_info_csv[:,4]
data_nH = data_info_csv[:,5].astype(float)
data_redshift = data_info_csv[:,6].astype(float)
data_kT = data_info_csv[:,7].astype(float)
data_CvrFract = data_info_csv[:,8].astype(float)

n_sources = data_dir_names.shape[0]

data_load_array = np.array(data_dir_names, dtype=object)

# print(data_dir_names)

for i in range(n_sources):
    data_load_array[i] = data_path + '/' + data_dir_names[i] + '/' + data_dir_names[i] + 'pc.pi'

col_labels = np.array([
    'Source',
    'H_0 [km / (Mpc s)]',
    'Ω_m',
    'Ω_Λ',
    'Redshift [z]',
    'D_L [Mpc]',
    'Count Rate (2.0-10.0 keV) [s]',
    'Count Rate +/-',
    'Kcorr',
    'Observed Flux (2.0-10.0 keV) [erg/cm^2/s]',
    '-1 Sigma',
    '+1 Sigma',
    'Luminosity (2.0-10.0 keV) [erg/s]',
    '-1 Sigma',
    '+1 Sigma'
]) 

temp_array = np.zeros((n_sources, len(col_labels)))

results_DF = pds.DataFrame(data = temp_array, columns=col_labels)

'''
Load data into sherpa environment

Notes:
    > The error below is expected as the RMF and AMF grid starts with 0
       https://cxc.cfa.harvard.edu/sherpa/bugs/

        "The minimum ENERG_LO in the RMF 'x' was 0 and has 
        been replaced by 1e-10 warnings.warn(wmsg)"
'''
logging.basicConfig(level=logging.INFO, filename="fit.log", filemode="w")
logger = logging.getLogger("sherpa")
logger.setLevel(logging.WARN)

sh.clean() # clean sherpa env

# set fitting/error statistic to use, in this case chi2Gehrels is used due 
# to the possibility of 0 counts in a bin
sh.set_stat('chi2gehrels')
# set fitting optimization method to neldermead since fitting parameters 
# are suspected to be correlated
sh.set_method('neldermead')

# for_index = np.array(0,n_sources)
for_index = np.array([
    0,
    1
])

for i in for_index:
    data_row = np.copy(col_labels)

    c_d_id = data_dir_names[i] # temp save current data id (c_d_id)
    
    model_type = data_model[i]
    # model_type = 'gp'

    nH = data_nH[i]
    z = data_redshift[i]
    kT = data_kT[i]
    CvrFract = data_CvrFract[i]

    sh.load_pha(c_d_id, data_load_array[i]) # load data

    # regroup data
    sh.ungroup(c_d_id)
    sh.notice_id(c_d_id)
    sh.notice_id(c_d_id, data_E_filter_lo[i], data_E_filter_hi[i])

    dmask = sh.get_data(c_d_id).mask
    sh.group_counts(c_d_id, data_binning_value[i], tabStops=~dmask)

    # setup source & bkg model
    if (model_type == 'gp'):
        gAbsorber_powerlaw(c_d_id, i, nH)
        
    elif (model_type == 'gzp'):
        gAbsorber_zAbsorber_powerlaw(c_d_id, i, nH, z)
        
    elif (model_type == 'gzep'):
        gAbsorber_zAbsorber_emissionSpec_powerlaw(c_d_id, i, nH, z, kT)
    
    elif (model_type == 'gpp'):
        gAbsorber_pAbsorber_powerlaw(c_d_id, i, nH, z, CvrFract)

    else:
        sys.exit('Model type + ' + model_type + ', is undefined!!')

    bkg_model = f'abs{i} * powlaw1d.bkgp{i}'
    sh.set_bkg_model(c_d_id, bkg_model)

    # run fit on current source
    sh.fit(c_d_id)
    fit_res = sh.get_fit_results()

    print('\nFit Results:')
    print(fit_res.format(),'\n')

    sh.covar(c_d_id)
    covar_results = sh.get_covar_results()
    covar_matrix = covar_results.extra_output

    print('\nCovariance (parameter confidence interval estimation) Results:')
    print(covar_results.format(), '\n')



    '''
    Calculate the observed and rest-frame fluxes
    '''
    # calculate the median flux, and errors via random sampling of the
    # model parameters (within their ranges) along a normal distribution.
    # This assumes the parameters are correlated in some way
    sample_flux_res = sh.sample_flux(
        lo = obslo, hi = obshi, id = c_d_id,
        scales = covar_matrix, correlated=True, num=1000
        )

    # this returns a 3 element array containing the median, upper and lower 
    # 1-sigma limits from the calculated flux distribution
    obs_flux_res = sample_flux_res[0]

    # calculate the kcorrection that converts from the observed energy band
    # to the rest-frame energy band, at a specified redshift.
    # Effectively this function is doing the following calculation
    # kcorr = rest_flux / obs_flux
    kcorr = sh.calc_kcorr(z, obslo, obshi, id=c_d_id)

    r_f_flux_res = obs_flux_res * kcorr

    print('\nObs flux : [val, +1 sigma, -1 sigma] [erg/cm^2/s]')
    print(obs_flux_res)
    print('\nRest-frame flux : [val, +1 sigma, -1 sigma] [erg/cm^2/s]')
    print(r_f_flux_res)



    '''
    Calculate the observed and rest-frame luminosities
    '''
    D_l = Planck15.luminosity_distance(z)

    units = au.erg / au.cm / au.cm / au.s

    l_f_obs = 4 * np.pi * D_l * D_l * obs_flux_res * units
    l_f_r_f = 4 * np.pi * D_l * D_l * r_f_flux_res * units

    l_f_obs.to(au.erg/au.s)
    l_f_r_f.to(au.erg/au.s)

    print('\nObs luminosity : [val, +1 sigma, -1 sigma] [erg/s]')
    print(l_f_obs.value)
    print('\nRest-frame luminosity : [val, +1 sigma, -1 sigma] [erg/s]')
    print(l_f_obs.value)
    print()



    '''
    Calculate count rate
    '''
    c_rate = 1
    c_rate_err = 1


    
    '''
    Saving data to output csv
    '''
    results_DF.iat[i, 0] = c_d_id
    results_DF.iat[i, 1] = Planck15.H0.value
    results_DF.iat[i, 2] = Planck15.Om0
    results_DF.iat[i, 3] = Planck15.Ode0
    results_DF.iat[i, 4] = z
    results_DF.iat[i, 5] = D_l.value
    results_DF.iat[i, 6] = c_rate
    results_DF.iat[i, 7] = c_rate_err
    results_DF.iat[i, 8] = kcorr
    results_DF.iat[i, 9] = r_f_flux_res[0]
    results_DF.iat[i, 10] = r_f_flux_res[1]
    results_DF.iat[i, 11] = r_f_flux_res[2]
    results_DF.iat[i, 12] = l_f_r_f[0].value
    results_DF.iat[i, 13] = l_f_r_f[1].value
    results_DF.iat[i, 14] = l_f_r_f[2].value

results_DF.to_csv('fit_results.csv')
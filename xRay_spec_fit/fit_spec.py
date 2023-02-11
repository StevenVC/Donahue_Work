import numpy as np
import matplotlib.pyplot as plt
import csv
import configparser
import sys
import logging

import astropy.units as apu
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

obslo = conf_file['FLUX']['obslo']
obshi = conf_file['FLUX']['obshi']


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
logger.setLevel(logging.ERROR)

sh.clean()
sh.set_stat('chi2xspecvar') # set fitting/error statistic to use

temp = 2

for i in range(temp):
    c_d_id = data_dir_names[i] # temp save current data id (c_d_id)
    
    model_type = data_model[i]

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

    bkg_model =  'abs' + str(i) + ' * powlaw1d.bkgp' + str(i)
    sh.set_bkg_model(c_d_id, bkg_model)

    # run fit on current source
    sh.simulfit(c_d_id)

    print('\n1')
    sh.covar(c_d_id)
    covar_results = sh.get_covar_results()
    covar_matrix = covar_results.extra_output
    # covar_err = covar_results.parmaxes
    print('\n2')

    # print(covar_err)
    # covar_err = np.array(covar_err, dtype=float)
    # print(type(covar_err[0]))

    # calculate the median flux, and errors via random sampling of the
    # model parameters (within their ranges) along a normal distribution.
    # This assumes the parameters are correlated in some way
    sample_flux_res = sh.sample_flux(
        lo = obslo, hi = obshi, id = c_d_id, 
        scales = covar_matrix, correlated=True, num=1000
        )

    # this returns a 3 element array containing the median, upper and lower 
    # 1-sigma limits
    obs_flux_res = sample_flux_res[0]

    # calculate the kcorrection that converts from the observed energy band
    # to the rest-frame energy band, at a specified redshift.
    # Effectively this function is doing the following calculation
    # kcorr = rest_flux / obs_flux
    kcorr = sh.calc_kcorr(z, obslo, obshi, id=c_d_id)
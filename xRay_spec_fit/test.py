import astropy.units as au
import numpy as np
from astropy.cosmology import Planck15
import pandas as pds

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

temp_array = np.zeros((1, len(col_labels)))

results_DF = pds.DataFrame(data = temp_array, columns=col_labels)

a = 1 * au.erg/au.s

b = np.array([
    2,
    5,
    7
])

c = a*b

i = 0

results_DF.iat[i, 1] = Planck15.H0.value
results_DF.iat[i, 2] = Planck15.Om0
results_DF.iat[i, 3] = Planck15.Ode0

results_DF.iat[i, 5] = Planck15.luminosity_distance(0.99).value

results_DF.iat[i, 12] = b[0]
results_DF.iat[i, 13] = b[1]
results_DF.iat[i, 14] = b[2]


print(len(col_labels))
print(results_DF)
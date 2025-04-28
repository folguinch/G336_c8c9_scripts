"""Fix the header of a FERIA model PV map.

The native BUNIT of FERIA is `JY/BEAM`, which is not recognized by astropy
(`Jy/beam`).

This script also returns the peak value which can be used in plots to scale the
flux to that of the observations.
"""
from astropy.io import fits
import numpy as np

from common_paths import RESULTS

# Specify pv map
pvmap = RESULTS / 'G336.01-0.82/c8c9/CH3OH/feria_models_PA-55'
pvmap = pvmap / 'feria_model_updated.fits_PV-PA-55deg-CentRA0.0Dec0.0.fits'

# Open data
img = fits.open(pvmap)[0]
img.header['BUNIT'] = 'Jy/beam'
img.writeto(pvmap, overwrite=True)

#Peak value
peak = np.nanmax(img.data)
print('Peak value: ', peak)

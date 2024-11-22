"""Deprojection of continuum map coordinates.

Following the steps in [Persson et al.
2016](https://ui.adsabs.harvard.edu/abs/2016A%26A...590A..33P/abstract)
"""
from pathlib import Path

from astropy.io import fits
from casatasks import tclean
from casatools import ms
import numpy as np

# Values
final_data = Path('/data/share/dihca2/combined_projects/final_data')
vis = ('/data/share/dihca2/combined_projects/'
       'uvdata/G336.018-0.827.c8c9.selfcal.cvel.cont_avg.ms')
vis = Path(vis)
spws = ['0', '1', '2', '3']
PA = np.radians(125.)
incl = np.radians(65.)

# Open ms

# Work per spw
msfiles = []
for spw in spws:
    # Convert to FITS per spw
    fitsname = vis.with_suffix(f'.spw{spw}.uvfits')
    mstool = ms()
    if not fitsname.exists():
        mstool.open(f'{vis}')
        mstool.tofits(f'{fitsname}', spw=spw)
        mstool.close()

    # Reproject
    data = fits.open(fitsname)
    ruv = np.sqrt(data[0].data['UU']**2 + data[0].data['VV']**2)
    gamma = np.arctan2(data[0].data['VV'], data[0].data['UU']) - PA
    data[0].data['UU'] = ruv * np.cos(gamma) * np.cos(incl)
    data[0].data['VV'] = ruv * np.sin(gamma)
    vis_deproj = fitsname.with_suffix(f'.deproj.uvfits')
    data.writeto(vis_deproj)

    # Convert back to MS
    msfile = vis_deproj.with_suffix(f'.ms')
    mstool.fromfits(msfile=f'{msfile}', fitsfile=f'{vis_deproj}')
    mstool.close()
    msfiles.append(msfile)

# Clean
imagename = final_data / vis.with_suffix('.deproj').name
tclean(vis=[f'{msfile}' for msfile in msfiles],
       imagename=f'{imagename}',
       imsize=10800,
       cell='0.004arcsec',
       specmode='mfs',
       outframe='LSRK',
       gridder='standard',
       deconvolver='hogbom',
       interactive=False,
       weighting='briggs',
       robust=0.5,
       niter=0,
       threshold='0.01mJy')

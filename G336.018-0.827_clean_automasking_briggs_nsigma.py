"""
Antennae stats: 20%ile=1437.7m
                25%ile=1732.4m
                30%ile=2043.9m
                75%ile=6216.9m
                90%ile=8491.2m
"""
from pathlib import Path

from casatasks import tclean, exportfits
from astropy.io import fits
from astropy.stats import mad_std

imsize = 14400
cell = '0.004arcsec'
nsigma = 1.25
vis = ('/data/share/dihca2/combined_projects/goco/uvdata/'
       'G336.018-0.827.c8c9.selfcal.cvel.cont_avg.ms')
tclean_args = {'vis': vis,
               'imsize': imsize,
               'cell': cell,
               'nsigma': 1.25,
               'specmode': 'mfs',
               'outframe': 'LSRK',
               'gridder': 'standard',
               'deconvolver': 'mtmfs',
               'nterms': 2,
               'weighting': 'briggs',
               #'robust': 0,
               'niter': 100000,
               'usemask': 'auto-multithresh',
               'noisethreshold': 5,
               'sidelobethreshold': 2.5,
               'lownoisethreshold': 1.5,
               'minbeamfrac': 0.3,
               'negativethreshold': 0,
               'fastnoise': True,
               'pbcor': True,
               }
imagedir = '/data/share/dihca2/combined_projects/auto_continuum_nsigma/'

for robust in [0.5, 2.0]:
    imagebase = ('G336.018-0.827.c8c9.selfcal.cvel.cont_avg'
                 f".{tclean_args['deconvolver']}"
                 f".{tclean_args['weighting']}"
                 f".robust{robust}")

    # Clean down to threshold
    tclean_args['imagename'] = imagedir + imagebase
    imagename = Path(tclean_args['imagename'] + '.image.tt0')
    if not imagename.exists():
        tclean(robust=robust, **tclean_args)

    # Export FITS
    print('Exporting images')
    fitsimage = Path(tclean_args['imagename'] + '.image.tt0.fits')
    exportfits(imagename=f'{imagename}', fitsimage=f'{fitsimage}',
               overwrite=True)
    exportfits(imagename=tclean_args['imagename'] + '.image.tt0.pbcor',
               fitsimage=tclean_args['imagename'] + '.image.tt0.pbcor.fits',
               overwrite=True)
    exportfits(imagename=tclean_args['imagename'] + '.alpha',
               fitsimage=tclean_args['imagename'] + '.alpha.fits',
               overwrite=True)


"""Continuum cleaning.

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
vis = '../uvdata/G336.018-0.827.c8c9.selfcal.cvel.cont_avg.ms'
tclean_args = {'vis': vis,
               'imsize': imsize,
               'cell': cell,
               'nsigma': 1.25,
               'specmode': 'mfs',
               'outframe': 'LSRK',
               'gridder': 'standard',
               'deconvolver': 'hogbom',
               'weighting': 'briggs',
               'robust': 0.5,
               'niter': 1000000,
               'usemask': 'auto-multithresh',
               'noisethreshold': 5,
               'sidelobethreshold': 2.5,
               'lownoisethreshold': 1.5,
               'minbeamfrac': 0.3,
               'negativethreshold': 0,
               'fastnoise': True,
               'pbcor': True,
               }
imagedir = '../data/continuum/'
imagebase = ('G336.018-0.827.c8c9.selfcal.cvel.cont_avg'
             f".{tclean_args['deconvolver']}"
             f".{tclean_args['weighting']}")
if tclean_args['weighting'] == 'briggs':
    imagebase += f".robust{tclean_args['robust']}"

# Clean down to threshold
tclean_args['imagename'] = imagedir + imagebase
imagename = Path(tclean_args['imagename'] + '.image')
if not imagename.exists():
    imagename.parent.mkdir(parents=True, exist_ok=True)
    tclean(**tclean_args)

# Export FITS
print('Exporting images')
to_export = ['.image', '.residual', '.image.pbcor']
for suffix in to_export:
    print('Exporting: ', suffix)
    fitsimage = imagename.with_suffix(f'{suffix}.fits')
    exportfits(imagename=f'{imagename.with_suffix(suffix)}', fitsimage=f'{fitsimage}',
           overwrite=True)
    
    img = fits.open(fitsimage)[0].data
    rms = mad_std(img, ignore_nan=True) * 1E3
    print(f'{suffix[1:].title()} rms: ', rms, 'mJy/beam')

# Old version
# Dirty
#tclean_args['imagename'] = imagedir + imagebase + '.dirty'
#tclean(**tclean_args)
#exportfits(imagename=tclean_args['imagename'] + '.image',
#           fitsimage=tclean_args['imagename'] + '.image.fits')
#
## Get rms
#img = fits.open(tclean_args['imagename'] + '.image.fits')[0].data
#rms = mad_std(img, ignore_nan=True) * 1E3
#print(f'Initial rms: {rms} mJy')
#threshold = 2 * nsigma * rms
#
## Clean down to threshold
#tclean_args['imagename'] = imagedir + imagebase
#tclean_args['niter'] = 1000000
#tclean_args['threshold'] = f'{threshold}mJy'
#tclean(**tclean_args)
#exportfits(imagename=tclean_args['imagename'] + '.image',
#           fitsimage=tclean_args['imagename'] + '.run1.image.fits')
#
## Improve rms
#img = fits.open(tclean_args['imagename'] + '.run1.image.fits')[0].data
#rms = mad_std(img, ignore_nan=True) * 1E3
#print(f'First run rms: {rms} mJy')
#threshold = nsigma * rms
#
## Final run
#print(f'Final threshold: {threshold} mJy')
#tclean_args['threshold'] = f'{threshold}mJy'
#tclean_args['pbcor'] = True
#tclean_args['calcpsf'] = False
#tclean_args['calcres'] = False
#tclean(**tclean_args)
#exportfits(imagename=tclean_args['imagename'] + '.image',
#           fitsimage=tclean_args['imagename'] + '.final.image.fits')

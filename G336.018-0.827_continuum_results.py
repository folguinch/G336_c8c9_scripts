"""Measure fluxes from continuum."""
from pathlib import Path

from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.stats import mad_std
from astropy.wcs import WCS
from regions import CircleSkyRegion
from casatasks import imstat, imfit
import astropy.units as u
import numpy as np

from common_paths import REGIONS, CONTINUUM, RESULTS

# Directories
imagename = ('G336.018-0.827.c8c9.selfcal.cvel.cont_avg.hogbom.'
             'briggs.robust0.5.image.pbcor')
imagename = CONTINUUM / imagename

# Peak position
img_nonpbcor = fits.open(imagename.with_suffix('.fits'))[0]
img = fits.open(imagename.with_suffix('.pbcor.fits'))[0]
wcs = WCS(img, naxis=2)
ind = np.nanargmax(np.squeeze(img.data))
ypeak, xpeak = np.unravel_index(ind, np.squeeze(img.data).shape)
position = SkyCoord.from_pixel(xpeak, ypeak, wcs)
print('Peak coordinate: ', position.to_string(style='hmsdms'))
print('Peak coordinate: ', position.to_string(style='decimal', precision=10))

# Source region for Gaussian fit
region = REGIONS / 'central_source.crtf'
if not region.exists():
    #center = SkyCoord('16:35:09.2585735600', '-48:46:47.6579993356',
    #                  unit=(u.hourangle, u.deg), frame='icrs')
    radius = Angle(0.05, u.arcsec)
    region = CircleSkyRegion(position, radius)
    region.write(REGIONS / 'central_source.crtf', format='crtf',
                 overwrite=True)
stats = imstat(imagename=f'{imagename}', region=f'{region}')
print('Region: ', region)
print('Flux in region: ', stats['flux'][0]*1e3, 'mJy')
rms = mad_std(img.data, ignore_nan=True)
rms_nonpbcor = mad_std(img_nonpbcor.data, ignore_nan=True)
print('Image pbcor rms (uJy): ', rms * 1e6)
print('Image non-pbcor rms (uJy): ', rms_nonpbcor * 1e6)

# Gaussian fit
results = RESULTS / 'G336.01-0.82/c8c9/continuum_results'
results.mkdir(parents=True, exist_ok=True)
kwargs = {
    'residual': f'{results}/gaussian_fit.residual',
    'model': f'{results}/gaussian_fit.model',
    'logfile': f'{results}/gaussian_fit.log',
    # Initial estimates from CARTA fit for the same region
    'estimates': f'{REGIONS}/gaussian_fit_estimates.dat',
    'summary': f'{results}/gaussian_fit_summary.txt',
    'overwrite': True,
}
fit = imfit(imagename=f'{imagename}', region=f'{region}', **kwargs)

# Measure streamer fluxes
for region in REGIONS.glob('*stream_region.crtf'):
    stats = imstat(imagename=f'{imagename}', region=f'{region}')
    print('Region: ', region)
    print('Flux in region: ', stats['flux'][0]*1e3, 'mJy')

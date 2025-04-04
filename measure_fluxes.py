"""Measure fluxes from continuum."""
from pathlib import Path

from casatasks import imstat

from common_paths import REGIONS, CONTINUUM

imagename = ('G336.018-0.827.c8c9.selfcal.cvel.cont_avg.hogbom.'
             'briggs.robust0.5.image.pbcor')
imagename = CONTINUUM / imagename

for region in REGIONS.glob('*stream_region.crtf'):
    stats = imstat(imagename=f'{imagename}', region=f'{region}')
    print('Region: ', region)
    print('Flux in region: ', stats['flux'][0]*1e3, 'mJy')

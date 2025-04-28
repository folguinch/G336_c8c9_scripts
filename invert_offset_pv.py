"""Invert the Offset axis of the pv map published in Olguin et al. 2023 to
match the change in slit P.A. (from 125deg to -55deg)."""
from astropy.io import fits

# Files
pvmap = ('../data/olguin_2023/'
         'G336.01-0.82_rotation.CH3OH_spw0.'
         'ra248.78859_dec-48.77991.PA125.fits')
outpvmap = ('../data/olguin_2023/'
            'G336.01-0.82_rotation.CH3OH_spw0.'
            'ra248.78859_dec-48.77991.PA-55.fits')

# Open images
img = fits.open(pvmap)[0]
img.data = img.data[:,::-1]
img.writeto(outpvmap)

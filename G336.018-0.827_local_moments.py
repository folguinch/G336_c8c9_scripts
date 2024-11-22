from pathlib import Path

from astropy.io import fits
from spectral_cube import SpectralCube
import astropy.units as u
import numpy as np

# Some definitions
datadir = Path('../results/G336.01-0.82/c8c9/CH3OH/')
#datadir = results / 'binary_project/G336.01-0.82/moments/c8c9/CH3OH'
cubename = datadir / ('CH3OH_18_3_15_-17_4_14_A_vt_0_b6_c8c9_spw0_520_590'
                      '_robust2_multiscale_width40_nsig3.subcube.fits')
#momentname = datadir / ('CH3OH_18_3_15_-17_4_14_A_vt_0_b6_c8c9_spw0_520_590'
#                        '_width40_nsig3.subcube.moment0.local.fits')
maskname = datadir / ('CH3OH_18_3_15_-17_4_14_A_vt_0_b6_c8c9_spw0_520_590'
                      '_robust2_multiscale_width40_nsig3.subcube.moment0.'
                      'local.mask.fits')
line_freq = 233.795666 * u.GHz
rms = 1.3 * u.mJy/u.beam
nsigma = 3
vlsr = -47.2 #* u.km/u.s
vel_range = 5 #km/s
momentname = datadir / ('CH3OH_18_3_15_-17_4_14_A_vt_0_b6_c8c9_spw0_520_590'
                        f'_robust2_multiscale_width{vel_range*2}_nsig3.'
                        'subcube.moment0.local.fits')

# Open cube
cube = SpectralCube.read(cubename)
cube = cube.with_spectral_unit(u.km/u.s,
                               velocity_convention='radio',
                               rest_value=line_freq)

# Mask data
data = np.ma.masked_invalid(cube.filled_data[:].value)

# Peak intensity in velocity
vel = cube.spectral_axis.to(u.km/u.s).value
index_max = np.nanargmax(data, axis=0)
peak_vel_map = vel[index_max]
#peak_vel_map = np.ma.masked_where(np.all(data_ma.mask, axis=0), vel_map)

# Velocity cube
vel_cube = np.zeros(cube.shape).transpose()
vel_cube = vel_cube + vel
vel_cube = vel_cube.transpose()

# Shift the velocity and create a mask
vel_cube = vel_cube - peak_vel_map
range_mask = (vel_cube >= -vel_range) & (vel_cube <= vel_range)
range_mask[data.mask] = False
if maskname is not None:
    hdu = fits.PrimaryHDU(range_mask.astype(int), header=cube.header)
    hdu.writeto(maskname)

# Mask cube and create moment 0 map
if rms is not None:
    cube = cube.with_mask(cube > nsigma*rms)
cube = cube.with_mask(range_mask)
moment = cube.moment(order=0)
moment.write(momentname)


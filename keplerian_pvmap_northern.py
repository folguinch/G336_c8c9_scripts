"""Generate projected Keplerian veolocities for PV."""
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.io import fits
from line_little_helper.pvmap_extractor import get_spectral_slab
from matplotlib import font_manager
from matplotlib.ticker import AutoMinorLocator
from pvextractor import extract_pv_slice
from pvextractor import Path as PVPath
from regions import Regions
from spectral_cube import SpectralCube
from tile_plotter.utils import get_extent
import astropy.units as u
import astropy.constants as ct
import matplotlib.pyplot as plt
import numpy as np

from common_paths import DATA, CONFIGS, RESULTS, FIGURES

# Import fonts
font_dirs = [Path('./fonts'), Path('/usr/share/fonts/OTF')]
font_dirs = filter(lambda x: x.exists(), font_dirs)
print('Font directories:', font_dirs)
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)

for font_file in font_files:
    font_manager.fontManager.addfont(font_file)
fonts = font_manager.get_font_names()

# Plot styles
styles = ['science_double_col', 'inferno']
if 'Myriad Pro' in fonts:
    print('Using Myriad style')
    styles.append('science_myriad')

# Directories
cube = DATA / 'cubes/G336.018-0.827_spw0_CH3OH_robust2_multiscale.image.fits'
pvmap_file = (RESULTS / 'G336.01-0.82/c8c9/pvmaps' /
              'G336_north_streamer_inverted_CH3OH_spw0.fits')
region = CONFIGS / 'pvmaps/regions/G336_north_streamer_inverted.crtf'
figures = FIGURES / 'G336.01-0.82/c8c9/paper'

# Source parameters
source = SkyCoord('16h35m09.2587s -48d46m47.657s', frame='icrs')
mass = 10 * u.M_sun
incl = 65 * u.deg
pa = 125 * u.deg
rcb = 250 * u.au
distance = 3.1 * u.kpc
vlsr = -47.2 * u.km / u.s
rot = -1
line_freq = 233.795666 * u.GHz

# Cube
cube = SpectralCube.read(cube)
wcs = cube.wcs.sub(2)
spacing = 1
width = 0.04 * u.arcsec
freq_slab = 28 * u.MHz

# Slab
if not pvmap_file.exists():
    cube = cube.with_spectral_unit(u.km/u.s,
                                   velocity_convention='radio',
                                   rest_value=line_freq)
    cube = get_spectral_slab(cube, line_freq, freq_slab, vlsr=vlsr)

# Read region
reg = Regions.read(region, format='crtf').pop()

# Path
pv_path = PVPath(reg.vertices, width=width)
x, y = pv_path.sample_points(spacing, wcs=wcs)
positions = SkyCoord.from_pixel(x, y, wcs=wcs)
proj_d = source.separation(positions).to(u.arcsec)
pa_point = source.position_angle(positions).to(u.deg)

# De-project distance
# Verified using https://ui.adsabs.harvard.edu/abs/2012A%26A...547A..84T/abstract
# Note that definition of the azimuthal angle phi is inverted with respect to
# that paper so the velocity directions are inverted here
phi = 180 * u.deg + pa - pa_point
deproj_d = proj_d * np.sqrt((np.sin(phi)/np.cos(incl))**2 + 
                            np.cos(phi)**2)
deproj_d = deproj_d.to(u.arcsec).value * distance.to(u.pc).value * u.au
proj_d = proj_d.to(u.arcsec).value * distance.to(u.pc).value * u.au

# Keplerian velocity
vel_kp = rot * np.sqrt(ct.G * mass / deproj_d).to(u.km / u.s)
proj_vel_kp = vel_kp * np.sin(incl) * np.cos(phi)

# Infall + Keplerian
vel_inf = np.sqrt(2 * ct.G * mass / deproj_d).to(u.km / u.s)
proj_vel_kp_inf = (vel_kp * np.sin(incl) * np.cos(phi) +
                   vel_inf * np.sin(incl) * np.sin(phi))

# IRE
vel_rot = rot * np.sqrt(2 * ct.G * mass * rcb) / deproj_d
vel_inf = np.sqrt(2 * ct.G * mass * (deproj_d - rcb)) / deproj_d
proj_vel_ire = (vel_rot * np.sin(incl) * np.cos(phi) +
                vel_inf * np.sin(incl) * np.sin(phi))
proj_vel_ire[deproj_d < rcb] = proj_vel_kp[deproj_d < rcb]

# PV map
if pvmap_file.exists():
    print('Opening pv file')
    pv_map = fits.open(pvmap_file)[0]
else:
    print('Calculating pv from cube')
    pv_map = extract_pv_slice(cube, pv_path)
extent = get_extent(pv_map)
extent = [extent[0].to(u.arcsec).value * distance.to(u.pc).value,
          extent[1].to(u.arcsec).value * distance.to(u.pc).value,
          extent[2].to(u.km/u.s).value - vlsr.value,
          extent[3].to(u.km/u.s).value - vlsr.value]
dist = np.linspace(extent[0], extent[1], len(proj_vel_kp))

# Second axis functions
def forward(x):
    return np.interp(x, dist, deproj_d.value)

def inverse(x):
    return np.interp(x, deproj_d.value, dist)

# Plot
fontsize = 8
ratio = 11 / 5
width = 5.67
height = 2.9
plt.style.use(styles)
fig = plt.figure(figsize=(width, height))
ax = fig.add_axes((0.08, 0.12, 0.91, 0.77))
ax.imshow(pv_map.data, extent=extent, origin='lower')
ax.set_aspect('auto')
ax.set_xlabel('Distance along streamer path (au)')
ax.set_ylabel('Velocity (km/s)')
ax.set_ylim(10, -15)
ax.set_xlim(0, 1500)
ax.axvline(inverse(500), linewidth=1.5, color='#1ae6ed')
ax.annotate('$R_c$', (0.35, 0.9), color='#1ae6ed', xytext=(0.35, 0.9),
            xycoords='axes fraction', size=fontsize)
ax.annotate('Inner', (0.27, 0.1), color='#1ae6ed', xytext=(0.27, 0.1),
            xycoords='axes fraction', size=fontsize)
ax.annotate('Outer', (0.35, 0.1), color='#1ae6ed', xytext=(0.35, 0.1),
            xycoords='axes fraction', size=fontsize)
line1, = ax.plot(dist, proj_vel_kp.to(u.km/u.s).value, ls='-', color='#6adb02',
                 lw=1.5)
line2, = ax.plot(dist, proj_vel_kp_inf.to(u.km/u.s).value, ls=':',
                 color='#6adb02', lw=1.5)
line3, = ax.plot(dist, proj_vel_ire.to(u.km/u.s).value, ls='--',
                 color='#6adb02', lw=1.5)
ax.legend([line1, line2, line3],
          ['Keplerian rotation', 'Keplerian + free-fall', 'IRE'],
          fontsize=fontsize,
          labelcolor='#6adb02')
ax.tick_params(length=4, color='w', top=False)
secax = ax.secondary_xaxis('top', functions=(forward, inverse))
secax.xaxis.set_minor_locator(AutoMinorLocator())
secax.set_xlabel('Deprojected radial distance to source (au)')
secax.tick_params(length=4, color='w')
secax.tick_params(length=2.5, color='w', which='minor')
fig.savefig(figures / 'keplerian_pv_northern.png')
fig.savefig(figures / 'keplerian_pv_northern.pdf')

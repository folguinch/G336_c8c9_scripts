"""Plot projected velocities over red shifted streamer PV."""
from pathlib import Path

from astropy.table import QTable
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
              'G336_south_streamer_inverted_CH3OH_spw0_red.fits')
region = CONFIGS / 'pvmaps/regions/G336_south_streamer_inverted.crtf'
figures = FIGURES / 'G336.01-0.82/c8c9/papers'
results = RESULTS / 'G336.01-0.82/c8c9/CH3OH'
outer_streamer = (results /
                  'south_best/south_best_incl65_pa125/regions/'
                  'stream_south_rmin500_r02500_rc500_theta080_phi0280_vr02.ecsv')

# PV parameters
distance = 3.1 * u.kpc
vlsr = -47.2 * u.km / u.s
line_freq = 233.795666 * u.GHz
width = 0.04 * u.arcsec
freq_slab = 28 * u.MHz
spacing = 1

# For model
source = SkyCoord('16h35m09.2585s -48d46m47.661s', frame='icrs')
incl = 65 * u.deg
pa = -55 * u.deg

# Path and pv
cube = SpectralCube.read(cube)
wcs = cube.wcs.sub(2)
reg = Regions.read(region, format='crtf').pop()
pv_path = PVPath(reg.vertices, width=width)
if not pvmap_file.exists():
    print('Calculating pv from cube')
    cube = cube.with_spectral_unit(u.km/u.s,
                                   velocity_convention='radio',
                                   rest_value=line_freq)
    cube = get_spectral_slab(cube, line_freq, freq_slab, vlsr=vlsr)

    # PV map
    pv_map = extract_pv_slice(cube, pv_path)
else:
    print('Opening pv file')
    pv_map = fits.open(pvmap_file)[0]

# Deprojected distance
x, y = pv_path.sample_points(spacing, wcs=wcs)
positions = SkyCoord.from_pixel(x, y, wcs=wcs)
proj_d = source.separation(positions).to(u.arcsec)
pa_point = source.position_angle(positions).to(u.deg)
phi = pa_point - pa
deproj_d = proj_d * np.sqrt((np.sin(phi)/np.cos(incl))**2 + 
                            np.cos(phi)**2)
deproj_d = deproj_d.to(u.arcsec).value * distance.to(u.pc).value * u.au
#proj_d = proj_d.to(u.arcsec).value * distance.to(u.pc).value * u.au

extent = get_extent(pv_map)
extent = [extent[0].to(u.arcsec).value * distance.to(u.pc).value,
          extent[1].to(u.arcsec).value * distance.to(u.pc).value,
          extent[2].to(u.km/u.s).value - vlsr.value,
          extent[3].to(u.km/u.s).value - vlsr.value]
dist = np.linspace(extent[0], extent[1], len(deproj_d))

# Streamers
outer_streamer = QTable.read(outer_streamer, format='ascii.ecsv')
outer_streamer_positions = SkyCoord(outer_streamer['ra'],
                                    outer_streamer['dec'],
                                    frame='icrs')
outer_streamer_proj_d = source.separation(outer_streamer_positions).to(u.arcsec)
outer_streamer_pa_point = source.position_angle(outer_streamer_positions).to(u.deg)
outer_streamer_phi = outer_streamer_pa_point - pa
outer_streamer_deproj_d = (outer_streamer_proj_d *
                           np.sqrt((np.sin(outer_streamer_phi)/np.cos(incl))**2 + 
                                   np.cos(outer_streamer_phi)**2))
outer_streamer_deproj_d = (outer_streamer_deproj_d.to(u.arcsec).value *
                           distance.to(u.pc).value * u.au)

# Second axis functions
def forward(x):
    return np.interp(x, dist, deproj_d.value)
def inverse(x):
    return np.interp(x, deproj_d.value, dist)

# Plot
fontsize = 8
width = 5.67
height = 2.9
plt.style.use(styles)
fig = plt.figure(figsize=(width, height))
ax = fig.add_axes((0.07, 0.12, 0.92, 0.84))
ax.imshow(pv_map.data, extent=extent, origin='lower')
ax.set_aspect('auto')
ax.set_xlabel('Distance along streamer path (au)')
ax.set_ylabel('Velocity (km/s)')
ax.set_ylim(-2, 10)
ax.set_xlim(0, 1750)
line1, = ax.plot(inverse(outer_streamer_deproj_d.value),
                 outer_streamer['vlsr'].to(u.km/u.s).value,
                 ls='-', color='#baa800', lw=1.5)
ax.legend([line1], ['Streamline model'], fontsize=fontsize, labelcolor='#baa800')
ax.tick_params(length=4, color='w', top=True)
# Doesn't work well because the slice crosses the centrifugal radius twice
#ax.tick_params(length=4, color='w', top=False)
#secax = ax.secondary_xaxis('top', functions=(forward, inverse))
#secax.xaxis.set_minor_locator(AutoMinorLocator())
#secax.set_xlabel('Deprojected radial distance to source (au)')
#secax.tick_params(length=4, color='w')
#secax.tick_params(length=2.5, color='w', which='minor')
fig.savefig(figures / 'southern_streamer_pv.png')
fig.savefig(figures / 'southern_streamer_pv.pdf')

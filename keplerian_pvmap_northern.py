"""Plot projected velocities over blue shifted streamer PV."""
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
              'G336_north_streamer_inverted_CH3OH_spw0_blue.fits')
region = CONFIGS / 'pvmaps/regions/G336_north_streamer_inverted.crtf'
figures = FIGURES / 'G336.01-0.82/c8c9/papers'
results = RESULTS / 'G336.01-0.82/c8c9/CH3OH'
outer_streamer = (results /
                  'north_outer_best/north_outer_best_incl65_pa125/regions/'
                  'stream_north_rmin500_r02500_rc500_theta080_phi055_vr00.ecsv')
inner_streamer = (results /
                  'north_inner_best/north_inner_best_incl65_pa125/regions/'
                  'stream_north_rmin200_r0500_rc200_theta089_phi0145_vr00.1.ecsv')

# Source parameters
source = SkyCoord('16h35m09.2585s -48d46m47.661s', frame='icrs')
mass = 10 * u.M_sun
incl = 65 * u.deg
#pa = 125 * u.deg
pa = -55 * u.deg
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
# Velocity deprojection verified using Tang et al. (2012)
# https://ui.adsabs.harvard.edu/abs/2012A%26A...547A..84T/abstract
# Note that definition of the azimuthal angle phi is inverted with respect to
# that paper so the velocity directions are inverted here.
# Uncomment/comment the values below to change between the conventions.
# Published results use the Tang et al. convention as it is the same as FERIA.
#phi = 180 * u.deg + pa - pa_point
phi = pa_point - pa
deproj_d = proj_d * np.sqrt((np.sin(phi)/np.cos(incl))**2 + 
                            np.cos(phi)**2)
deproj_d = deproj_d.to(u.arcsec).value * distance.to(u.pc).value * u.au
proj_d = proj_d.to(u.arcsec).value * distance.to(u.pc).value * u.au

# Keplerian velocity
vel_kp = rot * np.sqrt(ct.G * mass / deproj_d).to(u.km / u.s)
#proj_vel_kp = vel_kp * np.sin(incl) * np.cos(phi)
# Tang et al. (2012) check
proj_vel_kp = vel_kp * np.sin(incl) * np.cos(phi)

# Infall + Keplerian
#vel_inf = np.sqrt(2 * ct.G * mass / deproj_d).to(u.km / u.s)
#proj_vel_kp_inf = (vel_kp * np.sin(incl) * np.cos(phi) +
#                   vel_inf * np.sin(incl) * np.sin(phi))
# Tang et al. (2012) check
vel_inf = -np.sqrt(2 * ct.G * mass / deproj_d).to(u.km / u.s)
proj_vel_kp_inf = (vel_kp * np.sin(incl) * np.cos(phi) +
                   vel_inf * np.sin(incl) * np.sin(phi))

# IRE
vel_rot = rot * np.sqrt(2 * ct.G * mass * rcb) / deproj_d
#vel_inf = np.sqrt(2 * ct.G * mass * (deproj_d - rcb)) / deproj_d
#proj_vel_ire = (vel_rot * np.sin(incl) * np.cos(phi) +
#                vel_inf * np.sin(incl) * np.sin(phi))
# Tang et al. (2012) check
vel_inf = -np.sqrt(2 * ct.G * mass * (deproj_d - rcb)) / deproj_d
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

# Streamers
outer_streamer = QTable.read(outer_streamer, format='ascii.ecsv')
outer_streamer_positions = SkyCoord(outer_streamer['ra'],
                                    outer_streamer['dec'],
                                    frame='icrs')
outer_streamer_proj_d = source.separation(outer_streamer_positions).to(u.arcsec)
outer_streamer_pa_point = source.position_angle(outer_streamer_positions).to(u.deg)
#outer_streamer_phi = 180 * u.deg + pa - outer_streamer_pa_point
outer_streamer_phi = outer_streamer_pa_point - pa
outer_streamer_deproj_d = (outer_streamer_proj_d *
                           np.sqrt((np.sin(outer_streamer_phi)/np.cos(incl))**2 + 
                                   np.cos(outer_streamer_phi)**2))
outer_streamer_deproj_d = (outer_streamer_deproj_d.to(u.arcsec).value *
                           distance.to(u.pc).value * u.au)
inner_streamer = QTable.read(inner_streamer, format='ascii.ecsv')
inner_streamer_positions = SkyCoord(inner_streamer['ra'],
                                    inner_streamer['dec'],
                                    frame='icrs')
inner_streamer_proj_d = source.separation(inner_streamer_positions).to(u.arcsec)
inner_streamer_pa_point = source.position_angle(inner_streamer_positions).to(u.deg)
#inner_streamer_phi = 180 * u.deg + pa - inner_streamer_pa_point
inner_streamer_phi = inner_streamer_pa_point - pa
inner_streamer_deproj_d = (inner_streamer_proj_d *
                           np.sqrt((np.sin(inner_streamer_phi)/np.cos(incl))**2 + 
                                   np.cos(inner_streamer_phi)**2))
inner_streamer_deproj_d = (inner_streamer_deproj_d.to(u.arcsec).value *
                           distance.to(u.pc).value * u.au)

# Plot
fontsize = 8
width = 5.67
height = 2.9
plt.style.use(styles)
fig = plt.figure(figsize=(width, height))
ax = fig.add_axes((0.08, 0.12, 0.91, 0.77))
ax.imshow(pv_map.data, extent=extent, origin='lower')
ax.set_aspect('auto')
ax.set_xlabel('Distance along streamer path (au)')
ax.set_ylabel('Velocity (km/s)')
ax.set_ylim(4, -12)
ax.set_xlim(0, 1500)
ax.axvline(inverse(500), linewidth=1.5, color='#66ccff')
ax.annotate('$R_c$', (0.35, 0.9), color='#66ccff', xytext=(0.35, 0.9),
            xycoords='axes fraction', size=fontsize)
ax.annotate('Inner', (0.27, 0.1), color='#66ccff', xytext=(0.27, 0.1),
            xycoords='axes fraction', size=fontsize)
ax.annotate('Outer', (0.35, 0.1), color='#66ccff', xytext=(0.35, 0.1),
            xycoords='axes fraction', size=fontsize)
line1, = ax.plot(dist, proj_vel_kp.to(u.km/u.s).value, ls='-', color='#009900',
                 lw=1.5)
line2, = ax.plot(dist, proj_vel_kp_inf.to(u.km/u.s).value, ls=':',
                 color='#009900', lw=1.5)
line3, = ax.plot(dist, proj_vel_ire.to(u.km/u.s).value, ls='--',
                 color='#009900', lw=1.5)
line4, = ax.plot(inverse(outer_streamer_deproj_d.value),
                 outer_streamer['vlsr'].to(u.km/u.s).value,
                 ls=(0, (3, 1, 1, 1)), color='#3366ff', lw=1.5)
line5, = ax.plot(inverse(inner_streamer_deproj_d.value),
                 inner_streamer['vlsr'].to(u.km/u.s).value,
                 ls=(0, (3, 1, 1, 1, 1, 1)), color='#3366ff', lw=1.5)
ax.legend([line1, line2, line3, line4, line5],
          ['Keplerian rotation', 'Keplerian + free-fall', 'IRE',
           'Outer streamline', 'Inner streamline'],
          fontsize=fontsize,
          labelcolor=['#009900', '#009900', '#009900', '#3366ff', '#3366ff'])
ax.tick_params(length=4, color='w', top=False)
secax = ax.secondary_xaxis('top', functions=(forward, inverse))
secax.xaxis.set_minor_locator(AutoMinorLocator())
secax.set_xlabel('Deprojected radial distance to source (au)')
secax.tick_params(length=4, color='w')
secax.tick_params(length=2.5, color='w', which='minor')
fig.savefig(figures / 'keplerian_pv_northern.png')
fig.savefig(figures / 'keplerian_pv_northern.pdf')

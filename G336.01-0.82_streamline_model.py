"""Calculate streamer models."""
from typing import Dict, Optional, Tuple
from itertools import product
from pathlib import Path
from dataclasses import dataclass, field

from astropy.coordinates import SkyCoord, ICRS
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import rc
from matplotlib.collections import LineCollection
from regions import Regions, PolygonSkyRegion
from scipy import stats
from tile_plotter.multi_plotter import OTFMultiPlotter
#import aplpy
import astropy.constants as ct
import astropy.units as u
import numpy as np
import numpy.typing as npt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import pickle
import velocity_tools.coordinate_offsets as c_offset
import velocity_tools.stream_lines as SL

@dataclass
class FitPars:
    """Keep track of parameters used for fitting."""
    rmin: npt.ArrayLike
    r0: npt.ArrayLike
    rc: npt.ArrayLike
    theta0: npt.ArrayLike
    phi0: npt.ArrayLike
    v_r0: npt.ArrayLike

    def __iter__(self):
        for vals in product(*self._as_tuple()):
            name = FitPars.generate_name(vals)
            yield name, vals

    def _as_tuple(self):
        return (self.rmin, self.r0, self.rc, self.theta0, self.phi0, self.v_r0)

    @staticmethod
    def generate_name(vals: Tuple, sep: str = '_'):
        names = ('rmin', 'r0', 'rc', 'theta0', 'phi0', 'vr0')
        fname = []
        for name, val in zip(names, vals):
            fname.append(f'{name}{val.value:n}')

        return sep.join(fname)

@dataclass
class ModelPars:
    """Keep track of model parameters."""
    position: SkyCoord
    distance: u.Quantity
    v_lsr: u.Quantity
    Mstar: u.Quantity
    # Source angles
    inc: u.Quantity
    PA_ang: u.Quantity
    components: Dict[str, FitPars] = field(default_factory=dict)
    ranges: Dict[str, tuple] = field(default_factory=dict)

    def __iter__(self):
        for item in self.components.items():
            yield item

@dataclass
class ObsData:
    """Store observed data."""
    position: SkyCoord
    moment1: Path
    continuum: Path
    components: Dict[str, Regions] = field(default_factory=dict)
    spines: Dict[str, Regions] = field(default_factory=dict)
    mom1_map: Optional[fits.PrimaryHDU] = None
    cont_map: Optional[fits.PrimaryHDU] = None
    mom1_wcs: Optional[WCS] = None
    cont_wcs: Optional[WCS] = None

    def __post_init__(self):
        if self.mom1_map is None:
            self.mom1_map = fits.open(self.moment1)[0]
            self.mom1_wcs = WCS(self.mom1_map, naxis=2)
        if self.cont_map is None:
            self.cont_map = fits.open(self.continuum)[0]
            self.cont_wcs = WCS(self.cont_map, naxis=2)

def get_vc_r(image_file, region_file, position, distance):
    """
    Returns the centroid velocity and projected separation in the sky.

    `r_proj` is in `u.au` and `V_los` is in `u.km/u.s`
    """
    # load region file and WCS structures
    regions = Regions.read(region_file, format='crtf')
    wcs_Vc = WCS(f'{image_file}')
    #
    hd_Vc = fits.getheader(image_file)
    results = c_offset.generate_offsets(hd_Vc, position.ra,
                                        position.dec, pa_angle=0*u.deg,
                                        inclination=0*u.deg)
    rad_au = results.r * distance.to(u.pc)
    rad_au = rad_au.to(u.au, equivalencies=u.dimensionless_angles())
    # Vc_all =
    #
    mask_Vc = (regions[0].to_pixel(wcs_Vc)).to_mask()
    Vc_cutout = mask_Vc.cutout(fits.getdata(image_file))
    rad_cutout = mask_Vc.cutout(rad_au)
    #
    gd = mask_Vc.data == 1
    v_los = Vc_cutout[gd] * u.km/u.s
    r_proj = rad_cutout[gd]

    return r_proj, v_los

#def setup_plot(fig, center, radius=1.*u.arcsec, label_col='black',
#               star_col='red', distance=None, marker=None):
#    """
#    Setup of plots, since they will show all the same format.
#    """
#    fig.set_system_latex(False)
#    fig.ticks.set_color(label_col)
#    fig.recenter(*center, radius=radius.to(u.deg).value)
#    fig.set_nan_color('0.9')
#    fig.add_beam(color=label_col)
#    if distance is not None:
#        ang_size = (500 / distance.to(u.pc).value) * u.arcsec
#        fig.add_scalebar(ang_size, label='500 au', color=label_col)
#    if marker is not None:
#        fig.show_markers(*marker, marker='*', s=60, layer='star',
#                         edgecolor=star_col, facecolor=label_col, zorder=31)
#    fig.tick_labels.set_xformat('hh:mm:ss.ss')
#    fig.tick_labels.set_yformat('dd:mm:ss.s')
#    fig.ticks.set_length(7)
#    fig.axis_labels.set_xtext(r'Right Ascension (J2000)')
#    fig.axis_labels.set_ytext(r'Declination (J2000)')
def plot_stream_map(obsdata, streamer, streamer_vel, filename, v_lsr):
    # Generate plot
    loc = (0, 0)
    map_plot = OTFMultiPlotter(nrows='1',
                               right='0.8',
                               left='1.3',
                               vertical_cbar='true',
                               styles='maps vik',
                               label_xpad='0.45',
                               label_ypad='-0.7',
                               vcbarpos='0')
    handler = map_plot.gen_handler(
        loc,
        'moment',
        projection=obsdata.mom1_wcs,
        include_cbar=True,
        name='Velocity',
        unit='km/s',
        ticks_color='k',
        yticks_fmt='dd:mm:ss.s',
        xticks_fmt='hh:mm:ss.ss',
    )
    handler.plot_map(obsdata.mom1_map,
                     position=obsdata.position,
                     radius=0.6*u.arcsec,
                     zorder=1,
                     shift_data=-v_lsr,
                     )
    handler.plot_contours(obsdata.cont_map,
                          ignore_units=True,
                          stretch='log',
                          linewidths=1,
                          colors='#7f7f7f',
                          zorder=2,
                          )

    # Plot streamer
    #points = np.array([streamer.ra.deg, streamer.dec.deg]).T.reshape(-1, 1, 2)
    #segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #lc = LineCollection(segments, cmap='vik',
    #                    norm=handler.vscale.get_normalization(),
    #                    transform=handler.get_transform(),
    #                    zorder=4)
    ##lc.set_array(streamer_vel)
    #lc.set(lw=6, array=streamer_vel[:-1].value, alpha=1,
    #       #path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()],
    #       )
    ##print(lc.properties(), lc.get_fc())
    #line = handler.ax.add_collection(lc)
    #handler.plot(fil.ra.deg, fil.dec.deg, transform=handler.get_transform(),
    #             lw=2, ls='-', )
    handler.scatter(streamer.ra, streamer.dec, c=streamer_vel.value,
                    cmap='vik', norm=handler.vscale.get_normalization(),
                    zorder=5)
    handler.plot(streamer.ra, streamer.dec, linestyle='-', linewidth=10,
                 solid_capstyle='round',
                 color='k', zorder=4, transform=handler.get_transform())

    # Configuration
    map_plot.apply_config(loc, handler, 'moment')
    handler.plot_cbar(map_plot.fig,
                      map_plot.axes[loc].cborientation)
    map_plot.savefig(filename)

def convert_into_mili(file_name):
    """
    It converts a file into one rescaled by 1e3.
    This is useful to convert between Jy -> mJy or m/s into km/s
    for plotting purposes (e.g. to use with aplpy).

    Usage:
    fig = aplpy.FITSFigure(convert_into_mili(file_in_Jy), figsize=(4,4))
    fig.show_colorscale(vmin=0, vmax=160, cmap='inferno')
    fig.add_colorbar()
    fig.colorbar.set_axis_label_text(r'Integrated Intensity (mJy beam$^{-1}$ km s$^{-1}$)')

    :param file_name: string with filename to process
    :return: hdu
    """
    data, hd = fits.getdata(file_name, header=True)
    return fits.PrimaryHDU(data=data*1e3, header=hd)

def fit_streamer(data: ObsData, model: ModelPars, outdir: Path):
    # Iterate components
    for name, component in model:
        print(f'Fitting component: {name}')
        # Compute observed KDE
        #r_proj, v_los = get_vc_r(data.moment1, data.components[name],
        #                         data.position, model.distance)
        #xmin, xmax, velmin, velmax = model.ranges[name]
        #xx, vv = np.mgrid[xmin:xmax:100j, velmin:velmax:100j]
        #positions = np.vstack([xx.ravel(), vv.ravel()])
        #
        #gd_vlos = np.isfinite(r_proj * v_los)
        #values = np.vstack([r_proj[gd_vlos].value, v_los[gd_vlos].value])
        #kernel = stats.gaussian_kde(values)
        #kernel = np.reshape(kernel(positions).T, xx.shape)
        #kernel /= kernel.max()

        # Iterate models
        for fname, vals in component:
            # Output region
            reg_name = outdir / 'regions' / f'stream_{name}_{fname}.crtf'
            if reg_name.exists():
                print(f'Skipping {fname}: already calculated')
                continue

            # Parameters
            rmin, r0, rc, theta0, phi0, v_r0 = vals
            omega0 = np.sqrt(rc * ct.G * model.Mstar) / r0**2
            omega0 = omega0.to(1 / u.s)
            if np.abs(rmin - r0) < 100*u.au:
                r_step = 0.1
            else:
                r_step = 10

            # Stream model
            (dra, _, ddec), (_, vel, _) = SL.xyz_stream(
                mass=model.Mstar,
                r0=r0,
                theta0=theta0,
                phi0=phi0,
                omega=omega0,
                v_r0=v_r0,
                inc=model.inc,
                pa=model.PA_ang,
                rmin=rmin,
                step=r_step,
            )
            d_sky_au = np.sqrt(dra**2 + ddec**2)

            # Stream line into arcsec
            dra = -dra.value / model.distance.to(u.pc).value * u.arcsec
            ddec = ddec.value / model.distance.to(u.pc).value * u.arcsec
            if not np.any(np.isfinite(vel)):
                print(f'Skipping {fname}: no finite solutions')
                continue

            # Save as region
            fil = SkyCoord(dra, ddec,
                           frame=model.position.skyoffset_frame())
            fil = fil.transform_to(ICRS)
            fil_table = fil.to_table()
            fil_table.add_column(vel, name='vlsr')
            fil_table.write(reg_name.with_suffix('.ecsv'), overwrite=True)
            fil_reg = PolygonSkyRegion(vertices=fil)
            print(f'Saving region: {reg_name}')
            fil_reg.write(str(reg_name), format='crtf', overwrite=True)

            # Plot map
            fig_name = outdir / f'stream_{name}_{fname}.png'
            plot_stream_map(data, fil, vel, fig_name, model.v_lsr)

if __name__ == '__main__':
    basedir = Path('/data/share/dihca2/combined_projects')
    results = basedir / 'results/G336.01-0.82/c8c9/CH3OH'
    regions = basedir / 'scripts/configs/pvmaps/regions'
    #figures = results / 'streamer_south_models_incl65_cb500_rout2500_v02'
    figures = results / 'refined_inner_streamer' / 'streamer_north_models_incl65_cb475_rout500_v0'
    moment0 = results / 'CH3OH_18_3_15_-17_4_14_A_vt_0_b6_c8c9_spw0_520_590_width40_nsig3.subcube.moment0.fits'
    moment1 = results / 'CH3OH_18_3_15_-17_4_14_A_vt_0_b6_c8c9_spw0_520_590_width40_nsig3.subcube.moment1.fits'
    continuum = (basedir / 'manual_selfcal' /
                 'G336.018-0.827.4.cont.image.fits')

    # Streamers
    region_files = [regions / 'north_streamer_5sigma_poly.crtf',
                    regions / 'south_streamer_5sigma_poly.crtf']

    # Load data

    # Source properties
    position = SkyCoord('16h35m9.26085s -48d46m47.65854s', frame=ICRS)
    distance = 3.1 * u.kpc
    v_lsr = -47.2 * u.km/u.s
    obs_data = ObsData(
        position,
        moment1,
        continuum,
        components={'north': (regions /
                              'G336.01-0.82_north_stream_ch3oh_poly.crtf'),
                    'south': (regions /
                              'G336.01-0.82_south_stream_ch3oh_poly.crtf')},
        spines={'north': regions / 'G336.01-0.82_north_stream_ch3oh.crtf',
                'south': regions / 'G336.01-0.82_south_stream_ch3oh.crtf'}
    )

    # Model parameters
    # Large scale refinement
    #north_pars = FitPars(
    #    np.array([500]) * u.au,
    #    np.array([2500]) * u.au,
    #    np.array([500]) * u.au,
    #    np.array([75, 80, 85]) * u.deg,
    #    np.array([50, 55, 60, 65, 70]) * u.deg,
    #    np.array([0]) * u.km/u.s,
    #)
    # Small scale refinement
    north_pars = FitPars(
        np.array([475]) * u.au,
        np.array([500]) * u.au,
        np.array([475]) * u.au,
        #np.array([85, 86, 87, 88, 89, 90]) * u.deg,
        #np.array([150, 155, 160, 165, 170, 175, 180, 185, 190, 195]) * u.deg,
        # For best
        np.array([89]) * u.deg,
        np.array([155]) * u.deg,
        np.array([0]) * u.km/u.s,
    )
    #south_pars = FitPars(
    #    np.array([500]) * u.au,
    #    #np.array([1000, 1500, 2000, 2500, 3000]) * u.au,
    #    np.array([2000, 2500, 3000, 3500, 4000]) * u.au,
    #    np.array([500]) * u.au,
    #    #np.array([40, 50, 60, 70]) * u.deg,
    #    np.array([50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]) * u.deg,
    #    #np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 110, 120]) * u.deg,
    #    #np.array([255, 260, 265, 270, 275]) * u.deg,
    #    #np.array([265, 270, 275, 280]) * u.deg,
    #    np.array([270, 275, 280, 285, 290, 295, 300]) * u.deg,
    #    #np.array([0, 2]) * u.km/u.s,
    #    np.array([2]) * u.km/u.s,
    #)
    model_pars = ModelPars(
        position,
        distance,
        v_lsr,
        10 * u.Msun,
        #5 * u.Msun,
        (90 - 65) * u.deg,
        (125 + 90) * u.deg,
        #components={'north': north_pars, 'south': south_pars},
        components={'north': north_pars},
        #components={'south': south_pars},
        ranges={'north': (0, 2000, -53, -44),
                'south': (0, 2000, -47, -40)}
    )

    # Iterate over model parameters
    fit_streamer(obs_data, model_pars, figures)

    # Create a reference coordinate system
    #center = SkyCoord(ra, dec, frame='icrs')
    #refframe = center.skyoffset_frame()


    #stream_model = {'ra': fil.ra.value*u.deg, 'dec': fil.dec.value*u.deg,
    #                'd_sky_au': d_sky_au, 'vlsr': v_lsr + vy1}
    #with open(stream_pickle, 'wb') as f:
    #    pickle.dump(stream_model, f)
    #
    #KDE_vel_rad = {'radius': xx, 'v_lsr': yy, 'dens': zz}
    #with open(vlsr_rad_kde_pickle, 'wb') as f:
    #    pickle.dump(KDE_vel_rad, f)


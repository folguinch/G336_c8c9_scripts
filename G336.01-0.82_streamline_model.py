"""Calculate streamer models.

Some of the functions were taken from the original implementation of
velocity_tools. Part of the code is a work in progress to fit a streamer, but
the engine to calculate streamer models is working.
"""
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

from common_paths import RESULTS, CONFIGS, CONTINUUM

@dataclass
class FitPars:
    """Keep track of parameters used for fitting."""
    rmin: npt.ArrayLike
    r0: npt.ArrayLike
    rc: npt.ArrayLike
    theta0: npt.ArrayLike
    phi0: npt.ArrayLike
    v_r0: npt.ArrayLike
    tied_r: bool = True

    def __iter__(self):
        for vals in product(*self._as_tuple()):
            if self.tied_r and vals[0] != vals[2]:
                print('Skipping model with rmin=', vals[0], ' and rc=', vals[2])
                continue
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

def fit_streamer(data: ObsData, model: ModelPars, outdir: Path):
    """Fit a streamer (work in progress).
    """
    # Iterate components
    for name, component in model:
        print(f'Fitting component: {name}')
        # Iterate models
        for fname, vals in component:
            # Output region
            reg_name = outdir / 'regions' / f'stream_{name}_{fname}.crtf'
            reg_name.parent.mkdir(parents=True, exist_ok=True)
            if reg_name.exists():
                print(f'Skipping {fname}: already calculated')
                continue

            # Parameters
            rmin, r0, rc, theta0, phi0, v_r0 = vals
            #if rmin != rc:
            #    print('Skipping model with rmin=', rmin, ' and rc=', rc)
            #    continue
            omega0 = np.sqrt(rc * ct.G * model.Mstar) / r0**2
            omega0 = omega0.to(1 / u.s)
            if np.abs(rmin - r0) < 100*u.au:
                r_step = 0.1 * u.au
            else:
                r_step = 10 * u.au
            r_step = 0.1 * u.au

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
                deltar=r_step,
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
    # Directories
    results = RESULTS / 'G336.01-0.82/c8c9/CH3OH'
    moment1 = (results /
               ('CH3OH_18_3_15_-17_4_14_A_vt_0_b6_c8c9_spw0_520_590_robust2_'
                'multiscale_width40_nsig3.subcube.moment1.fits'))
    continuum = (CONTINUUM /
                 ('G336.018-0.827.c8c9.selfcal.cvel.cont_avg.hogbom.'
                  'briggs.robust0.5.image.fits'))
    # In FERIA standard
    disk_pa_orig = 125 * u.deg
    disk_incl_orig = 65 * u.deg # 0 face-on
    # In velocity_tools standard
    # The rotation direction is set by the outflow axis. For clockwise this
    # should point south.
    disk_pa = 90 * u.deg + disk_pa_orig # rotation axis 
    disk_incl = 90 * u.deg - disk_incl_orig # 0 edge-on

    # Fit pars
    # For parameter exploration
    #streamers = ['north_outer', 'north_inner']
    # For best models
    streamers = ['north_outer_best',
                 'north_inner_best',
                 'south_best']
    streamers = ['north_outer_best']
    stream_params = {
        'north_outer': FitPars(
            np.array([400, 450, 500, 550, 600]) * u.au,
            np.array([2500]) * u.au,
            np.array([400, 450, 500, 550, 600]) * u.au,
            np.array([75, 80, 85]) * u.deg,
            np.array([50, 55, 60, 65, 70]) * u.deg,
            np.array([0]) * u.km/u.s,
        ),
        'north_inner': FitPars(
            # Following central source radius and pv map 
            np.array([60, 200]) * u.au,
            np.array([500]) * u.au,
            # Following central source radius and pv map
            np.array([60, 200]) * u.au,
            np.array([85, 88, 89]) * u.deg,
            # First round
            #np.array([150, 155, 160, 165, 170, 175, 180, 185, 190, 195]) * u.deg,
            # Second round
            #np.array([140, 145, 150, 155, 160, 165]) * u.deg,
            np.array([140, 145, 150, 155]) * u.deg,
            np.array([0, 0.1, 0.25, 0.5, 0.75]) * u.km/u.s,
        ),
        #'north_inner_untied': FitPars(
        #    #np.array([200, 300, 400, 450, 470, 475, 480]) * u.au,
        #    np.array([200]) * u.au,
        #    np.array([500]) * u.au,
        #    np.array([60]) * u.au,
        #    np.array([85, 88, 89]) * u.deg,
        #    # First round
        #    #np.array([150, 155, 160, 165, 170, 175, 180, 185, 190, 195]) * u.deg,
        #    # Second round
        #    #np.array([140, 145, 150, 155, 160, 165]) * u.deg,
        #    np.array([145, 155]) * u.deg,
        #    np.array([0, 0.5, 1]) * u.km/u.s,
        #    tied_r=False,
        #),
        'north_outer_best': FitPars(
            np.array([500]) * u.au,
            np.array([2500]) * u.au,
            np.array([500]) * u.au,
            np.array([80]) * u.deg,
            np.array([55]) * u.deg,
            np.array([0]) * u.km/u.s,
        ),
        'north_inner_best': FitPars(
            np.array([200]) * u.au,
            np.array([500]) * u.au,
            np.array([200]) * u.au,
            np.array([89]) * u.deg,
            np.array([145]) * u.deg,
            np.array([0.1]) * u.km/u.s,
        ),
        'south_best': FitPars(
            np.array([500]) * u.au,
            np.array([2500]) * u.au,
            np.array([500]) * u.au,
            np.array([80]) * u.deg,
            np.array([280]) * u.deg,
            np.array([2]) * u.km/u.s,
        ),
    }

    # Source properties
    position = SkyCoord('16h35m9.2585s -48d46m47.661s', frame=ICRS)
    distance = 3.1 * u.kpc
    v_lsr = -47.2 * u.km/u.s
    obs_data = ObsData(
        position,
        moment1,
        continuum,
    )

    # Run for each streamer
    for streamer in streamers:
        if 'north' in streamer:
            components = {'north': stream_params[streamer]}
        else:
            components = {'south': stream_params[streamer]}
        model_pars = ModelPars(
            position,
            distance,
            v_lsr,
            10 * u.Msun,
            disk_incl,
            disk_pa,
            components=components,
            ranges={'north': (0, 2000, -53, -44),
                    'south': (0, 2000, -47, -40)}
        )

        # Iterate over model parameters
        outdir = results / streamer
        outdir = outdir / (f'{streamer}_incl{disk_incl_orig.value:g}'
                           f'_pa{disk_pa_orig.value:g}')
        outdir.mkdir(parents=True, exist_ok=True)
        fit_streamer(obs_data, model_pars, outdir)


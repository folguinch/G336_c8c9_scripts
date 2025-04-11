"""Export streamer model in 3-D cartesian space.

This can be used to plot the model in 3-D. Not used in the paper.
"""
from pathlib import Path

import astropy.units as u
import astropy.constants as ct
import numpy as np
import velocity_tools.stream_lines as SL
from astropy.coordinates import SkyCoord

# Constants
outdir = Path(('/data/share/dihca2/combined_projects/results/G336.01-0.82/'
               'c8c9/CH3OH/streamers_best_positions'))
mstar = 10 * u.Msun

# Blue streamer
# Large scale
rmin = rc = 500 * u.au
r0 = 2500 * u.au
theta = 80 * u.deg
phi = 55 * u.deg
vr0 = 0 * u.km / u.s

# Omega
omega0 = np.sqrt(rc * ct.G * mstar) / r0**2
omega0 = omega0.to(1 / u.s)

# Coords
name = 'north_large'
fname = (f'rmin{int(rmin.value)}_'
         f'r0{int(r0.value)}_'
         f'rc{int(rc.value)}_'
         f'theta0{int(theta.value)}_'
         f'phi0{int(phi.value)}_'
         f'vr0{int(vr0.value)}')
tab_name = outdir / f'stream_positions_{name}_{fname}.ecsv'
(x, y, z), *_ = SL.xyz_stream(mass=mstar, r0=r0, theta0=theta,
                              phi0=phi, omega=omega0, v_r0=vr0,
                              inc=0*u.deg, pa=0*u.deg, rmin=rmin)
fil = SkyCoord(x, y, z, representation_type='cartesian')
fil_table = fil.to_table()
fil_table.write(tab_name, overwrite=True)

# Small scale
rmin = rc = 475 * u.au
r0 = 500 * u.au
theta = 89 * u.deg
phi = 155 * u.deg
vr0 = 0 * u.km / u.s

# Omega
omega0 = np.sqrt(rc * ct.G * mstar) / r0**2
omega0 = omega0.to(1 / u.s)

# Coords
name = 'north_small'
fname = (f'rmin{int(rmin.value)}_'
         f'r0{int(r0.value)}_'
         f'rc{int(rc.value)}_'
         f'theta0{int(theta.value)}_'
         f'phi0{int(phi.value)}_'
         f'vr0{int(vr0.value)}')
tab_name = outdir / f'stream_positions_{name}_{fname}.ecsv'
(x, y, z), *_ = SL.xyz_stream(mass=mstar, r0=r0, theta0=theta,
                              phi0=phi, omega=omega0, v_r0=vr0,
                              inc=0*u.deg, pa=0*u.deg, rmin=rmin)
fil = SkyCoord(x, y, z, representation_type='cartesian')
fil_table = fil.to_table()
fil_table.write(tab_name, overwrite=True)

# Red streamer
rmin = rc = 500 * u.au
r0 = 2500 * u.au
theta = 80 * u.deg
phi = 280 * u.deg
vr0 = 2 * u.km / u.s

# Omega
omega0 = np.sqrt(rc * ct.G * mstar) / r0**2
omega0 = omega0.to(1 / u.s)

# Coords
name = 'south'
fname = (f'rmin{int(rmin.value)}_'
         f'r0{int(r0.value)}_'
         f'rc{int(rc.value)}_'
         f'theta0{int(theta.value)}_'
         f'phi0{int(phi.value)}_'
         f'vr0{int(vr0.value)}')
tab_name = outdir / f'stream_positions_{name}_{fname}.ecsv'
(x, y, z), *_ = SL.xyz_stream(mass=mstar, r0=r0, theta0=theta,
                              phi0=phi, omega=omega0, v_r0=vr0,
                              inc=0*u.deg, pa=0*u.deg, rmin=rmin)
fil = SkyCoord(x, y, z, representation_type='cartesian')
fil_table = fil.to_table()
fil_table.write(tab_name, overwrite=True)

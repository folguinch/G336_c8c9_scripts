import astropy.constants as ct
import astropy.units as u
import numpy as np
import velocity_tools.stream_lines as SL

# Parameters for 24/07/01
r0 = 2500 * u.au
rc = 500 * u.au
mass = 10 * u.M_sun
theta0 = 80 * u.deg
v_r0 = 0 * u.km/u.s


# Calculate spherical vel
r = np.array([rc.value]) * rc.unit
omega = np.sqrt(rc * ct.G * mass) / r0**2
omega = omega.to(1 / u.s)
theta = SL.stream_line(r, mass=mass, r0=r0, theta0=theta0, omega=omega,
                       v_r0=v_r0)
v_r, v_theta, v_phi = SL.stream_line_vel(r, theta, mass=mass, r0=r0,
                                         theta0=theta0, omega=omega, v_r0=v_r0)

print(f'Radial velocity at rc = {v_r}')

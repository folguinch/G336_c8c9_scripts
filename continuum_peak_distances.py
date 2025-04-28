"""Determine the distance between the streamers peak to the central source."""
from astropy.coordinates import SkyCoord
import astropy.units as u

# Positions from Gaussian fit for ALMA1 and CARTA for the streamers
peak = SkyCoord('16h35m09.2585s -48d46m47.661s', frame='icrs')
blue_peak = SkyCoord('16h35m09.2482s -48d46m47.588s', frame='icrs')
red_peak = SkyCoord('16h35m09.2693s -48d46m47.708s', frame='icrs')

# Distances
dist_blue = peak.separation(blue_peak).to(u.arcsec)
dist_red = peak.separation(red_peak).to(u.arcsec)
print('Distance to blue streamer: ', dist_blue)
print('Distance to red streamer: ', dist_red)

"""Berkhoff et al. (1982) elliptic shoal on a 1:50 slope rotated by 20 degrees, as set up in the NHWAVE example
of this case."""
import numpy as np


def bathymetry(x, y):
    theta = np.deg2rad(20.)
    xp = x*np.cos(theta)-y*np.sin(theta)
    yp = x*np.sin(theta)+y*np.cos(theta)
    plane = np.maximum(.07, .45-.02*np.maximum(0., xp+5.84))
    inside = (xp/3)**2+(yp/4)**2 <= 1
    mound = np.where(inside, -.3+.5*np.sqrt(np.maximum(0., 1-(xp/3.75)**2-(yp/5)**2)), 0.)
    return plane-mound

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import UnivariateSpline
import pickler.shelf as shelf
from DBS.dbgrabber import dbsPull

h = 0.6777



def orbital_path(frame_set, orbital_args):
    '''
    creates coords mapping to a circular_path about a target coord
    Coords are initial generated in the orbital planes frame of reference, then transformed in to world coords
    '''
    thetas = frame_set * orbital_args.w + orbital_args.phi
    #coords in plane basis
    plane_xs = orbital_args.r * np.cos(thetas)
    plane_ys = orbital_args.r * np.sin(thetas)
    plane_zs = orbital_args.r *  0  *  thetas
    plane_coords = np.asarray([plane_xs, plane_ys, plane_zs]).T
    # transform in to world coords
    world_coords = coord_transform(orbital_args.basis[0], orbital_args.basis[1], orbital_args.basis[2], orbital_args.centre, plane_coords, inv=False, homog=False, kieren=False).T
    return world_coords


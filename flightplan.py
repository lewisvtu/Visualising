#!/bin/env python2.7
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import pickler.shelf as shelf
import math
from DBS.dbgrabber import dbsPull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import UnivariateSpline
import timeit


h = 0.6777
# SQL = """
#     SELECT
#         DES.GalaxyID,
#         PROG.SnapNum,
#         PROG.Mass,
#         PROG.CentreOfPotential_x,
#         PROG.CentreOfPotential_y,
#         PROG.CentreOfPotential_z,
#         PROG.Redshift
#     FROM
#         RefL0100N1504_Subhalo as PROG with(forceseek),
#         RefL0100N1504_Subhalo as DES,
#         RefL0100N1504_Aperture as AP
#     WHERE
#         DES.SnapNum = 28 and
#         DES.MassType_Star > 1.0e9 and
#         DES.MassType_DM > 5.0e10 and
#         PROG.GalaxyID between DES.GalaxyID and DES.TopLeafID and
#         AP.ApertureSize = 30 and
#         AP.GalaxyID = DES.GalaxyID and
#         AP.Mass_Star > 1.0e9
#     ORDER BY
#         PROG.GalaxyID,
#         PROG.SnapNum
# """

# # Grabs new data from db based on sql. If file name already exists, it loads that data instead

# filename = "FollowProgs2.p"

# raw_dbs = dbsPull(SQL, filename)

# shelf.push(raw_dbs, "followup2")

# dbs_data = shelf.pull("followup2")

# interesting_ids = {

#     13660659: 28,
#     13793733: 28,
#     13722615: 28,
#     20440704: 28,
#     17891603: 28,
#     14784533: 28

# }
# gals = np.asarray([list(gal)[3:] for gal in dbs_data if gal[0] in interesting_ids.keys() and gal[1] == interesting_ids[gal[0]]])
# # Makes all galaxies have the form [sf,x,y,z] sorted by sf
# gals[:,-1] = 1.0 / (1.0 + gals[:,-1])
# gals = gals[np.argsort(gals[:,3])]
# gals = gals[:,[3,0,1,2]]

class path():
    def __init__(self, no_of_frames, collection):
        self.collection = collection
        self.frames = np.arange(no_of_frames)
        self.coord_spline = self.gen_coord_spline()
        self.coords = self.coord_spline(self.frames)
        self.look_at_dirs = self.get_look_ats()
        self.tangent_dirs = self.gen_tangent_vectors()
        self.final_dirs = self.gen_final_basis()

    def gen_coord_spline(self):
        coords = []
        frames = []
        for galaxy, path_function, frame_set, path_args in self.collection:
            coords = coords + list(path_function(frame_set, path_args))
            frames = frames + list(frame_set)
        coords = np.asarray(coords)
        frames = np.asarray(frames)
        bundle = np.c_[frames, coords]
        spl = spline3D(bundle)
        return spl



    def get_look_ats(self):
        frames = []
        look_at_points = []
        for t_gal, path_f, fs, path_args in self.collection:
            frames = frames + list(fs)
            look_at_dirs = t_gal[1:] - path_f(fs, path_args)
            look_at_dirs = look_at_dirs / np.linalg.norm(look_at_dirs, axis=1)[:,None]
            look_at_points = look_at_points + list(look_at_dirs)
        frames = np.asarray(frames)
        look_at_points = np.asarray(look_at_points)
        return spline3D(np.c_[frames, look_at_points], k=2)
        
    def gen_tangent_vectors(self):
        frame_nos = self.frames
        path_function = self.coord_spline
        look_at_dirs = self.look_at_dirs
        derivs = np.zeros((len(frame_nos), 3))
        d_frame = 0.01
        for index in range(len(frame_nos)):
            frame_no = frame_nos[index]
            derivs[index] = path_function(frame_no + d_frame/2) - path_function(frame_no - d_frame/2)
        #einsum computes the vector dots for each pair in the array
        basis_1 = derivs - np.transpose(np.einsum("ij,ij->i", derivs, look_at_dirs) * np.transpose(look_at_dirs))
        basis_1 = basis_1 / np.linalg.norm(basis_1, axis=1)[:,None]
        return basis_1

    def gen_final_basis(self):
        basis_2 = np.cross(self.look_at_dirs, self.tangent_dirs)
        basis_2 = basis_2 / np.linalg.norm(basis_2, axis=1)[:,None]
        return basis_2

def circular_path(frame_nos, args):
    '''
    Given a frame number, returns x,y,z coords for the camera at that frame

    args:
        frame_nos: the frame number/ array of frame numbers to compute the positions of
        args: list of form []
    '''
    decay = 0.5
    frame_nos = np.asarray(frame_nos)
    target_coords = np.asarray(args[0])
    rad = args[1]
    orbits = args[2]
    frames = args[3]
    dir = args[4]

    z_factor = args[5]
    ang_int = orbits * 2*np.pi / frames
    x_coords = target_coords[0] + dir * rad * np.sin(frame_nos * ang_int)
    y_coords = target_coords[1] + dir * rad * np.cos(frame_nos * ang_int)
    z_coords = target_coords[2] + dir * (rad* np.sin(frame_nos * ang_int)) * z_factor

    return np.transpose(np.asarray([x_coords, y_coords, z_coords]))

def straight_path(frame_nos, args):
    '''
    Gens points along a straight line between two coords
    Args:
        Frame_nos: Single/list of frame/s to generate points for
        args: [start_coord, end_coord, frame_length_of_path]
    Returns:
        coords: [[x,y,z], ...] list of coord/s on the line for the given frames
            '''
    frame_nos = np.asarray([frame_nos])
    frames = args[2]
    start_coords = args[0]
    end_coords = args[1]
    c_coords = np.transpose(np.asarray([start_coords]*len(frame_nos))) + np.transpose(np.asarray([end_coords - start_coords]*len(frame_nos))) * (frame_nos/frames)
    return np.transpose(c_coords)

def get_scalefactors(start_sf, end_sf, frames):
    '''
    returns list of #frames scale factors scaled uniformly in log10 between start_sf and end_sf
    '''
    array_log_sf = np.linspace(np.log10(start_sf), np.log10(end_sf), frames)
    array_sf = np.power(10, array_log_sf)
    return array_sf[::-1]


class spline3D():
    '''
    class for 3d splines.

    '''
    def __init__(self, bundle, k=3):
        '''
        Creates the spline3D object
        Args:
            bundle: [fks, xks, yks, zks] for pairs of know frames, xs, ys, and zs
            k: integer order of the fitting function for numpy splines. Values accepted [1,5]
        Returns:
            spline3D object
        '''
        fks, xks, yks, zks = np.transpose(bundle)
        self.x_spline = UnivariateSpline(fks, xks, k=k)
        self.y_spline = UnivariateSpline(fks, yks, k=k)
        self.z_spline = UnivariateSpline(fks, zks, k=k)

    def __call__(self, fs, args=None):
        '''
        Generates new points for given frames
        Args:
            fs: New frames to generate points for
        Returns:
            [[x,y,z], ...] list of new coords in order of given frames
        '''
        xs = self.x_spline(fs)
        ys = self.y_spline(fs)
        zs = self.z_spline(fs)
        return np.transpose(np.asarray([xs, ys, zs]))
    
'''
galaxy : [frames, path_function, path_args]
'''
# collection = np.asarray([
#     [gals[0], circular_path, np.arange(20, dtype=float), [gals[0,1:], 5.0, 1, 20, -1, 0.5]],
#     #[gals[1], circular_path, np.arange(20, dtype=float) + 40, [gals[1,1:], 5.0, 1, 20, 1, -0.5]],
#     [gals[2], circular_path, np.arange(20, dtype=float) + 40, [gals[2,1:], 5.0, 1, 20, 1, 2.5]],
#     [gals[3], circular_path, np.arange(20, dtype=float) + 120, [gals[3,1:], 5.0, 1, 20, 1, 0.75]],
#     #[gals[4], circular_path, np.arange(20, dtype=float) + 160, [gals[4,1:], 5.0, 1, 20, 1, -0.5]],
#     #[gals[5], circular_path, np.arange(20, dtype=float) + 200, [gals[5,1:], 5.0, 1, 20, 1, 1.5]]
# ])
# interested = [0,2,3]
# # collection = np.asarray([
# #     [gals[0], straight_path, np.arange(50), [gals[0,1:] + [3,3,-10], gals[0,1:] + [3,3, 10], 50.0]]
# # ])
# everything = path(140, collection)
# frames = everything.frames
# sfs = get_scalefactors(1,1,len(frames))
# xs, ys, zs = np.transpose(everything.coords)
# v3xs, v3ys, v3zs = np.transpose(everything.look_at_dirs)
# v1xs, v1ys, v1zs = np.transpose(everything.tangent_dirs)
# v2xs, v2ys, v2zs = np.transpose(everything.final_dirs)
# fig = plt.figure()
# ax = fig.add_subplot(111, projection="3d")
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("z")
# ax.plot(xs, ys, zs)
# for ele in interested:
#     ax.scatter(gals[ele,1], gals[ele,2], gals[ele,3])

# ax.quiver(xs,ys,zs, v1xs, v1ys, v1zs, color="#682860", pivot="tail")
# ax.quiver(xs,ys,zs, v2xs, v2ys, v2zs, color="#000000", pivot="tail")
# ax.quiver(xs,ys,zs, v3xs, v3ys, v3zs, color="#FF0000", pivot="tail")
# plt.show()


# setspace = np.transpose(np.asarray([frames, sfs, xs, ys, zs, v1xs, v1ys, v1zs, v2xs, v2ys, v2zs, v3xs, v3ys, v3zs]))
# #print setspace
# np.savetxt("straight_galbig.txt", setspace, fmt="%i %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f" )
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
from timeit import default_timer as timer


h = 0.6777
SQL = """
    SELECT
        DES.GalaxyID,
        PROG.SnapNum,
        PROG.MassType_Star,
        (PROG.CentreOfPotential_x * %f) as x,
        (PROG.CentreOfPotential_y * %f) as y,
        (PROG.CentreOfPotential_z * %f) as z,
        PROG.Redshift
    FROM
        RefL0100N1504_Subhalo as PROG with(forceseek),
        RefL0100N1504_Subhalo as DES,
        RefL0100N1504_Aperture as AP
    WHERE
        DES.SnapNum = 28 and
        DES.MassType_Star > 1.0e10 and
        DES.MassType_DM > 1.0e11 and
        PROG.GalaxyID between DES.GalaxyID and DES.TopLeafID and
        AP.ApertureSize = 30 and
        AP.GalaxyID = DES.GalaxyID and
        AP.Mass_Star > 1.0e9
    ORDER BY
        PROG.GalaxyID,
        PROG.SnapNum
""" % (h,h,h)

# # Grabs new data from db based on sql. If file name already exists, it loads that data instead

filename = "scaledDB.p"

# dbs_data = dbsPull(SQL, filename)
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

# shelf.push(gals, filename)

class Path():
    def __init__(self, no_of_frames, collection):
        self.collection = collection
        self.frames = np.arange(no_of_frames, dtype=int)
        self.coord_spline = self.gen_coord_spline()
        self.coords = self.coord_spline(self.frames)
        # self.basis_z = self.get_to_targets()
        # self.basis_x = self.orthonormalise(self.gen_tangent_vectors(), self.basis_z)
        # self.basis_y = self.cross_basis(self.basis_z, self.basis_x)
        self.basis_z = self.gen_tangent_vectors()
        self.basis_y = self.orthonormalise(np.asarray([[0,0,1]]*no_of_frames), self.basis_z)
        self.basis_x = self.cross_basis(self.basis_y, self.basis_z)

    def gen_coord_spline(self):
        coords = []
        frames = []
        for galaxy, path_function, frame_set, path_args in self.collection:
            coords = coords + list(path_function(frame_set, path_args))
            frames = frames + list(frame_set)
        coords = np.asarray(coords)
        frames = np.asarray(frames)
        bundle = np.c_[frames, coords]
        spl = Spline3D(bundle)
        return spl

    def get_to_targets(self):
        look_at_dirs = np.zeros((len(self.frames),3))
        calced_frames = []
        for galaxy, path_function, frame_set, path_args in self.collection:
            target_coords = np.asarray(galaxy[1:])
            path_coords = self.coord_spline(frame_set)
            #print target_coords, path_coords
            look_vectors = target_coords - path_coords
            look_at_dirs[frame_set] = look_vectors / np.linalg.norm(look_vectors, axis=1)[:,None]
            calced_frames = calced_frames + list(frame_set)
        look_at_dirs = np.c_[self.frames, look_at_dirs]
        interp_frames = [frame[0] for frame in look_at_dirs if np.linalg.norm(frame[1:]) == 0]
        print interp_frames
        return look_at_dirs[:,1:]
        
    def gen_tangent_vectors(self):
        frame_nos = self.frames
        path_function = self.coord_spline
        derivs = np.zeros((len(frame_nos), 3))
        d_frame = 0.01
        for index in range(len(frame_nos)):
            frame_no = frame_nos[index]
            derivs[index] = path_function(frame_no + d_frame/2) - path_function(frame_no - d_frame/2)
        derivs = derivs / np.linalg.norm(derivs, axis=1)[:,None]
        return derivs

    def orthonormalise(self, new_vects, basis_1):
        #einsum computes the vector dots for each pair in the array
        basis_2 = new_vects - np.transpose(np.einsum("ij,ij->i", new_vects, basis_1) * np.transpose(basis_1))
        basis_2 = basis_2 / np.linalg.norm(basis_2, axis=1)[:,None]
        return basis_2

    def cross_basis(self, basis_1, basis_2):
        basis_3 = np.cross(basis_1, basis_2)
        basis_3 = basis_3 / np.linalg.norm(basis_3, axis=1)[:,None]
        return basis_3


def circular_path(frame_nos, circ_args):
    '''
    Given a frame number, returns x,y,z coords for the camera at that frame

    args:
        frame_nos: the frame number/ array of frame numbers to compute the positions of
        circ_args: a SpiralArgs object defining the path parameters
    '''
    #Working in std spherical polars s.t. z=rcos(theta')
    thetas = (frame_nos * circ_args.theta_dot) + circ_args.theta_off
    phis   = (frame_nos * circ_args.phi_dot)   + circ_args.phi_off
    radii  = (frame_nos * circ_args.rad_dot)   + circ_args.rad_off
    print thetas, phis, radii
    x_coords = circ_args.centre[0] + radii * np.sin(thetas) * np.cos(phis)
    y_coords = circ_args.centre[1] + radii * np.sin(thetas) * np.sin(phis)
    z_coords = circ_args.centre[2] + radii * np.cos(thetas)
    #print x_coords, y_coords, z_coords
    return np.transpose(np.asarray([x_coords, y_coords, z_coords]))

class SpiralArgs():
    '''
    hold params for spiral paths
    Args:
        centre: the centre of the path of the form [x,y,z]
        polar_vels: polar velocities, r_dot in Mpsc/h, theta_dot/ phi_dot in orbits per frame
    '''
    def __init__(self, centre, polar_vels, coord_offsets):
        self.centre    = centre
        self.rad_dot   = polar_vels[0]
        self.theta_dot = polar_vels[1] * 2 * np.pi
        self.phi_dot   = polar_vels[2] * 2 * np.pi
        self.rad_off   = coord_offsets[0]
        self.theta_off = coord_offsets[1]
        self.phi_off   = coord_offsets[2]



def straight_path(frame_nos, args):
    '''
    Gens points along a straight line between two coords
    Args:
        Frame_nos: Single/list of frame/s to generate points for
        args: [start_coord, end_coord, frame_length_of_path]
    Returns:
        coords: [[x,y,z], ...] list of coord/s on the line for the given frames
    '''
    assert len(args) == 3
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

class Spline3D:
    '''
    class for 3d splines.

    '''
    def __init__(self, bundle, k=3):
        '''
        Creates the S object
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


if __name__ == "__main__":
    start = timer()
    gals = shelf.pull(filename)    
    '''
    galaxy : [frames, path_function, path_args]
        circle path_args: [centre_pos, radius, #orbits, #frames, direction, z_scale, phi]
    '''
    collection = np.asarray([
        [gals[0], circular_path, np.arange(20, dtype=int), SpiralArgs(gals[0,1:], [0.,1/20.,1/20], [5.,0.,0.,])]
        # [gals[2], circular_path, np.arange(20, dtype=int) + 15, SpiralArgs(gals[2,1:], [1.,1.,1.], [0.,0.,0.,]).set_vels(20., 0.75)],
        # [gals[3], circular_path, np.arange(10, dtype=int) + 45, SpiralArgs(gals[3,1:], [1.,1.,1.], [0.,0.,0.,]).set_vels(10, 0.5)],
    ])
    interested = [0]
    # collection = np.asarray([
    #     [gals[0], straight_path, np.arange(50), [gals[0,1:] + [3,3,-10], gals[0,1:] + [3,3, 10], 50.0]]
    # ])
    everything = Path(20, collection)
    end = timer()
    
    print "Time taken: %f" % (end-start)
    frames = everything.frames
    sfs = get_scalefactors(1,1,len(frames))
    xs, ys, zs = np.transpose(everything.coords)
    v3xs, v3ys, v3zs = np.transpose(everything.basis_z)
    v1xs, v1ys, v1zs = np.transpose(everything.basis_x)
    v2xs, v2ys, v2zs = np.transpose(everything.basis_y)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.plot(xs, ys, zs)
    for ele in interested:
        ax.scatter(gals[ele,1], gals[ele,2], gals[ele,3])

    ax.quiver(xs,ys,zs, v1xs, v1ys, v1zs, color="#682860", pivot="tail")
    ax.quiver(xs,ys,zs, v2xs, v2ys, v2zs, color="#000000", pivot="tail")
    ax.quiver(xs,ys,zs, v3xs, v3ys, v3zs, color="#FF0000", pivot="tail")
    plt.show()


    setspace = np.transpose(np.asarray([frames, sfs, xs, ys, zs, v1xs, v1ys, v1zs, v2xs, v2ys, v2zs, v3xs, v3ys, v3zs]))
    #print setspace
    np.savetxt("tangential_splines.txt", setspace, fmt="%i %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f" )



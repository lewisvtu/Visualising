from __future__ import division

import matplotlib.pyplot as plt
import numpy as np

import utils
from utils import coord_transform, cross_basis, orthonormalise, Spline3D


class SplinePath(object):
    '''Spline paths that stitch together general paths'''
    def __init__(self, knowns, k=3):
        '''Args:
            frame_set [reals]: set of frame numbers associated with the coordinates, used to
                generate the spline
            knowns [array]: set of known values, [frame,x,y,z]
            k [0,1,2,3,4]: spline smoothing factor, defined by scipy'''
        self.spline = Spline3D(knowns, k)

    def __call__(self, frame_set):
        '''get coords for new frame value/s
        Args:
            frames [real/s]: frame number/s to generate coords for
        Returns:
            coords [array]: coords of new frames [x,y,z]'''
        return self.spline(frame_set)

class OrbitalPath(object):
    '''Creates coords mapping to a circular_path about a target coord
    Coords are initial generated in the orbital planes frame of reference,
    then transformed in to world coords'''
    def __init__(self, centre, plane_normal, rad_vel, rpf,
                 rad_off, rev_off, helix_vel=0, helix_off=0):
        '''sets up the arguments for the orbital path
        Args:
            centre [array]: centre of the orbital path, [x,y,z]
            plane_normal [array]: vector [x,y,z] in the normal to the orbital plane
            rad [real]: radius of orbital path
            rpf [real]: angular velocity in revolutions per frame
            rev_off [real]: angular offset of the start of the orbital path, in units of
                revolution'''
        self.centre = centre
        self.norm = plane_normal / np.linalg.norm(plane_normal)
        self.ang_vel = 2*np.pi*rpf
        self.ang_off = 2*np.pi*rev_off
        self.rad_vel = rad_vel
        self.rad_off = rad_off
        self.hel_vel = helix_vel
        self.hel_off = helix_off
        self.basis = self.get_plane_basis(self.norm)


    def get_plane_basis(self, plane_normal):
        '''Generates a basis for the plane to be used when transforming in to world coords
        Args:
            plane_normal [array]: Normal to the orbital plane
        Returns:
            basis [array]: Array of basis vectors for the orbital reference frame
                [x_basis, y_basis, z_basis] for basis [bx,by,bz] for bi real'''
        temp_vect = np.array([1., 0., 0.])
        if abs(np.dot(temp_vect, plane_normal)) == 1.:
            #If normal in x dir, form basis from y instead
            temp_vect = np.asarray([0., 1., 0.])
        basis_1 = temp_vect / np.linalg.norm(temp_vect)
        basis_1 -= np.dot(basis_1, plane_normal) * basis_1
        basis_1 /= np.linalg.norm(basis_1)
        basis_2 = np.cross(basis_1, plane_normal)
        return np.asarray([basis_1, basis_2, plane_normal])

    def __call__(self, frame_set):
        thetas = frame_set * self.ang_vel + self.ang_off
        radii = frame_set * self.rad_vel + self.rad_off
        #coords in obital axis frame
        frame_xs = radii     * np.cos(thetas)
        frame_ys = radii     * np.sin(thetas)
        frame_zs = frame_set * self.hel_vel + self.hel_off
        frame_coords = np.asarray([frame_xs, frame_ys, frame_zs])
        #transform in to world coords
        world_coords = coord_transform(self.basis[0], self.basis[1], self.basis[2],
                                       self.centre, frame_coords, inv=False, homog=False, tran=True)
        return np.asarray(world_coords)

def vector_derivs(frame_set, path_function, d_frame=0.01):
    '''Calculates the vector derivatives/ tangents to the path.
    Args:
        frame_set [real/s]: list of frame numbers to generate tangents at
        path_function [PathObject]: callable object defining the path. Takes any number of
            frames as an arg, returning the pos at those frames
        d_frame (optional) [real]: The dx used to find path difference
    Returns:
        derivs [array]: array of tangent vectors, [dx,dy,dz] normalised'''
    derivs = np.zeros((len(frame_set), 3))
    for index in range(len(frame_set)):
        frame_no = np.asarray([frame_set[index]])
        derivs[index] = path_function(frame_no + d_frame/2) - path_function(frame_no - d_frame/2)
    derivs /= np.linalg.norm(derivs, axis=1)[:, None]
    return derivs

def look_at_vectors(path_coords, target_coords):
    '''calculates the look at vectors for given target coords and camera path coords
    Args:
        path_coords [array]: array of coords of the camera [x,y,z]
        target_coords [array]: array of coords to look at [x,y,z]
    Returns:
        target_coords [array]: array of look at vectors'''
    look_dirs = target_coords - path_coords
    look_dirs /= np.linalg.norm(look_dirs, axis=1)[:, None]
    return look_dirs


if __name__ == "__main__":#
    print "Actually started running -_- z z z"
    no_frames = 30
    gal1_coords = [11.2204,16.5994,12.0005]
    #gal1_coords = np.asarray([0., 0., 0.])
    test_coords = [
        [0, 0, 0],
        [10, 0, 0]
    ]

    path = OrbitalPath(gal1_coords, [0, 0, 1], 0, 0.25/30, 5, 0, 0, 0)
    frames = np.arange(no_frames)
    look_pos = np.tile(gal1_coords, (no_frames,1))
    path_coords = path(frames)
    basis_3 = look_at_vectors(path_coords, look_pos)
    tangents = vector_derivs(frames, path)
    basis_1 = orthonormalise(tangents, basis_3)
    basis_2 = cross_basis(basis_3, basis_1)

    # frames = np.arange(60)
    # look_pos = np.asarray([[0,0,0]]*30 + [[10,0,0]]*30)
    # first_path = OrbitalPath(test_coords[0], [1,0,0], -5/20, 2/30, 5, 0, -3/20, 0)
    # sec_path = OrbitalPath(test_coords[1], [-1,0,0], 5/20, 2/30, 5, 0, -3/20, )

    # first_set = np.c_[frames[:20], first_path(frames[:20])]
    # sec_set = np.c_[frames[-20:], sec_path(frames[-20:])]
    # tot_set = np.r_[first_set, sec_set]
    # path = sec_path
    # path_coords = path(frames)
    # print path_coords

    # basis_3 = look_at_vectors(path_coords, look_pos)
    # tangents = vector_derivs(frames, path)
    # basis_1 = orthonormalise(tangents, basis_3)
    # basis_2 = cross_basis(basis_3, basis_1)
    #Plotting bits
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.scatter(gal1_coords[0], gal1_coords[1], gal1_coords[2])
    ax.quiver(path_coords[:, 0], path_coords[:, 1], path_coords[:, 2],
              basis_1[:,0], basis_1[:, 1], basis_1[:, 2], pivot="tail", color="#FF0000")
    ax.quiver(path_coords[:, 0], path_coords[:, 1] ,path_coords[:, 2],
              basis_2[:, 0], basis_2[:, 1], basis_2[:, 2], pivot="tail", color="#00FF00")
    ax.quiver(path_coords[:, 0],path_coords[:, 1],path_coords[:, 2],
              basis_3[:, 0],basis_3[:, 1],basis_3[:, 2], pivot="tail", color="#0000FF")
    plt.show()
    #make file
    sfs = utils.get_scalefactors(0.44,0.6,no_frames)
    utils.gen_flight_file(frames, sfs, path_coords, np.asarray([basis_1, basis_2,basis_3]), "Orbit_through_les_time.txt")

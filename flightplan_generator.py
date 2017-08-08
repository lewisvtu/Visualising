from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import UnivariateSpline
from utils import coord_transform, orthonormalise, cross_basis

h = 0.6777



def orbital_path_function(frame_set, orbital_args):
    '''Creates coords mapping to a circular_path about a target coord
    Coords are initial generated in the orbital planes frame of reference,
    then transformed in to world coords
    
    Args:
        frame_set [real/s]: Frame numbers for which to create points
        orbital_args [OrbitalArgs]: argument class defining the orbital path
    Returns:
        world_coords [array]: list of coords for the path, [x,y,z]'''
    thetas = frame_set * orbital_args.w + orbital_args.phi
    #coords in plane basis
    plane_xs = orbital_args.r * np.cos(thetas)
    plane_ys = orbital_args.r * np.sin(thetas)
    plane_zs = orbital_args.r *  0  *  thetas
    plane_coords = np.asarray([plane_xs, plane_ys, plane_zs]).T
    # transform in to world coords
    world_coords = coord_transform(orbital_args.basis[0], orbital_args.basis[1], orbital_args.basis[2], orbital_args.centre, plane_coords, inv=False, homog=False).T
    return world_coords


class OrbitalArgs:
    '''holds args for orbital paths'''
    def __init__(self, centre, plane_normal, rad, rpf, rev_off):
        '''sets up the arguments for the orbital path function
        Args:
            centre [array]: centre of the orbital path, [x,y,z]
            plane_normal [array]: vector [x,y,z] in the normal to the orbital plane
            rad [real]: radius of orbital path
            rpf [real]: angular velocity in revolutions per frame
            rev_off [real]: angular offset of the start of the orbital path, in units of revolution'''
        self.centre = centre
        self.norm   = plane_normal / np.linalg.norm(plane_normal)
        self.w      = 2*np.pi*rpf
        self.phi    = 2*np.pi*rev_off
        self.r      = rad
        self.basis  = self.get_plane_basis(self.norm)

    def get_plane_basis(self, plane_normal):
        '''Generates a basis for the plane to be used when transforming in to world coords
        Args:
            plane_normal [array]: Normal to the orbital plane
        Returns:
            basis [array]: Array of basis vectors for the orbital reference frame [x_basis, y_basis, z_basis]
                for basis [bx,by,bz] for bi real'''
        temp_vect = np.asarray([1,0,0])
        if np.dot(temp_vect, plane_normal) == 1.:
            #If normal in x dir, form basis from y instead
            temp_vect = np.asarray([0,1,0])
        basis_1  = temp_vect / np.linalg.norm(temp_vect)
        basis_1 -= np.dot(basis_1, plane_normal) * basis_1
        basis_1 /= np.linalg.norm(basis_1)
        basis_2  = np.cross(basis_1, plane_normal)
        return np.asarray([basis_1, basis_2, plane_normal])

def vector_derivs(frame_set, path_function, path_args, d_frame = 0.01):
    '''Calculated the vector derivatives/ tangents to the path.
    Args:
        frame_set [real/s]: list of frame numbers to generate tangents at
        path_function [function]: function defining the path. Must take frame number as first arg,
            followed by path argument object
        path_args [arg_object]: object containing arguments defining the path
        d_frame (optional) [real]: The dx used to find path difference
    Returns:
        derivs [array]: array of tangent vectors, [dx,dy,dz] normalised'''
    derivs = np.zeros((len(frame_set), 3))
    for index in range(len(frame_set)):
        frame_no = np.asarray([frame_set[index]])
        derivs[index] = path_function(frame_no + d_frame/2, path_args) - path_function(frame_no - d_frame/2, path_args)
    derivs /= np.linalg.norm(derivs, axis=1)[:,None]
    return derivs

def look_at_vectors(path_coords, target_coords):
    look_dirs  = target_coords - path_coords
    look_dirs  /= np.linalg.norm(look_dirs, axis=1)[:,None]
    return look_dirs

if __name__ == "__main__":
    args = OrbitalArgs([0,0,0], [0,1,1], 5, 1/15, 0)
    frames = np.arange(15)
    temp_set = np.asarray([[0,0,1]]*15)
    path_coords = orbital_path_function(frames, args)
    basis_1 = look_at_vectors(path_coords, temp_set)
    tangents = vector_derivs(frames, orbital_path_function, args)
    basis_2 = orthonormalise(tangents, basis_1)
    basis_3 = cross_basis(basis_1, basis_2)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.quiver(path_coords[:,0],path_coords[:,1],path_coords[:,2], basis_1[:,0],basis_1[:,1],basis_1[:,2], pivot="tail", color="#FF0000")
    ax.quiver(path_coords[:,0],path_coords[:,1],path_coords[:,2], basis_2[:,0],basis_2[:,1],basis_2[:,2], pivot="tail", color="#00FF00")
    ax.quiver(path_coords[:,0],path_coords[:,1],path_coords[:,2], basis_3[:,0],basis_3[:,1],basis_3[:,2], pivot="tail", color="#0000FF")
    plt.show()
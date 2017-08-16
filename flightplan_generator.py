from __future__ import division

import matplotlib.pyplot as plt
import numpy as np

import utils
from utils import coord_transform, cross_basis, orthonormalise, Spline3D, Interp3D


class SplinePath(object):
    '''Spline paths that stitch together general paths'''
    def __init__(self, knowns, k=3):
        '''Args:
            frame_set [reals]: set of frame numbers associated with the coordinates, used to
                generate the spline
            knowns [array]: set of known values, [frame,x,y,z]
            k [0,1,2,3,4]: spline smoothing factor, defined by scipy'''
        self.spline = Spline3D(knowns, k)

    def __call__(self, frame_set, mdf=False):
        '''get coords for new frame value/s
        Args:
            frames [real/s]: frame number/s to generate coords for
        Returns:
            coords [array]: coords of new frames [x,y,z]'''
        if mdf:
            frame_set = frame_set[:,0]
        coords = self.spline(frame_set)
        # print "spline path call"
        # print frame_set, "\n ----------\n", coords
        return coords

class OrbitalPath(object):
    '''Creates coords mapping to a circular_path about a target coord
    Coords are initial generated in the orbital planes frame of reference,
    then transformed in to world coords'''
    def __init__(self, centre_bundle, nx,ny,nz, rad_vel, rpf,
                 rad_off, rev_off, helix_vel=0, helix_off=0):
        '''sets up the arguments for the orbital path
        Args:
            centre [array]: centre of the orbital path, [x,y,z]
            plane_normal [array]: vector [x,y,z] in the normal to the orbital plane
            rad [real]: radius of orbital path
            rpf [real]: angular velocity in revolutions per frame
            rev_off [real]: angular offset of the start of the orbital path, in units of
                revolution'''
        plane_normal = np.asarray([nx,ny,nz])
        self.norm = plane_normal / np.linalg.norm(plane_normal)
        self.ang_vel = 2*np.pi*rpf
        self.ang_off = 2*np.pi*rev_off
        self.rad_vel = rad_vel
        self.rad_off = rad_off
        self.hel_vel = helix_vel
        self.hel_off = helix_off
        self.basis = self.get_plane_basis(self.norm)
        self.centre_func = Interp3D(centre_bundle)


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

    def __call__(self, frame_set, mdf=False):
        #print frame_set
        if mdf:
            frame_set = frame_set[:,0]
        frame_set = np.asarray(frame_set)
        thetas = frame_set * self.ang_vel + self.ang_off
        radii = frame_set * self.rad_vel + self.rad_off
        #coords in obital axis frame
        frame_xs = radii     * np.cos(thetas)
        frame_ys = radii     * np.sin(thetas)
        frame_zs = frame_set * self.hel_vel + self.hel_off
        frame_coords = np.asarray([frame_set, frame_xs, frame_ys, frame_zs]).T
        #print frame_coords
        #transform in to world coords
        world_coords = []
        for bundle in np.atleast_2d(frame_coords):
            centre = self.centre_func(bundle[0])
            frame_coord = bundle[1:]
            frame_coord.shape = (3,1)
            centre.shape = 3
            world_coords.append(coord_transform(self.basis[0], self.basis[1], self.basis[2],
                                       centre, frame_coord, inv=False, homog=False, tran=True))
        world_coords = np.asarray(world_coords)
        world_coords.shape = (len(world_coords), 3)
        #print "obital path call"
        return world_coords

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
        frame_no = frame_set[index]
        derivs[index] = path_function(frame_no + d_frame/2) - path_function(frame_no - d_frame/2)
    derivs /= np.linalg.norm(derivs, axis=1)[:, None]
    return derivs

def look_at_vectors(path_coords, target_coords, weights):
    prim_look = target_coords[:, :3] - path_coords
    sec_look =  target_coords[:, 3:] - path_coords
    #print prim_look, "\n----------\n", sec_look, "\n-----------\n"
    tot_look = prim_look * weights[:,np.newaxis] + sec_look * (1 - weights[:,np.newaxis])
    tot_look = tot_look.T / np.linalg.norm(tot_look, axis=1)[:None]
    return tot_look.T

def gen_centre_bundle(frames, coords, const=True):
    if const:
        rep = [list(coords)]*len(frames)
        return np.c_[frames, rep]

class CombinedPath(object):
    def __init__(self, func_domain):
        self.func_domain = func_domain
        self.fix_domain_holes()
        #print self.func_domain[:,1]

    def fix_domain_holes(self):
        new_func_domain = list(np.copy(self.func_domain))
        for index in range(len(self.func_domain)-1):
            current_end = self.func_domain[index][1]
            next_start = self.func_domain[index+1][0]
            if next_start - current_end > 0:
                #if gap, make new spline
                before_func = self.func_domain[index][2]
                before_frames = np.arange(current_end - 1, current_end + 1)
                after_func = self.func_domain[index + 1][2]
                after_frames = np.arange(next_start, next_start+2)
                before_bundle = np.c_[before_frames, before_func(before_frames)]
                after_bundle = np.c_[after_frames, after_func(after_frames)]
                tot_bundle = np.r_[before_bundle, after_bundle]
                spl_path = SplinePath(tot_bundle)
                spl_func_dom = np.asarray([current_end-1, next_start+1, spl_path])
                #print self.func_domain, "\n---------------\n", spl_func_dom
                new_func_domain.append(list(spl_func_dom))
        self.func_domain = np.asarray(new_func_domain)

    def __call__(self, frames):
        conditions = []
        frames = np.atleast_1d(frames)
        for dom_s, dom_e, func in self.func_domain:
            dom_mask = (frames >= dom_s) * (frames < dom_e)
            conditions.append(dom_mask)
        frames = np.asarray([list(frames)]*3, dtype="f8").T
        out = np.piecewise(frames, conditions, self.func_domain[:,2], mdf=True)
        return out

def gen_look_bundle(t_data, no_frames):
    look_bundle = np.zeros((no_frames, 6))
    tmp = np.asarray([[no_frames,no_frames,0.,0.,0.]])
    print t_data, "\n---------\n", tmp
    t_data = np.r_[t_data, tmp]
    for index in range(len(t_data) - 1):
        cds, cde, ctgx, ctgy, ctgz = t_data[index]
        ctg = np.asarray([ctgx,ctgy,ctgz])
        nds, nde, ntgx, ntgy, ntgz = t_data[index+1]
        ntg = np.asarray([ntgx,ntgy,ntgz])
        print ctg,ntg
        look_bundle[cds:nds] = np.r_[ctg,ntg]
    return look_bundle


if __name__ == "__main__":
    print "Actually started running -_- z z z"
    #gal1_coords = [11.2204,16.5994,12.0005]
    gal1_coords = np.asarray([0., 0., 0.])
    test_coords = np.asarray([
        [3.44255, 14.94486, 2.22045],
        [9.17897, 13.82553, 10.63093],
        [1,1,1],
        [5,5,5]
    ])

    centre_coords = test_coords
    inp_data = np.asarray([
    #   [domain], coords at centre of montion                               ,rotaxis,rv,av  ,ro,ao,hv,ho]
        [-1.,10., centre_coords[0,0], centre_coords[0,1], centre_coords[0,2], 1,2,1, 0, 0.5/10, 5, 1/2, 0, 0],
        [20.,31., centre_coords[1,0], centre_coords[1,1], centre_coords[1,2], 1,0,1, 0, 0.5/10, 5, 1/4, 0, 0],
        [45.,61., centre_coords[2,0], centre_coords[2,1], centre_coords[2,2], -0.5,0,1, 0, 0.5/10, 5, 0, 0, 0]
        #[70.,76., centre_coords[3,0], centre_coords[3,1], centre_coords[3,2], -0.5,0,1, 0, 0.5/10, 5, 0, 0, 0]
    ])

    targ_data = np.asarray([
        [0.,10., test_coords[0,0], test_coords[0,1], test_coords[0,2]],
        [20.,30., test_coords[1,0], test_coords[1,1], test_coords[1,2]],
        [45.,60., test_coords[2,0], test_coords[2,1], test_coords[2,2]]
    ])

    dom = inp_data[:,:2]
    no_frames_each = dom[:,1] - dom[:,0] + 1
    frames = [np.arange(no_frames) + start for no_frames, start in zip(no_frames_each, dom[:,0])]
    centre_bundles = [gen_centre_bundle(frame_set, central_coord) for frame_set, central_coord in zip(frames, inp_data[:,2:5])]
    path_functions = [OrbitalPath(centre_bundle, *args) for centre_bundle, args in zip(centre_bundles, inp_data[:, 5:])]
    
    dom_path_pair = np.c_[dom, path_functions]

    path = CombinedPath(dom_path_pair)
    frames = np.arange(60, dtype="f8")

    look_pos = gen_look_bundle(targ_data,60)
    print look_pos
    weights = np.asarray([1] * 10 + list(np.linspace(1,0,10)) + [1] * 10 + [1] * 10 + list(np.linspace(1,0,10)) + [1] * 10 )
    #print look_pos, "\n----------\n", weights, "\n-----------\n"

    print "computing path coords -------"
    path_coords = path(frames)

    print "computing cam basis -------"
    basis_3 = look_at_vectors(path_coords, look_pos, weights)
    tangents = vector_derivs(frames, path, d_frame=2.)
    basis_1 = orthonormalise(tangents, basis_3)
    basis_2 = cross_basis(basis_3, basis_1)

    #Plotting bits
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.plot(path_coords[:, 0], path_coords[:, 1], path_coords[:, 2])
    ax.scatter(path_coords[:, 0], path_coords[:, 1], path_coords[:, 2])
    ax.scatter(test_coords[:, 0], test_coords[:, 1], test_coords[:, 2])
    ax.quiver(path_coords[:, 0], path_coords[:, 1], path_coords[:, 2],
              basis_1[:, 0], basis_1[:, 1], basis_1[:, 2], pivot="tail", color="#FF0000")
    ax.quiver(path_coords[:, 0], path_coords[:, 1] ,path_coords[:, 2],
              basis_2[:, 0], basis_2[:, 1], basis_2[:, 2], pivot="tail", color="#00FF00")
    ax.quiver(path_coords[:, 0],path_coords[:, 1],path_coords[:, 2],
              basis_3[:, 0], basis_3[:, 1], basis_3[:, 2], pivot="tail", color="#0000FF")
    for x,y,z,f in np.c_[path_coords, frames]:
        ax.text(x,y,z,f)
    plt.show()
    #make file
    # sfs = utils.get_scalefactors(0.5,0.6,no_frames)
    # utils.gen_flight_file(frames, sfs, path_coords, np.asarray([basis_1, basis_2,basis_3]), "Paths\combined_path.txt")


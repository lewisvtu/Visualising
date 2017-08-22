from __future__ import division

import matplotlib.pyplot as plt
import numpy as np

import utils
from utils import coord_transform, cross_basis, orthonormalise, Spline3D, Interp3D


class SplinePath(object):
    '''Spline paths that stitch together general paths'''
    def __init__(self, knowns, k=3):
        '''Args:
            knowns [array]: numpy array of containing, for each coordinate on the known line,
                [f,x,y,z] for frame number f, and coord [x,y,z]
            k [real < 5]: smoothing factor defined by scipy'''
        self.spline = Spline3D(knowns, k)

    def __call__(self, frame_set, mdf=False):
        '''generate the coordinates for given frame numbers
        Args:
            frame_set [array/real]: frame number/s to calculate the new coords for
            mdf (optional) [bool]: if the frames given are multidimentional, take only the
                first column. Blame numpy.piecewise for requiring the same array sizes for
                input and output >:|
        Returns:
            coords [array]: numpy array of coords [x,y,z] for new frame numbers'''
        if mdf:
            frame_set = frame_set[:,0]
        coords = self.spline(frame_set)
        # print "spline path call"
        # print frame_set, "\n ----------\n", coords
        return coords

class OrbitalPath(object):
    '''Generates paths central to a set of points, typically orbits, spirals of helical paths'''
    def __init__(self, centre_bundle, nx,ny,nz, rad_vel, rpf,
                 rad_off, rev_off, helix_vel, helix_off):
        '''sets up the arguments for the orbital path
        Args:
            Centre_bundle [array]: array of frames and the centre positions at those points, used to interpolate
                the galaxy positions when moving between snapshots
            nx,ny,nz [reals]: the orbital axis vector in world coord basis, orients the orbital path
            rad_vel [real]: radial velocity, 0 for circles, negative spiral inwards, positive spiral out, in units
                cMpc/h per frame
            rpf [real]: angular velocity, 0 for lines, any other valuer results in circular path, in units revolution per frame
            rad_off [real]: starting radius for orbit in units cMpc/h
            rev_off [real]: starting orbital angle, used to shift the patyh around an orbit, in units revolution
            helix_vel [real]: velocity in the direction of orbital axis, used to create helical paths. in units cMpc/h per frame
            helix_off [real]: displacement from the centre in the orbital axis direction, in units cMpc/h'''
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
        self.init_frame = centre_bundle[0][0]


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
        frame_set = np.asarray(frame_set) - self.init_frame
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
        if np.linalg.norm(derivs[index]) != 0:
            pass
        elif index == 0:
            derivs[index] = np.asarray([0,0,1])
        else:
            derivs[index] = derivs[index - 1]
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
    #print no_frames
    look_bundle = np.zeros((no_frames, 7))
    #print look_bundle
    tmp = np.asarray([[no_frames,no_frames,0.,0.,0.]])
    #print t_data, "\n---------\n", tmp
    t_data = np.r_[t_data, tmp]
    for index in range(len(t_data) - 1):
        cds, cde, ctgx, ctgy, ctgz = t_data[index]
        cds, cde = int(cds), int(cde)
        ctg = np.asarray([ctgx,ctgy,ctgz])
        nds, nde, ntgx, ntgy, ntgz = t_data[index+1]
        nds, nde = int(nds), int(nde)
        ntg = np.asarray([ntgx,ntgy,ntgz])
        #print ctg,ntg
        look_bundle[cds:nds, :6] = np.r_[ctg,ntg]
        look_bundle[cds:cde, 6] = 1
        #print cde, nds
        look_bundle[cde:nds, 6] = np.logspace(0,-2,nds-cde)
        #print look_bundle[:, 6]
    return look_bundle

def create_flight_path(inp_data, mult_h):
    #Gen old style target array
    h= 0.6777
    if mult_h:
        inp_data[:,[4,5,6,10,12,14,15]] = inp_data[:,[4,5,6,10,12,14,15]] * h
    targ_data = np.copy(inp_data[:, [0,1,4,5,6]])
    targ_data[0,0] = targ_data[0,0] + 1
    targ_data[-1, 1] = targ_data[-1,1] - 1


    no_frames = int(inp_data[-1,1] - 1)
    print "No of frames: %s" % no_frames

    #Split input data into chunks
    dom = inp_data[:,:2]
    no_frames_each = dom[:,1] - dom[:,0] + 1
    #Get frames for each galaxy target
    targ_frames = [np.arange(no_of_frames) + start for no_of_frames, start in zip(no_frames_each, dom[:,0])]
    #Specify the centre of each orbital path
    centre_bundles = [gen_centre_bundle(frame_set, central_coord) for frame_set, central_coord in zip(targ_frames, inp_data[:,4:7])]
    #gen the path functions for each sub path and match to their resp frame domains
    path_functions = [OrbitalPath(centre_bundle, *args) for centre_bundle, args in zip(centre_bundles, inp_data[:, 7:])]
    dom_path_pair = np.c_[dom, path_functions]
    #Gen a combined, piecewise path for the motion
    path = CombinedPath(dom_path_pair)

    frames = np.arange(no_frames, dtype="f8")
    look_pos = gen_look_bundle(targ_data,no_frames)
    weights = look_pos[:,-1]
    look_pos = look_pos[:,:-1]
    #print look_pos, "\n----------\n", weights, "\n-----------\n"
    targ_sfs = [np.linspace(ssf, esf, no_of_frames) for ssf, esf, no_of_frames in zip(inp_data[:,2], inp_data[:,3], no_frames_each)]
    targ_sfs = np.concatenate(targ_sfs).ravel()
    targ_frames = np.concatenate(targ_frames).ravel()
    sfs = utils.get_scalefactors(targ_sfs, targ_frames, frames)

    print "computing path coords -------"
    path_coords = path(frames)

    print "computing cam basis -------"
    basis_3 = look_at_vectors(path_coords, look_pos, weights)
    tangents = vector_derivs(frames, path, d_frame=2.)
    basis_1 = orthonormalise(tangents, basis_3)
    basis_2 = cross_basis(basis_3, basis_1)
    utils.gen_flight_file(frames, sfs, path_coords, np.asarray([basis_1, basis_2,basis_3]), "galaxy_tour2_.txt")
    return True


if __name__ == "__main__":
    print "Actually started running -_- z z z"

    inp_data = np.asarray([
    #   [domain    , sf,  coords at centre of montion  ,rotaxis, rv,   av,ro,ao,hv,ho]
        [-1.,  0., .25, .25, 7., 11.,  10.22,  0,0,1,  0, 0, 0, 0, 0, 0],
        [100.,   300., .29, .32, 12.2787, 19.071,  17.22,  0,0,1,  0, 1/200, 1, 0, 0, 0],
        [400.,  600.,  .475, .525,  8.434,  9.601,   4.25,  1,1,1,  0, 1/200, 2, 0, 0, 0],
        [700., 900.,  1., 1., 16.557,  24.49, 17.708, -1,1,0,  0, 1/200, 3, 0, 0, 0],
        [999.,   1000., 1., 1., 7., 11.,  10.22,  0,0,1,  0, 0, 0, 0, 0, 0]

    ])


    create_flight_path(inp_data, True)

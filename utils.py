import matplotlib.pyplot as plt
import numpy as np
import pickler.shelf as shelf
from DBS.dbgrabber import dbsPull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import pi  
from scipy.interpolate import UnivariateSpline
def get_scalefactors(start_sf, end_sf, frames):
    '''
    returns list of #frames scale factors scaled uniformly in log10 between start_sf and end_sf
    '''
    array_log_sf = np.linspace(np.log10(start_sf), np.log10(end_sf), frames)
    array_sf = np.power(10, array_log_sf)
    return array_sf

def gen_flight_file(frames, sfs, coords, basis_vects, fname):
    setspace = np.asarray([frames, sfs, coords[:,0]       , coords[:,1]       , coords[:,2],
										basis_vects[0,:,0], basis_vects[0,:,1], basis_vects[0,:,2],
										basis_vects[1,:,0], basis_vects[1,:,1], basis_vects[1,:,2],
										basis_vects[2,:,0], basis_vects[2,:,1], basis_vects[2,:,2]])
    #print setspace
    np.savetxt(fname, setspace.T, fmt="%i %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f", header="RefL0100N1504" )


class Spline3D:
    '''
    class for 3d splines.

    '''
    def __init__(self, bundle, k=3):
        '''
        Creates the Spline object
        '''
        fks, xks, yks, zks = np.transpose(bundle)
        self.x_spline = UnivariateSpline(fks, xks, k=k)
        self.y_spline = UnivariateSpline(fks, yks, k=k)
        self.z_spline = UnivariateSpline(fks, zks, k=k)

    def __call__(self, fs):
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


class Data():
    def __init__(self, sql_query, f_name):
        self.dbs_data = dbsPull(sql_query, f_name)
    
    def galaxies_from_ids(self, interesting_ids):
        '''
        gives useful galaxy data for interesting ids and snapshots, in a useful form

        Args:
            interesting_ids [dictionary]: Maps keys of galaxyID to values of snapshots
        
        Returns:
            interesting_gals [array]: [sf,x,y,z] for each galaxy in interesting_gals
        '''
        interesting_gals = np.asarray([list(gal)[3:] for gal in dbs_data if gal[0] in interesting_ids.keys() and gal[1] == interesting_ids[gal[0]]])
        interesting_gals[:,-1] = 1.0 / (1.0 + interesting_gals[:,-1])
        gals = gals[np.argsort(gals[:,3])]
        gals = gals[:,[3,0,1,2]]
        return interesting_gals


def orthonormalise(new_vects, basis_1):
    basis_2 = new_vects - (np.einsum("ij,ij->i", new_vects, basis_1) * basis_1.T).T
    basis_2 = basis_2 / np.linalg.norm(basis_2, axis=1)[:,None]
    return basis_2

def cross_basis(basis_1, basis_2):
    basis_3 = np.cross(basis_1, basis_2)
    basis_3 = basis_3 / np.linalg.norm(basis_3, axis=1)[:,None]
    return basis_3    

def coord_transform(x_basis, y_basis, z_basis, cam_position, particles, inv=True, homog=True, tran=False):
	if tran:
		particles = particles.T
	coords_none_trans = np.transpose(np.c_[particles, np.ones(len(particles))])
	M_world_camera = np.array([

		(x_basis[0], y_basis[0], z_basis[0], cam_position[0]),
		(x_basis[1], y_basis[1], z_basis[1], cam_position[1]),
		(x_basis[2], y_basis[2], z_basis[2], cam_position[2]),
		(0, 0, 0, 1)

		])
	if inv:
		M_world_camera = np.linalg.inv(M_world_camera)
	coords_in_cam = np.dot(M_world_camera, coords_none_trans)
	if not homog:
		coords_in_cam = coords_in_cam[:-1,:]
	if tran:
		coords_in_cam = coords_in_cam.T
	return coords_in_cam

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def perspective_transfomation(coords_in_cam, region):


	#coords_in_cam = coord_transform(x_basis, y_basis, z_basis, cam_position, particles)

	fov = np.pi / 4
	d               = 1./(np.tan(fov/2.))
	aspect_ratio    = np.true_divide(region[0], region[1])
	near            = region[2] * 0.0001
	far             = region[2] * 1.0000

	M_projection = np.array([
		(d*(1./aspect_ratio),0,0,0),
		(0,d,0,0),
		(0,0,(-near-far)/(near-far), (2*near*far)/(near-far)),
		(0,0,1,0)
	])



	perpec_in_cam = np.dot(M_projection, coords_in_cam)
	coords = perpec_in_cam
	coords      = coords.T
	coords[:,0] = coords[:,0]/coords[:,3]
	coords[:,1] = coords[:,1]/coords[:,3]
	coords[:,2] = coords[:,2]/coords[:,3]
	coords[:,3] = coords[:,3]/coords[:,3]

	#clips all galaxies that are not in your field of view
	coords_alt = np.asarray([np.append(coord, index) for index, coord  in enumerate(coords) if (abs(coord[0]) <= 1) and (abs(coord[1]) <= 1) and abs(coord[2]) <= 2

		])
	#print coords_alt 

	# w_s         = region[0]
	# h_s         = region[1]
	# s_x         = 0.0
	# s_y         = 0.0

	# M_viewport  = np.array([(w_s/2.,0,0,s_x+w_s/2.),
	# 						(0,h_s/2.,0,s_y+h_s/2.),
	# 						(0,0,(far-near)/2.,(near+far)/2.),
	# 						(0,0,0,0)])
	# coords      = np.dot(M_viewport, coords_alt.T)



	return coords_alt


def find_snapnums(scale_factor):

	Expansion_F_snaps = np.array([0.05, 0.06, 0.09, 0.10, 0.11, 0.12, 0.14, 0.15, 0.17,
	                     0.18, 0.20, 0.22, 0.25, 0.29, 0.31, 0.33, 0.37, 0.40,
	                     0.44, 0.50, 0.54, 0.58, 0.62, 0.67, 0.73, 0.79,0.85,
	                     0.91, 1.00])

	nearest_snapnumber = (np.abs(Expansion_F_snaps - scale_factor)).argmin()
	if Expansion_F_snaps[nearest_snapnumber] == scale_factor:
		beforeSnap = nearest_snapnumber
		afterSnap = nearest_snapnumber

	elif Expansion_F_snaps[nearest_snapnumber] > scale_factor:
		beforeSnap = nearest_snapnumber-1
		afterSnap = nearest_snapnumber

	elif Expansion_F_snaps[nearest_snapnumber] < scale_factor:
		beforeSnap = nearest_snapnumber
		afterSnap = nearest_snapnumber+1

	return [beforeSnap, afterSnap]	


def galaxy_interpolation(scale_factor,dbs_data):

	all_raw_data = np.asarray(dbs_data)
	sideSnaps = find_snapnums(scale_factor)
	beforeSnap, afterSnap = sideSnaps[0], sideSnaps[1]




# xyz = perspective_transfomation(x_basis, y_basis, z_basis, cam_position, particles)


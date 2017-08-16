import matplotlib.pyplot as plt
import numpy as np
import pickler.shelf as shelf
from DBS.dbgrabber import dbsPull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import pi  
from scipy.interpolate import UnivariateSpline, interp1d, spline

Expansion_F_snaps = np.array([0.05, 0.06, 0.09, 0.10, 0.11, 0.12, 0.14, 0.15, 0.17,
					 0.18, 0.20, 0.22, 0.25, 0.29, 0.31, 0.33, 0.37, 0.40,
					 0.44, 0.50, 0.54, 0.58, 0.62, 0.67, 0.73, 0.79,0.85,
					 0.91, 1.00])
def plot_from_file(f_name):
	fs,sfs,xs,ys,zs,v1xs,v1ys,v1zs,v2xs,v2ys,v2zs,v3xs,v3ys,v3zs = np.loadtxt(f_name, unpack=True)
	#Plotting bits
	fig = plt.figure()
	ax = fig.add_subplot(111, projection="3d")
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.set_zlabel("z")
	path_coords = np.asarray([xs,ys,zs]).T
	basis_1 = np.asarray([v1xs,v1ys,v1zs]).T
	basis_2 = np.asarray([v2xs,v2ys,v2zs]).T
	basis_3 = np.asarray([v3xs,v3ys,v3zs]).T
	ax.plot(path_coords[:, 0], path_coords[:, 1], path_coords[:, 2])
	ax.scatter(path_coords[:, 0], path_coords[:, 1], path_coords[:, 2])
	ax.quiver(path_coords[:, 0], path_coords[:, 1], path_coords[:, 2],
				basis_1[:, 0], basis_1[:, 1], basis_1[:, 2], pivot="tail", color="#FF0000")
	ax.quiver(path_coords[:, 0], path_coords[:, 1] ,path_coords[:, 2],
				basis_2[:, 0], basis_2[:, 1], basis_2[:, 2], pivot="tail", color="#00FF00")
	ax.quiver(path_coords[:, 0],path_coords[:, 1],path_coords[:, 2],
				basis_3[:, 0], basis_3[:, 1], basis_3[:, 2], pivot="tail", color="#0000FF")
	for x,y,z,f in np.c_[path_coords, fs]:
		ax.text(x,y,z,f)
	plt.show()


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

class Interp3D(object):
	def __init__(self, bundle):
		fks, xks, yks, zks = np.transpose(bundle)
		self.x_interp = interp1d(fks, xks, fill_value="extrapolate")
		self.y_interp = interp1d(fks, yks, fill_value="extrapolate")
		self.z_interp = interp1d(fks, zks, fill_value="extrapolate")

	def __call__(self, frames):
		xs = self.x_interp(frames)
		ys = self.y_interp(frames)
		zs = self.z_interp(frames)
		return np.asarray([xs, ys, zs]).T

class Spline3D:
    '''
    class for 3d splines.

    '''
    def __init__(self, bundle, k=3):
		'''
		Creates the Spline object
		'''
		fks, xks, yks, zks = bundle.T
		self.x_spline = UnivariateSpline(fks, xks, k=k)
		self.y_spline = UnivariateSpline(fks, yks, k=k)
		self.z_spline = UnivariateSpline(fks, zks, k=k)
		# self.fks = fks
		# self.xks = xks
		# self.yks = yks
		# self.zks = zks

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
		# xs = spline(self.fks, self.xks, fs)
		# ys = spline(self.fks, self.yks, fs)
		# zs = spline(self.fks, self.zks, fs)
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

	#to clip aything out of the region on z 



	perpec_in_cam = np.dot(M_projection, coords_in_cam)
	coords = perpec_in_cam
	coords      = coords.T
	coords[:,0] = coords[:,0]/coords[:,3]
	coords[:,1] = coords[:,1]/coords[:,3]
	coords[:,2] = coords[:,2]/coords[:,3]
	coords[:,3] = coords[:,3]/coords[:,3]

	#clips all galaxies that are not in your field of view
	coords_alt = np.asarray([np.append(coord, index) for index, coord  in enumerate(coords) if (abs(coord[0]) <= 1) and (abs(coord[1]) <= 1) and abs(coord[2]) <= 1.
		])
	# mask = np.where(np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.logical_and(
	# 		coords[:,0] <= 1, coords[:,0] >=-1),
	# 		coords[:,1]<=1), coords[:,1]>=-1),
	# 		coords[:,2]<=1), coords[:,2]>=-1))
	# coords_alt = coords[mask]
	# print coords_alt.T 
	# print coords_alt2
	# w_s         = region[0]
	# h_s         = region[1]
	# s_x         = 0.0
	# s_y         = 0.0

	# M_viewport  = np.array([(w_s/2.,0,0,s_x+w_s/2.),
	# 						(0,h_s/2.,0,s_y+h_s/2.),
	# 						(0,0,(far-near)/2.,(near+far)/2.),
	# 						(0,0,0,0)])
	# coords      = np.dot(M_viewport, coords_alt.T)
	# print coords


	#return M_viewport
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


def orderGals_all(dbs_data, snapshot_num):

	snaps = {}
	for galID in gals.keys()[:]:
		gal = gals[galID]

		for galsnap in gal:
			if galsnap[1] not in snaps.keys():
				snaps[galsnap[1]] = [galsnap]
			else: 
				snaps[galsnap[1]].append(galsnap)

	return snaps



def galsTree(dbs_data):
	#creates a dictionary of all the galaxies and all there snapshots 
	gals = {}
	for ele in dbs_data:
	  galID = ele[0]
	  if galID in gals.keys():
		  gals[galID].append([ele[i] for i in range(0, len(ele))])
	  else:
		  gals[galID] = [[ele[i] for i in range(0, len(ele))]]

	return gals  



def gal_interpolation(scale_factor, dbs_data):

	sideSnaps = find_snapnums(scale_factor)
	beforeSnap, afterSnap = sideSnaps[0], sideSnaps[1]

	#galaxies of the before snapshots
	# mask_B = np.where(dbs_data['SnapNum'] == beforeSnap) 
	mask_B = np.where(np.logical_and(dbs_data['SnapNum'] == beforeSnap, dbs_data['MassType_DM'] >= 1e11))
	beforeGals = dbs_data[mask_B] 

	mask_After = np.where(np.in1d( dbs_data['ID'], beforeGals['DesID']))
	afterGals = dbs_data[mask_After]

	if beforeSnap == afterSnap:
		interpGals = np.asarray([ beforeGals['ID'],beforeGals['SnapNum'],beforeGals['MassType_DM'],(beforeGals['x']),
								(beforeGals['y']),(beforeGals['z']),beforeGals['Redshift']]).T
	else:
		delta_x = afterGals['x'] - beforeGals['x']
		delta_y = afterGals['y'] - beforeGals['y']
		delta_z = afterGals['z'] - beforeGals['z']
		delta_coords = np.array([delta_x,delta_y,delta_z]).T
		delta_time =  scale_factor - Expansion_F_snaps[beforeSnap]
		fracTime = delta_time / (Expansion_F_snaps[afterSnap] - Expansion_F_snaps[beforeSnap])
		fracCoords = delta_coords * fracTime

		interpGals = np.asarray([ beforeGals['ID'],beforeGals['SnapNum'],beforeGals['MassType_DM'],(beforeGals['x']+fracCoords[:,0]),
								(beforeGals['y']+fracCoords[:,1]),(beforeGals['z']+fracCoords[:,2]),beforeGals['Redshift']]).T

	return interpGals


def galaxy_interpolation(scale_factor, dbs_data, snaps):
	#gals is all snaps dictionary

	all_raw_data = np.asarray(dbs_data)
	sideSnaps = find_snapnums(scale_factor)
	beforeSnap, afterSnap = sideSnaps[0], sideSnaps[1]

	beforeGals = np.asarray(snaps[beforeSnap])
	afterGals = np.asarray(snaps[afterSnap])


	if beforeSnap == afterSnap:
		return beforeGals

	else:	

		#romoves any galaxies that aren't in both snapshots 
		afterGalsS = np.asarray([after_gal for after_gal in afterGals if after_gal[0] in beforeGals[:,0]])
		beforeGalsS = np.asarray([before_gal for before_gal in beforeGals if before_gal[0] in afterGals[:,0]])

		#sort them by the id's so that they are in the same order
		afterGalsS = afterGalsS[afterGalsS[:,0].argsort()]
		beforeGalsS = beforeGalsS[beforeGalsS[:,0].argsort()]
		#print beforeGalsS
		delta_coords = afterGalsS[:,3:6] - beforeGalsS[:,3:6]
		delta_time =  scale_factor - Expansion_F_snaps[beforeSnap]

		fracTime = delta_time / (Expansion_F_snaps[afterSnap] - Expansion_F_snaps[beforeSnap])
		fracCoords = delta_coords * fracTime

		beforeGalsS[:,3:6] = beforeGalsS[:,3:6] + fracCoords

		InterpGals = beforeGalsS

		return InterpGals 



def get_centre(basis_vectors, cam_position, region):
	""" Return centre as seen from the camera. """
	bv = basis_vectors
	M = np.array([(bv[0][0], bv[1][0], bv[2][0], cam_position[0]),
				  (bv[0][1], bv[1][1], bv[2][1], cam_position[1]),
				  (bv[0][2], bv[1][2], bv[2][2], cam_position[2]),
				  (0,0,0,1)])
	return np.dot(M, np.array([0,0,region[2]/2.,1]))[:-1]


def periodic_wrap(pos, boxsize, centre=None):
	"""
	Wrap the coordinates in pos to the periodic copy nearest centre
	assuming box size boxsize. Coordinates may be more than one
	boxsize away from centre.

	Method is to shift coords such that centre is at 0.5*boxsize,
	use numpy.mod to wrap coordinates into range 0-boxsize, then shift
	back.

	If centre is None just wrap into range 0-boxsize.
	"""
	if centre is None:
		return np.mod(pos, boxsize)
	else:
		return np.mod(pos-centre+0.5*boxsize, boxsize)+centre-0.5*boxsize

# xyz = perspective_transfomation(x_basis, y_basis, z_basis, cam_position, particles)


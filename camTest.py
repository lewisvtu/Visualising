#rough test 
import matplotlib.pyplot as plt
import numpy as np
import pickler.shelf as shelf
import math
from DBS.dbgrabber import dbsPull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import pi  


def coord_transform(x_basis, y_basis, z_basis, cam_position, particles, inv=True, homog=True):

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
		return coords_in_cam[:-1,:]
	return coords_in_cam



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
	coords_alt = np.asarray([np.append(coord, index) for index, coord  in enumerate(coords) if (abs(coord[0]) <= 1) and (abs(coord[1]) <= 1) and abs(coord[2]) <= 1])
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



# xyz = perspective_transfomation(x_basis, y_basis, z_basis, cam_position, particles)



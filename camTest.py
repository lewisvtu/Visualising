#rough test 
import flightplan as fp
import matplotlib.pyplot as plt
import numpy as np
import pickler.shelf as shelf
import math
from DBS.dbgrabber import dbsPull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import pi  

def eulerAnglesToRotationMatrix(theta):


	c1, s1 = np.cos(theta[0]), np.sin(theta[0])
	c2, s2 = np.cos(theta[1]), np.sin(theta[1])
	c3, s3 = np.cos(theta[2]), np.sin(theta[2])

	R_x = np.matrix('{} {} {}; {} {} {}; {} {} {}'.format(1.0, 0.0, 0.0, 0.0, c1, -s1, 0.0, s1, c1))
	R_y = np.matrix('{} {} {}; {} {} {}; {} {} {}'.format(c2, 0.0, s2, 0.0, 1.0, 0.0, -s2, 0.0, c2))
	R_z = np.matrix('{} {} {}; {} {} {}; {} {} {}'.format(c3, -s3, 0.0, s3, c3, 0.0, 0.0, 0.0, 1.))

	R = R_z * R_y * R_x
	return R
 

frame, ts, xs, ys, zs, b1,b2,b3,b4,b5,b6,b7,b8,b9 = np.loadtxt("cam_test_9994243.txt", unpack=True)
thetas = np.asarray([[pi/4,pi/4,pi/4]] * len(xs))
#print thetas



z_basis = np.transpose(np.array([b7, b8, b9]))[0]
#print z_basis
y_basis = np.transpose(np.array([b4, b5, b6]))[0]
x_basis = np.transpose(np.array([b1, b2, b3]))[0]
cam_position = np.transpose(np.array([xs, ys, zs]))[0]


# print z_basis
# #z is the look direction therefore we can get the yaw and pitch from them 
# pitch = np.arcsin(z_basis[:,1])
# yaw = np.arctan2(z_basis[:,0], z_basis[:,2])
# roll = [0.0]*len(z_basis)


# print pitch
# print yaw

Rotation_Marticies = [eulerAnglesToRotationMatrix(theta) for theta in thetas]





# thetas = np.concatenate((roll, pitch, yaw), axis=1)
# #to get the roll we need to look at the 

# #to now transform the objects positions to the camera frame of refernce 

# cam_frame_coords = eulerAnglesToRotationMatrix(thetas) * new_coords
# print cam_frame_coords



# gals = fp.gals

# no_of_frames = 20
# frame_array = np.arange(no_of_frames)
# circle_args = [gals[0], 5.0, 1.0, no_of_frames]
# xs, ys, zs = np.transpose(fp.circular_path(frame_array, circle_args))
# v1xs, v1ys, v1zs, v2xs, v2ys, v2zs, v3xs, v3ys, v3zs = np.transpose(fp.cam_vectors(frame_array, gals[0], fp.circular_path, circle_args))


# fig = plt.figure()
# ax = fig.add_subplot(111, projection="3d")
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("z")

# z_basis2 = np.asarray([np.dot(Rotation_Martix, z_b) for Rotation_Martix, z_b in zip(Rotation_Marticies, z_basis)])
# print z_basis2

# vectTest = np.asarray([0., 0., 1.0])
# vect45 = np.dot(Rotation_Marticies[0], vectTest)
# vect45_45 = np.dot(Rotation_Marticies[1], vectTest)
# vect45_45_45 = np.dot(Rotation_Marticies[2], vectTest)

# ax.quiver(xs, ys, zs ,z_basis[:,0], z_basis[:,1], z_basis[:,2], color="r", pivot="tail")
# ax.quiver(xs, ys, zs ,z_basis2[:,0,0], z_basis2[:,0,1], z_basis2[:,0,2], color="r", pivot="tail")

# # ax.quiver(0.,0.,0.,vectTest[0], vectTest[1], vectTest[2], color="#682860", pivot="tail")
# # ax.quiver(0.,0.,0.,vect45[:,0], vect45[:,1], vect45[:,2], color="k", pivot="tail")
# # ax.quiver(0.,0.,0.,vect45_45[:,0], vect45_45[:,1], vect45_45[:,2], color="c", pivot="tail")
# # ax.quiver(0.,0.,0.,vect45_45_45[:,0], vect45_45_45[:,1], vect45_45_45[:,2], color="c", pivot="tail")
# plt.show()


x, y, z = [70.], [70.], [10.]

particles = np.transpose(np.array([x, y, z]))
#print particles

def perspective_transfomation(x_basis, y_basis, z_basis, cam_position, particles):
	coords_none_trans = np.transpose(np.c_[particles, np.ones(len(particles))])
	M_world_camera = np.linalg.inv(np.array([

		(x_basis[0], y_basis[0], z_basis[0], cam_position[0]),
		(x_basis[1], y_basis[1], z_basis[1], cam_position[1]),
		(x_basis[2], y_basis[2], z_basis[2], cam_position[2]),
		(0, 0, 0, 1)

		]))
	print "these are the coordiantes of the non tranformed "
	print coords_none_trans
	coords_in_cam = np.dot(M_world_camera, coords_none_trans)
	return coords_in_cam 


xyz = perspective_transfomation(x_basis, y_basis, z_basis, cam_position, particles)

print "transformed coordinates of the stufffffffffffffffffffffffffffffffffffffffff"
print xyz


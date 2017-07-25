#!/bin/env python2.7
import matplotlib.pyplot as plt
import numpy as np
import pickler.shelf as shelf
import math
from DBS.dbgrabber import dbsPull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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

dbs_data = shelf.pull("followup2")

interesting_ids = {
    13660659: 20,
    13660659: 21,
    13660659: 22,
    13660659: 23,
    13660659: 28

}
gals = np.asarray([list(gal)[3:] for gal in dbs_data if gal[0] in interesting_ids.keys() and gal[1] == interesting_ids[gal[0]]])
# Makes all galaxies have the form [sf,x,y,z] sorted by sf
gals[:,-1] = 1.0 / (1.0 + gals[:,-1])
gals = gals[np.argsort(gals[:,3])]
gals = gals[:,[3,0,1,2]]


def circular_path(frame_nos, args):
    '''
    Given a frame number, returns x,y,z coords for the camera at that frame
    
    args:
        frame_nos: the frame number/ array of frame numbers to compute the positions of
        args: list of form []
    '''
    frame_nos = np.asarray(frame_nos)
    target_coords = np.asarray(args[0])
    rad = args[1]
    orbits = args[2]
    frames = args[3]
    ang_int = orbits * 2*np.pi / frames
    x_coords = target_coords[1] + rad * np.sin(frame_nos * ang_int)
    y_coords = target_coords[2] + rad * np.cos(frame_nos * ang_int)
    z_coords = target_coords[3] + rad * np.sin(frame_nos * ang_int) 
    return np.transpose(np.asarray([x_coords, y_coords, z_coords]))

def cam_vectors(frame_nos, target_coords, path_function, args):
    look_at_dirs = target_coords[1:] - path_function(frame_nos, args)
    look_at_dirs = look_at_dirs / np.linalg.norm(look_at_dirs, axis=1)[:,None]
    derivs = np.zeros((len(frame_nos), 3))
    d_frame = 0.01
    for index in range(len(frame_nos)):
        frame_no = frame_nos[index]
        derivs[index] = path_function(frame_no + d_frame/2, args) - path_function(frame_no - d_frame/2, args)
    #einsum computes the vector dots for each pair in the array
    basis_1 = derivs - np.transpose(np.einsum("ij,ij->i", derivs, look_at_dirs) * np.transpose(look_at_dirs))
    basis_1 = basis_1 / np.linalg.norm(basis_1, axis=1)[:,None]
    basis_2 = np.cross(look_at_dirs, basis_1)
    basis_2 = basis_2 / np.linalg.norm(basis_2, axis=1)[:,None]
    return np.concatenate((basis_1, basis_2, look_at_dirs), axis=1)

def straight_path(frame_nos, args):
    frames = args[2]
    start_coords = args[0]
    end_coords = args[1]
    x_coords = np.linspace(start_coords[0], end_coords[0], frames)
    y_coords = np.linspace(start_coords[1], end_coords[1], frames)
    z_coords = np.linspace(start_coords[2], end_coords[2], frames)
    return np.transpose(np.asarray([x_coords, y_coords, z_coords]))    

def get_scalefactors(start_sf, end_sf, frames):
    array_log_sf = np.linspace(np.log10(start_logsf), np.log10(end_logsf), frames)
    array_sf = np.power(10, array_log_sf)
    return array_sf[::-1]

no_of_frames = 20
frame_array = np.arange(no_of_frames)
circle_args = [gals[0], 5.0, 1.0, no_of_frames]
straight_args = [gals[0,1:] - [0,5,5], gals[0,1:] + [0,5,5], no_of_frames]
xs, ys, zs = np.transpose(straight_path(frame_array, straight_args))
v1xs, v1ys, v1zs, v2xs, v2ys, v2zs, v3xs, v3ys, v3zs = np.transpose(cam_vectors(frame_array, gals[0], straight_path, straight_args))


fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
# ax.plot(xs,ys,zs)
# ax.scatter(gals[:,0], gals[:,1], gals[:,2], marker="o", s=200.0, c="#682860")
# for galaxy in nearby_gal_datas:
#     ax.scatter(galaxy[0], galaxy[1], galaxy[2], marker="o", s=50.0)
ax.quiver(xs,ys,zs, v1xs, v1ys, v1zs, color="#682860", pivot="tail")
ax.quiver(xs,ys,zs, v2xs, v2ys, v2zs, color="#000000", pivot="tail")
ax.quiver(xs,ys,zs, v3xs, v3ys, v3zs, color="#FF0000", pivot="tail")
plt.show()

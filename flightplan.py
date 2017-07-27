#!/bin/env python2.7
import matplotlib.pyplot as plt
import numpy as np
import pickler.shelf as shelf
import math
from DBS.dbgrabber import dbsPull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import spline

h = 0.6777
SQL = """
    SELECT
        DES.GalaxyID,
        PROG.SnapNum,
        PROG.Mass,
        PROG.CentreOfPotential_x,
        PROG.CentreOfPotential_y,
        PROG.CentreOfPotential_z,
        PROG.Redshift
    FROM
        RefL0100N1504_Subhalo as PROG with(forceseek),
        RefL0100N1504_Subhalo as DES,
        RefL0100N1504_Aperture as AP
    WHERE
        DES.SnapNum = 28 and
        DES.MassType_Star > 1.0e9 and
        DES.MassType_DM > 5.0e10 and
        PROG.GalaxyID between DES.GalaxyID and DES.TopLeafID and
        AP.ApertureSize = 30 and
        AP.GalaxyID = DES.GalaxyID and
        AP.Mass_Star > 1.0e9
    ORDER BY
        PROG.GalaxyID,
        PROG.SnapNum
"""

# Grabs new data from db based on sql. If file name already exists, it loads that data instead

filename = "FollowProgs2.p"

raw_dbs = dbsPull(SQL, filename)

shelf.push(raw_dbs, "followup2")

dbs_data = shelf.pull("followup2")

interesting_ids = {
<<<<<<< HEAD
    # 13660659: 28,
    # 13793733: 28,
    # 13722615: 28
    14784533: 28
=======
    13660659: 28,
    13793733: 28,
    13722615: 28
>>>>>>> c89193123abcc2538bbd3ac0baf4a435f01e8a4e

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
    decay = 0.5
    frame_nos = np.asarray(frame_nos)
    target_coords = np.asarray(args[0])
    rad = args[1]
    orbits = args[2]
    frames = args[3]
    dir = args[4]
    ang_int = orbits * 2*np.pi / frames
    x_coords = target_coords[0] + dir * rad * np.sin(frame_nos * ang_int)
    y_coords = target_coords[1] + dir * rad * np.cos(frame_nos * ang_int)
    z_coords = target_coords[2] + dir * (rad* np.sin(frame_nos * ang_int)) * 0.5
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
    frame_nos = np.asarray([frame_nos])
    frames = args[2]
    start_coords = args[0]
    end_coords = args[1]
    c_coords = np.transpose(np.asarray([start_coords]*len(frame_nos))) + np.transpose(np.asarray([end_coords - start_coords]*len(frame_nos))) * (frame_nos/frames)
    return np.transpose(c_coords)

def get_scalefactors(start_sf, end_sf, frames):
    array_log_sf = np.linspace(np.log10(start_logsf), np.log10(end_logsf), frames)
    array_sf = np.power(10, array_log_sf)
    return array_sf[::-1]

def gen_file(no_of_frames, target_gal, path_function, path_args, file_name):
    frame_array = np.arange(no_of_frames, dtype=float)
    sf_array = np.asarray([1.0]*no_of_frames)
    xs, ys, zs = np.transpose(path_function(frame_array, path_args))
    v1xs, v1ys, v1zs, v2xs, v2ys, v2zs, v3xs, v3ys, v3zs = np.transpose(cam_vectors(frame_array, target_gal, path_function, path_args))
    #v1xs, v1ys, v1zs, v2xs, v2ys, v2zs, v3xs, v3ys, v3zs = [1]*no_of_frames,[0]*no_of_frames,[0]*no_of_frames,[0]*no_of_frames,[1]*no_of_frames,[0]*no_of_frames,[0]*no_of_frames,[0]*no_of_frames,[1]*no_of_frames
    linepoints = np.transpose(np.asarray([frame_array, sf_array, xs*h, ys*h, zs*h, v1xs, v1ys, v1zs, v2xs, v2ys, v2zs, v3xs, v3ys, v3zs]))
    np.savetxt(file_name, linepoints,fmt='%i %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f', header='RefL0100N1504',comments='#')
    print "Saved to file: " + file_name

def draw_graph(file_name, target_gal):
    fs, sfs, xs, ys, zs, v1xs, v1ys, v1zs, v2xs, v2ys, v2zs, v3xs, v3ys, v3zs = np.transpose(np.loadtxt(file_name))
    fig = plt.figure()
    xs, ys, zs = xs/h, ys/h, zs/h
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    # ax.plot(xs,ys,zs)
    #ax.scatter(target_gal[0], target_gal[1], target_gal[2], marker="o", s=200.0, c="#682860")
    # for galaxy in nearby_gal_datas:
    #     ax.scatter(galaxy[0], galaxy[1], galaxy[2], marker="o", s=50.0)
    ax.quiver(xs,ys,zs, v1xs, v1ys, v1zs, color="#682860", pivot="tail")
    ax.quiver(xs,ys,zs, v2xs, v2ys, v2zs, color="#000000", pivot="tail")
    ax.quiver(xs,ys,zs, v3xs, v3ys, v3zs, color="#FF0000", pivot="tail")
    plt.show()


first_frames = np.arange(10, dtype=float)
sec_frames = np.arange(10, dtype=float) + 20
third_frames = np.arange(10, dtype=float) + 40
first_coords = circular_path(first_frames, [gals[0,1:], 5.0, 1, 10, -1])
sec_coords = circular_path(sec_frames, [gals[1,1:], 5.0, 1, 10, 1])
third_coords = circular_path(third_frames, [gals[2,1:], 5.0, 1, 10, 1])

fks = np.concatenate((first_frames, sec_frames, third_frames))
xks = np.concatenate((first_coords[:,0], sec_coords[:,0], third_coords[:,0]))
yks = np.concatenate((first_coords[:,1], sec_coords[:,1], third_coords[:,1]))
zks = np.concatenate((first_coords[:,2], sec_coords[:,2], third_coords[:,2]))

frame_array = np.arange(50)
xs = spline(fks, xks, frame_array)
ys = spline(fks, yks, frame_array)
zs = spline(fks, zks, frame_array)
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.plot(xs, ys, zs)

plt.show()


# no_of_frames = 50
# circle_args = [gals[0,1:], 5.0, 1.0, no_of_frames]
# straight_args = [gals[0,1:] + [5,0,-5], gals[0,1:] + [5,0,5], no_of_frames]
# file_name = "C_path_r=5_o=1_f=50.txt"
# gen_file(no_of_frames, gals[0], circular_path, circle_args, file_name)
# draw_graph(file_name, gals[0,1:])

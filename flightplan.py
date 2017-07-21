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


part_snap_galaxies = [list(gal) for gal in dbs_data if gal[1] == 19]
#Grabs all snapshot data for a particular galaxy
part_gal = [list(gal) for gal in dbs_data if gal[0] == 9994243][::-1]
#Cuts the data down to just x,y,z,rs
imp_data = [gal[3:] for gal in part_gal]
interesting_ids = {
    9994243: 19,
   10777540: 18, 
   10351144: 20
}
gals = np.asarray([list(gal)[3:] for gal in dbs_data if gal[0] in interesting_ids.keys() and gal[1] == interesting_ids[gal[0]]])
gals = gals[np.argsort(gals[:,3])]
def path(frame_no):
    '''
    Given a frame number, returns an x,y,z and sf coord for the camera at that frame
    '''
    rad = 5
    orbits = 0.5
    frames = 100
    target_galaxy = gals[0]
    ang_int = orbits * 2*np.pi / frames
    x_coord = target_galaxy[0] + np.sin(frame_no * ang_int)
    y_coord = target_galaxy[1] + np.sin(frame_no * ang_int)
    z_coord = target_galaxy[2] + np.sin(frame_no * ang_int)
    rs_coord = target_galaxy[3]
    return numpy.asarray([x_coord, y_coord, z_coord, rs_coord])


def nearby(centre_gal, other_gal):
    '''
    takes an important galaxy and another, less important galaxy and returns if it is considered close, within 5mpc
    inputs:
        centre_gal: important galaxy we are centering on
        other_gal: another galaxy that may be considered close
    returns:
        True/ False: if the other galaxy is within/ not within 5mpc
    '''
    if abs(centre_gal[0] - other_gal[0]) < 5.0 and abs(centre_gal[1] - other_gal[1]) < 5.0 and abs(centre_gal[2] - other_gal[2]) < 5.0:
    
        return True
    else:
        return False




nearby_gal_datas = [galaxy[3:] for galaxy in part_snap_galaxies if nearby(gal, galaxy[3:])]


def get_scalefactors(start_rs, end_rs, frames):
    start_sf = 1/(1.0 + start_rs)
    end_sf = 1/(1.0 + end_rs)
    start_logsf = np.log10(start_sf)
    end_logsf = np.log10(end_sf)
    array_log_sf = np.linspace(start_logsf, end_logsf, frames)
    array_sf = np.power(10, array_log_sf)
    return array_sf[::-1]






# np.savetxt("curvedmanygals.txt", coords, fmt="%i %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f", header="RefL0100N1504")

# fns, ts, xs, ys, zs, v1xs, v1ys, v1zs, v2xs, v2ys, v2zs, v3xs, v3ys, v3zs = np.transpose(coords) 


# fig = plt.figure()
# ax = fig.add_subplot(111, projection="3d")

# ax.plot(xs,ys,zs)
# ax.scatter(gals[:,0], gals[:,1], gals[:,2], marker="o", s=200.0, c="#682860")
# for galaxy in nearby_gal_datas:
#     ax.scatter(galaxy[0], galaxy[1], galaxy[2], marker="o", s=50.0)
# ax.quiver(xs,ys,zs, v1xs, v1ys, v1zs, color="#682860", pivot="tail")
# ax.quiver(xs,ys,zs, v2xs, v2ys, v2zs, color="#000000", pivot="tail")
# ax.quiver(xs,ys,zs, v3xs, v3ys, v3zs, color="#FF0000", pivot="tail")
# plt.show()
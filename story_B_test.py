import matplotlib.pyplot as plt
import numpy as np 
import pickler.shelf as shelf
import math
from DBS.dbgrabber import dbsPull
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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

filename = "FollowProgs19.p"

raw_dbs = dbsPull(SQL, filename)

shelf.push(raw_dbs, "followup19")

dbs_data = shelf.pull("followup19")

# r = 5
# ts = 0.25

# timeframes = 100

# def gen_surrounding_points(ip, rad, time_space):
#     offsets = [
#         np.array([rad,0,0,-time_space]),
#         np.array([rad,rad,0,0]),
#         np.array([0,rad,0,time_space])
#     ]

#     points = []
#     for ofs in offsets:
#         newpoint = ip + ofs
#         points.append(list(newpoint))
#     return points

# ips = [
#     np.array([10,30,50,4]),
#     np.array([30,80,60,7]),
#     np.array([20,40,30,9]),
#     np.array([60,30,40,15]),
#     np.array([40,50,10,25])
# ]

# def gen_functional_points(ips):
#     func_points = []
#     for poi in ips:
#         func_points = func_points + gen_surrounding_points(poi, r, ts)
#     return func_points

part_snap_galaxies = [list(gal) for gal in dbs_data if gal[1] == 19]

#Grabs all snapshot data for a particular galaxy
part_gal = [list(gal) for gal in dbs_data if gal[0] == 9994243][::-1]
#Cuts the data down to just x,y,z,rs
imp_data = [gal[3:] for gal in part_gal]



def nearby(centre_gal, other_gal):
    if abs(centre_gal[0] - other_gal[0]) < 10.0 and abs(centre_gal[1] - other_gal[1]) < 10.0 and abs(centre_gal[2] - other_gal[2]) < 10.0:
        return True
    else:
        return False

def gen_line(sv, ev, frames, centre_gal):
    '''
    input:
        sv, ev : vectors of form [x,y,z,redshift], start and end of line resp.
        frames: no of frames between two coords
        centre_gal:the coordinates of the centre of the galaxy of interest
    returns:
        sfs, xs, ys, zs b1,b2,b3,b4,b5,b6,b7,b8,b9,: A list of linearly interpolated x, y, z coords and 
            scalefactors spaced uniformly in log10(sf) and the basis vectors for the cammera, pointing towards the object

    ''' 
    start_sf= 1.0/(1.0+sv[3])
    end_sf = 1.0/(1.0+ev[3])
    start_logsf = np.log10(start_sf)
    end_logsf = np.log10(end_sf)
    array_log_sf = np.linspace(start_logsf, end_logsf, frames)
    array_sf = np.power(10, array_log_sf)

    frameNo = np.asarray(range(frames-1,-1,-1))
    Path_x, Path_y, Path_z = np.linspace(sv[0], ev[0], frames), np.linspace(sv[1], ev[1], frames), np.linspace(sv[2], ev[2], frames) 

    cen_position = np.zeros([frames, 3])
    cen_position[:] = (centre_gal[:3])

    b1 = (cen_position[:,0] - Path_x) 
    b2 = (cen_position[:,1] - Path_y) 
    b3 = (cen_position[:,2] - Path_z) 

    vect_Z = np.transpose(np.asarray([b1,b2,b3]))
    zvec = vect_Z

    xvec = np.cross([0.,1.,0.], zvec)

    # Find vector perpendicular to image x axis and view direction -
    # this will be the image y axis
    yvec = np.cross(zvec, xvec)

    # Normalize values to return
    v_x = xvec / np.linalg.norm(xvec,axis=1)[:,None]
    v_y = yvec / np.linalg.norm(yvec,axis=1)[:,None]
    v_z = zvec / np.linalg.norm(zvec,axis=1)[:,None]

    vectors = np.empty((frames, 3, 3), dtype='f8')

    direc_Of_travel = np.array([np.ediff1d(Path_x), np.ediff1d(Path_y), np.ediff1d(Path_z)])
    direc_Of_travel = direc_Of_travel / np.linalg.norm(direc_Of_travel, axis=1)


#Added things!!!

    print direc_Of_travel
    vectors[:,0,:] = v_x
    vectors[:,1,:] = v_y
    vectors[:,2,:] = v_z
    #setting the up direction as the positive y 
    # b4 = [0]*frames
    # b5 = [1]*frames
    # b6 = [0]*frames

    return np.array([frameNo,array_sf,Path_x ,Path_y , Path_z, vectors[:,0,0], vectors[:,0,1], vectors[:,0,2], vectors[:,1,0], vectors[:,1,1], vectors[:,1,2], vectors[:,2,0], vectors[:,2,1], vectors[:,2,2]])

####[1]*frames, [0]*frames, [0]*frames, [0]*frames, [1]*frames, [0]*frames, [0]*frames, [0]*frames, [1]*frames

gal = np.asarray(imp_data[17]) # take a specific snapshot data
# construct a line next to the position of a galaxy in that snapshot
sv = gal + np.array([5,5,-10,1.0])
ev = gal + np.array([5,5,10,1.0])

nearby_gal_datas = [galaxy[3:] for galaxy in part_snap_galaxies if nearby(gal, galaxy[3:])]

linePoints = np.transpose(gen_line(sv,ev,100, gal)) 
linePoints = linePoints[::-1]





#saves the paths of the flight plan to a text file
#np.savetxt("timeStoped.txt", linePoints,fmt='%i %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f %0.5f', header='RefL0100N1504',comments='#')




#the epansion factore that corresponds to the snapshot number. in order form snapshot 0 to 28,
#ie z=20 to z=0
# Expansion_F_snaps = [0.05, 0.06, 0.09, 0.10, 0.11, 0.12, 0.14, 0.15, 0.17,
#                      0.18, 0.20, 0.22, 0.25, 0.29, 0.31, 0.33, 0.37, 0.40,
#                      0.44, 0.50, 0.54, 0.58, 0.62, 0.67, 0.73, 0.79,0.85,
#                      0.91, 1.00]



# #load again to draw on graph
frame, ts, xs, ys, zs, b1,b2,b3,b4,b5,b6,b7,b8,b9 = np.loadtxt("timeStoped.txt", unpack=True)
path_of_cam =  np.loadtxt("timeStoped.txt")

# first_sf = path_of_cam[0,1]
# last_sf = path_of_cam[-1,1]
# print last_sf
# idx = (np.abs(Expansion_F_snaps - first_sf)).argmin()
# idxL = (np.abs(Expansion_F_snaps - last_sf)).argmin()
# print idx,idxL
# numberOfSnaps = abs(idx - idxL)
# print numberOfSnaps

# Inclusive_snapshots = Expansion_F_snaps[idx:idxL+1]
# print Inclusive_snapshots
# #to find the posistion of the line at the point which has the closest scale factor 
# # closest_point_path = 

# Line_plots_shown = []
# for Inclusive_snap in Inclusive_snapshots:
    

#     index = (np.abs(ts - Inclusive_snap)).argmin()
#     Line_plots_shown.append([xs[index], ys[index], zs[index]])

# print Line_plots_shown



#plots a 3d graph of the line and the galaxies near the object of interest at snap 19
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

ax.plot(xs,ys,zs)
ax.scatter(gal[0], gal[1], gal[2], marker="o", s=200.0, c="#682860")

for galaxy in nearby_gal_datas:
    ax.scatter(galaxy[0], galaxy[1], galaxy[2], marker="o", s=50.0)

ax.quiver(xs, ys, zs, b1,b2,b3,pivot="tail",color='r')
ax.quiver(xs, ys, zs, b4,b5,b6,pivot="tail",color='k')
ax.quiver(xs, ys, zs, b7,b8,b9,pivot="tail",color='c')

plt.show()







# #to first plot the x,y,z for our test. manually retrieved data
# Line_points_snaphots15_28 = [[83.14082, 83.72607, 16.10832],[83.14082, 83.72607, 14.08812],[83.14082, 83.72607, 12.67398],[83.14082, 83.72607, 11.05782],
#                              [83.14082, 83.72607, 8.63358],[83.14082, 83.72607, 7.21943],[83.14082, 83.72607, 6.00731],[83.14082, 83.72607, 4.79519],
#                               [83.14082, 83.72607, 3.38105],[83.14082, 83.72607, 1.76489],[83.14082, 83.72607, 0.35075],[83.14082, 83.72607, -1.06339],
#                               [83.14082, 83.72607, -2.27552],[83.14082, 83.72607, -2.27552]]


# print points
viewing_distance = 5.0
txt_name = "timeStoped"



# def story_board(dbs_data, viewing_distance, txt_name,Line_points_snaphots15_28):

#     data = dbs_data

#     gals = {}
#     for ele in data:
#       galID = ele[0]
#       if galID in gals.keys():
#           gals[galID].append([ele[i] for i in range(0, len(ele))])
#       else:
#           gals[galID] = [[ele[i] for i in range(0, len(ele))]]

#     #to creat a dictionary in order to extract the positions of the glaxies at each snapshot
#     snaps = {}
#     for galID in gals.keys()[:]:
#         gal = gals[galID]

#         for galsnap in gal:
#             if galsnap[1] not in snaps.keys():
#                 snaps[galsnap[1]] = [galsnap]
#             else: 
#                 snaps[galsnap[1]].append(galsnap)


#     fig = plt.figure()
#     for i in range(14):

#         All_snaps28 = snaps[15]
#         xyz_glas = []
#         too_close_behind = []
#         too_close_infornt = []


#         for k in range(len(All_snaps28)):
#             xyz_glas.append([All_snaps28[k][3], All_snaps28[k][4], All_snaps28[k][5], All_snaps28[k][0]])

    

#         center_x = Line_points_snaphots15_28[i][0]
#         center_y = Line_points_snaphots15_28[i][1]
#         center_z = Line_points_snaphots15_28[i][2]




#         for j in range(len(xyz_glas)):

#             if (abs((xyz_glas[j][0] - center_x)) <=viewing_distance and abs((xyz_glas[j][1] - center_y)) <=viewing_distance
#                 and abs((xyz_glas[j][2] - center_z)) <=viewing_distance):
                
#                 too_close_behind.append(xyz_glas[j])



#         nparr1 = np.asarray(too_close_behind)

#         labels = nparr1[:,3]
#         #print nparr1
#         plt.scatter(nparr1[:,0],nparr1[:,1],marker='x', c='b', s=50.) 
        
#         for g, txt in enumerate(labels):
#             plt.annotate(txt, (nparr1[g][0],nparr1[g][1]))
#         #plt.scatter(nparr2[:,0],nparr2[:,1],marker='x', c='g', s=25.) #infornt

#         plt.scatter(center_x,center_y,marker='o',s=200.0, c ='r')
#         plt.ylim(center_y - viewing_distance, center_y + viewing_distance)
#         plt.xlim(center_x - viewing_distance, center_x + viewing_distance)
#         plt.savefig(txt_name + str(i+15))
#         plt.clf()






def story_board(dbs_data, viewing_distance, txt_name, path_file):

    #the epansion factore that corresponds to the snapshot number. in order form snapshot 0 to 28,
    #ie z=20 to z=0
    Expansion_F_snaps = [0.05, 0.06, 0.09, 0.10, 0.11, 0.12, 0.14, 0.15, 0.17,
                         0.18, 0.20, 0.22, 0.25, 0.29, 0.31, 0.33, 0.37, 0.40,
                         0.44, 0.50, 0.54, 0.58, 0.62, 0.67, 0.73, 0.79,0.85,
                         0.91, 1.00]

    path_of_camera = np.loadtxt(path_file)
    frame, ts, xs, ys, zs, b1,b2,b3,b4,b5,b6,b7,b8,b9 = np.loadtxt(path_file, unpack=True)


    #to find the snapshot which the path starts at from the scale factor, ie the index 
    first_sf = path_of_camera[0,1]
    last_sf = path_of_camera[-1,1]
    idx = (np.abs(Expansion_F_snaps - first_sf)).argmin()
    idxL = (np.abs(Expansion_F_snaps - last_sf)).argmin()

    numberOfSnaps = abs(idx - idxL)

    #to find the scale factor at the snapshot you are at and to find the nearest x,y,z from the flight path
    Inclusive_snapshots = Expansion_F_snaps[idx:idxL+1]
    Line_plots_shown = []

    for Inclusive_snap in Inclusive_snapshots:

        index = (np.abs(ts - Inclusive_snap)).argmin()
        Line_plots_shown.append([xs[index], ys[index], zs[index]])
    
    print Line_plots_shown


    data = dbs_data

    #creates a dictionary of all the galaxies and all there snapshots 
    gals = {}
    for ele in data:
      galID = ele[0]
      if galID in gals.keys():
          gals[galID].append([ele[i] for i in range(0, len(ele))])
      else:
          gals[galID] = [[ele[i] for i in range(0, len(ele))]]

    #to creat a dictionary in order to extract the positions of the glaxies at each snapshot
    snaps = {}
    for galID in gals.keys()[:]:
        gal = gals[galID]

        for galsnap in gal:
            if galsnap[1] not in snaps.keys():
                snaps[galsnap[1]] = [galsnap]
            else: 
                snaps[galsnap[1]].append(galsnap)


    fig = plt.figure()
    for i in range(numberOfSnaps+1):

        All_snaps28 = snaps[idx+i]
        xyz_glas = []
        too_close_behind = []
        too_close_infornt = []


        for k in range(len(All_snaps28)):
            xyz_glas.append([All_snaps28[k][3], All_snaps28[k][4], All_snaps28[k][5], All_snaps28[k][0]])


        center_x = Line_plots_shown[i][0]
        center_y = Line_plots_shown[i][1]
        center_z = Line_plots_shown[i][2]




        for j in range(len(xyz_glas)):

            if (abs((xyz_glas[j][0] - center_x)) <=viewing_distance and abs((xyz_glas[j][1] - center_y)) <=viewing_distance
                and abs((xyz_glas[j][2] - center_z)) <=viewing_distance):
                
                too_close_behind.append(xyz_glas[j])



        nparr1 = np.asarray(too_close_behind)

        labels = nparr1[:,3]
        #print nparr1
        plt.scatter(nparr1[:,0],nparr1[:,1],marker='x', c='b', s=50.) 
        
        for g, txt in enumerate(labels):
            plt.annotate(txt, (nparr1[g][0],nparr1[g][1]))

        sf_time = Inclusive_snapshots[i]

        plt.title("Scale Factor: " + str(sf_time), loc='left', size=16)
        plt.scatter(center_x,center_y,marker='o',s=200.0, c ='r')
        plt.ylim(center_y - viewing_distance, center_y + viewing_distance)
        plt.xlim(center_x - viewing_distance, center_x + viewing_distance)
        plt.savefig(txt_name + str(i+15))
        plt.clf()

#story_board(dbs_data, viewing_distance, txt_name, "line_pointed_9994243.txt")


'''ddd ffffffffffffffffff rrrrrrrrrrrrrrrrr yyyyyyytreddddddddddddddddd
dddddddddddddd ggggggggggggggggggggg kkkkkkkkkkk shsjnsuger gwuiehrgiunsdfkfjg weiuuiggndfjgp 
sduhgueh guhghthi erugsoiensnfigsibg nsjdfgfpijg ;isrvs'''
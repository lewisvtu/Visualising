#rough test 
import numpy as np

# def story_board(dbs_data, viewing_distance, txt_name, path_file):

# 	#the epansion factore that corresponds to the snapshot number. in order form snapshot 0 to 28,
# 	#ie z=20 to z=0
# 	Expansion_F_snaps = [0.05, 0.06, 0.09, 0.10, 0.11, 0.12, 0.14, 0.15, 0.17,
# 	                     0.18, 0.20, 0.22, 0.25, 0.29, 0.31, 0.33, 0.37, 0.40,
# 	                     0.44, 0.50, 0.54, 0.58, 0.62, 0.67, 0.73, 0.79,0.85,
# 	                     0.91, 1.00]

# 	path_of_camera = np.loadtxt(path_file)
# 	frame, ts, xs, ys, zs, b1,b2,b3,b4,b5,b6,b7,b8,b9 = np.loadtxt(path_file, unpack=True)


# 	#to find the snapshot which the path starts at from the scale factor, ie the index 
# 	first_sf = path_of_camera[0,1]
# 	last_sf = path_of_camera[-1,1]
# 	idx = (np.abs(Expansion_F_snaps - first_sf)).argmin()
# 	idxL = (np.abs(Expansion_F_snaps - last_sf)).argmin()

# 	numberOfSnaps = abs(idx - idxL)

# 	#to find the scale factor at the snapshot you are at and to find the nearest x,y,z from the flight path
# 	Inclusive_snapshots = Expansion_F_snaps[idx:idxL+1]

# 	for Inclusive_snap in Inclusive_snapshots:

# 		Line_plots_shown = []
# 		index = (np.abs(ts - Inclusive_snap)).argmin()
# 		Line_plots_shown.append([xs[index], ys[index], zs[index]])



#     data = dbs_data

#     #creates a dictionary of all the galaxies and all there snapshots 
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
#     for i in range(numberOfSnaps+1):

#         All_snaps28 = snaps[idx+i]
#         xyz_glas = []
#         too_close_behind = []
#         too_close_infornt = []


#         for k in range(len(All_snaps28)):
#             xyz_glas.append([All_snaps28[k][3], All_snaps28[k][4], All_snaps28[k][5], All_snaps28[k][0]])


#         center_x = Line_plots_shown[i][0]
#         center_y = Line_plots_shown[i][1]
#         center_z = Line_plots_shown[i][2]




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

# def cam_basis(Path_x,Path_y,Path_z):




# 	#Path_x, Path_y, Path_z = np.linspace(sv[0], ev[0], frames), np.linspace(sv[1], ev[1], frames), np.linspace(sv[2], ev[2], frames) 

# 	cen_position = np.zeros([frames, 3])
# 	cen_position[:] = (centre_gal[:3])

# 	b1 = (cen_position[:,0] - Path_x) 
# 	b2 = (cen_position[:,1] - Path_y) 
# 	b3 = (cen_position[:,2] - Path_z) 

# 	vect_Z = np.transpose(np.asarray([b1,b2,b3]))
# 	zvec = vect_Z

# 	direc_Of_travel = np.array([np.ediff1d(Path_x), np.ediff1d(Path_y), np.ediff1d(Path_z)])
# 	direc_Of_travel = direc_Of_travel / np.linalg.norm(direc_Of_travel)


 





# 	xvec = np.cross([0.,1.,0.], zvec)

# 	# Find vector perpendicular to image x axis and view direction -
# 	# this will be the image y axis
# 	yvec = np.cross(zvec, xvec)

# 	# Normalize values to return
# 	v_x = xvec / np.linalg.norm(xvec,axis=1)[:,None]
# 	v_y = yvec / np.linalg.norm(yvec,axis=1)[:,None]
# 	v_z = zvec / np.linalg.norm(zvec,axis=1)[:,None]

# 	vectors = np.empty((frames, 3, 3), dtype='f8')

# 	vectors[:,0,:] = v_x
# 	vectors[:,1,:] = v_y
# 	vectors[:,2,:] = v_z	

# def gen_line(sv, ev, frames, centre_gal):
#     '''
#     input:
#         sv, ev : vectors of form [x,y,z,redshift], start and end of line resp.
#         frames: no of frames between two coords
#         centre_gal:the coordinates of the centre of the galaxy of interest
#     returns:
#         sfs, xs, ys, zs b1,b2,b3,b4,b5,b6,b7,b8,b9,: A list of linearly interpolated x, y, z coords and 
#             scalefactors spaced uniformly in log10(sf) and the basis vectors for the cammera, pointing towards the object

#     ''' 
#     start_sf= 1.0/(1.0+sv[3])
#     end_sf = 1.0/(1.0+ev[3])
#     start_logsf = np.log10(start_sf)
#     end_logsf = np.log10(end_sf)
#     array_log_sf = np.linspace(start_logsf, end_logsf, frames)
#     array_sf = np.power(10, array_log_sf)

#     frameNo = np.asarray(range(frames-1,-1,-1))
#     Path_x, Path_y, Path_z = np.linspace(sv[0], ev[0], frames), np.linspace(sv[1], ev[1], frames), np.linspace(sv[2], ev[2], frames) 

#     cen_position = np.zeros([frames, 3])
#     cen_position[:] = (centre_gal[:3])

#     b1 = (cen_position[:,0] - Path_x) 
#     b2 = (cen_position[:,1] - Path_y) 
#     b3 = (cen_position[:,2] - Path_z) 

#     vect_Z = np.transpose(np.asarray([b1,b2,b3]))
#     zvec = vect_Z

#     xvec = np.cross([0.,1.,0.], zvec)

#     # Find vector perpendicular to image x axis and view direction -
#     # this will be the image y axis
#     yvec = np.cross(zvec, xvec)

#     # Normalize values to return
#     v_x = xvec / np.linalg.norm(xvec,axis=1)[:,None]
#     v_y = yvec / np.linalg.norm(yvec,axis=1)[:,None]
#     v_z = zvec / np.linalg.norm(zvec,axis=1)[:,None]

#     vectors = np.empty((frames, 3, 3), dtype='f8')

#     vectors[:,0,:] = v_x
#     vectors[:,1,:] = v_y
#     vectors[:,2,:] = v_z
#     #setting the up direction as the positive y 
#     # b4 = [0]*frames
#     # b5 = [1]*frames
#     # b6 = [0]*frames

#     return np.array([frameNo,array_sf,Path_x ,Path_y , Path_z, vectors[:,0,0], vectors[:,0,1], vectors[:,0,2], vectors[:,1,0], vectors[:,1,1], vectors[:,1,2], vectors[:,2,0], vectors[:,2,1], vectors[:,2,2]])



#to rotate the positions of the galaxies for the 2d projection to work 
#this will use yaw, pitch and roll 

#the amtrix of rotation, a:about z, b:about y and c:about x axis all counterclockwise

# Calculates Rotation Matrix given euler angles.
def eulerAnglesToRotationMatrix(thetas) :
     
    R_x = np.array([[[1]*len(thetas),         [0]*len(thetas),        [0]*len(thetas)                   ],
                    [[0]*len(thetas),         np.cos(thetas[:,0]), -np.sin(thetas[:,0]) ],
                    [[0]*len(thetas),         np.sin(thetas[:,0]), np.cos(thetas[0])  ]
                    ])
         
         
                     
    R_y = np.array([[np.cos(thetas[1]),    [0]*len(thetas),      np.sin(thetas[1])  ],
                    [[0]*len(thetas),       [1]*len(thetas),      [0]*len(thetas)     ],
                    [-np.sin(thetas[1]),   [0]*len(thetas),      np.cos(thetas[1])  ]
                    ])
                 
    R_z = np.array([[np.cos(thetas[2]),    -np.sin(thetas[2]),    [0]*len(thetas)],
                    [np.sin(thetas[2]),    np.cos(thetas[2]),     [0]*len(thetas)],
                    [[0]*len(thetas),       [0]*len(thetas),        [1]*len(thetas)]
                    ])
                     
                     
    R = R_z * R_y * R_z
 
    return R 

print eulerAnglesToRotationMatrix(np.array([[1.,2.,3.],[1.1,2.2,3.3],[4.,5.6,6.]]))
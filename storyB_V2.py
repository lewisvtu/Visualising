import matplotlib.pyplot as plt
import numpy as np 
import pickler.shelf as shelf
from DBS.dbgrabber import dbsPull
import matplotlib.pyplot as plt
from camTest import perspective_transfomation
from camTest import coord_transform
from timeit import default_timer as timer
strt = timer()
h = 0.6777
region = [10., 10., 10.]

#SQL to grab from the database 
SQL = """
    SELECT
        DES.GalaxyID,
        PROG.SnapNum,
        PROG.MassType_Star,
        (PROG.CentreOfPotential_x * %0.5f) as x,
        (PROG.CentreOfPotential_y * %0.5f) as y,
        (PROG.CentreOfPotential_z * %0.5f) as z,
        PROG.Redshift
    FROM
        RefL0100N1504_Subhalo as PROG with(forceseek),
        RefL0100N1504_Subhalo as DES,
        RefL0100N1504_Aperture as AP
    WHERE
        DES.SnapNum = 28 and
        DES.MassType_Star > 1.0e10 and
        DES.MassType_DM > 1.0e11 and
        PROG.GalaxyID between DES.GalaxyID and DES.TopLeafID and
        AP.ApertureSize = 30 and
        AP.GalaxyID = DES.GalaxyID and
        AP.Mass_Star > 1.0e10
    ORDER BY
        PROG.GalaxyID,
        PROG.SnapNum
""" % (h,h,h)

txt_name = "ClippedAttempt2_"
filename = "scaled_DBS.p"

raw_dbs = dbsPull(SQL, filename)
shelf.push(raw_dbs, "scaled_DBS")
dbs_data = shelf.pull("scaled_DBS")



def orderGals(dbs_data, snapshot_num):
	#creates a dictionary of all the galaxies and all there snapshots 
    gals = {}
    for ele in dbs_data:
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

    return snaps[snapshot_num]



def story_board(txt_name, path_file, All_galaxies):

	frame, ts, xs, ys, zs, b1,b2,b3,b4,b5,b6,b7,b8,b9 = np.loadtxt(path_file, unpack=True)
	z_basis = np.transpose(np.asarray([b7, b8, b9]))
	y_basis = np.transpose(np.asarray([b4, b5, b6]))
	x_basis = np.transpose(np.asarray([b1, b2, b3]))
	line_coords = np.transpose(np.asarray([xs, ys, zs]))


	fig = plt.figure()
	#to loop over every frame which will be plot on the story board
	for i in range(len(frame)):

		#to find the new camera position on the path
		cam_position = [xs[i], ys[i], zs[i]]
		x_bas = x_basis[i]
		y_bas = y_basis[i]
		z_bas = z_basis[i]


		galaxies_trans = coord_transform(x_bas, y_bas, z_bas, cam_position, All_galaxies)

		galaxies_afterT = np.transpose(galaxies_trans)


		galaxies_to_plot = perspective_transfomation(galaxies_trans, region)




		perspec = []

		indexList = np.asarray(galaxies_to_plot[:,4], dtype=int)

		galaxZs = galaxies_afterT[indexList][:,2]
		print galaxZs
		print cam_position



		perspec = 10000./galaxZs**2




		plt.scatter(galaxies_to_plot[:,0],galaxies_to_plot[:,1],marker='o', s=perspec, c='c')
		# plt.ylim( 0, region[1])
		# plt.xlim( - region[0]/2., region[0]/2.)
		plt.ylim( -1., 1.)
		plt.xlim( - 1., 1.)
		plt.savefig(txt_name + str(i+1))
		plt.clf()




All_galaxiesDATA = np.asarray(orderGals(dbs_data,28))
All_galaxiesXYZ = All_galaxiesDATA[:,[3,4,5]]

story_board( txt_name, "gla200_0_orbit.txt", All_galaxiesXYZ)
stp = timer()
print "time taken: %f" %(stp - strt) 

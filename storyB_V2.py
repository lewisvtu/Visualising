import matplotlib.pyplot as plt
import numpy as np 
import pickler.shelf as shelf
from DBS.dbgrabber import dbsPull
import matplotlib.pyplot as plt
from utils import perspective_transfomation
from utils import coord_transform
from timeit import default_timer as timer
import utils
from scipy.misc import imread

h = 0.6777
region = [15., 15., 15.]

Expansion_F_snaps = np.array([0.05, 0.06, 0.09, 0.10, 0.11, 0.12, 0.14, 0.15, 0.17,
                     0.18, 0.20, 0.22, 0.25, 0.29, 0.31, 0.33, 0.37, 0.40,
                     0.44, 0.50, 0.54, 0.58, 0.62, 0.67, 0.73, 0.79,0.85,
                     0.91, 1.00])


#SQL to grab from the database 
SQL = """
    SELECT
    	PROG.GalaxyID as ID,
    	PROG.DescendantID as DesID,
        DES.GalaxyID,
        PROG.SnapNum,
        PROG.MassType_DM,
        (PROG.CentreOfPotential_x * %0.5f) as x,
        (PROG.CentreOfPotential_y * %0.5f) as y,
        (PROG.CentreOfPotential_z * %0.5f) as z,
        PROG.Redshift
    FROM
        RefL0025N0376_Subhalo as PROG with(forceseek),
        RefL0025N0376_Subhalo as DES

    WHERE
        DES.SnapNum = 28 and
        DES.MassType_DM > 1.0e11 and
        PROG.GalaxyID between DES.GalaxyID and DES.LastProgID 

    ORDER BY
        PROG.GalaxyID,
        PROG.SnapNum
""" % (h,h,h)
#        PROG.MassType_DM > 1.0e11 and
txt_name = "no_images_new_"
filename = "allProgs11_DBS.p"
boxsize = 25 * h
dbs_data = dbsPull(SQL, filename)

strt = timer()






#print dbs_data['ID'], dbs_data['DesID']
#assert False 



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



def orderGals(gals, snapshot_num):          

    #to creat a dictionary in order to extract the positions of the glaxies at each snapshot
    snaps = {}
    for galID in gals.keys()[:]:
        gal = gals[galID]

        for galsnap in gal:
            if galsnap[1] not in snaps.keys():
                snaps[galsnap[1]] = [galsnap]
            else: 
                snaps[galsnap[1]].append(galsnap)

    return snaps



def story_board(txt_name, path_file):

	''' This function returns '''

	frame, ts, xs, ys, zs, b1,b2,b3,b4,b5,b6,b7,b8,b9 = np.loadtxt(path_file, unpack=True)
	z_basis = np.transpose(np.asarray([b7, b8, b9]))
	y_basis = np.transpose(np.asarray([b4, b5, b6]))
	x_basis = np.transpose(np.asarray([b1, b2, b3]))
	line_coords = np.transpose(np.asarray([xs, ys, zs]))



	fig = plt.figure()
	#to loop over every frame which will be plot on the story board
	for i in range(len(frame)):
		print "creating image: " + str(i)
		#to find the new camera position on the path
		cam_position = [xs[i], ys[i], zs[i]]
		x_bas = x_basis[i]
		y_bas = y_basis[i]
		z_bas = z_basis[i]
		scale_factor = ts[i]

		#wrap coordinates image data get center image manipulations preiodic wrap
		centre = utils.get_centre([x_bas,y_bas,z_bas], cam_position, region)

		# #times the box size by 0.6777 

		All_galaxies = utils.gal_interpolation(scale_factor, dbs_data)
		All_galaxies[:,[3,4,5]] = utils.periodic_wrap(All_galaxies[:,[3,4,5]], boxsize, centre)


		#print All_galaxies
		galaxies_trans = coord_transform(x_bas, y_bas, z_bas, cam_position, All_galaxies[:,[3,4,5]])

		galaxies_afterT = np.transpose(galaxies_trans)
		#print galaxies_afterT

		galaxies_to_plot = perspective_transfomation(galaxies_trans, region)

		#this is to catch any frames that are of zero lenght ie no galaxies in view 
		if len(galaxies_to_plot) >= 1:

			perspec = []

			indexList = np.asarray(galaxies_to_plot[:,4], dtype=int)

			galaxZs = galaxies_afterT[indexList][:,2]
			galZsMass = All_galaxies[indexList][:,2]



			perspec = 1./galaxZs**3
			perspec *= (galZsMass)**0.43
			perspec.shape = (1, len(perspec))



			galaxies_to_plot = np.asarray(sorted(np.concatenate((galaxies_to_plot, perspec.T), axis=1), key=lambda coords: -coords[2]))
			#img = imread("gas_%06i.png"%(i))
			plt.scatter(galaxies_to_plot[:,0],galaxies_to_plot[:,1],marker='o', s=galaxies_to_plot[:,5], c='#7E317B',edgecolors='k')
			plt.ylim( 0, region[1])
			plt.xlim( - region[0]/2., region[0]/2.)
			plt.ylim( -1., 1.)
			plt.xlim( - 1., 1.)
			#plt.imshow(img,extent=[-1.,1.,-1.,1.], aspect='auto')
			plt.savefig(txt_name + str(i))
			plt.clf()

		else:
			plt.scatter(0.0,0.0, s=0.0)
			plt.ylim( -1., 1.)
			plt.xlim( - 1., 1.)
			plt.savefig(txt_name + str(i))
			plt.clf()



# gals = galsTree(dbs_data)
# snaps = orderGals(gals, 1)

#story_board( txt_name, "Orbit_through_les_time.txt")
stp = timer()

print "time taken: %f" %(stp - strt) 










def runTkinter():

	h = 0.6777
	region = [15., 15., 15.]

	Expansion_F_snaps = np.array([0.05, 0.06, 0.09, 0.10, 0.11, 0.12, 0.14, 0.15, 0.17,
	                 0.18, 0.20, 0.22, 0.25, 0.29, 0.31, 0.33, 0.37, 0.40,
	                 0.44, 0.50, 0.54, 0.58, 0.62, 0.67, 0.73, 0.79,0.85,
	                 0.91, 1.00])

	#SQL to grab from the database 
	SQL = """
	SELECT
		PROG.GalaxyID as ID,
		PROG.DescendantID as DesID,
	    DES.GalaxyID,
	    PROG.SnapNum,
	    PROG.MassType_DM,
	    (PROG.CentreOfPotential_x * %0.5f) as x,
	    (PROG.CentreOfPotential_y * %0.5f) as y,
	    (PROG.CentreOfPotential_z * %0.5f) as z,
	    PROG.Redshift
	FROM
	    RefL0025N0376_Subhalo as PROG with(forceseek),
	    RefL0025N0376_Subhalo as DES

	WHERE
	    DES.SnapNum = 28 and
	    DES.MassType_DM > 1.0e11 and
	    PROG.GalaxyID between DES.GalaxyID and DES.LastProgID 

	ORDER BY
	    PROG.GalaxyID,
	    PROG.SnapNum
	""" % (h,h,h)

	
	txt_name = "allProgsTest_"
	filename = "allProgs8_DBS.p"

	dbs_data = dbsPull(SQL, filename)
	story_board( txt_name, "Orbit_through_les_time.txt")


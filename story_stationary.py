import matplotlib.pyplot as plt
import numpy as np 
import pickler.shelf as shelf
import math
from DBS.dbgrabber import dbsPull
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from camTest import perspective_transfomation
h = 0.6777

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

region = [10., 10., 10.] 
viewing_distance = region[2]
txt_name = "orbitWithCliped_"
filename = "FollowProgs19.p"
raw_dbs = dbsPull(SQL, filename)
shelf.push(raw_dbs, "followup19")
dbs_data = shelf.pull("followup19")





def story_board(dbs_data, viewing_distance, txt_name, path_file):

    frame, ts, xs, ys, zs, b1,b2,b3,b4,b5,b6,b7,b8,b9 = np.loadtxt(path_file, unpack=True)
    z_basis_s = np.transpose(np.asarray([b7, b8, b9]))
    y_basis_s = np.transpose(np.asarray([b4, b5, b6]))
    x_basis_s = np.transpose(np.asarray([b1, b2, b3]))


    Line_plots_shown = []
    z_basis = []
    y_basis = []
    x_basis = []

    for test in range(len(xs)):
    	Line_plots_shown.append([xs[test], ys[test], zs[test]])
    	z_basis.append(z_basis_s[test])
    	y_basis.append(y_basis_s[test])
    	x_basis.append(x_basis_s[test])

    #print Line_plots_shown

    x_basis = np.asarray(x_basis)
    y_basis = np.asarray(y_basis)
    z_basis = np.asarray(z_basis)
    #print x_basis, y_basis, z_basis
    # #print Line_plots_shown

    data = dbs_data
    #print data

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
    #print "snaps "
    print snaps[28]
    
    fig = plt.figure()

    for i in range(len(frame)):

        All_snaps28 = snaps[28]
        xyz_glas = []
        too_close_behind = []


        for k in range(len(All_snaps28)):
            xyz_glas.append([All_snaps28[k][3], All_snaps28[k][4], All_snaps28[k][5], All_snaps28[k][0], All_snaps28[k][2]])
        print "xyzglas ---------------------------------------"
        print xyz_glas
        assert False
        center_x = Line_plots_shown[i][0]
        center_y = Line_plots_shown[i][1]
        center_z = Line_plots_shown[i][2]

        cam_position = np.transpose(np.array([center_x, center_y, center_z]))
        x_bas = x_basis[i]
        y_bas = y_basis[i]
        z_bas = z_basis[i]




        for j in range(len(xyz_glas)):

            if (abs((xyz_glas[j][0] - center_x)) <=viewing_distance and abs((xyz_glas[j][1] - center_y)) <=viewing_distance
                and abs((xyz_glas[j][2] - center_z)) <=viewing_distance):
                
                particles = np.transpose( np.array([[xyz_glas[j][0]], [xyz_glas[j][1]], [xyz_glas[j][2]]]) )

                arr = perspective_transfomation(x_bas, y_bas, z_bas, cam_position, particles, region)
                new_arr = arr.flatten()[:3]

                values_for_plot = np.r_[new_arr, xyz_glas[j][3]]
                values_for_plot = np.r_[values_for_plot, xyz_glas[j][4]]
                too_close_behind.append(values_for_plot)


        positive_plot = []
        for galaxy in too_close_behind:
            #print galaxy
            if galaxy[2] >= 0.0:
                positive_plot.append(galaxy)







        nparr1 = np.asarray(positive_plot)





        #print "nparr1 nnnnnnnnnnnnnnnnnnnnnnn"
        #print nparr1

        if len(nparr1) >= 1.0:
            labels = nparr1[:,3]
            perspec_size = 100.#10000 / nparr1[:,2]**2
            colors = []

            for mass in nparr1[:,4]:
                if mass <= 1.0e10:
                    colors.append("g")
                elif mass > 1.0e10 and mass <= 1.0e11:
                    colors.append("#FF8C00")
                elif mass > 1.0e11:
                    colors.append("r")    

            print colors

            plt.scatter(nparr1[:,0],nparr1[:,1],marker='o', s=perspec_size, c=colors) 
            

            for g, txt in enumerate(labels):
                plt.annotate(txt, (nparr1[g][0],nparr1[g][1]))

            plt.ylim( - viewing_distance, viewing_distance)
            plt.xlim( - viewing_distance, viewing_distance)
            plt.savefig(txt_name + str(i+1))
            plt.clf()

story_board(dbs_data, viewing_distance, txt_name, "gla200_0_orbit.txt")
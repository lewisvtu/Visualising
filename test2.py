from eagleSqlTools import connect, execute_query
import matplotlib.pyplot as plt
import numpy as np 
import pickle 
from mpl_toolkits.mplot3d import Axes3D

USERID = "llight"
PASSWORD = "CX2vy392"

CON = connect(USERID, PASSWORD)
PICKFILE = "1pcGals.p"


def dbsQuery(sql):
    return execute_query(CON,sql)

#to import the data

# SQL = """
# SELECT 
#  	DES.GalaxyID, 
#  	PROG.SnapNum,
#  	PROG.Mass,   
#  	PROG.CentreOfPotential_x, 
#   	PROG.CentreOfPotential_y, 
#  	PROG.CentreOfPotential_z,
# 	PROG.GalaxyID as ProgID                        
#  FROM  
#  	RefL0100N1504_Subhalo as PROG with(forceseek),       
#  	RefL0100N1504_Subhalo as DES, 
#  	RefL0100N1504_Aperture as AP           
#  WHERE  
#  	DES.SnapNum = 28                       
#  	and DES.MassType_Star between 1.0e10 and 6e10                              
#     and DES.MassType_DM between 5.0e11 and 2.0e12                           
#  	and DES.RandomNumber < 0.01 
#  	and PROG.GalaxyID between DES.GalaxyID and DES.LastProgID 
#  	and AP.ApertureSize = 30 
#  	and AP.GalaxyID = DES.GalaxyID 
#  	and AP.Mass_Star > 1.0e9  
#  ORDER BY  
 	 
#  	PROG.GalaxyID,
#  	PROG.SnapNum
#     """
# data = dbsQuery(SQL)

# gals = {}
# for ele in data:
# 	galID = ele[0]
# 	if galID in gals.keys():
# 		gals[galID].append([ele[i] for i in range(1, len(ele))])
# 	else:
# 		gals[galID] = [[ele[i] for i in range(1, len(ele))]]

# with open(PICKFILE, "wb") as pfile:
# 	pickle.dump(gals, pfile)

# with open(PICKFILE, "rb") as pfile:
# 	gals2 = pickle.load(pfile)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# gal1 = gals2[9744960]
# print gal1
# with open("first_gal_progs.p", "wb") as pfile:
# 	pickle.dump(gal1, pfile)

with open("first_gal_progs.p", "rb") as pfile:
	first_gal_progs = pickle.load(pfile)
#print first_gal_progs

for sg in first_gal_progs:
	print sg[0], sg[5]
#first_gal_progs.sort(key=lambda snap: snap[])

# for galID in gals2.keys()[:1]:
# 	print galID
# 	gal = gals2[galID]
# 	xs = []
# 	ys = []
# 	zs = []
	
# 	for snap in gal:

# 		mass = snap[1]
# 		x,y,z = snap[2], snap[3], snap[4]
# 		mark_s = mass/1e9
# 		ax.scatter(x,y,z,marker='o', c='#7E317B')
# 		xs.append(x)
# 		ys.append(y)
# 		zs.append(z)

# 	ax.plot(xs,ys,zs)



plt.show()
import matplotlib.pyplot as plt
import numpy as np
import pickler.shelf as shelf
import math
from DBS.dbgrabber import dbsPull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


SQL = """
    SELECT DISTINCT
        DES.GalaxyID,
        PROG.SnapNum,
        PROG.Mass,
        PROG.CentreOfPotential_x,
        PROG.CentreOfPotential_y,
        PROG.CentreOfPotential_z,
        PROG.Redshift,
        DES.GroupNumber,
        DES.SubGroupNumber
    FROM
        RefL0100N1504_Subhalo as PROG with(forceseek),
        RefL0100N1504_Subhalo as DES,
        RefL0100N1504_Aperture as AP
    WHERE
        PROG.SnapNum = 28 and
        PROG.GalaxyID between DES.GalaxyID and DES.TopLeafID and
        AP.GalaxyID = DES.GalaxyID and
        DES.GroupNumber IN (200) and
        DES.SubGroupNumber = 0

"""

# Grabs new data from db based on sql. If file name already exists, it loads that data instead

filename = "aBigThings.p"

# raw_dbs = dbsPull(SQL, filename)
# print raw_dbs
fs, sfs, xs, ys, zs, v1xs, v1ys, v1zs, v2xs, v2ys, v2zs, v3xs, v3ys, v3zs = np.loadtxt("orbit150.txt", unpack=True)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.plot(xs, ys, zs)


ax.quiver(xs,ys,zs, v1xs, v1ys, v1zs, color="#682860", pivot="tail")
ax.quiver(xs,ys,zs, v2xs, v2ys, v2zs, color="#000000", pivot="tail")
ax.quiver(xs,ys,zs, v3xs, v3ys, v3zs, color="#FF0000", pivot="tail")
plt.show()
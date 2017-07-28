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
        DES.GroupNumber IN (200, 175, 150, 125, 100, 75) and
        DES.SubGroupNumber = 0

"""

# Grabs new data from db based on sql. If file name already exists, it loads that data instead

filename = "6BigThings.p"

raw_dbs = dbsPull(SQL, filename)
print raw_dbs

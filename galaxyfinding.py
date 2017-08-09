import matplotlib.pyplot as plt
import numpy as np
import pickler.shelf as shelf
import math
from DBS.dbgrabber import dbsPull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
h=0.6777

SQL = """
    SELECT DISTINCT
        DES.GalaxyID,
        PROG.SnapNum,
        PROG.Mass,
        (PROG.CentreOfPotential_x * %0.5f) as x,
        (PROG.CentreOfPotential_y * %0.5f) as y,
        (PROG.CentreOfPotential_z * %0.5f) as z,
        PROG.Redshift,
        DES.GroupNumber,
        DES.SubGroupNumber
    FROM
        RefL0025N0376_Subhalo as PROG with(forceseek),
        RefL0025N0376_Subhalo as DES,
        RefL0025N0376_Aperture as AP
    WHERE
        PROG.SnapNum = 28 and
        PROG.GalaxyID between DES.GalaxyID and DES.TopLeafID and
        AP.GalaxyID = DES.GalaxyID and
        DES.GroupNumber IN (1) and
        DES.SubGroupNumber = 0

""" % (h,h,h)

# Grabs new data from db based on sql. If file name already exists, it loads that data instead

filename = "aBigThing.p"

raw_dbs = dbsPull(SQL, filename)
print raw_dbs

from DBS.dbgrabber import dbsPull
import numpy as np
h=0.6777

SQL = """
    SELECT DISTINCT
        DES.GalaxyID as id,
        PROG.SnapNum as sn,
        PROG.Mass as ma,
        (PROG.CentreOfPotential_x * %0.5f) as x,
        (PROG.CentreOfPotential_y * %0.5f) as y,
        (PROG.CentreOfPotential_z * %0.5f) as z,
        PROG.Redshift as rs,
        DES.GroupNumber as gn,
        DES.SubGroupNumber as sgn
    FROM
        RefL0025N0376_Subhalo as PROG with(forceseek),
        RefL0025N0376_Subhalo as DES,
        RefL0025N0376_Aperture as AP
    WHERE
        PROG.SnapNum = 28 and
        PROG.GalaxyID between DES.GalaxyID and DES.TopLeafID and
        AP.GalaxyID = DES.GalaxyID and
        DES.GroupNumber IN (2,4) and
        DES.SubGroupNumber = 0

""" % (h,h,h)

# Grabs new data from db based on sql. If file name already exists, it loads that data instead

filename = "BigThing.p"

raw_dbs = dbsPull(SQL, filename, make=False, search=False)

print raw_dbs

import matplotlib.pyplot as plt
import numpy as np 
import pickler.shelf as shelf
import math
from DBS.dbgrabber import dbsPull
import scipy
TIMESTEPS = 100

SQL = """
    SELECT  DISTINCT
        PROG.DescendantID,
        PROG.Snapno,
        PROG.CentreOfPotential_x,
        PROG.CentreOfPotential_y,
        PROG.CentreOfPotential_z,
        PROG.Redshift
    FROM
        RefL0100N1504_Subhalo as PROG,
        RefL0100N1504_Aperture as AP
    WHERE
        PROG.Mass > 5E10 and
        PROG.RandomNumber < 0.002 and
        AP.ApertureSize = 30 and
        PROG.Spurious = 0 and
        PROG.CentreOfPotential_x BETWEEN 25 and 75 and
        PROG.CentreOfPotential_y BETWEEN 25 and 75 and
        PROG.CentreOfPotential_z BETWEEN 25 and 75
"""
filename = "Interesting1.p"

raw_dbs = dbsPull(SQL, filename)
print raw_dbs
interesting_gals = {}
for galsnap in raw_dbs:
    if galsnap[1] not in interesting_gals.keys():
        interesting_gals[galsnap[1]] = galsnap

shelf.push(interesting_gals, "IntGalsPSnap")

interesting_gals = shelf.pull("IntGalsPSnap")

def get_camera_pos(t):
    sgal = interesting_gals[int(math.floor(t))]
    egal = interesting_gals[int(math.ceil(t))]
    dt, st = math.modf(t)
    sx,sy,sz = sgal[2], sgal[3], sgal[4]
    ex,ey,ez = egal[2], egal[3], egal[4]
    dx, dy, dz = ex-sx, ey-sy, ez-sz

    return sx + dt*dx, sy + dt*dy, sz + dt*dz
    
print get_camera_pos(16.4563)

sx, sy, sz = 50.0,50.0,50.0
tx, ty, tz = 60.0,60.0,60.0


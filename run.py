import utils
import flightplan_generator
import numpy as np


h = 0.6777
SQL = """
    SELECT
        DES.GalaxyID,
        PROG.SnapNum,
        PROG.MassType_Star,
        (PROG.CentreOfPotential_x * %f) as x,
        (PROG.CentreOfPotential_y * %f) as y,
        (PROG.CentreOfPotential_z * %f) as z,
        PROG.Redshift
    FROM
        RefL00250N0376_Subhalo as PROG with(forceseek),
        RefL00250N0376_Subhalo as DES,
        RefL00250N0376_Aperture as AP
    WHERE
        DES.SnapNum = 28 and
        DES.MassType_Star > 1.0e10 and
        DES.MassType_DM > 1.0e11 and
        PROG.GalaxyID between DES.GalaxyID and DES.TopLeafID and
        AP.ApertureSize = 30 and
        AP.GalaxyID = DES.GalaxyID and
        AP.Mass_Star > 1.0e9
    ORDER BY
        PROG.GalaxyID,
        PROG.SnapNum
""" % (h,h,h)

dbs_data = utils.Data(SQL, "new_sim")

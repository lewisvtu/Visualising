import eagleSqlTools as dbt
import pickle
import os

USERID = "llight"
PASSWORD = "CX2vy392"

def dbsQuery(sql):
    con = dbt.connect(USERID, PASSWORD)
    return dbt.execute_query(con,sql)

def dbsPull(sql, fname):
    if not fname[-2:] == ".p":
        fname = fname + ".p"
    cwd = os.getcwd()
    path = os.path.join("records", fname)
    path = os.path.join("DBS", path)
    path = os.path.join(cwd, path)
    try:
        with open(path, "rb") as pfile:
            print "Pulling from records"
            return pickle.load(pfile)
    except:
        f = open(path, "w")
        f.close()
        print "Pulling from DB"
        data = dbsQuery(sql)
        with open(path, "wb") as pfile:
            pickle.dump(data, pfile)
        return data

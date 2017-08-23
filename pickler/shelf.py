import pickle
import os

def getpath(fname):
    if not fname[-2:] == ".p":
        fname = fname + ".p"
    cwd = os.getcwd()
    path = os.path.join("shelves", fname)
    path = os.path.join("pickler", path)
    path = os.path.join(cwd, path)
    return path

def pull(fname):
    path = getpath(fname)
    try:
        with open(path, "rb") as pfile:
            return pickle.load(pfile)
    except:
        return False

def push(data, fname):
    path = getpath(fname)
    try:
        with open(path, "wb") as pfile:
            pickle.dump(pfile)
    except:
        f = open(path, "w")
        f.close()
        with open(path, "wb") as pfile:
            pickle.dump(data, pfile)
import cPickle as pickle
import numpy as np
import scipy
def div_fact(num,denom):
        try:
                return round(float(num)/denom,3)
        except:
                return 0

def readPickle(fname):
    with open (fname) as fin:
        has=pickle.load(fin)
    return has

def writePickle (fname,has):
    with open(fname, "w") as fin:
        pickle.dump(has, fin)

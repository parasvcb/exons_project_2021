import cPickle as pickle
import numpy as np
import gzip,sys
import scipy,os
import constructing_data.Classes_exons as Classes_exons
sys.modules['Classes_exons'] = Classes_exons


def div_fact(num,denom):
        try:
                return round(float(num)/denom,3)
        except:
                return 0

def readFileSepslashN(fname,lengthofLine=0):
    if ".gz" in fname:
        with gzip.open(fname) as fin:
            dat = [i for i in fin.read().split("\n") if len(i) > lengthofLine]
    else:
        with open(fname) as fin:
            dat = [i for i in fin.read().split("\n") if len(i) > lengthofLine]
    return dat

def stats(lis):
    lis=[i for i in lis if i>0]
    if len(lis):
        string="\t".join(map(str,[len(lis),sum(lis),np.max(lis),np.min(lis),round(np.mean(lis),3),round(np.median(lis),3),scipy.stats.mode(lis)[0][0],scipy.stats.mode(lis)[1][0],round(np.std(lis),3)]))
    else:
        string="\t".join(map(str,[len(lis),sum(lis),0,0,0,0,0,0,0]))
    return string

def chunks_based_on_cores(l, n):
    div_fac = int(ceil(float(len(l))/n))
    """Yield successive n-sized chunks_based_on_cores from l."""
    for i in range(0, n+1):
        yield l[i*div_fac: (i+1)*div_fac]

def chunks_based_on_element_size(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

def readPickle(fname):
    with open (fname) as fin:
        has=pickle.load(fin)
    return has

def writePickle (fname,has):
    with open(fname, "w") as fin:
        pickle.dump(has, fin)
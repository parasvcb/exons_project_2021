
import numpy as np

def CDF_String(lis,columns):
    string=''
    temp_binst=np.arange(0,1.001,0.001,dtype=float)
    total_t=len(lis)
    string=",".join(columns)+"\n"
    for i in temp_binst:
        valst=0
        for j in lis:
            if j<=i:
                valst+=1
        try:
        	val_round=round(float(valst)/total_t,3)
        except:
                val_round=0
        string+='%s,%s,%s\n'%(i,valst,val_round)
    return string
'''
usage :
function will take in the list or array normalized from 0 to 1 in 'lis' parameter and columns as array of three values
like ['positionaffected','total','frequency'] as headers to csv file
will return the csv string writable to file
'''

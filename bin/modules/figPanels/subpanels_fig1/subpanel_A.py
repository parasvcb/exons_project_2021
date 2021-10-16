import figPanels.modules_common as cm
import os
def panel1_protein_length(has,res_dir):
    all_tup=[]
    has_bins={(0,100):0,(101,250):0,(251,500):0,(501,750):0,(751,1000):0,(1001,3000):0,(3001,7000):0,(7001,15000):0,(15001,37000):0}
    for i in has:
        #print (i,'key')
        
        j=has[i].PI
        all_tup+=[(has[i].ID,j.seqlen)]
        for key in has_bins:
            if key[0] <= j.seqlen <= key[1]:
                has_bins[key] += 1
                break
    values = has_bins.values()
    su = sum(values)
    keys = has_bins.keys()
    with open(os.path.join(res_dir,"F1_protein_length_variation_pi_histogramOrganized.csv"), "w") as fin:
        fin.write("Length_range,Total,Frequency\n")
        keys.sort()
        for i in keys:
            fin.write("%s,%s,%s\n" %("-".join(map(str,i)), has_bins[i], cm.div_fact(has_bins[i],su)))
    with open (os.path.join(res_dir,'F1_protein_length_variation_pi_raw'),'w') as fout:
        fout.write('Gene\tTranscript(RefIsf)\tRefIsfLength\n')
        for i in has:
            fout.write('%s\t%s\t%s\n'%(i,has[i].PI.ID,has[i].PI.seqlen))
#1st aspect: The above module takes the whole proteome and the tries to get the distribution of there principal isoforms

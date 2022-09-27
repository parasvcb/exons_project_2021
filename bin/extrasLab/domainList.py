import pandas as pd
import numpy as np
import re, scipy.stats, os
import cPickle as pickle
import sys
import figPanels.modules_analysis as am

print (len(sys.argv))
if len(sys.argv)!=3:
    print ('Please type 1. object 2. outputdir+append')
    sys.exit()

prog,source_gene_object,results_file =sys.argv
additionalFilter=False

def isoform_giver(gene_ob,fout):
    for gene in gene_ob:    
        for trans in gene_ob[gene].transcripts:
            domainTranscript=trans.pfam_list
            if domainTranscript:
                relev_doms=[]
                for i in domainTranscript:
                    fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(i[0][0],i[0][1],i[1],i[3],i[2],gene,trans.ID))
                    if i[2]>=0.5:
                            relev_doms+=[i]
                    '''
                    pfam_list looks like follwing
                    [((1, 64), 'PF02295.17', 0.96, 'z-alpha', 'Domain'), ((209, 274), 'PF00035.26', 0.99, 'dsrm', 'Domain'), ((320, 385), 'PF00035.26', 0.99, 'dsrm', 'Domain'), ((432, 497), 'PF00035.26', 0.99, 'dsrm', 'Domain'), ((591, 920), 'PF02137.18', 1.0, 'A_deamin', 'Family')]
                    ''' 

def loadpickle(fname):
    with open(fname,'rb') as fin:
        has=pickle.load(fin)
    return has

fout=open(results_file,'w')
fout.write('start\tend\tdomainID\domainNAME\tcoverageFromModel\tgene\ttranscript\n')
human= loadpickle(source_gene_object)
isoform_giver(human,fout)
fout.close()
#python bin/extrasLab/domainList.py output/derived_data/results/objectsave_9606_0.6.pick output/domainlist.tsv

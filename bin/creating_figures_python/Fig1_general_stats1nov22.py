'''
Before creation of object Builder type 2, the object needs pre filtering for the crrect depiction of the panels having keywords as AGF in exon nomenclature 
'''

import sys
import re, scipy, os
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import common.general_modules as cm

import figPanels.subpanels_fig1.subpanel_A as panelA
import figPanels.subpanels_fig1.subpanel_B as panelB
import figPanels.subpanels_fig1.subpanel_C as panelC
import figPanels.subpanels_fig1.subpanel_D as panelD
import figPanels.subpanels_fig1.subpanel_E as panelE
import figPanels.subpanels_fig1.subpanel_F as panelF

from constructing_data import Classes_exons

if len(sys.argv)!=3:
        print ("Please give 1 object 2 results dir")
        sys.exit()
prog,source_gene_object,results_dir_csv = sys.argv

humanGeneObject=cm.readPickle(source_gene_object)


def runnerWithGeneCondition(secDir,lengthThreshold=3000, UniqueTransCount=2,UniqueCodingExonCount=2): 
    if not os.path.isdir(secDir):
        os.makedirs(secDir)  
    filewriter=open(os.path.join(secDir,'F1_GeneralStatsGenePILength.log'),'w')
    CONDITION_GENES=panelB.transcripts_type(humanGeneObject,filewriter,secDir,lengthThreshold,UniqueTransCount, UniqueCodingExonCount)
    panelA.panel1_protein_length(humanGeneObject,secDir)
    panelC.exon_counting(humanGeneObject,secDir,CONDITION_GENES,filewriter)
    panelD.alternate_WEF(humanGeneObject,CONDITION_GENES,secDir,filewriter)
    panelE.exon_length(humanGeneObject,secDir, CONDITION_GENES,filewriter)
    panelF.change_in_length(humanGeneObject,secDir,CONDITION_GENES,filewriter)

    return CONDITION_GENES,filewriter

secDir2_2=os.path.join(results_dir_csv,'F1_Condition_%sPilength_%sIsf_%sExCount'%(3000,2,2))
secDir4_4=os.path.join(results_dir_csv,'F1_Condition_%sPilength_%sIsf_%sExCount'%(3000,4,4))
CONDITION_GENES_2_2, filewriter_2_2=runnerWithGeneCondition(secDir2_2,lengthThreshold=3000, UniqueTransCount=2,UniqueCodingExonCount=2)
CONDITION_GENES_4_4, filewriter_4_4=runnerWithGeneCondition(secDir4_4,lengthThreshold=3000, UniqueTransCount=4,UniqueCodingExonCount=4)

if not os.path.isfile(os.path.join(secDir2_2,"condition_genes.pick")):
    cm.writePickle(os.path.join(secDir2_2,"condition_genes.pick"),CONDITION_GENES_2_2)
if not os.path.isfile(os.path.join(secDir4_4,"condition_genes.pick")):
    cm.writePickle(os.path.join(secDir4_4,"condition_genes.pick"),CONDITION_GENES_4_4)

def glob_mem(has,cond,fout):
    has_prop={'Mem':0,'Glob':0,'Inter':0}
    for gene in has:
        if gene in cond:
            val= has[gene].localization[0]
            if val == 1:
                key="Mem"
            elif val == 0:
                key = "Glob"
            else:
                key = 'Inter'
            has_prop[key]+=1

    fout.write("\n\nIn the selected dataset follwing are the types of genes \n%s\n"%has_prop)

glob_mem(humanGeneObject,CONDITION_GENES_2_2,filewriter_2_2)
glob_mem(humanGeneObject,CONDITION_GENES_4_4,filewriter_4_4)

def principal_isoform(human,fout, cond):
    N=0;X=0
    largest=0;notlargest=0
    Ng=0;Xg=0
    largestG=0;notlargestG=0
    for gene in human:
        if human[gene].PI.ID[0]=='N':
            N+=1
        else:
            X+=1
        ma=max([i.seqlen for i in human[gene].transcripts])
        if human[gene].PI.seqlen==ma:
            largest+=1
        else:
            notlargest+=1
        if gene in cond:
            if human[gene].PI.ID[0]=='N':
                Ng+=1
            else:
                Xg+=1
            if human[gene].PI.seqlen==ma:
                largestG+=1
            else:
                notlargestG+=1
    fout.write("\nIN total set, %s genes have PI as reviewed and %s as model"%(cm.div_fact(N,(N+X)), cm.div_fact(float(X)/(N+X),2)))
    fout.write("\nIN analysed set, %s genes have PI as reviewed and %s as model"%(cm.div_fact(Ng,(Ng+Xg)),cm.div_fact(float(Xg)/(Ng+Xg),2)))
    fout.write("\nIN total set, %s genes have PI as longest"%(cm.div_fact(largest,(largest+notlargest))))
    fout.write("\nIN analysed set, %s genes have PI as longest\n"%(cm.div_fact(largestG,(largestG+notlargestG))))

principal_isoform(humanGeneObject,filewriter_2_2, CONDITION_GENES_2_2)
principal_isoform(humanGeneObject,filewriter_4_4, CONDITION_GENES_4_4)

# THREE CONDITION TO PUT IN PLACE, 
# AtLEAST TWO PROTEIN CODING ISOFORMS
# LENGTH <3000
# CODING_EXONS_MUTS BE atleast TWO always

filewriter_2_2.close()
filewriter_4_4.close()
'''
for removing the bias in analysis of laternate exon, choose only those AE which are prsent in PI
ADD ONE FUNCTION TO GIVE PI
Parent should be callable from Exon itself

python bin/Fig1_general_stats.py results/objectsave_9606_0.6.pick results/fig1/
'''
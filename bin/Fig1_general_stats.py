import sys
import re, scipy, os
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import modules_common as cm

import subpanels_fig1.subpanel_A as panelA
import subpanels_fig1.subpanel_B as panelB
import subpanels_fig1.subpanel_C as panelC
import subpanels_fig1.subpanel_D as panelD
import subpanels_fig1.subpanel_E as panelE

if len(sys.argv)!=3:
        print "Please give 1 object 2 results dir"
        sys.exit()
prog,source_gene_object,results_dir_csv = sys.argv

filewriter=open(results_dir_csv+'F1_General_stats1.log','w')
humanGeneObject=cm.readPickle(source_gene_object)

panelA.panel1_protein_length(humanGeneObject,results_dir_csv)
CONDITION_GENES=panelB.transcripts_type(humanGeneObject,filewriter,results_dir_csv)
panelC.exon_counting(humanGeneObject,results_dir_csv,CONDITION_GENES,filewriter)
panelD.alternate_WEF(humanGeneObject,CONDITION_GENES,results_dir_csv,filewriter)
panelE.exon_length(humanGeneObject,results_dir_csv, CONDITION_GENES,filewriter)
filewriter.close()
sys.exit()
if not os.path.isfile(results_dir_csv+"new_condition_genes.pick"):
    cm.writePickle(results_dir_csv+"new_condition_genes.pick",CONDITION_GENES)

def glob_mem(has,CONDITION_GENES,fout):
    has_prop={'Mem':0,'Glob':0,'Inter':0}
    for gene in has:
        if gene in CONDITION_GENES:
            val= has[gene].localization[0]
            if val == 1:
                key="Mem"
            elif val == 0:
                key = "Glob"
            else:
                key = 'Inter'
            has_prop[key]+=1

    fout.write("\n\nIn the selected dataset follwing are the types of genes \n%s\n"%has_prop)

glob_mem(has,CONDITION_GENES,filewriter)
            

def principal_isoform(human,fout):
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
        if gene in CONDITION_GENES:
            if human[gene].PI.ID[0]=='N':
                Ng+=1
            else:
                Xg+=1
            if human[gene].PI.seqlen==ma:
                largestG+=1
            else:
                notlargestG+=1
    fout.write("\nIN total set, %s genes have PI as reviewed and %s as model"%(div_fact(N,(N+X)), div_fact(float(X)/(N+X),2)))
    fout.write("\nIN analysed set, %s genes have PI as reviewed and %s as model"%(div_fact(Ng,(Ng+Xg)),div_fact(float(Xg)/(Ng+Xg),2)))
    fout.write("\nIN total set, %s genes have PI as longest"%(div_fact(largest,(largest+notlargest))))
    fout.write("\nIN analysed set, %s genes have PI as longest\n"%(div_fact(largestG,(largestG+notlargestG))))

principal_isoform(has,filewriter)

# THREE CONDITION TO PUT IN PLACE, 
# AtLEAST TWO PROTEIN CODING ISOFORMS
# LENGTH <3000
# CODING_EXONS_MUTS BE atleast TWO always


'''
for removing the bias in analysis of laternate exon, choose only those AE which are prsent in PI
ADD ONE FUNCTION TO GIVE PI
Parent should be callable from Exon itself
'''
python bin/Fig1_general_stats.py results/objectsave_9606_0.6.pick results/fig1/
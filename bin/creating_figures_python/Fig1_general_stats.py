'''
Before creation of object Builder type 2, the object needs pre filtering for the crrect depiction of the panels having keywords as AGF in exon nomenclature 
DONE above

sets the file names in advance
    creates the conditional objects
    

'''

import sys
import os
import matplotlib.pyplot as plt
import common.general_modules as cm
from progress.bar import Bar
import figPanels.modules_analysis as ca

import figPanels.subpanels_fig1.subpanel_A as panelA
import figPanels.subpanels_fig1.subpanel_B as panelB
import figPanels.subpanels_fig1.subpanel_C as panelC
# import figPanels.subpanels_fig1.subpanel_D as panelD
import figPanels.subpanels_fig1.subpanel_D_advancedExonSelection_dontuse as panelD
import figPanels.subpanels_fig1.subpanel_E as panelE
import figPanels.subpanels_fig1.subpanel_F as panelF

from constructing_data import Classes_exons

if len(sys.argv)!=5:
        print ("Please give 1 object 2 org 2 alnFile 3 results dir")
        sys.exit()
prog,source_gene_object, org, alnFile, results_dir_csv = sys.argv

humanGeneObject=cm.readPickle(source_gene_object)

def alignmentCallerRead(alnFile):
    preDoneHas = {}
    preDoneHas1 = {}
    if os.path.isfile(alnFile):
        with open (alnFile) as fin:
            # gene   utmrd0  utmrd2  type	exon1	exon2   aln1	aln2	aln3	statsLis
            # 103   D   A   default exonStr
            # 103   D   G   default exonStr
            # 104   D   G   NA
            dat=[i for i in fin.read().split('\n') if len(i)>2]
            if len(dat)>1:
                for i in dat[1:]:
                    ele = i.split('\t')
                    if len(ele)>=4:
                        gene, utmrd0, utmrd2, type = ele[:4]
                        if gene not in preDoneHas:
                            preDoneHas[gene]={}
                        if (utmrd0, utmrd2) not in preDoneHas[gene]:
                            preDoneHas[gene][(utmrd0, utmrd2)]=[[],{}]
                        if len(ele) == 5 and ele[3] == 'default':
                            gene, utmrd0, utmrd2, type, placeholderTags = ele
                            placeholderList = list(map(int,placeholderTags.split(',')))
                            preDoneHas[gene][(utmrd0, utmrd2)][0]=placeholderList
                        elif len(ele) == 10:
                            # 2       T       G       assign  5:0:0   2       VKFRVVSMDENFHPLNEL      -------MDENFHPLNEL             |||||||||||      18,11,0.611,1.0,0.611,1.0
                            gene, utmrd0, utmrd2, type, exon1, exon2, aln1, aln2, aln3, statLis = i.split('\t')
                            seq1Length, seq2Length, identity1, identity2, cov1, cov2 = list(map(float, statLis.split(',')))
                            preDoneHas[gene][(utmrd0, utmrd2)][1][(exon1,exon2)]=[aln1,aln2,aln3, [seq1Length, seq2Length, identity1, identity2, cov1, cov2]]
        preDoneHas1 = {int(i): preDoneHas[i] for i in preDoneHas}
    return preDoneHas1

def alnRecapturer(ob,fname):
    hasDone = alignmentCallerRead(fname)
    if len(hasDone)!=len(ob):
        bar = Bar('Processing genes for alignment MAIN():', max=len(ob))
        fhandler = open (fname,'w')
        fhandler.write('gene\tutmrd0\tutmrd2\ttype\texon1\texon2\taln1\taln2\taln3\tstatsLis\n')
        for gene in ob:
            exonMatrix = ca.positionalExonMatrix_forExCharacterization(ob[gene])
            ca.giveExonSelectionAdvanced(exonMatrix, ob[gene].exons, gene, fhandler)
            bar.next()
        bar.finish()
        fhandler.close()
        hasDone = alignmentCallerRead(fname)
    #print (hasDone[79955])
    #sys.exit()
    return hasDone
    #sys.exit()

def runnerWithGeneCondition(secDir, org, lengthThreshold=3000, UniqueTransCount=2, UniqueCodingExonCount=2 ): 
    if not os.path.isdir(secDir):
        os.makedirs(secDir)  
    filewriter=open(os.path.join(secDir,'F1_GeneralStatsGenePILength.log'),'w')
    CONDITION_GENES=panelB.transcripts_type(humanGeneObject,filewriter,secDir,lengthThreshold,UniqueTransCount, UniqueCodingExonCount)
    hasAln = alnRecapturer(humanGeneObject, alnFile)
    panelA.panel1_protein_length(humanGeneObject,secDir)
    panelC.exon_counting(humanGeneObject, org, UniqueTransCount, hasAln, secDir,CONDITION_GENES,filewriter)
    panelD.alternate_WEF(humanGeneObject,CONDITION_GENES,secDir,filewriter)
    panelE.exon_length(humanGeneObject,hasAln, secDir, CONDITION_GENES,filewriter)
    panelF.change_in_length(humanGeneObject,secDir,CONDITION_GENES,filewriter)

    return CONDITION_GENES,filewriter,hasAln

secDir2_2=os.path.join(results_dir_csv,'F1_Condition_%sPilength_%sIsf_%sExCount'%(3000,2,2))
secDir4_4=os.path.join(results_dir_csv,'F1_Condition_%sPilength_%sIsf_%sExCount'%(3000,4,4))
CONDITION_GENES_2_2, filewriter_2_2, hasAln2 =runnerWithGeneCondition(secDir2_2, org, lengthThreshold=3000, UniqueTransCount=2,UniqueCodingExonCount=2)
CONDITION_GENES_4_4, filewriter_4_4, hasAln4 =runnerWithGeneCondition(secDir4_4, org, lengthThreshold=3000, UniqueTransCount=4,UniqueCodingExonCount=4)

if not os.path.isfile(os.path.join(secDir2_2,"condition_genes.pick")) or 1:
    cm.writePickle(os.path.join(secDir2_2,"condition_genes.pick"),CONDITION_GENES_2_2)
if not os.path.isfile(os.path.join(secDir4_4,"condition_genes.pick")) or 1:
    cm.writePickle(os.path.join(secDir4_4,"condition_genes.pick"),CONDITION_GENES_4_4)


def getFrameChGenes(cond,has,fname):
    no=[]
    assign =[]
    frame = []
    for gene in has:
        if gene in cond:
            tno=[]
            for tupkey in has[gene]:
                if not has[gene][tupkey][0]:
                    tno+=[tupkey]
                else:
                    for exonvar in has[gene][tupkey][1]:
                        aln1,aln2,aln3, statlis = has[gene][tupkey][1][exonvar]
                        seq1Length, seq2Length, identity1, identity2, cov1, cov2 = statlis
                        if identity1 == 1 or identity2 ==1:
                            assign += [gene]
                        else:
                            frame += [gene]
            if len(tno)==len(has[gene]):
                no+=[gene]
    assign = set(assign)
    frame = set(frame)
    assignANDframe = assign | frame
    assignBUTNOframe= assign - frame
    frameBUTNOassign= frame - assign
    assignORframe= frame | assign
    hasWrite = {
        'noAaChangeGeneList':no,
        'L_assignANDframeINTERSECTINGGeneList':assignANDframe,
        'L_assignBUTNOframeGeneList':assignBUTNOframe,
        'L_frameBUTNOassignGeneList':frameBUTNOassign,
        'L_assignORframeALLGeneList':assignORframe
        }
    for i in hasWrite:
        with open (fname+i,'w') as fout:
            for j in hasWrite[i]:
                fout.write('%s\n'%j)

getFrameChGenes(CONDITION_GENES_2_2,hasAln2,os.path.join(results_dir_csv,'_2X2_') )
getFrameChGenes(CONDITION_GENES_4_4,hasAln2,os.path.join(results_dir_csv,'_4X4_') )

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
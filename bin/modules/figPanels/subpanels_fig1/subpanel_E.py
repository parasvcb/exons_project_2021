import re,os
# import figPanels.modules_common as cm
import common.general_modules as cm
import figPanels.modules_analysis as ca
from progress.bar import Bar

def constitutive_alternate_with_freq(geneob, gene, hasAln):
    CodingStrict = []
    CodingMajorlyConstitutive = []
    CodingMajorlyWaachange = []
    CodingMajorlyWaachangeAssign = []
    CodingMajorlyWaachangeFrame = [] 
    CodingAlternate = []
    CodingAlternateWithSS = []
    CodingConstitutive = []
    CodingConstitutiveWaachange = []
    CodingConstitutiveWaachangeAssign = []
    CodingConstitutiveWaachangeFrame = []
    CodingAltWaachange = []
    CodingAltWaachangeAssign = []
    CodingAltWaachangeFrame = []

    exonMatrix = ca.positionalExonMatrix_forExCharacterization(geneob)
            #{1: ['U', 'F', [2, 0, 2], 0, 1, 0, 0], 2: ['U', 'A', [0, 0, 0], 0, 1, 0, 0], 3: ['T', 'A', [0, 0, 0], 0, 0, 0, 0], 4: ['T', 'F', [0, 1, 0], 0, 0, 0, 0]}

    temp1, temp2, temp3, temp4, temp5 = ca.giveExonSelectionBasic(exonMatrix,"T")
    pos_CodingStrict = temp1
    pos_CodingMajorlyConstitutive = temp2
    pos_CodingAlternate = temp3
    pos_CodingAlternateWithSS = temp4
    pos_CodingConstitutive = temp5

    def AssignerFramerdefault (tupkey,hasFramechAssign):
        if hasFramechAssign[gene][tupkey][0]:
            l1 = hasFramechAssign[gene][tupkey][0]
            if hasFramechAssign[gene][tupkey][1]:
                assign = []
                frame = []
                for exonAndVar in hasFramechAssign[gene][tupkey][1]:
                    aln1,aln2,aln3, statlis = hasFramechAssign[gene][tupkey][1][exonAndVar]
                    seq1Length, seq2Length, identity1, identity2, cov1, cov2 = statlis
                    if identity1 == 1 or identity2 ==1:
                        assign += [int(exonAndVar[0].split(':')[0])]
                    else:
                        frame += [int(exonAndVar[0].split(':')[0])]
                l2 = list(set(assign))
                l3 = list(set(frame))
            else:
                l2 = []
                l3 = []
        else:
            l1= []
            l2= []
            l3= []
        return l1, l2, l3
    
    pos_CodingConstitutiveWaachange, pos_CodingConstitutiveWaachangeAssign, pos_CodingConstitutiveWaachangeFrame = AssignerFramerdefault(('T','G'), hasAln)
    pos_CodingAltWaachange, pos_CodingAltWaachangeAssign, pos_CodingAltWaachangeFrame = AssignerFramerdefault(('T','A'), hasAln)
    pos_CodingMajorlyWaachange, pos_CodingMajorlyWaachangeAssign, pos_CodingMajorlyWaachangeFrame = AssignerFramerdefault(('T','F'), hasAln)
    
    # temp1, temp2, temp3 = ca.giveExonSelectionAdvanced(exonMatrix, 'T', geneob.exons, 'G', gene)
    # pos_CodingConstitutiveWaachange = temp1 
    # pos_CodingConstitutiveWaachangeAssign = temp2
    # pos_CodingConstitutiveWaachangeFrame = temp3
    
    # temp1, temp2, temp3 = ca.giveExonSelectionAdvanced(exonMatrix, 'T', geneob.exons, 'A', gene)
    # pos_CodingAltWaachange = temp1 
    # pos_CodingAltWaachangeAssign = temp2
    # pos_CodingAltWaachangeFrame = temp3 

    for pos in pos_CodingStrict:
        tlis = []
        for i in geneob.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingStrict += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingMajorlyConstitutive:
        tlis = []
        for i in geneob.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingMajorlyConstitutive += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingAlternate:
        tlis = []
        for i in geneob.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingAlternate += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingAlternateWithSS:
        tlis = []
        for i in geneob.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingAlternateWithSS += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingConstitutive:
        tlis = []
        for i in geneob.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingConstitutive += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingConstitutiveWaachange:
        tlis = []
        for i in geneob.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingConstitutiveWaachange += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingConstitutiveWaachangeAssign:
        tlis = []
        for i in geneob.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingConstitutiveWaachangeAssign += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingConstitutiveWaachangeFrame:
        tlis = []
        for i in geneob.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingConstitutiveWaachangeFrame += [round(float(sum(tlis))/len(tlis),3)]


    for pos in pos_CodingAltWaachange:
        tlis = []
        for i in geneob.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingAltWaachange += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingAltWaachangeAssign:
        tlis = []
        for i in geneob.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingAltWaachangeAssign += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingAltWaachangeFrame:
        tlis = []
        for i in geneob.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingAltWaachangeFrame += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingMajorlyWaachange:
        tlis = []
        for i in geneob.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingMajorlyWaachange += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingMajorlyWaachangeAssign:
        tlis = []
        for i in geneob.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingMajorlyWaachangeAssign += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingMajorlyWaachangeFrame:
        tlis = []
        for i in geneob.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingMajorlyWaachangeFrame += [round(float(sum(tlis))/len(tlis),3)]
    
    
    # i will only be considering positions which are T, though child forms may adapt different forms, lets consider them also
    # categories would be 
    # T without distinction, takle average of all the forms other than R for such events 
    # T with F cases, will hacve splice site changes also, lste cosnider them also, but only the coding form, not the UTR form
    # T with A but no ss, 
    # T with A but ss
    # T with G
    # T with G but change in aa
    
    return CodingStrict, CodingMajorlyConstitutive, CodingAlternate, CodingAlternateWithSS, CodingConstitutive, CodingConstitutiveWaachange, CodingConstitutiveWaachangeFrame, CodingConstitutiveWaachangeAssign, CodingAltWaachange, CodingAltWaachangeFrame, CodingAltWaachangeAssign, CodingMajorlyWaachange, CodingMajorlyWaachangeFrame, CodingMajorlyWaachangeAssign 

def exon_length(has,hasAln, res_dir,genes_cond,output_stats):
    A_CodingStrict = []
    A_CodingMajorlyConstitutive = []
    A_CodingAlternate = []
    A_CodingAlternateWithSS = []
    A_CodingConstitutive = []
    A_CodingConstitutiveWaachange = []
    A_CodingConstitutiveWaachangeFrame =[]
    A_CodingConstitutiveWaachangeAssign = []
    A_CodingAltWaachange = []
    A_CodingAltWaachangeFrame =[]
    A_CodingAltWaachangeAssign =[]
    A_CodingMajorlyWaachange = []
    A_CodingMajorlyWaachangeFrame =[]
    A_CodingMajorlyWaachangeAssign =[]
    bar = Bar('Processing genes E:', max=len(genes_cond))
    for gene in has:
        if gene in genes_cond:
            CodingStrict, CodingMajorlyConstitutive, CodingAlternate, CodingAlternateWithSS, CodingConstitutive, CodingConstitutiveWaachange, CodingConstitutiveWaachangeFrame, CodingConstitutiveWaachangeAssign, CodingAltWaachange, CodingAltWaachangeFrame, CodingAltWaachangeAssign, CodingMajorlyWaachange, CodingMajorlyWaachangeFrame, CodingMajorlyWaachangeAssign = constitutive_alternate_with_freq(has[gene], gene, hasAln)
            A_CodingStrict += CodingStrict
            A_CodingMajorlyConstitutive += CodingMajorlyConstitutive
            A_CodingAlternate += CodingAlternate
            A_CodingAlternateWithSS += CodingAlternateWithSS
            A_CodingConstitutive += CodingConstitutive
            A_CodingConstitutiveWaachange += CodingConstitutiveWaachange
            A_CodingConstitutiveWaachangeFrame += CodingConstitutiveWaachangeFrame
            A_CodingConstitutiveWaachangeAssign += CodingConstitutiveWaachangeAssign
            A_CodingAltWaachange += CodingAltWaachange
            A_CodingAltWaachangeFrame += CodingAltWaachangeFrame
            A_CodingAltWaachangeAssign += CodingAltWaachangeAssign
            A_CodingMajorlyWaachange += CodingMajorlyWaachange
            A_CodingMajorlyWaachangeFrame += CodingMajorlyWaachangeFrame
            A_CodingMajorlyWaachangeAssign += CodingMajorlyWaachangeAssign
            bar.next()
    bar.finish()

    tabularLength = "Category\tVar_count\tSum\tMax\tMin\tMean\tMedian\tMode\tMode_count\tStddev\n"
    tabularLength += "0_CodingStrict\t%s\n"%cm.stats(A_CodingStrict)
    tabularLength += "1_CodingAlternate\t%s\n"%cm.stats(A_CodingAlternate)
    tabularLength += "1_CodingAlternateWithSS\t%s\n"%cm.stats(A_CodingAlternateWithSS)
    tabularLength += "1_CodingAltWaachange\t%s\n"%cm.stats(A_CodingAltWaachange)
    tabularLength += "1_CodingAltWaachangeAssign\t%s\n"%cm.stats(A_CodingAltWaachangeAssign)
    tabularLength += "1_CodingAltWaachangeFrame\t%s\n"%cm.stats(A_CodingAltWaachangeFrame)
    
    tabularLength += "2_CodingMajorlyConstitutive\t%s\n"%cm.stats(A_CodingMajorlyConstitutive)
    tabularLength += "2_CodingMajorlyWaachange\t%s\n"%cm.stats(A_CodingMajorlyWaachange)
    tabularLength += "2_CodingMajorlyWaachangeAssign\t%s\n"%cm.stats(A_CodingMajorlyWaachangeAssign)
    tabularLength += "2_CodingMajorlyWaachangeFrame\t%s\n"%cm.stats(A_CodingMajorlyWaachangeFrame)
    
    tabularLength += "3_CodingConstitutive\t%s\n"%cm.stats(A_CodingConstitutive)
    tabularLength += "3_CodingConstitutiveWaachange\t%s\n"%cm.stats(A_CodingConstitutiveWaachange)
    tabularLength += "3_CodingConstitutiveWaachangeAssign\t%s\n"%cm.stats(A_CodingConstitutiveWaachangeAssign)
    tabularLength += "3_CodingConstitutiveWaachangeFrame\t%s\n"%cm.stats(A_CodingConstitutiveWaachangeFrame)
    
    with open (os.path.join(res_dir,"tabular_ExonLengthData.csv"),'w') as outf:
        outf.write("%s"%tabularLength)

    output_stats.write("\n\n=============>\n")
    output_stats.write("Length Comparison of constituive and alternate exons\n%s"%tabularLength)

    output_stats.write("\nCSV_File:%s"%res_dir+"length_distribution_exons.csv")

    with open(os.path.join(res_dir,'panelE_Length_distribution_exonsRAW.csv'),"w") as fin:
        fin.write("Length,Category\n")
        for i in A_CodingStrict:
            fin.write("%s,CodingStrict\n"%(i))
        for i in A_CodingAlternate:
            fin.write("%s,CodingAlternate\n"%(i))
        for i in A_CodingAlternateWithSS:
            fin.write("%s,CodingAlternateWithSS\n"%(i))
        for i in A_CodingMajorlyConstitutive:
            fin.write("%s,CodingMajorlyConstitutive\n"%(i))
        for i in A_CodingMajorlyWaachange:
            fin.write("%s,CodingMajorlyWaachange\n"%(i))
        for i in A_CodingMajorlyWaachangeFrame:
            fin.write("%s,CodingMajorlyWaachangeFrame\n"%(i))
        for i in A_CodingMajorlyWaachangeAssign:
            fin.write("%s,CodingMajorlyWaachangeAssign\n"%(i))
        for i in A_CodingConstitutive:
            fin.write("%s,CodingConstitutive\n"%(i))
        for i in A_CodingConstitutiveWaachange:
            fin.write("%s,CodingConstitutiveWaachange\n"%(i))
        for i in A_CodingConstitutiveWaachangeFrame:
            fin.write("%s,CodingConstitutiveWaachangeFrame\n"%(i))
        for i in A_CodingConstitutiveWaachangeAssign:
            fin.write("%s,CodingConstitutiveWaachangeAssign\n"%(i))
        for i in A_CodingAltWaachange:
            fin.write("%s,CodingAltWaachange\n"%(i))
        for i in A_CodingAltWaachangeFrame:
            fin.write("%s,CodingAltWaachangeFrame\n"%(i))
        for i in A_CodingAltWaachangeAssign:
            fin.write("%s,CodingAltWaachangeAssign\n"%(i))
        

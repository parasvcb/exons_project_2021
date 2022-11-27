import re,os
# import figPanels.modules_common as cm
import common.general_modules as cm
import figPanels.modules_analysis as ca

def constitutive_alternate_with_freq(gene):
    CodingStrict = []
    CodingMajorlyConstitutive = []
    CodingAlternate = []
    CodingAlternateWithSS = []
    CodingConstitutive = []
    CodingConstitutiveWaachange = []

    exonMatrix = ca.positionalExonMatrix_forExCharacterization(gene)
            #{1: ['U', 'F', [2, 0, 2], 0, 1, 0, 0], 2: ['U', 'A', [0, 0, 0], 0, 1, 0, 0], 3: ['T', 'A', [0, 0, 0], 0, 0, 0, 0], 4: ['T', 'F', [0, 1, 0], 0, 0, 0, 0]}

    pos_CodingStrict=[i for i in exonMatrix if exonMatrix[i][0]=='T']
    pos_CodingMajorlyConstitutive=[i for i in exonMatrix if exonMatrix[i][0]=='T' and exonMatrix[i][1]=='F']   
    pos_CodingAlternate=[i for i in exonMatrix if exonMatrix[i][0]=='T' and exonMatrix[i][1]=='A' and sum(exonMatrix[i][2]) ==0]
    pos_CodingAlternateWithSS=[i for i in exonMatrix if exonMatrix[i][0]=='T' and exonMatrix[i][1]=='A' and sum(exonMatrix[i][2])>0]  
    pos_CodingConstitutive=[i for i in exonMatrix if exonMatrix[i][0]=='T' and exonMatrix[i][1]=='G']
    pos_CodingConstitutiveWaachange=[i for i in exonMatrix if exonMatrix[i][0]=='T' and exonMatrix[i][1]=='G' and exonMatrix[i][3]>0]

    for pos in pos_CodingStrict:
        tlis = []
        for i in gene.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingStrict += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingMajorlyConstitutive:
        tlis = []
        for i in gene.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingMajorlyConstitutive += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingAlternate:
        tlis = []
        for i in gene.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingAlternate += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingAlternateWithSS:
        tlis = []
        for i in gene.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingAlternateWithSS += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingConstitutive:
        tlis = []
        for i in gene.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingConstitutive += [round(float(sum(tlis))/len(tlis),3)]
    
    for pos in pos_CodingConstitutiveWaachange:
        tlis = []
        for i in gene.exons:
            if i.ID[0]!='R' and i.ID[0] == 'T' and int(i.ID.split('.')[3])  == pos and i.length >0:
                tlis+=[i.length]
        if tlis:
            CodingConstitutiveWaachange += [round(float(sum(tlis))/len(tlis),3)]
    
    # i will only be considering positions which are T, though child forms may adapt different forms, lets consider them also
    # categories would be 
    # T without distinction, takle average of all the forms other than R for such events 
    # T with F cases, will hacve splice site changes also, lste cosnider them also, but only the coding form, not the UTR form
    # T with A but no ss, 
    # T with A but ss
    # T with G
    # T with G but change in aa
    
    return CodingStrict, CodingMajorlyConstitutive, CodingAlternate, CodingAlternateWithSS, CodingConstitutive, CodingConstitutiveWaachange

def exon_length(has,res_dir,genes_cond,output_stats):
    A_CodingStrict = []
    A_CodingMajorlyConstitutive = []
    A_CodingAlternate = []
    A_CodingAlternateWithSS = []
    A_CodingConstitutive = []
    A_CodingConstitutiveWaachange = []

    for gene in has:
        if gene in genes_cond:
            CodingStrict, CodingMajorlyConstitutive, CodingAlternate, CodingAlternateWithSS, CodingConstitutive, CodingConstitutiveWaachange = constitutive_alternate_with_freq(has[gene])
            A_CodingStrict += CodingStrict
            A_CodingMajorlyConstitutive += CodingMajorlyConstitutive
            A_CodingAlternate += CodingAlternate
            A_CodingAlternateWithSS += CodingAlternateWithSS
            A_CodingConstitutive += CodingConstitutive
            A_CodingConstitutiveWaachange += CodingConstitutiveWaachange

    tabularLength = "Category\tVar_count\tSum\tMax\tMin\tMean\tMedian\tMode\tMode_count\tStddev\n"
    tabularLength += "0_CodingStrict\t%s\n"%cm.stats(A_CodingStrict)
    tabularLength += "1_CodingAlternate\t%s\n"%cm.stats(A_CodingAlternate)
    tabularLength += "1_CodingAlternateWithSS\t%s\n"%cm.stats(A_CodingAlternateWithSS)
    tabularLength += "2_CodingMajorlyConstitutive\t%s\n"%cm.stats(A_CodingMajorlyConstitutive)
    tabularLength += "3_CodingConstitutive\t%s\n"%cm.stats(A_CodingConstitutive)
    tabularLength += "3_CodingConstitutiveWaachange\t%s\n"%cm.stats(A_CodingConstitutiveWaachange)

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
        for i in A_CodingConstitutive:
            fin.write("%s,CodingConstitutive\n"%(i))
        for i in A_CodingConstitutiveWaachange:
            fin.write("%s,CodingConstitutiveWaachange\n"%(i))
        
        

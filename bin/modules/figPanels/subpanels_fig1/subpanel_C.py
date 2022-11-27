
'''
This block will read the genes in object, will filter them on the basis of atleast two protein coding different siforms
atleast two protein coding exons and length of pricinpal soform less than 3000 aa
when that criteria is met, it will ask to change the parameter valeu from unique to varition 
that then will calculate the listed cases on the basis of the unique genomic coods as once instance 
and in latter to count them on the basis of the variation also in those cases

correspsonding data is written into the files also and dispaleys on screen with genral statustucs
'''
import re,os,sys
# import figPanels.modules_common as cm
import common.general_modules as cm
import figPanels.modules_analysis as ca

def writer_and_displayer(res_dir,filename,fout,hasPacked):
    fout.write(filename)
    fout.write("\n=================\nUnique exons count based on genomic coordinates\n")

    tabularData = "Category\tGeneCount\tTotal\tMax\tMin\tMean\tMedian\tMode\tMode_count\tStddev\n"
    tabularData += "Total_exons\t%s\n"%cm.stats(hasPacked['tot_exons'])
    tabularData += "Total_exons_ncb&R\t%s\n"%cm.stats(hasPacked['tot_exonsWncbAndRet'])

    tabularData += "UTRStrict\t%s\n"%cm.stats(hasPacked['UTRStrict'])
    tabularData += "UTRAlternate\t%s\n"%cm.stats(hasPacked['UTRAlternate'])
    tabularData += "UTRAlternateWithSS\t%s\n"%cm.stats(hasPacked['UTRAlternateWithSS'])
    tabularData += "UTRMajorlyConstitutive\t%s\n"%cm.stats(hasPacked['UTRMajorlyConstitutive'])
    tabularData += "UTRConstitutive\t%s\n"%cm.stats(hasPacked['UTRConstitutive'])
    
    tabularData += "CodingStrict\t%s\n"%cm.stats(hasPacked['CodingStrict'])
    tabularData += "CodingAlternate\t%s\n"%cm.stats(hasPacked['CodingAlternate'])
    tabularData += "CodingAltWaachange\t%s\n"%cm.stats(hasPacked['CodingAltWaachange'])
    tabularData += "CodingAlternateWithSS\t%s\n"%cm.stats(hasPacked['CodingAlternateWithSS'])
    tabularData += "CodingMajorlyConstitutive\t%s\n"%cm.stats(hasPacked['CodingMajorlyConstitutive'])
    tabularData += "CodingConstitutive\t%s\n"%cm.stats(hasPacked['CodingConstitutive'])
    tabularData += "CodingConstitutiveWaachange\t%s\n"%cm.stats(hasPacked['CodingConstitutiveWaachange'])
    
    tabularData += "DualStrict\t%s\n"%cm.stats(hasPacked['DualStrict'])
    tabularData += "DualAlternate\t%s\n"%cm.stats(hasPacked['DualAlternate'])
    tabularData += "DualAltWaachange\t%s\n"%cm.stats(hasPacked['DualAltWaachange'])
    tabularData += "DualAlternateWithSS\t%s\n"%cm.stats(hasPacked['DualAlternateWithSS'])
    tabularData += "DualMajorlyConstitutive\t%s\n"%cm.stats(hasPacked['DualMajorlyConstitutive'])
    tabularData += "DualConstitutive\t%s\n"%cm.stats(hasPacked['DualConstitutive'])
    tabularData += "DualConstitutiveWaachange\t%s\n"%cm.stats(hasPacked['DualConstitutiveWaachange'])
    
    tabularData += "Intron_retained_coding\t%s\n"%cm.stats(hasPacked['ir_retained_coding'])
    tabularData += "Intron_retained_UTR\t%s\n"%cm.stats(hasPacked['ir_retained_UTR'])
    tabularData += "StrictMcases(1ntLength)\t%s\n"%cm.stats(hasPacked['strict_M'])
    tabularData += "Strict_aa_removed\t%s\n"%cm.stats(hasPacked['strict_aa_removed'])

    fout.write( "%s\n"%tabularData)

    with open (os.path.join(res_dir,"tabular_ExonsData.csv"),'w') as outf:
        outf.write("%s"%tabularData)
    
    fout.write("\nCSV_file:%s"%res_dir+filename+".csv")
    
    with open(os.path.join(res_dir,filename+".csv"),"w") as fin:
        fin.write("Exons,Category,Type\n")
        for i in hasPacked['tot_exons']:
            fin.write("%s,TotalExons,Extras\n"%(i))
        for i in hasPacked['tot_exonsWncbAndRet']:
            fin.write("%s,TotalExonsWncb&r,Extras\n"%(i))
        for i in hasPacked['ir_retained_UTR']:
            fin.write("%s,IR_UTR,Extras\n"%(i))
        for i in hasPacked['ir_retained_coding']:
            fin.write("%s,IR_coding,Extras\n"%(i))
        for i in hasPacked['strict_aa_removed']:
            fin.write("%s,DeprivedAA,Extras\n"%(i))
        for i in hasPacked['strict_M']:
            fin.write("%s,MCases,Extras\n"%(i))
        
        for i in hasPacked['UTRStrict']:
            fin.write("%s,Strict,UTR\n"%(i))
        for i in hasPacked['UTRAlternate']:
            fin.write("%s,Alternate,UTR\n"%(i))
        for i in hasPacked['UTRAlternateWithSS']:
            fin.write("%s,AltWithSSchange,UTR\n"%(i))
        for i in hasPacked['UTRMajorlyConstitutive']:
            fin.write("%s,MajorlyConst,UTR\n"%(i))
        for i in hasPacked['UTRConstitutive']:
            fin.write("%s,Const,UTR\n"%(i))

        for i in hasPacked['CodingStrict']:
            fin.write("%s,Strict,Coding\n"%(i))
        for i in hasPacked['CodingAlternate']:
            fin.write("%s,Alternate,Coding\n"%(i))
        for i in hasPacked['CodingAltWaachange']:
            fin.write("%s,AlternateWaaChange,Coding\n"%(i))
        for i in hasPacked['CodingAlternateWithSS']:
            fin.write("%s,AltWithSSchange,Coding\n"%(i))
        for i in hasPacked['CodingMajorlyConstitutive']:
            fin.write("%s,MajorlyConst,Coding\n"%(i))
        for i in hasPacked['CodingConstitutive']:
            fin.write("%s,Const,Coding\n"%(i))
        for i in hasPacked['CodingConstitutiveWaachange']:
            fin.write("%s,ConstWaaChange,Coding\n"%(i))
        
        for i in hasPacked['DualStrict']:
            fin.write("%s,Strict,Dual\n"%(i))
        for i in hasPacked['DualAlternate']:
            fin.write("%s,Alternate,Dual\n"%(i))
        for i in hasPacked['DualAltWaachange']:
            fin.write("%s,AlternateWaaChange,Dual\n"%(i))
        for i in hasPacked['DualAlternateWithSS']:
            fin.write("%s,AltWithSSchange,Dual\n"%(i))
        for i in hasPacked['DualMajorlyConstitutive']:
            fin.write("%s,MajorlyConst,Dual\n"%(i))
        for i in hasPacked['DualConstitutive']:
            fin.write("%s,Const,Dual\n"%(i))
        for i in hasPacked['DualConstitutiveWaachange']:
            fin.write("%s,ConstWaaChange,Dual\n"%(i))

def exon_counting(has,res_dir,condition,fout ):
    # newly added categories ################
    # 04 -11 -22
    tot_exonsWncbAndRet = []
    DualConstitutiveWaachange = []
    DualAltWaachange = []
    CodingConstitutiveWaachange = []
    CodingAltWaachange = []
    #########################################
    tot_exons=[]
    UTRStrict=[]
    UTRAlternate=[]
    UTRAlternateWithSS=[]
    UTRMajorlyConstitutive=[]
    UTRConstitutive=[]
    '''
    strict is one without segregation and subtypes, 
        alternate is the one with A,
            will have splice sites 
        constitutive with G
        majorly consective is F
    '''
    CodingStrict = []
    CodingAlternate = []
    CodingAlternateWithSS=[]
    CodingMajorlyConstitutive=[]
    CodingConstitutive=[]

    DualStrict=[]
    DualAlternate=[]
    DualAlternateWithSS=[]
    DualMajorlyConstitutive=[]
    DualConstitutive=[]

    strict_M=[]
    strict_aa_removed=[]    
    ir_retained_coding=[]
    ir_retained_UTR=[]

    PARAMETER=""
    #PARAMETER="Variation"
    filename="Variation_allowed" if PARAMETER == "Variation" else "Unique_Coods"
    output_stats = open("General_stats_fig1.log", "a")
    # A,B,C,D,E=[0,0,0,0,0]
    for gene in has:
        if gene in condition:
        #if gene==319701:
            '''
            dividing in three categories
            1st [Strict UTR regions, Coding regions, Dual regions, M cases].
                2nd child elements, [constitutive, majorly constitutive, alternate pure cases, alternate which undergoes splice sites change and if retention starts from this exon]
                ## still cond=fused about this event, retention from 1 to 3 , starts from 1 but ends at three, shall i gave retention to exons 1-3 or only consider 1 to be starting exon and hence being retained ?
                lets stick with the former approach.
                consider aa change frequency along with 0-1 
            '''
            exonMatrix = ca.positionalExonMatrix_forExCharacterization(has[gene])
            #{1: ['U', 'F', [2, 0, 2], 0, 1, 0, 0], 2: ['U', 'A', [0, 0, 0], 0, 1, 0, 0], 3: ['T', 'A', [0, 0, 0], 0, 0, 0, 0], 4: ['T', 'F', [0, 1, 0], 0, 0, 0, 0]}

            tot_exons+=[len(exonMatrix)]
            tot_exonsWncbAndRet+=[len(exonMatrix)+sum([exonMatrix[i][2][0]+exonMatrix[i][2][1]+exonMatrix[i][2][2]+exonMatrix[i][4]+exonMatrix[i][5] for i in exonMatrix])]

            UTRStrict+=[len([i for i in exonMatrix if exonMatrix[i][0]=='U'])]
            UTRMajorlyConstitutive+=[len([i for i in exonMatrix if exonMatrix[i][0]=='U' and exonMatrix[i][1]=='F'])]   
            UTRAlternate+=[len([i for i in exonMatrix if exonMatrix[i][0]=='U' and exonMatrix[i][1]=='A' and sum(exonMatrix[i][2]) ==0])]
            UTRAlternateWithSS+=[len([i for i in exonMatrix if exonMatrix[i][0]=='U' and exonMatrix[i][1]=='A' and sum(exonMatrix[i][2])>0])]  
            UTRConstitutive+=[len([i for i in exonMatrix if exonMatrix[i][0]=='U' and exonMatrix[i][1]=='G'])]
            
            # A+=len([i for i in exonMatrix if exonMatrix[i][0]=='U'])
            # B+=len([i for i in exonMatrix if exonMatrix[i][0]=='U' and exonMatrix[i][1]=='F'])
            # C+=len([i for i in exonMatrix if exonMatrix[i][0]=='U' and exonMatrix[i][1]=='A' and sum(exonMatrix[i][2]) ==0])
            # D+=len([i for i in exonMatrix if exonMatrix[i][0]=='U' and exonMatrix[i][1]=='A' and sum(exonMatrix[i][2])>0])
            # E+=len([i for i in exonMatrix if exonMatrix[i][0]=='U' and exonMatrix[i][1]=='G'])
            # if UTRStrict[-1]!=(UTRMajorlyConstitutive[-1]+UTRAlternate[-1]+UTRAlternateWithSS[-1]+UTRConstitutive[-1]):
            #     print (gene)
            #     print (exonMatrix)
            #     print (UTRStrict[-1],UTRMajorlyConstitutive[-1],UTRAlternate[-1],UTRAlternateWithSS[-1],UTRConstitutive[-1])
            #     print '\n'
            #     sys.exit()


            #print (UTRMajorlyConstitutive)
            DualStrict+=[len([i for i in exonMatrix if exonMatrix[i][0]=='D'])]
            DualMajorlyConstitutive+=[len([i for i in exonMatrix if exonMatrix[i][0]=='D' and exonMatrix[i][1]=='F'])]   
            DualAlternate+=[len([i for i in exonMatrix if exonMatrix[i][0]=='D' and exonMatrix[i][1]=='A' and sum(exonMatrix[i][2]) ==0])]
            DualAlternateWithSS+=[len([i for i in exonMatrix if exonMatrix[i][0]=='D' and exonMatrix[i][1]=='A' and sum(exonMatrix[i][2])>0])]  
            DualConstitutive += [len([i for i in exonMatrix if exonMatrix[i][0]=='D' and exonMatrix[i][1]=='G'])]
            DualConstitutiveWaachange += [len([i for i in exonMatrix if exonMatrix[i][0]=='D' and exonMatrix[i][1]=='G' and exonMatrix[i][3]>0])]
            DualAltWaachange += [len([i for i in exonMatrix if exonMatrix[i][0]=='D' and exonMatrix[i][1]=='A' and sum(exonMatrix[i][2])==0 and exonMatrix[i][3]>0])]

            CodingStrict+=[len([i for i in exonMatrix if exonMatrix[i][0]=='T'])]
            CodingMajorlyConstitutive+=[len([i for i in exonMatrix if exonMatrix[i][0]=='T' and exonMatrix[i][1]=='F'])]   
            CodingAlternate+=[len([i for i in exonMatrix if exonMatrix[i][0]=='T' and exonMatrix[i][1]=='A' and sum(exonMatrix[i][2]) ==0])]
            CodingAlternateWithSS+=[len([i for i in exonMatrix if exonMatrix[i][0]=='T' and exonMatrix[i][1]=='A' and sum(exonMatrix[i][2])>0])]  
            CodingConstitutive+=[len([i for i in exonMatrix if exonMatrix[i][0]=='T' and exonMatrix[i][1]=='G'])]
            CodingConstitutiveWaachange+=[len([i for i in exonMatrix if exonMatrix[i][0]=='T' and exonMatrix[i][1]=='G' and exonMatrix[i][3]>0])]
            # print (CodingConstitutiveWaachange)
            CodingAltWaachange += [len([i for i in exonMatrix if exonMatrix[i][0]=='T' and exonMatrix[i][1]=='A' and sum(exonMatrix[i][2])==0 and exonMatrix[i][3]>0])]
            
            ir_retained_coding+=[len([i for i in exonMatrix if exonMatrix[i][0] if exonMatrix[i][5]])]
            ir_retained_UTR+=[len([i for i in exonMatrix if exonMatrix[i][0] if exonMatrix[i][4]])]
            # remember we are calculating here the non coding and coding variants per position, not how oftent their variations exists, means if exon 5 has retention veent starting, we would say that this position has been affected, not calculating how many times retention events occurs in psition 5 with poymorphs          

            strict_M+=[len([i for i in exonMatrix if exonMatrix[i][0]=='M'])]
            strict_aa_removed+=[len([i for i in exonMatrix if exonMatrix[i][6]==1])]
     
            # they needs to be added to the top, what if thats empty ?
            # ans: empty will be taken into cionsideration while plotting but stats count wont do that :)
            # excluding their D's from above as they may likely partcipate in the coding fraction atkeast once. 
            # define category strict alternate as the one with no change in the splice sites howsoever, and rest as different, for this different (AE with change in splice sites), count the instance of them once excluding their change in splice sites cases count.

    # print (A,B,C,D,E)
    hasPacked={'tot_exons':tot_exons, 'tot_exonsWncbAndRet':tot_exonsWncbAndRet, 'UTRStrict':UTRStrict,'UTRAlternate':UTRAlternate, 'UTRAlternateWithSS':UTRAlternateWithSS, 'UTRMajorlyConstitutive':UTRMajorlyConstitutive,'UTRConstitutive':UTRConstitutive, 'CodingStrict':CodingStrict, 'CodingAlternate':CodingAlternate,'CodingAltWaachange':CodingAltWaachange, 'CodingAlternateWithSS':CodingAlternateWithSS,'CodingMajorlyConstitutive':CodingMajorlyConstitutive, 'CodingConstitutive':CodingConstitutive,'CodingConstitutiveWaachange':CodingConstitutiveWaachange, 'DualStrict':DualStrict, 'DualAlternate':DualAlternate,'DualAltWaachange':DualAltWaachange, 'DualAlternateWithSS':DualAlternateWithSS, 'DualConstitutiveWaachange':DualConstitutiveWaachange, 'DualMajorlyConstitutive':DualMajorlyConstitutive,'DualConstitutive':DualConstitutive, 'strict_M':strict_M, 'strict_aa_removed':strict_aa_removed,'ir_retained_coding':ir_retained_coding, 'ir_retained_UTR':ir_retained_UTR}

    writer_and_displayer(res_dir,filename,fout,hasPacked)


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
from progress.bar import Bar
import time

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
    tabularData += "CodingAlternateWithSS\t%s\n"%cm.stats(hasPacked['CodingAlternateWithSS'])

    tabularData += "CodingAlternate\t%s\n"%cm.stats(hasPacked['CodingAlternate'])
    tabularData += "CodingAltWaachange\t%s\n"%cm.stats(hasPacked['CodingAltWaachange'])
    tabularData += "CodingAltWaachangeAssign\t%s\n"%cm.stats(hasPacked['CodingAltWaachangeAssign'])
    tabularData += "CodingAltWaachangeFrame\t%s\n"%cm.stats(hasPacked['CodingAltWaachangeFrame'])
    
    tabularData += "CodingMajorlyConstitutive\t%s\n"%cm.stats(hasPacked['CodingMajorlyConstitutive'])
    tabularData += "CodingMajorlyWaachange\t%s\n"%cm.stats(hasPacked['CodingMajorlyWaachange'])
    tabularData += "CodingMajorlyWaachangeAssign\t%s\n"%cm.stats(hasPacked['CodingMajorlyWaachangeAssign'])
    tabularData += "CodingMajorlyWaachangeFrame\t%s\n"%cm.stats(hasPacked['CodingMajorlyWaachangeFrame'])

    tabularData += "CodingConstitutive\t%s\n"%cm.stats(hasPacked['CodingConstitutive'])
    tabularData += "CodingConstitutiveWaachange\t%s\n"%cm.stats(hasPacked['CodingConstitutiveWaachange'])
    tabularData += "CodingConstitutiveWaachangeAssign\t%s\n"%cm.stats(hasPacked['CodingConstitutiveWaachangeAssign'])
    tabularData += "CodingConstitutiveWaachangeFrame\t%s\n"%cm.stats(hasPacked['CodingConstitutiveWaachangeFrame'])
    
    tabularData += "DualStrict\t%s\n"%cm.stats(hasPacked['DualStrict'])
    tabularData += "DualAlternateWithSS\t%s\n"%cm.stats(hasPacked['DualAlternateWithSS'])
    
    tabularData += "DualAlternate\t%s\n"%cm.stats(hasPacked['DualAlternate'])
    tabularData += "DualAltWaachange\t%s\n"%cm.stats(hasPacked['DualAltWaachange'])
    tabularData += "DualAltWaachangeAssign\t%s\n"%cm.stats(hasPacked['DualAltWaachangeAssign'])
    tabularData += "DualAltWaachangeFrame\t%s\n"%cm.stats(hasPacked['DualAltWaachangeFrame'])
    
    tabularData += "DualMajorlyConstitutive\t%s\n"%cm.stats(hasPacked['DualMajorlyConstitutive'])
    tabularData += "DualMajorlyWaachange\t%s\n"%cm.stats(hasPacked['DualMajorlyWaachange'])
    tabularData += "DualMajorlyWaachangeAssign\t%s\n"%cm.stats(hasPacked['DualMajorlyWaachangeAssign'])
    tabularData += "DualMajorlyWaachangeFrame\t%s\n"%cm.stats(hasPacked['DualMajorlyWaachangeFrame'])

    tabularData += "DualConstitutive\t%s\n"%cm.stats(hasPacked['DualConstitutive'])
    tabularData += "DualConstitutiveWaachange\t%s\n"%cm.stats(hasPacked['DualConstitutiveWaachange'])
    tabularData += "DualConstitutiveWaachangeAssign\t%s\n"%cm.stats(hasPacked['DualConstitutiveWaachangeAssign'])
    tabularData += "DualConstitutiveWaachangeFrame\t%s\n"%cm.stats(hasPacked['DualConstitutiveWaachangeFrame'])

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
        for i in hasPacked['CodingAltWaachangeAssign']:
            fin.write("%s,AlternateWaaChangeAssign,Coding\n"%(i))
        for i in hasPacked['CodingAltWaachangeFrame']:
            fin.write("%s,AlternateWaaChangeFrame,Coding\n"%(i))

            
        for i in hasPacked['CodingAlternateWithSS']:
            fin.write("%s,AltWithSSchange,Coding\n"%(i))
        for i in hasPacked['CodingMajorlyConstitutive']:
            fin.write("%s,MajorlyConst,Coding\n"%(i))
        for i in hasPacked['CodingMajorlyWaachange']:
            fin.write("%s,MajorlyConstWaaChange,Coding\n"%(i))
        for i in hasPacked['CodingMajorlyWaachangeAssign']:
            fin.write("%s,MajorlyConstWaaChangeAssign,Coding\n"%(i))
        for i in hasPacked['CodingMajorlyWaachangeFrame']:
            fin.write("%s,MajorlyConstWaaChangeFrame,Coding\n"%(i))

            
            
            

        for i in hasPacked['CodingConstitutive']:
            fin.write("%s,Const,Coding\n"%(i))
        for i in hasPacked['CodingConstitutiveWaachange']:
            fin.write("%s,ConstWaaChange,Coding\n"%(i))
        for i in hasPacked['CodingConstitutiveWaachangeAssign']:
            fin.write("%s,ConstWaaChangeAssign,Coding\n"%(i))
        for i in hasPacked['CodingConstitutiveWaachangeFrame']:
            fin.write("%s,ConstWaaChangeFrame,Coding\n"%(i))
            
            
        
        for i in hasPacked['DualStrict']:
            fin.write("%s,Strict,Dual\n"%(i))
        for i in hasPacked['DualAlternate']:
            fin.write("%s,Alternate,Dual\n"%(i))

        for i in hasPacked['DualAltWaachange']:
            fin.write("%s,AlternateWaaChange,Dual\n"%(i))
        for i in hasPacked['DualAltWaachangeAssign']:
            fin.write("%s,AlternateWaaChangeAssign,Dual\n"%(i))
        for i in hasPacked['DualAltWaachangeFrame']:
            fin.write("%s,AlternateWaaChangeFrame,Dual\n"%(i))
            
            
        for i in hasPacked['DualAlternateWithSS']:
            fin.write("%s,AltWithSSchange,Dual\n"%(i))
        for i in hasPacked['DualMajorlyConstitutive']:
            fin.write("%s,MajorlyConst,Dual\n"%(i))
        for i in hasPacked['DualMajorlyWaachange']:
            fin.write("%s,MajorlyWaaChange,Dual\n"%(i))
        for i in hasPacked['DualMajorlyWaachangeAssign']:
            fin.write("%s,MajorlyWaaChangeAssign,Dual\n"%(i))
        for i in hasPacked['DualMajorlyWaachangeFrame']:
            fin.write("%s,MajorlyWaaChangeFrame,Dual\n"%(i))
            
            
        for i in hasPacked['DualConstitutive']:
            fin.write("%s,Const,Dual\n"%(i))
        for i in hasPacked['DualConstitutiveWaachange']:
            fin.write("%s,ConstWaaChange,Dual\n"%(i))
        for i in hasPacked['DualConstitutiveWaachangeAssign']:
            fin.write("%s,ConstWaaChangeAssign,Dual\n"%(i))
        for i in hasPacked['DualConstitutiveWaachangeFrame']:
            fin.write("%s,ConstWaaChangeFrame,Dual\n"%(i))
        
def exon_counting(has, hasFramechAssign, res_dir,condition,fout ):
    # newly added categories ################
    # 04 -11 -22
    tot_exonsWncbAndRet = []
    
    DualConstitutiveWaachange = []
    DualConstitutiveWaachangeAssign =[]
    DualConstitutiveWaachangeFrame=[]

    DualMajorlyWaachange = []
    DualMajorlyWaachangeAssign = []
    DualMajorlyWaachangeFrame = []
    
    DualAltWaachange = []
    DualAltWaachangeAssign = []
    DualAltWaachangeFrame = []

    CodingConstitutiveWaachange = []
    CodingConstitutiveWaachangeAssign = []
    CodingConstitutiveWaachangeFrame = []
    
    CodingMajorlyWaachange = []
    CodingMajorlyWaachangeAssign = []
    CodingMajorlyWaachangeFrame = []
    
    CodingAltWaachange = []
    CodingAltWaachangeAssign = []
    CodingAltWaachangeFrame = []

    
    
    
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
    bar = Bar('Processing genes C:', max=len(condition))
    
    for gene in has:
        if gene in condition: # and gene == 11259:
            # st=time.time()
            # print ("**GENE:", gene)
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
            
            # exonMatrix[placeholder] = ['','',[0,0,0], 0,0,0, 0, [0,0,0]] #0, 1,    2   , 3,4,5, 6, 7
            # 0th ele, UTMRD tag, 
            # 1st ele, AGF tag, 
            # 2nd records children order of ncb, and 
            # 3rd occurrences of the aa change on this position, (should not be evaluated if values of the ncb counterparts is true
            # 4th coding intron retention events starting from this
            # 5th non coding intron retention events starting from this
            # 6th whether aa has been ever removed from this sequence (yes consider all the ncb variation also but not retention)
            # 7th whether ncb variations ever had change in aa case, like FMR1 gene does have, in that case it will be binary and value 1 means atleast 1 sduch case exists, 

            tot_exons+=[len(exonMatrix)]
            tot_exonsWncbAndRet+=[len(exonMatrix)+sum([exonMatrix[i][2][0]+exonMatrix[i][2][1]+exonMatrix[i][2][2]+exonMatrix[i][4]+exonMatrix[i][5] for i in exonMatrix])]

            temp1, temp2, temp3, temp4, temp5 = ca.giveExonSelectionBasic(exonMatrix,"U")
            UTRStrict += [len(temp1)]
            UTRMajorlyConstitutive += [len(temp2)]
            UTRAlternate += [len(temp3)]
            UTRAlternateWithSS += [len(temp4)]
            UTRConstitutive += [len(temp5)]

            temp1, temp2, temp3, temp4, temp5 = ca.giveExonSelectionBasic(exonMatrix,"D")
            DualStrict += [len(temp1)]
            DualMajorlyConstitutive += [len(temp2)]
            DualAlternate += [len(temp3)]
            DualAlternateWithSS += [len(temp4)]
            DualConstitutive += [len(temp5)]
            
            temp1, temp2, temp3, temp4, temp5 = ca.giveExonSelectionBasic(exonMatrix,"T")
            CodingStrict += [len(temp1)]
            CodingMajorlyConstitutive += [len(temp2)]
            CodingAlternate += [len(temp3)]
            CodingAlternateWithSS += [len(temp4)]
            CodingConstitutive += [len(temp5)]
            
            #
            # hasFramechAssign
            # preDoneHas[gene][(utmrd0, utmrd2)]=[[],{}]
            # gene, utmrd0, utmrd2, type, placeholderTags = ele
            # preDoneHas[gene][(utmrd0, utmrd2)][0]=placeholderList
            # preDoneHas[gene][(utmrd0, utmrd2)][1][(exon1,exon2)]=[aln1,aln2,aln3, [seq1Length, seq2Length, identity1, identity2, cov1, cov2]]
            # utmrd0 is T/D, utmrd2 is G/A/F, exon1 is #6:0:0, 6:n:1, exon 2 var, 2,3,4,... aln1 is parent, aln2 is child, aln3 is bars, identity1, identity2, cov1, cov2 in fractions
            # 100 100 100 100 is for assign automatically is length diff is 1/4th, and same is 0 0 0 0 for framechange , rest values will be globaldx assignment 

            
            def AssignerFramerdefault (hasFramechAssign, tupkey, l1, l2, l3):
                if hasFramechAssign[gene][tupkey][0]:
                    l1 +=  [len(hasFramechAssign[gene][tupkey][0])]
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
                        l2 += [len(set(assign))]
                        l3 += [len(set(frame))]
                    else:
                        l2 += [len([])]
                        l3 += [len([])]
                else:
                    l1+= [len([])]
                    l2+= [len([])]
                    l3+= [len([])]
                return l1, l2, l3
            
           
            DualConstitutiveWaachange, DualConstitutiveWaachangeAssign, DualConstitutiveWaachangeFrame = AssignerFramerdefault(hasFramechAssign, ('D','G'),DualConstitutiveWaachange, DualConstitutiveWaachangeAssign, DualConstitutiveWaachangeFrame)
            DualAltWaachange, DualAltWaachangeAssign, DualAltWaachangeFrame = AssignerFramerdefault(hasFramechAssign, ('D','A'),DualAltWaachange, DualAltWaachangeAssign, DualAltWaachangeFrame)
            DualMajorlyWaachange, DualMajorlyWaachangeAssign, DualMajorlyWaachangeFrame = AssignerFramerdefault(hasFramechAssign, ('D','F'),DualMajorlyWaachange, DualMajorlyWaachangeAssign, DualMajorlyWaachangeFrame)
            CodingConstitutiveWaachange, CodingConstitutiveWaachangeAssign, CodingConstitutiveWaachangeFrame = AssignerFramerdefault(hasFramechAssign, ('T','G'),CodingConstitutiveWaachange, CodingConstitutiveWaachangeAssign, CodingConstitutiveWaachangeFrame)
            CodingAltWaachange, CodingAltWaachangeAssign, CodingAltWaachangeFrame =AssignerFramerdefault(hasFramechAssign, ('T','A'),CodingAltWaachange, CodingAltWaachangeAssign, CodingAltWaachangeFrame)
            CodingMajorlyWaachange, CodingMajorlyWaachangeAssign, CodingMajorlyWaachangeFrame = AssignerFramerdefault(hasFramechAssign, ('T','F'),CodingMajorlyWaachange, CodingMajorlyWaachangeAssign, CodingMajorlyWaachangeFrame)

            ir_retained_coding+=[len([i for i in exonMatrix if exonMatrix[i][0] if exonMatrix[i][5]])]
            ir_retained_UTR+=[len([i for i in exonMatrix if exonMatrix[i][0] if exonMatrix[i][4]])]
            # remember we are calculating here the non coding and coding variants per position, not how oftent their variations exists, means if exon 5 has retention veent starting, we would say that this position has been affected, not calculating how many times retention events occurs in psition 5 with poymorphs          

            strict_M+=[len([i for i in exonMatrix if exonMatrix[i][0]=='M'])]
            strict_aa_removed+=[len([i for i in exonMatrix if exonMatrix[i][6]==1])]
            # ed = time.time()
            # spannedT= round((ed-st)/60,2)
            # print ("TIME:",spannedT)
            bar.next()
            # they needs to be added to the top, what if thats empty ?
            # ans: empty will be taken into cionsideration while plotting but stats count wont do that :)
            # excluding their D's from above as they may likely partcipate in the coding fraction atkeast once. 
            # define category strict alternate as the one with no change in the splice sites howsoever, and rest as different, for this different (AE with change in splice sites), count the instance of them once excluding their change in splice sites cases count.

    # print (A,B,C,D,E)
    bar.finish()
    hasPacked = {
        'tot_exons':tot_exons, 'tot_exonsWncbAndRet':tot_exonsWncbAndRet, 
        'UTRStrict':UTRStrict,'UTRAlternate':UTRAlternate, 'UTRAlternateWithSS':UTRAlternateWithSS, 'UTRMajorlyConstitutive':UTRMajorlyConstitutive,'UTRConstitutive':UTRConstitutive, 
        
        'CodingStrict':CodingStrict, 'CodingAlternate':CodingAlternate, 'CodingAlternateWithSS':CodingAlternateWithSS,'CodingMajorlyConstitutive':CodingMajorlyConstitutive,  'CodingConstitutive':CodingConstitutive, 
        'CodingAltWaachange':CodingAltWaachange, 'CodingAltWaachangeAssign':CodingAltWaachangeAssign, 'CodingAltWaachangeFrame':CodingAltWaachangeFrame,  
        'CodingMajorlyWaachange':CodingMajorlyWaachange, 'CodingMajorlyWaachangeAssign': CodingMajorlyWaachangeAssign, 'CodingMajorlyWaachangeFrame':CodingMajorlyWaachangeFrame,  
        'CodingConstitutiveWaachange':CodingConstitutiveWaachange, 'CodingConstitutiveWaachangeAssign': CodingConstitutiveWaachangeAssign, 'CodingConstitutiveWaachangeFrame':CodingConstitutiveWaachangeFrame,  

        'DualStrict':DualStrict, 'DualAlternate':DualAlternate, 'DualAlternateWithSS':DualAlternateWithSS, 'DualMajorlyConstitutive':DualMajorlyConstitutive,'DualConstitutive':DualConstitutive, 
        'DualAltWaachange':DualAltWaachange, 'DualAltWaachangeAssign':DualAltWaachangeAssign, 'DualAltWaachangeFrame':DualAltWaachangeFrame,
        'DualMajorlyWaachange':DualMajorlyWaachange,'DualMajorlyWaachangeAssign':DualMajorlyWaachangeAssign, 'DualMajorlyWaachangeFrame':DualMajorlyWaachangeFrame,
        'DualConstitutiveWaachange':DualConstitutiveWaachange,'DualConstitutiveWaachangeAssign':DualConstitutiveWaachangeAssign, 'DualConstitutiveWaachangeFrame':DualConstitutiveWaachangeFrame,   
         
        'strict_M':strict_M, 'strict_aa_removed':strict_aa_removed,
        'ir_retained_coding':ir_retained_coding, 'ir_retained_UTR':ir_retained_UTR
        }

    writer_and_displayer(res_dir,filename,fout,hasPacked)

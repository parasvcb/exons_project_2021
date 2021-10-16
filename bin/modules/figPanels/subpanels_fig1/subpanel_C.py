
'''
This block will read the genes in object, will filter them on the basis of atleast two protein coding different siforms
atleast two proetin coding exons and length of pricinpal soform less than 3000 aa
when that criteria is met, it will ask to change the parameter valeu from unique to varition 
that then will calculate the listed cases on the basis of the unique genomic coods as once instance 
and in latter to count them on the basis of the variation also in those cases

correspsonding data is written into the files also and dispaleys on screen with genral statustucs
'''
import re,os
import figPanels.modules_common as cm
import figPanels.modules_analysis as ca

def counter(exons,parameter):  
    #fout.write( len(exons),"inside"
    #fout.write( [i.ID for i in exons]
    if parameter == "Unique":
        uni_ex={}
        for i in exons:
            lett=".".join(i.ID.split(".")[2:])
            uni_ex[lett]=0
        #fout.write( len(uni_ex),"len"
        return len(uni_ex)
    else:
        return len(exons)


def writer_and_displayer(res_dir,filename,fout,hasPacked):
    fout.write(filename)
    fout.write("\n=================\nUnique exons count based on genomic coordinates\n")
    fout.write(
        "\tCategory\tVar_count\tSum\tMax\tMin\tMean\tMedian\tMode\tMode_count\tStddev\n")
    fout.write( "\n\tTotal_exons\t%s"%cm.stats(hasPacked['tot_exons']))

    fout.write( "\n\tUTRStrict\t%s"%cm.stats(hasPacked['UTRStrict']))
    fout.write( "\n\tUTRAlternate\t%s"%cm.stats(hasPacked['UTRAlternate']))
    fout.write( "\n\tUTRAlternateWithSS\t%s"%cm.stats(hasPacked['UTRAlternateWithSS']))
    fout.write( "\n\tUTRMajorlyConstitutive\t%s"%cm.stats(hasPacked['UTRMajorlyConstitutive']))
    fout.write( "\n\tUTRConstitutive\t%s"%cm.stats(hasPacked['UTRConstitutive']))
    
    fout.write( "\n\tCodingStrict\t%s"%cm.stats(hasPacked['CodingStrict']))
    fout.write( "\n\tCodingAlternate\t%s"%cm.stats(hasPacked['CodingAlternate']))
    fout.write( "\n\tCodingAlternateWithSS\t%s"%cm.stats(hasPacked['CodingAlternateWithSS']))
    fout.write( "\n\tCodingMajorlyConstitutive\t%s"%cm.stats(hasPacked['CodingMajorlyConstitutive']))
    fout.write( "\n\tCodingConstitutive\t%s"%cm.stats(hasPacked['CodingConstitutive']))
    
    fout.write( "\n\tDualStrict\t%s"%cm.stats(hasPacked['DualStrict']))
    fout.write( "\n\tDualAlternate\t%s"%cm.stats(hasPacked['DualAlternate']))
    fout.write( "\n\tDualAlternateWithSS\t%s"%cm.stats(hasPacked['DualAlternateWithSS']))
    fout.write( "\n\tDualMajorlyConstitutive\t%s"%cm.stats(hasPacked['DualMajorlyConstitutive']))
    fout.write( "\n\tDualConstitutive\t%s"%cm.stats(hasPacked['DualConstitutive']))
    
    fout.write( "\n\tIntron_retained_coding\t%s"%cm.stats(hasPacked['ir_retained_coding']))
    fout.write( "\n\tIntron_retained_UTR\t%s"%cm.stats(hasPacked['ir_retained_UTR']))
    fout.write( "\n\tStrictMcases(1ntLength)\t%s"%cm.stats(hasPacked['strict_M']))
    fout.write( "\n\tStrict_aa_removed\t%s"%cm.stats(hasPacked['strict_aa_removed']))
    fout.write("\nCSV_file:%s"%res_dir+filename+".csv")
    
    with open(os.path.join(res_dir,filename+".csv"),"w") as fin:

        fin.write("Exons,Category,Type\n")
        for i in hasPacked['tot_exons']:
            fin.write("%s,TotalExons,Extras\n"%(i))
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
        for i in hasPacked['CodingAlternateWithSS']:
            fin.write("%s,AltWithSSchange,Coding\n"%(i))
        for i in hasPacked['CodingMajorlyConstitutive']:
            fin.write("%s,MajorlyConst,Coding\n"%(i))
        for i in hasPacked['CodingConstitutive']:
            fin.write("%s,Const,Coding\n"%(i))
        
        for i in hasPacked['DualStrict']:
            fin.write("%s,Strict,Dual\n"%(i))
        for i in hasPacked['DualAlternate']:
            fin.write("%s,Alternate,Dual\n"%(i))
        for i in hasPacked['DualAlternateWithSS']:
            fin.write("%s,AltWithSSchange,Dual\n"%(i))
        for i in hasPacked['DualMajorlyConstitutive']:
            fin.write("%s,MajorlyConst,Dual\n"%(i))
        for i in hasPacked['DualConstitutive']:
            fin.write("%s,Const,Dual\n"%(i))

def exon_counting(has,res_dir,condition,fout ):
    '''
    There will be two sets now
    First: Unique coordinate wise
        1.Total number of all exons?
        2.total number of strict coding exons?
        3.Total strict UTR?
        4.Total strict M and 
        5.Total Strict aa change cases
        6.Total number of dual functioning exons
        7.Total number of aternate exons, coding 
        8.Total number of consectuive exons, coding
        9.Total changein splice sites
        10.Total Coding intron retention cases and non coding such cases 
    Second: Through in all the aa properties
    '''
    tot_exons=[]
    UTRStrict=[]
    UTRAlternate=[]
    UTRAlternateWithSS=[]
    UTRMajorlyConstitutive=[]
    UTRConstitutive=[]

    CodingStrict=[]
    CodingAlternate=[]
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

    PARAMETER="Unique"
    #PARAMETER="Variation"
    filename="Variation_allowed" if PARAMETER == "Variation" else "Unique_Coods"
    output_stats = open("General_stats_fig1.log", "a")
    for gene in has:
        if  gene in condition:
            tot_exons+=[counter([i for i in has[gene].exons],PARAMETER)]
            # dividing in three categories
            # 1st non coding regions
            # they will remain non codin in the entire duration,
            # their Alternate exons, Alt exons counted once with change insplice sites
            # constitutive exons, Alt exons counted once with chagne in splice sites when count gets eual to one
            
            
            # UTR_alternate+=[counter([i for i in has[gene].exons if re.match(r'^[UD]\.[-][2]+\.A',i.ID)],PARAMETER)]
            # UTR_consective+=[counter([i for i in has[gene].exons if re.match(r'^[UD]\.[-][2]+\.G',i.ID)],PARAMETER)]
            # majorly_ConstitutiveUTR+=[counter([i for i in has[gene].exons if re.match(r'^[UD]\.[-][2]\.F',i.ID)],PARAMETER)]
            
            UTRStrict+=[counter([i for i in has[gene].exons if re.match(r'^U',i.ID)],PARAMETER)]
            UTRConstitutive+=[counter([i for i in has[gene].exons if re.match(r'^[U]\.[-][2]+\.G',i.ID)],PARAMETER)]
            
            tempExAltUTR=[i for i in has[gene].exons if re.match(r'^[U]\.[-][2]+\.A',i.ID)]
            res1,res2=ca.split_cases(tempExAltUTR)
            #its output will give in first count of exons with no splice site change and 0 status, and then exons with chnage in splice sites, counted once for all variation
            UTRAlternate+=[res1]
            UTRAlternateWithSS+=[res2]
            
            tempExMajorUTR=[i for i in has[gene].exons if re.match(r'^[U]\.[-][2]\.F',i.ID)]
            res1,res2=ca.split_cases(tempExMajorUTR)
            UTRMajorlyConstitutive+=[res2]
            #they needs to be added to the top, what if thats empty ?
            #ans: empty will be taken into cionsideration while plotting but stats count wont do that :)
            # excluding their D's from above as they may likely partcipate in the coding fraction atkeast once. 
            # define category strict alternate as the one with no change in the splice sites howsoever, and rest as different, for this different (AE with change in splice sites), count the instance of them once excluding their change in splice sites cases count.

            # 2nd dual function exons, they may code for the coding as well as non coding
            DualStrict+=[counter([i for i in has[gene].exons if re.match(r'^D',i.ID)],PARAMETER)]
            DualConstitutive+=[counter([i for i in has[gene].exons if re.match(r'^D\.[-]?[1-n]+\.G',i.ID)],PARAMETER)]
            tempExAltDual=[i for i in has[gene].exons if re.match(r'^D\.[-]?[1-n]+\.A',i.ID)]
            res1,res2=ca.split_cases(tempExAltDual)
            DualAlternate+=[res1]
            DualAlternateWithSS+=[res2]
            tempExMajorDual=[i for i in has[gene].exons if re.match(r'^D\.[-]?[1-n]+\.F',i.ID)]
            res1,res2=ca.split_cases(tempExMajorDual)
            DualMajorlyConstitutive+=[res2]
            

            # 3rd will also be the same categories,
            CodingStrict+=[counter([i for i in has[gene].exons if re.match(r'^T',i.ID)],PARAMETER)]
            CodingConstitutive+=[counter([i for i in has[gene].exons if re.match(r'^T\.[1-n]+\.G',i.ID)],PARAMETER)]
            tempExAltCoding=[i for i in has[gene].exons if re.match(r'^T\.[1-n]+\.A',i.ID)]
            res1,res2=ca.split_cases(tempExAltCoding)
            CodingAlternate+=[res1]
            CodingAlternateWithSS+=[res2]
            tempExMajorCoding=[i for i in has[gene].exons if re.match(r'^T\.[1-n]+\.F',i.ID)]
            res1,res2=ca.split_cases(tempExMajorCoding)
            CodingMajorlyConstitutive+=[res2] 
            #previously it was T and D but now only D
           
           
            ir_retained_coding+=[counter([i for i in has[gene].exons if re.match(r'^R\:[1-n]',i.ID)],PARAMETER)]
            ir_retained_UTR+=[counter([i for i in has[gene].exons if re.match(r'^R\:[-][2]',i.ID)],PARAMETER)]            
            strict_M+=[counter([i for i in has[gene].exons if re.match(r'^M',i.ID)],PARAMETER)]
            strict_aa_removed+=[counter([i for i in has[gene].exons if re.match(r'^\w\.\-1',i.ID)],PARAMETER)]

    hasPacked={'tot_exons':tot_exons,'UTRStrict':UTRStrict,'UTRAlternate':UTRAlternate,'UTRAlternateWithSS':UTRAlternateWithSS,'UTRMajorlyConstitutive':UTRMajorlyConstitutive,'UTRConstitutive':UTRConstitutive,'CodingStrict':CodingStrict,'CodingAlternate':CodingAlternate,'CodingAlternateWithSS':CodingAlternateWithSS,'CodingMajorlyConstitutive':CodingMajorlyConstitutive,'CodingConstitutive':CodingConstitutive,'DualStrict':DualStrict,'DualAlternate':DualAlternate,'DualAlternateWithSS':DualAlternateWithSS,'DualMajorlyConstitutive':DualMajorlyConstitutive,'DualConstitutive':DualConstitutive,'strict_M':strict_M,'strict_aa_removed':strict_aa_removed,'ir_retained_coding':ir_retained_coding,'ir_retained_UTR':ir_retained_UTR}

    writer_and_displayer(res_dir,filename,fout,hasPacked)

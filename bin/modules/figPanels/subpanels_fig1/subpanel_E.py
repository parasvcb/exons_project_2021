import re,os
import figPanels.modules_common as cm
import figPanels.modules_analysis as ca

def constitutive_alternate_with_freq(gene):
    constitutive=[]
    Alternate=[]
    Majorly=[]
    for ex in gene.exons:
        if re.match(r'^[T]\.[1-n]+\.G',ex.ID):
            constitutive+=[ex.length]
        elif re.match(r'^[T]\.[1-n]+\.[AF]',ex.ID):
            specifier=ex.ID.split(".")[2]
            if specifier == 'F':
                Majorly+=[ex]
            else:
                Alternate+=[ex]

    res1,res2=ca.split_cases(Alternate,out='exons')
    #its output will give in first count of exons with no splice site change and 0 status, and then exons with chnage in splice sites, counted once for all variation
    CodingAlternate=res1
    # this will be list with single elements
    CodingAlternateWithSS=res2
    # this will be list of list, each list having possible combination of varitaions
    res1,res2=ca.split_cases(Majorly,out='exons')
    CodingMajorly=res1
    CodingMajorlyWithSS=res2

    exCommonAlt=[]
    exSpliceSiteChangeAlt=[]
    exSpliceSiteChangeMajorly=[]
    
    for i in CodingAlternate:
        exCommonAlt+=[i.length]
    for i in CodingAlternateWithSS:
        avglen = sum([ex.length for ex in i])/len(i)
        exSpliceSiteChangeAlt+=[avglen]
    for i in CodingMajorlyWithSS:
        avglen = sum([ex.length for ex in i])/len(i)
        exSpliceSiteChangeMajorly+=[avglen]

    return constitutive, exCommonAlt, exSpliceSiteChangeAlt, exSpliceSiteChangeMajorly


def exon_length(has,res_dir,genes_cond,output_stats):
    consecutive_exons=[]
    alternate_exons=[]
    alternate_exons_with_splicesites=[]
    majorly_exons=[]

    for gene in has:
        if gene in genes_cond:
            alt_ex={}
            consect,alt,altSS,majorlySS=constitutive_alternate_with_freq(has[gene])
            consecutive_exons+=consect
            alternate_exons+=alt
            alternate_exons_with_splicesites+=altSS
            majorly_exons+=majorlySS
    
    
    output_stats.write("\n\n=============>\n")
    output_stats.write("Length Comparison of constituive and alternate exons\n")
    output_stats.write("\n\tCategory\tVar_count\tSum\tMax\tMin\tMean\tMedian\tMode\tMode_count\tStddev")
    output_stats.write("\n\tconstitutive_exons\t%s"%cm.stats(consecutive_exons))
    output_stats.write("\n\talternate_exons\t%s"%cm.stats(alternate_exons))
    output_stats.write("\n\talternate_exons_SS(meanAlready)\t%s"%cm.stats(alternate_exons_with_splicesites))
    output_stats.write("\n\tmajorly_exons(meanAlready)\t%s"%cm.stats(majorly_exons))
    output_stats.write("\nCSV_File:%s"%res_dir+"length_distribution_exons.csv")

    with open(os.path.join(res_dir,'panelE_Length_distribution_exonsRAW.csv'),"w") as fin:
        fin.write("Length,Category\n")
        for i in consecutive_exons:
            fin.write("%s,Consecutive\n"%(i))
        for i in alternate_exons:
            fin.write("%s,Alternate-Coding\n"%(i))
        for i in alternate_exons_with_splicesites:
            fin.write("%s,Alternate-Coding_wSS\n"%(i))
        for i in majorly_exons:
            fin.write("%s,Majorly-Coding\n"%(i))   


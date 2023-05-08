
import sys,re, os
if len(sys.argv) != 4:
    print ("Please provide file that has newline spearated log files to parse for different oraganisms (General), output directory to store the CSV's")
    sys.exit()
prog, inpfile, fmode, outdir = sys.argv

if fmode =='2':
    fprefix = 'F1_Condition_3000Pilength_2Isf_2ExCount'
else:
    fprefix = 'F1_Condition_3000Pilength_4Isf_4ExCount'

def compilegeneTableTransandExonData (fname, string, organism):
    with open (fname) as fin:
        dat=[i for i in fin.read().split('\n') if len(i)>1]
    '''
    Different_exons_count	GeneCount	Frequency
	2	166	0.077
    '''
    for i in dat[1:]:
        categ, gcount, freq = i.split('\t')
        string += "%s\t%s\t%s\t%s\n"%(categ, organism, gcount, freq)
    return string

def compilegeneTableTransandExonDataStacked (fname, string, organism):
    with open (fname) as fin:
        dat=[i for i in fin.read().split('\n') if len(i)>1]
    '''
    ISFCount	Tag	GeneCount	Frequency
	2	166	0.077
    '''
    # Category\tTag\tOrganism\tValue
    for i in dat[1:]:
        categ, tag, count, freq = i.split('\t')
        string += "%s\t%s\t%s\t%s\t%s\n"%(categ, tag, organism, count, freq)
    return string

def compilegeneTableTransandExonDataContour (fname, string, organism):
    # Gene\tCodingExonCount\tNRIsfCount\tOrganism
    with open (fname) as fin:
        dat=[i for i in fin.read().split('\n') if len(i)>1]
    for i in dat[1:]:
        gene, exon, isf = i.split('\t')
        string += "%s\t%s\t%s\t%s\n"%(gene, exon, isf, organism)
    return string

def compileExonData (fname, string1,string2, organism, mode = "basic"):
    with open (fname) as fin:
        dat=[i for i in fin.read().split('\n') if len(i)>1]
    '''
    Category	GeneCount	Total	Max	Min	Mean	Median	Mode	Mode_count	Stddev
    Total_exons	2151	17271	41	2	8.029	7.0	4	255	4.977
    Total_exons_ncb&R	2151	21806	47	2	10.138	9.0	6	196	5.854

    # basic Category\tFacet\tOrganism\tValueMean\tValueStd\n
    # Category\tFacet\tOrganism\tValueFractionExons\tValueFractionGenes\n
    '''
    for i in dat[1:]:
        categ, gcount, Total, Max, Min, Mean, Median, Mode, M_count, Stddev = i.split('\t')
        # "UTRStrict","UTRAlternate","UTRAlternateWithSS","UTRMajorlyConstitutive","UTRConstitutive",
        # "CodingStrict","CodingAlternate","CodingAlternateWithSS","CodingMajorlyConstitutive","CodingConstitutive",
        # "DualStrict","DualAlternate",","DualAlternateWithSS","DualMajorlyConstitutive","DualConstitutive":
        if categ == "Total_exons":
            anchorGene = float (gcount)
            anchorExon = float (Total)
        facet = "Coding" if "Coding" == categ[:6] else "UTR" if "UTR" == categ[:3] else "Dual" if "Dual" == categ[:4] else False
        if facet and "aachange" not in categ:
            categ=categ[len(facet):]
            categ = "0TotalStrict" if categ =="Strict" else categ
            if mode == "basic":
                string1 += "%s\t%s\t%s\t%s\t%s\n" % (facet, categ, organism, Mean, Stddev)
            else:
                exFrac = round (int(Total)/anchorExon,3)
                geFrac = round (int(gcount)/anchorGene,3)
                string1 += "%s\t%s\t%s\t%s\t%s\n" % (facet, categ, organism, exFrac, geFrac)
        else:
            if categ[:5]!="Total":
            #"CodingAltWaachange", "CodingConstitutiveWaachange", "DualAltWaachange", "DualConstitutiveWaachange", "Intron_retained_coding","Intron_retained_UTR","StrictMcases(1ntLength)","Strict_aa_removed"
                if mode == "basic":
                    string2 += "%s\t%s\t%s\t%s\n" % (categ, organism, Mean, Stddev)
                else:
                    exFrac = round (int(Total)/anchorExon,3)
                    geFrac = round (int(gcount)/anchorGene,3)
                    string2 += "%s\t%s\t%s\t%s\n" % (categ, organism, exFrac, geFrac)
    return string1, string2

def compileExonLength (fname, string, organism):
    with open (fname) as fin:
        dat=[i for i in fin.read().split('\n') if len(i)>1]
    
    '''
    Category	Var_count	Sum	Max	Min	Mean	Median	Mode	Mode_count	Stddev
    CodingStrict	14373	1536770.024	2269.0	1.0	106.921	61.0	37.0	163	142.061
    CodingAlternate	3616	318294.833	2269.0	1.0	88.024	47.0	1.0	59	149.292
    # Category\tOrganism\tValueMean\tValueMedian\tValueStd\n
    '''

    for i in dat[1:]:
        categ, gcount, Total, Max, Min, Mean, Median, Mode, M_count, Stddev = i.split('\t')
        string += "%s\t%s\t%s\t%s\t%s\n"%(categ, organism, Mean, Median, Stddev)
    return string

def compileWEF (fname, string, organism):
    with open (fname) as fin:
        dat=[i for i in fin.read().split('\n') if len(i)>1]
    '''
    Inclusion_Range,Type,Total,Frequency
    0-0.2,AltOnly,7439,0.152
    0.201-0.4,AltOnly,6692,0.137
    0.401-0.6,AltOnly,9651,0.197
    # Category\tOrganism\trange\tFrequency\tTotal\n
    '''
    for i in dat[1:]:
        rangeInc, categ, Total, freq = i.split(',')
        string += "%s\t%s\t%s\t%s\t%s\n"%(categ, organism, rangeInc, freq, Total)
    return string

def compileWEF_raw (fname, string, organism):
    with open (fname) as fin:
        dat=[i for i in fin.read().split('\n') if len(i)>1]
    for i in dat[1:]:
        categ, value = i.split(',')
        string += "%s\t%s\t%s\n"%(categ, organism, value)
    return string

def compileTrans (fname, string, organism):
    with open (fname) as fin:
        dat=[i for i in fin.read().split('\n') if len(i)>1]
    '''
    Category	Var_count	Sum	Max	Min	Mean	Median	Mode	Mode_count	Stddev
    All_transcripts_statistics (including redundant, but passed criteria)	12976	98027	92	2	7.554	5.0	2	2263	7.422
    All_transcripts_statistics (only nonredundantPer gene)	12976	71250	91	2	5.491	4.0	2	3527	5.356
    Per_proteinISF_UTR_var (var name per protein UTR, category panel B UTR per protein) (for every unique protein coding isoform in a gene, how often UTR ragion varies while keeping protein sequqnce unchanged	71250	98027	41	1	1.376	1.0	1	58891	1.246
    Max_UTR_var_per_gene (for every unique protein coding isoform (if more than 1 listed in GT format), whats the maximum number of UTR varistions done on it, hence took the 1 represnetative isoform from gene undwerwent maximal UTR changes) 	12976	30505	41	1	2.351	2.0	1	6444	2.374
    Background_per_protein_UTR_var	1447	5296	29	2	3.66	3.0	2	654	2.751

    # basic Category\tOrganism\tValueMean\tValueStd\n
    '''
    for i in dat[1:]:
        # print (i, len(i.split('\t')))
        categ, gcount, Total, Max, Min, Mean, Median, Mode, M_count, Stddev = i.split('\t')
        if "Background" not in categ:
            categ = "0_Total_ISF" if "All_transcripts_statistics (including redundant" in categ else "1_NR_Isf" if "All_transcripts_statistics (only nonredundantPer gene" in categ else "2_UTR_per_NR_Isf" if "Per_proteinISF_UTR_var" in categ else "3_MaxUTRForIsf/gene"
            string += "%s\t%s\t%s\t%s\t%s\n"%(categ, organism, gcount, Mean, Stddev)
    return string


with open (inpfile) as fin:
    dat = fin.read().split('\n')
# grep Amongst the total UTR variants 5' 

hasTxid = {'10090': '3_Mouse', '9606':'4_Human', '7227':'1_Fly', '7955': '2_Fish', '6239':'0_Worm' }

stringA0_transData_supp = 'Category\tOrganism\tCount\tValueMean\tValueStd\n'
stringA1_GT_exonData_supp = 'Category\tOrganism\tCount\tValue\n'
stringA2_GT_transData_supp = 'Category\tOrganism\tCount\tValue\n'
stringA3_GT_exonData = 'Category\tTag\tOrganism\tCount\tValue\n'
stringA4_GT_transData = 'Category\tTag\tOrganism\tCount\tValue\n'
stringA5_GT_ContourExontransDataUniqueProteinCDS = 'Gene\tCodingExonCount\tNRIsfCount\tOrganism\n'


stringB1_exonData_meanStdCounts = 'Facet\tCategory\tOrganism\tValueMean\tValueStd\n'
stringB1_exonData_exonGeneFrac = 'Facet\tCategory\tOrganism\tValueFractionExons\tValueFractionGenes\n'

stringB2_exonData_meanStdCounts = 'Category\tOrganism\tValueMean\tValueStd\n'
stringB2_exonData_exonGeneFrac = 'Category\tOrganism\tValueFractionExons\tValueFractionGenes\n'

stringC_exonLength = 'Category\tOrganism\tValueMean\tValueMedian\tValueStd\n'

stringD_exonWEF = 'Category\tOrganism\tRange\tFrequency\tTotal\n'
stringE_exonWEF = 'Category\tOrganism\tValue\n'


for i in dat:
    if len(i) > 0:
        # extracting taxid from the filename address
        txid = re.search(r'/\d+[_]?[\d]?/derived_data', i).group()
        txid = re.sub(r'\/','',txid)
        txid = re.match (r'\d+', txid.split('_')[0]).group()
        org= hasTxid[txid]
        # print (txid, org)
        stringA1_GT_exonData_supp = compilegeneTableTransandExonData (os.path.join(i,fprefix, "tabular_GeneTableExonsDataProtCoding_IndependentCorrect.csv"), stringA1_GT_exonData_supp, org)
        stringA2_GT_transData_supp = compilegeneTableTransandExonData (os.path.join(i,fprefix, "tabular_GeneTableTransDataNR_IndependentCorrect.csv"), stringA2_GT_transData_supp, org)
        stringA3_GT_exonData = compilegeneTableTransandExonDataStacked (os.path.join(i,fprefix, "tabular_GeneTableExonsData.csv"), stringA3_GT_exonData, org)
        stringA4_GT_transData = compilegeneTableTransandExonDataStacked (os.path.join(i,fprefix, "tabular_GeneTableTransData.csv"), stringA4_GT_transData, org)
        stringA5_GT_ContourExontransDataUniqueProteinCDS = compilegeneTableTransandExonDataContour (os.path.join(i,fprefix, "raw_contourCodExNrIsf.csv"), stringA5_GT_ContourExontransDataUniqueProteinCDS, org)
        
        stringB1_exonData_meanStdCounts, stringB2_exonData_meanStdCounts = compileExonData (os.path.join(i,fprefix, "tabular_ExonsData.csv"), stringB1_exonData_meanStdCounts,stringB2_exonData_meanStdCounts, org)
        stringB1_exonData_exonGeneFrac, stringB2_exonData_exonGeneFrac = compileExonData (os.path.join(i,fprefix, "tabular_ExonsData.csv"), stringB1_exonData_exonGeneFrac, stringB2_exonData_exonGeneFrac, org, mode = "fraction")

        stringC_exonLength = compileExonLength (os.path.join(i,fprefix, "tabular_ExonLengthData.csv"), stringC_exonLength, org)

        stringD_exonWEF = compileWEF (os.path.join(i,fprefix, "Alternate_exons_WEF.csv"), stringD_exonWEF, org)
        stringE_exonWEF = compileWEF_raw (os.path.join(i,fprefix, "Alternate_exons_WEF_raw.csv"), stringE_exonWEF, org)

        stringA0_transData_supp = compileTrans (os.path.join(i,fprefix, "tabular_TransData.csv"), stringA0_transData_supp, org)


        
        
with open (os.path.join (outdir,"stringA0_transData_supp.tsv"),'w') as fout:
    fout.write ('%s'%(stringA0_transData_supp))

with open (os.path.join (outdir,"stringA1_GT_exonData_supp.tsv"),'w') as fout:
    fout.write ('%s'%(stringA1_GT_exonData_supp))

with open (os.path.join (outdir,"stringA2_GT_transData_supp.tsv"),'w') as fout:
    fout.write ('%s'%(stringA2_GT_transData_supp))

with open (os.path.join (outdir,"stringA3_GT_exonData.tsv"),'w') as fout:
    fout.write ('%s'%(stringA3_GT_exonData))

with open (os.path.join (outdir,"stringA4_GT_transData.tsv"),'w') as fout:
    fout.write ('%s'%(stringA4_GT_transData))

with open (os.path.join (outdir,"stringA5_GT_ContourExontransDataUniqueProteinCDS.tsv"),'w') as fout:
    fout.write ('%s'%(stringA5_GT_ContourExontransDataUniqueProteinCDS))



with open (os.path.join (outdir,"stringB1_exonData_meanStdCounts.tsv"),'w') as fout:
    fout.write ('%s'%(stringB1_exonData_meanStdCounts))

with open (os.path.join (outdir,"stringB1_exonData_exonGeneFrac.tsv"),'w') as fout:
    fout.write ('%s'%(stringB1_exonData_exonGeneFrac))

with open (os.path.join (outdir,"stringB2_exonData_meanStdCounts_additional.tsv"),'w') as fout:
    fout.write ('%s'%(stringB2_exonData_meanStdCounts))

with open (os.path.join (outdir,"stringB2_exonData_exonGeneFrac_additional.tsv"),'w') as fout:
    fout.write ('%s'%(stringB2_exonData_exonGeneFrac))

with open (os.path.join (outdir,"stringC_exonLength.tsv"),'w') as fout:
    fout.write ('%s'%(stringC_exonLength))

with open (os.path.join (outdir,"stringA0_transData_supp.tsv"),'w') as fout:
    fout.write ('%s'%(stringA0_transData_supp))

with open (os.path.join (outdir,"stringD_exonWEF.tsv"),'w') as fout:
    fout.write ('%s'%(stringD_exonWEF))

with open (os.path.join (outdir,"stringE_exonWEF.tsv"),'w') as fout:
    fout.write ('%s'%(stringE_exonWEF))





'''
tabular_GeneTableTransData
tabular_GeneTableExonsData
tabular_TransData
tabular_ExonsData
Alternate_exons_WEF.csv
'''
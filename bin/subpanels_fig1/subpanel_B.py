import modules_common as cm
import seaborn as sns
import scipy, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
def variations(geneob):
    var_gene_dual=[]
    retp1=[]
    retp2=[]
    def number_unique_exons(geneOrTrancript_ob):
        #will return the total number of unique non repetitive exons in the gene/tranciupt object 
        exons = [i for i in geneOrTrancript_ob.exons if i.length>0]
        uni_ex = {}
        for i in exons:
            #[UTMDR].[(-2)-n].[GAF].[1-n].[ncb0].[0-n]
            lett = ".".join(i.ID.split(".")[2:])
            uni_ex[lett] = 0

        return len(uni_ex)
    
    def internal_transition(trans):
        '''
        Well i dont know the purpose here
        What i can think is it traverses exons from left to right (and vice versa in next for loop) and start adding them to list till it enocunters first coding exon
        '''
        five_pi=[]
        three_pi=[]
        for i in trans.exons:
            if i.length==0:
                five_pi+=[i]
            else:
                five_pi+=[i]
                break
        for i in trans.exons[::-1]:
            if i.length==0:
                three_pi+=[i]
            else:
                three_pi+=[i]
                break
        return [five_pi,three_pi]
        
    ue=number_unique_exons(geneob)
    # as of now i am confused about the unique exons
    Main_has_clu={}
    #this will have unique aaseq of transcripts as key and lis of trancript objects as values
    Additional_exons_length_dual=[]
    # this will store tuple for every transcript (unique protein coding), which will carry information of the number of unique exons in trancfript and its length
    # exclude this from the analysis however 
    # -> below for loop will; do the analysis of the 
    for transcripts in geneob.transcripts:
        seqTranscriptaa=''.join(i.seq for i in transcripts.exons)
        if seqTranscriptaa not in Main_has_clu:
            Main_has_clu[seqTranscriptaa]=[]
            Additional_exons_length_dual += [
                    (number_unique_exons(transcripts), len(seqTranscriptaa))]
        Main_has_clu[seqTranscriptaa]+=[transcripts]
        if transcripts.PI:
            piseq=len(seqTranscriptaa)
    var_exons_pergene_dual = (len(Main_has_clu), ue)
    for clusters in Main_has_clu:
        # when the isoform remains same in cluster how ofetn do we see variation in the 5' non coding region and then the 3' region or in both cases
        BOTH=0;FIVE=0;THREE=0
        if len(Main_has_clu[clusters])>1:
            # checking if more than 1 isoform having same cluster
            refsides=internal_transition(Main_has_clu[clusters][0])
            # for cmment condition above, it gets list of exons from 5 and 3 prime, left to right and vice versa in second case till it encounters first coding exon (with a hope that some fraction of this may also code for the noncding region), for first transcript
            for j in Main_has_clu[clusters][1:]:
                transsides=internal_transition(j)
                #now check such list of exons for the another set of transcripts 
                if refsides[0]!=transsides[0] and refsides[1]!=transsides[1]:
                    BOTH+=1
                elif refsides[0]!=transsides[0]:
                    FIVE+=1
                elif refsides[1]!=transsides[1]:
                    THREE+=1
        retp2+=[(",".join([i.ID for i in Main_has_clu[clusters]]),len(Main_has_clu[clusters]),FIVE,THREE,BOTH)]
        #summ will alays be -1 of total as one is always hunted for thereferenec
    return ([piseq,ue,geneob.detail,Additional_exons_length_dual,var_exons_pergene_dual,[len(Main_has_clu[i]) for i in Main_has_clu],len(Main_has_clu)],retp2)
            #this lis will have UTR per prpteon and teher after niqie proteins



def freq_writer(lis):
    #getting list of total unique non redundant isoforms per gene as list elemmnt
    retstr=''
    for i in range(2, 11):
        retstr+="\t%s\t%s\t%s\n" %(i, lis.count(i), cm.div_fact(lis.count(i),len(lis)))
    ti = len([i for i in lis if i > 10])
    retstr+="\t11andMore\t%s\t%s\n" % (ti, cm.div_fact(ti,len(lis)))
    return retstr

def writer(total_transcripts, unique_protein_seq, per_protein_UTR, max_UTR_per_gene, res_dir):
    with open(res_dir+"Transcripts_and_UTR_per_gene.csv", "w") as fin:
        fin.write("Count,Category\n")
        for i in total_transcripts:
            fin.write("%s,All_transcripts\n" % i)
        for i in unique_protein_seq:
            fin.write("%s,Varying_Protein\n" % i)
        for i in per_protein_UTR:
            fin.write("%s,UTR_per_protein\n" % i)
        for i in max_UTR_per_gene:
            fin.write("%s,max_UTR_per_gene\n" % i)
def transcripts_type(has,fout,results_dir_csv):
    RAW_data=''
    RAW_data2=''
    total_transcripts = []
    background_UTR=[]
    unique_protein_seq = []
    CONDITION_GENES = {}
    per_protein_UTR = []
    max_UTR_per_protein = []
    genes_greater_than_3000len = []
    genes_gt1_isf_and_1_exon = []
    var_gene_dual = []
    exoncount_isf_length_allunique = []
    only1isf=0
    per_gene_exons=[]
    FIVETHREEBOTH=[0,0,0]
    for gene in has:
        info=variations(has[gene])

    #[pilen,ue,genename,Additional_exons_length_dual,var_exons_pergene_dual,total_var,unique_protein]
    #[commaspe_variants,len(Main_has_clu[clsters]),FIVE,THREE,BOTH)]
        if info[0][0] <= 3000: #length of PI to be studied shoud be <= 3000 amino acid
            if info[0][6] > 1:# numbero fo unique transcripts
                if info[0][1] > 1: # if more than 1 coding exon
                    per_gene_exons+=[info[0][1]]
                    total_transcripts += [sum(info[0][5])] # they are sum of total coding transcripts (redundant)
                    unique_protein_seq += [info[0][6]] # total nonredundant sequqnces
                    CONDITION_GENES[gene] = 1
                    var_gene_dual += [info[0][4]]
                    per_protein_UTR+=info[0][5]
                    # adding count of possible UTR variations for non redundant proteins, if more than 1 isoform is listed
                    max_UTR_per_protein += [max(info[0][5])]
                    for utr_variation in info[1]:
                        FIVETHREEBOTH[0]+=utr_variation[2]
                        FIVETHREEBOTH[1]+=utr_variation[3]
                        FIVETHREEBOTH[2]+=utr_variation[4]

                else:
                    genes_gt1_isf_and_1_exon += [(gene,
                                                  info[0][0], info[0][6], info[0][2])]
            else:
                if sum(info[0][5])>1:
                    background_UTR+=info[0][5]
                    #background UTR fraction only involves the cases when only 1 unique protein coding isoform is listed, this will be trated as the background fraction., this category should also have atleast 2 isoforms coding the same protein
                else:
                    only1isf+=1
    

        else:
            genes_greater_than_3000len += [(gene,info[0][0], info[0][2])]

        RAW_data+='%s\t%s\t%s\t%s\t%s\t%s\n'%(gene,info[0][2],info[0][0],sum(info[0][5]),info[0][6],info[0][1])
        for iteration in info[1]:
            RAW_data2+='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(gene,info[0][0],info[0][6],info[0][1],iteration[0],iteration[1],iteration[2],iteration[3],iteration[4])
        exoncount_isf_length_allunique+=info[0][3]
        
    special_genes=''
    for i in genes_greater_than_3000len:
        special_genes+="\t%s\n"%("\t".join(map(str, i)))
    
    linesep="============================>"
    array_to_write=[linesep]

    array_to_write+=["->Total_genes:%s" %len(has)]
    array_to_write+=["->Total genes with Pi length <=3000 are %s, and above are %s" % (len(has)-len(genes_greater_than_3000len), len(genes_greater_than_3000len))]
    array_to_write+=["->Their description is as follows\n\tGene\tLength\tName\n%s"% special_genes]
    array_to_write+=[linesep]
    
    array_to_write+=["-> Gene with PI length less than 3000 aa and only 1 listed isoform in genetable format (considering UTR variants but ignoring ncRNA and lncRNA) are : %s"% only1isf]

    array_to_write+=["-> Genes with more than 1 listed isoform but only single unique protein coding sequqnce are : %s"%len(background_UTR)]
    array_to_write+=["-> Genes with more than 1 unique protein coding isoform but only 1 exon are : %s"%len(genes_gt1_isf_and_1_exon)]
    array_to_write+=["-> Genes with PI len less than 3000 aa and atleast two different protein isoforms are :%s, amongst those, with atleast two coding exons are: %s" %(len(total_transcripts)+len(genes_gt1_isf_and_1_exon), len(CONDITION_GENES))]
    array_to_write+=[linesep]
    array_to_write+=["-> Number of transcripts in gene and frequency distribution of theirs and fraction of genes\n\tCDS_count\tGeneCount\tFrequency\n%s"% freq_writer(unique_protein_seq)]

    array_to_write+=["->Number of unique protein coding Exons in gene and fraction of genes\n\tDifferent_exons_count\tGeneCount\tFrequency\n%s"%freq_writer(per_gene_exons)]
    
    array_to_write+=["\n->Summarising\n\tCategory\tVar_count\tSum\tMax\tMin\tMean\tMedian\tMode\tMode_count\tStddev\n"]
    array_to_write+=["\tAll_transcripts_statistics (including redndant)\t%s"%cm.stats(total_transcripts)]
    array_to_write+=["\tAll_transcripts_statistics (only nonredundantPer gene)\t%s"%cm.stats(unique_protein_seq)]
    array_to_write+=["\tPer_proteinISF_UTR_var (var name per protein UTR, category panel B UTR per protein) (for every unique protein coding isoform in a gene, how often UTR ragion varies while keeping protein sequqnce unchanged\t%s"%cm.stats(per_protein_UTR)]
    array_to_write+=["\tMax_UTR_var_per_gene (for every unique protein coding isoform (if more than 1 listed in GT format), whats the maximum number of UTR varistions done on it, hence took the 1 represnetative isoform from gene undwerwent maximal UTR changes) \t%s"%cm.stats(max_UTR_per_protein)]
    array_to_write+=["\tBackground_per_protein_UTR_var\t%s"%cm.stats(background_UTR)]
    
    array_to_write+=["Amongst the total UTR variants 5' end got affected %s times, 3' end got affcted in :%s times and both ends got affected in %s times"%(cm.div_fact(FIVETHREEBOTH[0],sum(FIVETHREEBOTH)),cm.div_fact(FIVETHREEBOTH[1],sum(FIVETHREEBOTH)),
        cm.div_fact(FIVETHREEBOTH[2],sum(FIVETHREEBOTH)))]
    
    array_to_write+=["-->Wilcoxon_ranksums_between background per protein UTR var (when only 1 unique protein coding isoform is be listed along with UTR variations) and per protein Set (when more than 1 unique protein coding isoform is listed along with their UTR variations) var\t%s"%str(scipy.stats.ranksums(background_UTR, per_protein_UTR))]
    array_to_write+=["-->Wilcoxon_ranksums_between background per protein UTR var (when only 1 unique protein coding isoform is be listed (background_UTR) along with UTR variations) and MAx_UTR_per_gene (getting single variants undergo maximal changes in UTR from list of isoforms in gene) (var name max UTR per protein) in Set\t%s"%str(scipy.stats.ranksums(background_UTR, max_UTR_per_protein))]
    array_to_write+=["-->CorRelation (spearman) in (for every gene, number of nonred prptein coding isoforms it does have and list b as number of unique coding exon in that gene, var name var_gene_dual) : %s"%str(scipy.stats.spearmanr(a=[i[0] for i in var_gene_dual], b=[i[1] for i in var_gene_dual]))]
    array_to_write+=["-->CorRelation (spearman) in (for every unique protein coding isoform of a gene, get list a as, number of unique protein coding exons in it and then list b length of such sioform, var name exoncount_isf_length_allunique ) : %s"%str(scipy.stats.spearmanr(a=[i[0] for i in exoncount_isf_length_allunique],b=[i[1] for i in exoncount_isf_length_allunique]))]
    array_to_write+=["\ntag:All_transcripts ->total_transcripts"]
    array_to_write+=["tag:Varying_Protein ->unique_protein_seq"]
    array_to_write+=["tag:UTR_per_protein ->per_protein_UTR"]
    array_to_write+=["tag:max_UTR_per_gene ->max_UTR_per_gene"]


    fout.write("%s\n"%("\n".join(array_to_write)))

    
    with open (os.path.join(results_dir_csv,'genewiseinformation1.tab'),'w') as fin:
        fin.write("Gene\tGeneName\tLengthPI\tTotalVariants\tUniqueProteinVariants\tUniqueCodingExons\n")
        fin.write("%s"%RAW_data)
    with open (os.path.join(results_dir_csv,'genewiseinformation2_cluster_2.tab'),'w') as fin:
        #gene,pilen,cluster1var,totalUTR,5'3'both
        fin.write("Gene\tLengthPI\tUniqueProteinTranscripts\tCoding_exons_gene\tVariants(,)sep\tTotalUTR\t5'\t3'\tBoth\n")
        fin.write("%s"%RAW_data2)

    df = pd.DataFrame(var_gene_dual,columns=['Transcripts','Exons'])
    '''
    reg_plot=sns.jointplot(x=df["Exons"], y=df["Transcripts"], kind='reg',space=0,joint_kws={'line_kws':{'color':'cyan'}})
    reg_plot.ax_joint.legend_.remove()
    reg_plot.savefig(os.path.join(results_dir_csv,"Additional_transcripts_exons_varname(VAR_GENE_DUAL).png"))
    
    df1 = pd.DataFrame(exoncount_isf_length_allunique,columns=['Exons','PI_length'])
    reg1_plot=sns.jointplot(x=df1["Exons"], y=df1["PI_length"], kind='reg',space=0, color="g",joint_kws={'line_kws':{'color':'cyan'}})
    reg1_plot.ax_joint.legend_.remove()
    reg1_plot.savefig(os.path.join(results_dir_csv,"Additional_number_of_exons_and_protein_length_unique_dataset_all_varNAme(exoncount_isf_length_allunique).png"))
    '''
    writer(total_transcripts, unique_protein_seq, per_protein_UTR, max_UTR_per_protein, results_dir_csv)

    return CONDITION_GENES 


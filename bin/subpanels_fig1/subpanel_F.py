
def stats_new(lis):
    #print "1",len(lis),sum(lis),np.max(lis),np.min(lis),round(np.mean(lis),3),round(np.median(lis),3),scipy.stats.mode(lis)[0][0],scipy.stats.mode(lis)[1][0],round(np.std(lis),3)
  try:
    return "\t".join(map(str,[len(lis),sum(lis),np.max(lis),np.min(lis),round(np.mean(lis),3),round(np.median(lis),3),scipy.stats.mode(lis)[0][0],scipy.stats.mode(lis)[1][0],round(np.std(lis),3)]))
  except Exception as E:
    return 'False'
def change_in_length(has,res_dir,genes_cond,output_stats):
    change_XP=[]
    change_NP=[]
    change_all=[]
    for gene in has:
        if gene in genes_cond:
            has_seq={}
            pi_seq="".join([j.seq for j in has[gene].PI.exons])
            has_seq[pi_seq]=0
            pi_len=len(pi_seq)
            #print pi_seq
            #break
            for transcripts in has[gene].transcripts:
                seq="".join([j.seq for j in transcripts.exons])
                if seq not in has_seq:
                    has_seq[seq]=0
                    val = div_fact((pi_len - len(seq)),pi_len)
                    if val>90:
                        pass
                        '''
                        print has[gene].detail
                        print pi_len
                        print len(seq)
                        print val
                        '''
                    change_all += [val]
                    if transcripts.ID[0]=="X":
                        change_XP += [val]
                    else:
                        change_NP += [val]
    with open (res_dir+"csv/general/change_in_isoform_length.csv","w") as fin:
        fin.write("ChangeinLength,Type")
        for i in change_NP:
            fin.write("%s,NP\n"%i)
        for i in change_XP:
            fin.write("%s,XP\n"%i)
        for i in change_all:
            fin.write("%s,All\n"%i)
    #output_stats = open("General_stats_fig1.log", "a")
    output_stats.write("\n\n=============>\n")
    output_stats.write("Change in transcript length from principal isoform\n")
    output_stats.write("\n\tCategory\tVar_count\tSum\tMax\tMin\tMean\tMedian\tMode\tMode_count\tStddev")
    output_stats.write("\n\tChangeinNP\t%s"%stats_new(change_NP))
    output_stats.write("\n\tChangeinXP\t%s"%stats_new(change_XP))
    output_stats.write("\n\tChangeinall\t%s"%stats_new(change_all))
    output_stats.write("\nCSV_File:%s"%res_dir+"change_in_isoform_length.csv")
    output_stats.close()
change_in_length(has,results_dir_csv,genes_cond,filewriter)

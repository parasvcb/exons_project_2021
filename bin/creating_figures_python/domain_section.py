
import cPickle as pickle
import sys
if len(sys.argv)!=4:
	print ("Please type correct location of results dir having object stored as well as organism id  and pfam/cath as three different arguements")
	sys.exit()
results_dir = sys.argv[1]
organism_id = sys.argv[2]
source_gene_object = sys.argv[1]+'objectsave_%s.pick'%organism_id
CONST_val=sys.argv[3]

print CONST_val
import cPickle as pickle
import pandas as pd
import numpy as np
import re, scipy.stats, os
#CONST_val='pfam'

with open(source_gene_object) as fin:
    human=pickle.load(fin)
with open(results_dir+'new_condition_genes.pick') as fin:
    filterg=pickle.load(fin)

def div_fact(num,denom):
	try:
		return round(float(num)/denom,3)	
	except:
		return 0

fin=open(results_dir+'csv/domains/%s_numbers.txt'%CONST_val,'w')

print results_dir+'csv/domains/%s_numbers.txt'%CONST_val  
def isoform_giver(gene_ob):
    lis=gene_ob.PI.cath_list if CONST_val=='cath' else gene_ob.PI.pfam_list
    if lis:
        relev_exons=[i for i in gene_ob.PI.exons if i.length>0]
        if CONST_val=='cath':
            relev_doms=[]
            for i in lis:
                flag=1
                for j in i[2]:
                    if j<0.8:
                        flag=0
                if flag:
                    relev_doms+=[i]
        else:
            relev_doms=[i for i in gene_ob.PI.pfam_list if i[2]>=0.5]
            relev_doms=relev_doms if relev_doms else False
        return relev_exons,relev_doms
    else:
        return False,False


def quadrant_creator(lis,filename):
    cols = [1.01, 0.75, 0.50, 0.25, 0]
    rows = cols[::-1][:]
    dat = np.zeros((4, 4), dtype=int)
    for i in lis:
        ri = 0
        ci = 0
        for exi in range(0, len(cols)-1):
            if cols[exi+1] <= i[1] < cols[exi]:
                ci = exi
                break
        for dmi in range(0, len(rows)-1):
            if rows[dmi] <= i[0] < rows[dmi+1]:
                ri = dmi
                break
        dat[ci, ri] += 1

    def norma(z):
        return div_fact(z,np.ndarray.sum(dat))

    ndat = np.vectorize(norma)(dat)
    dfnorm = pd.DataFrame(
        ndat, columns=[1.00, 0.75, 0.50, 0.25][::-1], index=[1.00, 0.75, 0.50, 0.25])
    dfnorm.to_csv(filename)



'''
Ths function will give me the follwing information:
count the number fo domains
RAW: gene, PI, number of domains

will count the number of alt exons, and cont exons,
RAW: gene, PI, const, ALt

will count WEF of alt exons also, (wo domain knowledge)

if exon and domain intersects
prepare a template of hash with exids as keys(): fitst val list as fractiona nd sceond as domain id
if intersecton is 100%, add to the list
this iwll give us the number that out of total exons, how many code for nothin
ho wmany code for 100% and manya times

repea this with merged exons

so amongst total exons how nany do code for domains 100% of ties, and how many dont code for nothing
and how many codes fr spluit cases

split cases will have this categorization:
AA junction, GG junction and AG junction
took AG and prticeed
if not 100%, count the fractionof alt and const residue sion whole prtion and count them in middle 40%
prepare two sets, 

will add the exons to a list

'''
def files_writer(logic,lis,filename,columns='None'):
    if logic==1:
        with open(filename, "w") as fin:
            fin.write("Gene,Domain_fraction,Exon_fraction\n")
            for i in lis:
                fin.write("%s,%s,%s\n" % (i[0], i[1],i[2]))
    if logic ==2:
        with open(filename,"w") as fin:
            fin.write("%s\n"%(",".join(columns)))
            for i in lis:
                fin.write("%s\n"%(",".join(map(str,i))))
def update_hash(has,val,gene):
    #print has
    if val in has:
        has[val]+=[gene]
    else:
        has[val]=[gene]
    return has
def ret_sort(lis):
            return lis[1]
def set_retuner(logic,has):
    if logic==1:
        
        mat=[]
        mega_total=sum([len(i) for i in has.values()])
        
        for i in has:
            dom=i
            semi_total=len(has[i])
            tot_freq=div_fact(semi_total,mega_total)
            G_case=has[i].count('G')
            A_case=has[i].count('A')
            mat+=[(i,semi_total,tot_freq,G_case,A_case)]
        mat=sorted(mat,key=ret_sort,reverse=True)
        set_a=set([i[0] for i in mat])
        return mat,set_a
    else:
        
        mat=[]
        gene=[j for i in has.values() for j in i]
        for i in has:
            dom=i
            semi_total=len(has[dom])
            freq=div_fact(semi_total,len(gene))
            mat+=[(dom,semi_total,freq)]
        mat=sorted(mat,key=ret_sort,reverse=True)
        set_c=set(i[0] for i in mat)
        gene=set(gene)
        return mat,set_c,gene

def mini_cdf(lis1,lis2,filename,typi1,typi2):
    total_t1 = len(lis1)
    total_t2 = len(lis2)
    temp_binst = np.arange(0, 1.001, 0.001, dtype=float)
    with open (filename,"w") as fin:
        fin.write("Domain_fraction,Total,Freq,Group\n")
        for i in temp_binst:
            valst1 = 0
            valst2 = 0
            for j in lis1:
                if j <= i:
                    valst1 += 1
            for j in lis2:
                if j <= i:
                    valst2 += 1
            fin.write("%s,%s,%s,%s\n" %
                      (i, valst1, div_fact(valst1,total_t1),typi1))
            fin.write("%s,%s,%s,%s\n" %
                      (i, valst2, div_fact(valst2,total_t2),typi2))

def fraction_middle_whole(has,filename):
    total=len([j for i in has.values() for j in i])
    keys=['GG','AG','AA']
    strfil='Residuetypes\tTotal\tFreq\n'
    for i in keys:
        strfil+='%s\t%s\t%s\n'%(i,len(has[i]),div_fact(len(has[i]),total))
    A_freq=[i[0] for i in has['AG']]
    G_freq=[i[1] for i in has['AG']]
    mini_cdf(A_freq,G_freq,filename+"AG_CDF.csv","A","G")
    with open (filename+"_general.csv","w") as fin:
        fin.write("%s"%(re.sub(r'\t',',',strfil)))
    return strfil

def domain_ontology_file(mat,filename,columns,set_verified):
    with open (filename,'w') as fin:
        fin.write("%s\n"%"\t".join(columns))
        for i in mat:
            flag='True' if i[0] in set_verified else 'False'
            fin.write("%s\t%s\t%s\n"%(i[0],flag,"\t".join(map(str,i[1:]))))

        
def domain_exon_fraction(human, fout):
    gene_exons_fraction={}
    #this will give me the number of const exons, number of alt exopns, 
    #how many of them code for 100% of domains, and how mnay domains are completely
    #contained within the exons
    #how mnay code for null
    # if split domains, thi inofmration will be tased outr by domain whole and the whole middle
    gene_exon_wef={}
    domain_whole={'AG':[],'AA':[],'GG':[]}
    domain_middle={'AG':[],'AA':[],'GG':[]}
    domain_inf_middle={}
    domain_inf_whole={}
    dom_ex_frac = []
    dom_ex_frac_const = []
    total_dom = []
    total_intersec = 0
    total_aa = 0
    total_aa_dom = 0
    total_doms_background={}
    c = 0
    per_prot_dom=[]
    for gene in human:
        if gene in filterg:
            relev_ex,domlis=isoform_giver(human[gene])
            if domlis:
                total_dom += [(gene,len(domlis))]
                total_aa += sum([ex.length for ex in relev_ex])
                d_aa = 0
                gene_exon_wef[gene]={i.ID:i.WEF for i in relev_ex}
                gene_exons_fraction[gene]={i.ID:[] for i in relev_ex}
                for dom in domlis:
                    length={'A':0,'G':0}
                    length_m={'A':0,'G':0}
                    domain_constituents=''
                    domain_constituents_middle=''
                    flag_const=0
                    flag_const_m=0
                    span_temp=range(dom[0][0], dom[0][1]+1)
                    dspan = set(span_temp)
                    dspan_middle=set(span_temp[int(len(span_temp)*0.3):int(len(span_temp)*0.6)])
                    d_aa += len(dspan)
                    su = 0
                    for ex in relev_ex:
                        espan = set(range(su, su+ex.length))
                        su += ex.length
                        interori = espan & dspan
                        intermidd = espan & dspan_middle
                        if interori:
                            total_intersec += 1
                            dfrac = div_fact(len(interori),len(dspan))
                            
                            gene_exons_fraction[gene][ex.ID]+=[(dfrac, dom[1]+':'+dom[3])]
                            #push this up to see varying doimain definition for singl exns
                            efrac = div_fact(len(interori),len(espan))
                            dom_ex_frac += [(gene,dfrac, efrac)]
                            if ex.ID[0]!='R' and ex.ID.split('.')[2] == 'G':
                                 dom_ex_frac_const += [(gene,dfrac, efrac)]
                            if dfrac<1:
                                flag_const=1
                                key=ex.ID.split('.')[2] if ex.ID[0]!='R' else 'A'
                                key='G' if key == 'F' else key
                                domain_constituents+=key
                                length[key]+=len(interori)
                        
                        if intermidd:
                            dfrac_m=div_fact(len(intermidd),len(dspan_middle))
                            if dfrac_m<1:
                                flag_const_m=1
                                key_m=ex.ID.split('.')[2] if ex.ID[0]!='R' else 'A'
                                key_m='G' if key_m == 'F' else key_m
                                domain_constituents_middle+=key_m
                                length_m[key]+=len(intermidd)
                                
                    if flag_const:
                        constituents=re.sub(r'F','G',domain_constituents)
                        key="".join(list(set(constituents)))
                        key = key*2 if len(key)==1 else key
                        #print key, length
                        domain_whole[key]+=[(div_fact(length[key[0]],sum(length.values())), div_fact(length[key[1]],sum(length.values())))]
                        domain_inf_whole = update_hash(domain_inf_whole,dom[1]+':'+dom[3],gene)
                    if flag_const_m:
                        constituents_m=re.sub(r'F','G',domain_constituents_middle,gene)
                        #domain_middle["".join(list(set(constituents_m)))]+=1
                        key="".join(list(set(constituents_m)))
                        key = key*2 if len(key)==1 else key
                        domain_middle[key]+=[(div_fact(length_m[key[0]],sum(length_m.values())), div_fact(length_m[key[1]],sum(length_m.values())))]
                        domain_inf_middle = update_hash(domain_inf_middle,dom[1]+':'+dom[3],gene)
                    total_doms_background = update_hash(total_doms_background,dom[1]+':'+dom[3],gene)
                        
                                
                total_aa_dom += d_aa
                per_prot_dom+=[div_fact(total_aa_dom,total_aa)]
    #print results_dir+"csv/domains/%s/RAW_dual_fraction.csv"%CONST_val
    #print sum(total_dom), total_intersec
    # print "1111"
    '''
    files_writer: 1: Raw mode dual frac columns='None'
                  2: Raw mode single columns with columns
                  3:  
    '''
    dom100per={}
    gene_100=[]
    gene_mid=[]
    exons_count_comp_domains={'A':[],'G':[]}
    comp_domains_count={'A':0,'G':0}
    for gene in gene_exons_fraction:
        for exons in gene_exons_fraction[gene]:
            key=exons.split('.')[2] if exons[0]!='R' else 'A'
            key='G' if key == 'F' else key
            for spans in gene_exons_fraction[gene][exons]:
                if spans[0]==1:
                    if spans[1] in dom100per:
                        dom100per[spans[1]]+=key
                    else:
                        dom100per[spans[1]]=key
                    gene_100+=[gene]
                    exons_count_comp_domains[key]+=[(gene,exons)]
                    comp_domains_count[key]+=1
                else:
                    gene_mid+=[gene]
    exons_count_comp_domains={i:len(set(exons_count_comp_domains[i])) for i in exons_count_comp_domains}
                    
    gene_100=set(gene_100)
    gene_mid=set(gene_mid)
    mat_100,dom_100_set=set_retuner(1,dom100per)
    mat_midd,dom_middle_set,gene_middle=set_retuner(0,domain_inf_middle)
    mat_whole,dom_whole_set,gene_whole=set_retuner(0,domain_whole)
    mat_background,dom_all_b,gene_all_b=set_retuner(0,total_doms_background)
    conf_100=dom_100_set-dom_middle_set
    conf_middle=dom_middle_set-dom_100_set
    conf_NC_split=dom_whole_set-dom_middle_set-dom_100_set
    # print mat_whole
    domain_ontology_file(mat_100,results_dir+'csv/domains/%s/Domain_ontology_always_contained_cases.csv'%CONST_val,["Domain","Never_splitted","Total","Freq","C","G"],conf_100)
    domain_ontology_file(mat_whole,results_dir+'csv/domains/%s/Domain_ontology_gota_split.csv'%CONST_val,["Domain","Always_split_ony_NC","Total","Freq"],conf_NC_split)
    domain_ontology_file(mat_midd,results_dir+'csv/domains/%s/Domain_ontology_middle_cases.csv'%CONST_val                         ,["Domain","Always_middle_split","Total","Freq"],conf_middle)
    domain_ontology_file(mat_background,results_dir+'csv/domains/%s/Domain_ontology_background_cases.csv'%CONST_val                         ,["Domain","FALSE","Total","Freq"],set())
    str_mid=fraction_middle_whole(domain_middle,results_dir+"csv/domains/%s/Fraction_middledomain_40.csv"%CONST_val)
    str_whole=fraction_middle_whole(domain_whole,results_dir+"csv/domains/%s/Fraction_wholedomain.csv"%CONST_val)

    #prepare the histogranms and CDF in whle and in middle region between A and G 
    #and count the cases
    
                
    files_writer(1,dom_ex_frac,results_dir+"csv/domains/%s/RAW_dual_fraction_indi.csv"%CONST_val)
    files_writer(1,dom_ex_frac_const,results_dir+"csv/domains/%s/RAW_dual_fraction_indi_const.csv"                        %CONST_val)
    files_writer(2,total_dom,results_dir+'csv/domains/%s/RAW_Domain_Distribution.csv'%CONST_val,columns=['Gene','Domains'])
    quadrant_creator(dom_ex_frac, results_dir+'csv/domains/%s/quadrant.csv'%CONST_val)
    quadrant_creator(dom_ex_frac_const, results_dir+'csv/domains/%s/quadrant_const.csv'%CONST_val)
    
    dom_bins = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,'8+': 0}
    for i in total_dom:
        flag=0
        for j in range(1,8):
            if i[1]==j:
                dom_bins[j]+=1
                flag=1
                break
        if not flag:
            if i[1]>=8:
                dom_bins['8+']+=1
    
            
    val_file = ''
    with open(results_dir+'csv/domains/%s/Domain_Distribution.csv'%CONST_val, 'w') as fin:
        fin.write("Domains,Total,Freq\n")
        val_file += "Domains\tTotal\tFreq\n"
        vals = sum(dom_bins.values())
        keys = dom_bins.keys()
        keys.sort()
        keys.remove('8+')
        keys += ['8+']
        for i in keys:
            val_file += "%s\t%s\t%s\n" % (i,
                                          dom_bins[i], div_fact(dom_bins[i],vals))
            fin.write("%s,%s,%s\n" %
                      (i, dom_bins[i], div_fact(dom_bins[i],vals)))
    total_E=0
    alt_E=0
    const_E=0
    for gEne in gene_exon_wef:
        for ex in gene_exon_wef[gEne]:
            if ex[0]=='R' or ex.split(".")[2]=='A':
                alt_E+=1
            else:
                const_E+=1
            total_E+=1
    whole_dom_count=sum(i[1] for i in mat_100)
    primer_string = '''Domain_statistics\nTotal number of genes in condition set is:%s
    Genes that have domain predicted for principal isoform in that set is :%s
    Total domains in those genes = %s
    Amongst those, the domain distribution is present in plot and here
    %s\n
    Total_coverage = Total residues in domain/ Total residues in aaseq = %s/%s = %s\n
    
    Per_PI_Protein cov
    Mean:%s, median:%s stddev:%s

    Exon domain fraction dual quadrant in present in CSV and in plot
    
    Total number of exons in transcripts anaysed : %s
    Total_alt_exons : %s (%s)
    Total_const_exons : %s (%s)
    
    Amongst  the total number of domains
    Domains are wholly contained in exons in : %s times (%s freq)
    %s constitutive exons contain %s domains (%s) and 
    %s alternate exons contain %s domains (%s)
    This behaviopur is present in %s genes (COM)
    
    when domains are fractioned, %s genes have this phenomenon (MID)
    COM INTERSECT MID = %s
    
    when domans are being split by more than one exon
    following is the result
    
    when whole is affected:
    %s
    When middle is affcted: 
    %s
    ''' % (len(filterg), len(total_dom), sum([i[1] for i in total_dom]), val_file, total_aa_dom, total_aa,\
           div_fact(total_aa_dom,total_aa), np.mean(per_prot_dom),np.median(per_prot_dom),\
           np.std(per_prot_dom),\
          total_E,alt_E,div_fact(alt_E,total_E),const_E,div_fact(const_E,total_E),\
          whole_dom_count, div_fact(whole_dom_count,sum([i[1] for i in total_dom])),\
          exons_count_comp_domains['G'],comp_domains_count['G'], div_fact(comp_domains_count['G'],whole_dom_count),\
          exons_count_comp_domains['A'],comp_domains_count['A'], div_fact(comp_domains_count['A'],whole_dom_count),\
          len(gene_100),len(gene_mid),len(gene_mid&gene_100), str_whole,str_mid)

    fout.write('%s' % primer_string)
    #fout.close()
domain_exon_fraction(human=human, fout=fin)


# ### domain view story is presented:
# 

# In[18]:


'''exons_side_story'''
def CDFandHIST(has, filename_prefix,lis_a):
    #print 'fil',filename_prefix,fg
    #print has
    #fil_str=''
    def mini_cdf(lis, filename, typi):
        temp_binst = np.arange(0, 1.001, 0.001, dtype=float)
        total_t = len(lis)
        #print filename
        #print max(lis),min(lis)
        with open (filename,"w") as fin:
            fin.write("Domain_fraction,Total,Freq,Group\n")
            for i in temp_binst:
                valst = 0
                for j in lis:
                    if j <= i:
                        valst += 1
                fin.write("%s,%s,%s,%s\n" %
                          (i, valst, div_fact(valst,total_t),typi))

    def mini_raw(lis, filename):
        with open(filename, 'w') as fin:
            fin.write('Domain_Fraction\n')
            for i in lis:
                fin.write('%s\n' % i)
    
    
    has_a_freq={(0,0.25):{(0,0.25):0,(0.25,0.50):0,(0.50,0.75):0,(0.75,1.01):0},
               (0.25,0.50):{(0,0.25):0,(0.25,0.50):0,(0.50,0.75):0,(0.75,1.01):0},
               (0.50,0.75):{(0,0.25):0,(0.25,0.50):0,(0.50,0.75):0,(0.75,1.01):0},
               (0.75,1.0):{(0,0.25):0,(0.25,0.50):0,(0.50,0.75):0,(0.75,1.01):0},
               (1.0,1.1):{(0,0.25):0,(0.25,0.50):0,(0.50,0.75):0,(0.75,1.01):0}}
    has_a_freq_wo_dis_dom_part={i:0 for i in has_a_freq.keys() if i!=(1.0,1.1)}
    
    with open (filename_prefix+'RAW_alt_fraction_and_wef.csv','w') as fin:
        fin.write("Domain_fraction,WEF\n")
        for i in lis_a:
            fin.write("%s,%s\n"%(i[0],i[1]))
            for keys in has_a_freq:
                if keys[0]<=i[0]<keys[1]:
                    #has_a_freq[keys]['sum']+=1
                    for keysin in has_a_freq[keys]:
                        if keysin[0]<=i[1]<keysin[1]:
                            has_a_freq[keys][keysin]+=1
                            break
                    break
            for keysdist in has_a_freq_wo_dis_dom_part:
                    if keysdist[0]<=i[1]<keysdist[1]:
                        has_a_freq_wo_dis_dom_part[keysdist]+=1
                        break
    
    alt_ex_dom_part=''
    with open (filename_prefix+'ALt_exons_dom_participation_trend.csv','w') as fin:
        fin.write('WEF_range,Total,Freq\n')
        alt_ex_dom_part='WEF_range\tTotal\tFreq\n'
        keys=has_a_freq_wo_dis_dom_part.keys()
        total=sum(has_a_freq_wo_dis_dom_part.values())
        keys.sort()
        for i in keys:
            fin.write('%s,%s,%s\n'%("_".join(map(str,list(i))),has_a_freq_wo_dis_dom_part[i], div_fact(has_a_freq_wo_dis_dom_part[i],total)))
            alt_ex_dom_part+='%s\t%s\t%s\n'%("_".join(map(str,list(i))),has_a_freq_wo_dis_dom_part[i], div_fact(has_a_freq_wo_dis_dom_part[i],total))
    #print filename_prefix
    #print scipy.stats.spearmanr(a=[i[0] for i in lis_a], b=[i[1] for i in lis_a])
    #print np.corrcoef([i[0] for i in lis_a],[i[1] for i in lis_a])
                    
    with open (filename_prefix+'HIST_A_with_WEF.csv','w') as fin:
        min_str_a_ret='FractionDomain\tTotal\tFreq\t0-0.25\t0.25,0.50\t0.50-0.75\t0.75-1.01\n'
        fin.write('FractionDomain,Total,Freq,A,B,C,D\n')
        key = has_a_freq.keys()
        key.sort()
        mega_sum=sum([sum(has_a_freq[i].values()) for i in has_a_freq ])
        for i in key:
            semi_tot = sum(has_a_freq[i].values())
            subkeys=has_a_freq[i].keys()
            subkeys.sort()
            mini_str=''
            for j in subkeys:
                mini_str+=',%s'%str(div_fact(has_a_freq[i][j],semi_tot))
            
            fin.write("%s,%s,%s%s\n" % ('_'.join(map(str, list(i))), semi_tot,
                                             div_fact(semi_tot,mega_sum),mini_str))
            min_str_a_ret+="%s\t%s\t%s%s\n" % ('_'.join(map(str, list(i))), semi_tot,
                                             div_fact(semi_tot,mega_sum),re.sub(r",",r"\t",mini_str))
    
    
    mega = [j for i in has.values() for j in i]
    has_bins = {(0, 0.25): '', (0.25, 0.50): '',
                (0.50, 0.75): '', (0.75, 1.0): '',(1.0,1.1):''}
    # a=10
    #print has
    for i in has['A']:
        for key in has_bins:
            if key[0] <= i < key[1]:
                has_bins[key] += 'A'
    for i in has['F']:
        for key in has_bins:
            if key[0] <= i < key[1]:
                has_bins[key] += 'F'
    for i in has['G']:
        for key in has_bins:
            if key[0] <= i < key[1]:
                has_bins[key] += 'G'
                

    with open(filename_prefix+'_HIST_AFG.csv', 'w') as fin:
        min_str=''
        fin.write('FractionDomain,Total,Freq,A,F,G\n')
        key = has_bins.keys()
        key.sort()
        total = len(''.join(has_bins.values()))
        for i in key:
            
            semi_tot = len(has_bins[i])
            fin.write("%s,%s,%s,%s,%s,%s\n" % ('_'.join(map(str, list(i))), semi_tot,
                                             div_fact(semi_tot,total), div_fact(has_bins[i].count('A'),semi_tot),
                                             div_fact(has_bins[i].count('F'),semi_tot),
                                             div_fact(has_bins[i].count('G'),semi_tot)))
            min_str+="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%('_'.join(map(str, list(i))), semi_tot,
            div_fact(semi_tot,total),has_bins[i].count('A'),div_fact(has_bins[i].count('A'),semi_tot),
            has_bins[i].count('F'),div_fact(has_bins[i].count('F'),semi_tot),
            has_bins[i].count('G'),div_fact(has_bins[i].count('G'),semi_tot))
    fil_str='''%s section
    the following information is presented
    %s
    FractionDomain\tTotal\tFreq\tA\tAfreq\tF\tFreq\tG\tGfreq
    %s
    
    Alternate exons exclusive domain participation and WEF is shown following:
    %s
    
    Another perspective, all alterate exons that particiaptes in domains, and thier frequency distribution
    irrespetive of fraction covered but, they do intersects
    %s
    '''%(filename_prefix,filename_prefix+'_HIST_AFG.csv',min_str, min_str_a_ret, alt_ex_dom_part)
   
    mini_cdf(has['A'], filename_prefix+'_CDF_A.csv','A')
    mini_cdf(has['F'], filename_prefix+'_CDF_F.csv','F')
    mini_cdf(has['G'], filename_prefix+'_CDF_G.csv','G')
    mini_cdf(mega, filename_prefix+'_CDF_All.csv','All')
    mini_raw(has['A'], filename_prefix+'_RAW_A.csv')
    mini_raw(has['F'], filename_prefix+'_RAW_F.csv')
    mini_raw(has['G'], filename_prefix+'_RAW_G.csv')
    mini_raw(mega, filename_prefix+'_RAW_All.csv')
    with open (filename_prefix+'_CDF_groups.csv','w') as fin:
        with open (filename_prefix+'_CDF_All.csv') as f1:
            dat="\n".join([i for i in f1.read().split('\n') if len(i)>0])
        fin.write("%s\n"%dat)
        with open (filename_prefix+'_CDF_A.csv') as f1:
            dat="\n".join([i for i in f1.read().split('\n')[1:] if len(i)>0])
        fin.write("%s\n"%dat)
        with open (filename_prefix+'_CDF_F.csv') as f1:
            dat="\n".join([i for i in f1.read().split('\n')[1:] if len(i)>0])
        fin.write("%s\n"%dat)
        with open (filename_prefix+'_CDF_G.csv') as f1:
            dat="\n".join([i for i in f1.read().split('\n')[1:] if len(i)>0])
        fin.write("%s\n"%dat)
    os.remove(filename_prefix+'_CDF_All.csv')
    os.remove(filename_prefix+'_CDF_A.csv')
    os.remove(filename_prefix+'_CDF_F.csv')
    os.remove(filename_prefix+'_CDF_G.csv')
    return fil_str


# In[19]:


def np_to_csv(arr, filename, cols):
    fil_str=''
    values = arr[0]
    bins = arr[1]
    su = sum(values)
    with open(filename, "w") as fin:
        fin.write("%s\n" % (",".join(cols)))
        fil_str+='%s\n'%("\t".join(cols))
        for i in range(0, len(bins)-1):
            fin.write("%s,%s,%s\n" %
                      ("-".join(map(str, [round(bins[i], 2), round(bins[i+1], 2)])), values[i], div_fact(values[i],su)))
            fil_str+="%s\t%s\t%s\n"%("-".join(map(str, [round(bins[i], 2), round(bins[i+1], 2)])), values[i], div_fact(values[i],su))
    return fil_str


# In[20]:


def exons_view_ind(human,fout):
    exon_count_per_domain = []
    longest_fraction = {'A': [], 'G': [], 'F': []}
    longest_fraction_alt_f=[]
    exon_fraction = {'A': [], 'G': [], 'F': []}
    exon_fraction_alt_f=[]
    alt_exons_wef_background=[]
    for gene in human:
        if gene in filterg:
            relev_ex,domlis=isoform_giver(human[gene])
            if domlis:
                alt_exons_wef_background+=[ex.WEF for ex in relev_ex if ex.ID[0]!='R' and ex.ID.split('.')[2]=='A']
                for dom in domlis:
                    mat = []
                    dspan = set(range(dom[0][0], dom[0][1]))
                    su = 0
                    for ex in relev_ex:
                        espan = set(range(su, su+ex.length))
                        su += ex.length
                        inter = espan & dspan
                        if inter:
                            key = 'A' if ex.ID[0] == 'R' else ex.ID.split('.')[
                                2]
                            ex_frac = div_fact(len(inter),len(dspan))
                            mat += [(len(inter),ex.length,ex.WEF,ex_frac, key)]
                    exon_count_per_domain += [len(mat)]
                    mat.sort()
                    if mat:
                    	longest_fraction[mat[-1][-1]] += [mat[-1][-2]]
                    	if mat[-1][-1]=='A':
                        	longest_fraction_alt_f+=[(mat[-1][-2],mat[-1][2])]
                    	for i in range (0,len(mat)):
                        	exon_fraction[mat[i][-1]] += [mat[i][-2]]
                        	if mat[i][-1]=='A':
                            		exon_fraction_alt_f+=[(mat[i][-2],mat[i][2])]
    ind_frac_all=CDFandHIST(exon_fraction, results_dir +
               'csv/domains/%s/domain_fraction_per_exon'%CONST_val,exon_fraction_alt_f)
    lon_frac_ex=CDFandHIST(longest_fraction, results_dir +
               'csv/domains/%s/longest_exon_domain_fraction'%CONST_val,longest_fraction_alt_f)
    has_back={(0,0.25):0,(0.25,0.50):0,(0.50,0.75):0,(0.75,1.01):0}
    for i in alt_exons_wef_background:
        for keys in has_back:
            if keys[0]<=i<keys[1]:
                has_back[keys]+=1
                break
    alt_back=''
    with open (results_dir +
               'csv/domains/%s/Alt_exons_background_WEF_w_wo_domain'%CONST_val,'w') as fin:
        alt_back+='WEF_range\tTotal\tFreq\n'
        fin.write('WEF_range,Total,Freq\n')
        keys=has_back.keys()
        total=sum(has_back.values())
        keys.sort()
        for i in keys:
            fin.write('%s,%s,%s\n'%("_".join(map(str,list(i))),has_back[i],                                   div_fact(has_back[i],total)))
            alt_back+='%s,%s,%s\n'%("_".join(map(str,list(i))),has_back[i],                                   div_fact(has_back[i],total))
        
    np_arr1 = np.histogram(exon_count_per_domain, bins=[
                           1, 2, 3, 4, 5, 6, 10, 15, 20, 101])
        
    
    exon_count_dom=np_to_csv(np_arr1, results_dir +
               'csv/domains/%s/exons_numbers_codingforDom.csv'%CONST_val,
              cols=["NumberOfExons", "Total", "Freq"])
    fout.write("//\nINDIVIDUAL EXONS PERSPECTIVE\n exons count domain\n%s\n%s\n%s\n Alternate exons background\n%s"              %(exon_count_dom,ind_frac_all,lon_frac_ex,alt_back))
    #fout.close()
exons_view_ind(human,fin)
'''
information needed from this block, 
exon count value, once tab separted
exon fraction all tab sperated
exon fraction longest tab sperated
alt exons background trend
alt exon trend in domains participation
domains that are freqyuerntly contributed by aternate exons havingm ore tan 40%
'''


# In[21]:



def paired_exons_neutralizer(gene_ob,isoform_exons):
    const_exons = gene_ob.connst_togetherness_coding()
    const_exons_ref = []
    for i in const_exons:
        templis = []
        for j in i:
            templis += [j]
        const_exons_ref += [templis]
    const_exons = const_exons_ref[:]
    list_exons = isoform_exons
    # print list_exons
    indexons = {i: [i.length, i.ID,i.WEF, i.cons_flag] for i in list_exons}
    # has_ref = {i[0]: [i, sum([indexons[j] for j in i])] for i in const_exons}
    # print indexons,"indexons"
    has_ref = {}
    for i in const_exons:
        if i:
            has_ref[i[0]] = [i[1:]]
            lencons = 0
            consid = ''
            for j in i:
                lencons += indexons[j][0]
                consid += ":%s" % indexons[j][1]
            has_ref[i[0]] += [lencons, consid]
    # print has_ref
    merged_list = []
    for i in list_exons:
        if i not in has_ref:
            merged_list += [indexons[i]]
        else:
            tempMerge = [has_ref[i][1], has_ref[i][2],1, "G"]
            merged_list += [tempMerge]
            for j in has_ref[i][0]:
                list_exons.remove(j)
    '''
    each exon entry in merged_list has follwing information
    length,ID,cons_flag
    '''
    return merged_list


# In[22]:


'''exons_side_story'''
def exons_view_pair(human,fout):
    exon_count_per_domain=[]
    longest_fraction={'A':[],'G':[],'F':[]}
    exon_fraction={'A':[],'G':[],'F':[]}
    longest_fraction_alt_f=[]
    exon_fraction_alt_f=[]
    dom_ex_frac=[]
    dom_ex_frac_const=[]
    for gene in human:
        if gene in filterg:
            relev_ex,domlis=isoform_giver(human[gene])
            if domlis:
                lis_needed = paired_exons_neutralizer(human[gene],relev_ex)
                for dom in domlis:
                    mat=[]
                    dspan=set(range(dom[0][0],dom[0][1]))
                    su=0
                    for ex in lis_needed:
                        espan=set(range(su,su+ex[0]))
                        su+=ex[0]
                        inter = espan & dspan
                        if inter:
                            key=ex[3]
                            ex_frac=div_fact(len(inter),len(dspan))
                            mat+=[(len(inter),ex[0],ex[2],ex_frac,key)]
                            dfrac = div_fact(len(inter),len(dspan))
                            efrac = div_fact(len(inter),len(espan))
                            dom_ex_frac += [(dfrac, efrac)]
                            if ex[-1]=='G':
                                 dom_ex_frac_const += [(dfrac, efrac)]
                    exon_count_per_domain+=[len(mat)]
                    mat.sort()
                    if mat:
                      if mat[-1][-1] == 'A':
                         longest_fraction_alt_f+=[(mat[-1][-2],mat[-1][2])]
                      longest_fraction[mat[-1][-1]]+=[mat[-1][-2]]
                      for i in range(0,len(mat)):
                        if mat[i][-1]=='A':
                            exon_fraction_alt_f+=[(mat[i][-2],mat[i][2])]
                        exon_fraction[mat[i][-1]]+=[mat[i][-2]]
    quadrant_creator(dom_ex_frac, results_dir+'csv/domains/%s/quadrant_PAIRED.csv'%CONST_val)
    quadrant_creator(dom_ex_frac_const, results_dir+'csv/domains/%s/quadrant_PAIRED_const.csv'%CONST_val)
    ind_frac_all=CDFandHIST(exon_fraction,results_dir + 'csv/domains/%s/PAIRED_domain_fraction_per_exon'%CONST_val,exon_fraction_alt_f)
    lon_frac_ex=CDFandHIST(longest_fraction,results_dir + 'csv/domains/%s/PAIRED_longest_exon_domain_fraction'%CONST_val,longest_fraction_alt_f)
    np_arr1 = np.histogram(exon_count_per_domain, bins=[
                           1, 2, 3, 4, 5, 6, 10, 15, 20, 101])
    exon_count_dom=np_to_csv(np_arr1, results_dir +
               'csv/domains/%s/PAIRED_exons_numbers_codingforDom.csv'%CONST_val,
              cols=["NumberOfExons", "Total", "Freq"])
    fout.write("//\nPAIRED EXONS PERSPECTIVE\n exons count domain\n%s\n%s\n%s\n\n FINISH"              %(exon_count_dom,ind_frac_all,lon_frac_ex))
    fout.close()
exons_view_pair(human,fin)


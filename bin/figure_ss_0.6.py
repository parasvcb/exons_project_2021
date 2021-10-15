
import cPickle as pickle

import sys
if len(sys.argv)!=3:
	print "Please type correct location of results dir having object stored as well as organism id as two different arguements"
	sys.exit()
results_dir_csv = sys.argv[1]
organism_id = sys.argv[2]

source_gene_object = sys.argv[1]+'objectsave_%s_0.6.pick'%organism_id
#results_dir_csv = "/home/paras/project/protein_splicing/9606/results/"
fileswriter=open(results_dir_csv+"csv/ss0.6/ss_0.6.log","w")
window=4

with open(source_gene_object) as fin:
    has_gene = pickle.load(fin)
    print source_gene_object
with open (results_dir_csv+"new_condition_genes.pick") as fin:
    condition = pickle.load(fin)

def div_fact(num,denom):
	try:
		return round(float(num)/denom,3)	
	except:
		return 0

def ss_def(has, ex1, ex2, window):
    n = ex2[:window-1]
    c = ex1[-(window-1):]
    if len(set(n)) == 1:
        nr = list(set(n[:3]))[0]
    else:
        nr = "X"

    if len(set(c)) == 1:
        cr = list(set(c))[0]
    else:
        cr = "X"
    res = cr+nr
    if res not in has:
        has[res] = 1
    else:
        has[res] += 1
    return has

def csv_writer(has, filename,fout):
    key = has.keys()
    total = sum(has.values())
    key.sort()
    fout.write("\n%s"%filename)
    with open(filename, "w") as fin:
        fin.write("SSType,Total,Freq\n")
        fout.write("\nSSType\tTotal\tFreq")
        for i in key:
            fin.write("%s,%s,%s\n" % (i, has[i], div_fact(has[i],total)))
            fout.write("\n%s\t%s\t%s" % (i, has[i], div_fact(has[i],total)))

def sstype_ends_indi(gene_source,condition, filename, window,fout):
    cons_cons = {}
    alt_alt = {}
    cons_alt = {}
    all_junct = {}
    backgroud = {}
    irrelavant=0
    for gene in gene_source:
        flag_p = 0
        if gene in condition:
            tuples=[]
            PI=gene_source[gene].PI
            pi_seq=''.join([i.seq for i in PI.exons])
            ss_seq=''
            for i in PI.exons:
                if i.length>0:
                    ssres=i.out_secondseq(PI.ID)
                    if ssres:
                        ss_seq+=ssres
                    else:
                        pass
                        #print gene,PI.ID

            #ss_seq=''.join([i.out_secondseq(PI.ID) for i in PI.exons])
            if len(pi_seq)==len(ss_seq):
                for i in range (0 ,len(ss_seq)-5):
                    tuples+=[(ss_seq[i:i+6])]
                #print ss_seq
                #print tuples
                #return
            PI_representative=[]
            for i in range(0,len(PI.exons)-1):
                if PI.exons[i].length>9 and PI.exons[i+1].length>9:
                    PI_representative+=[(PI.exons[i],PI.exons[i+1])]
                else:
                    irrelavant+=1
            if PI_representative:
                for i in PI_representative:
                    i1seq=i[0].out_secondseq(PI.ID)
                    i2seq=i[1].out_secondseq(PI.ID)
                    if i[0].ID[0]!='R' and i[1].ID[0]!='R' and i1seq != False and i2seq != False:
                        c_af1=i[0].ID.split(".")[2]
                        c_af2=i[1].ID.split(".")[2]
                        if c_af1 + c_af2  in ['GG','FF','FG','GF']:
                            cons_cons = ss_def(cons_cons,  i1seq, i2seq, window)
                        elif c_af1 + c_af2  in ['AA']:
                            alt_alt = ss_def(alt_alt, i1seq, i2seq, window )
                        else:
                            cons_alt = ss_def(cons_alt, i1seq, i2seq, window )
                        all_junct = ss_def (all_junct, i1seq, i2seq, window)
            if tuples:
                for i in tuples:
                    backgroud = ss_def (backgroud, i[:3],i[3:],window)
    print irrelavant
    csv_writer(backgroud, filename+"background_window%s.csv" % window,fout)
    csv_writer(all_junct, filename+"alljunct_window%s.csv" % window,fout)
    csv_writer(cons_cons, filename+"ConsCons_window_%s.csv" % window,fout)
    csv_writer(alt_alt, filename+"AltAlt_window%s.csv" % window,fout)
    csv_writer(cons_alt, filename+"ConsAlt_window%s.csv" % window,fout)


sstype_ends_indi(has_gene, condition, results_dir_csv +'csv/ss0.6/_Individual_', window,fileswriter)


has_temp={'H':0,'E':0,'C':0}
for gene in has_gene:
    ss_lis=''
    for exons in has_gene[gene].PI.exons:
        ss_lis+=exons.out_secondseq(has_gene[gene].PI.ID) if exons.out_secondseq(has_gene[gene].PI.ID) else ''
    
    for res in ss_lis:
        has_temp[res]+=1
fileswriter.write("\n\nResidue\tTotal\tFreq\n")
keys =has_temp.keys()
tot=sum(has_temp.values())
for key in keys:
    fileswriter.write("%s\t%s\t%s\n"%(key,has_temp[key],div_fact(has_temp[key],tot)))
#fileswriter.close()


# In[17]:


def ss_def_customised( ex1, ex2, window):
    n = ex2[:window-1]
    c = ex1[-(window-1):]
    if len(set(n)) == 1:
        nr = list(set(n))[0]
    else:
        nr = "X"

    if len(set(c)) == 1:
        cr = list(set(c))[0]
    else:
        cr = "X"
    res = cr+nr
    return res
def key_searcher(box,key):
    box.sort()
    if key in box:
        return key
    else:
        small=min(key)
        big=max(key)
        for i in box:
            if i[0]==small or i[1] == big:
                return i
    
        
def conservation_junc_populator(repr_has, dup, trans,PI, window):
    texo=trans.exons[:]
    for i in range (0, len(texo)-1):
        if texo[i].length>9 and texo[i+1].length>9:
            i1seq=texo[i].out_secondseq(trans.ID)
            i2seq=texo[i+1].out_secondseq(trans.ID)
            if texo[i].ID[0]!='R' and texo[i+1].ID[0]!='R' and i1seq != False and i2seq != False:
                part1=int(texo[i].ID.split(".")[3])
                part2=int(texo[i+1].ID.split(".")[3])
                desired_key=key_searcher(dup.keys(),(part1,part2))
                key_refer=(texo[i].seq[-window:],texo[i+1].seq[:window],i1seq[-window:],i2seq[:window])
                if key_refer not in dup[desired_key]:
                    #print repr_has,1
                    #print dup,1
                    dup[desired_key]+=[key_refer]
                    #print repr_has,2
                    #print dup,2
                    #print (texo[i].ID,texo[i+1].ID,i1seq,i2seq,PI)
                    repr_has[desired_key]+=[(texo[i].ID,texo[i+1].ID,i1seq,i2seq,PI)]
                    
                #4 ele tuple of first two aaseq and last two sseq
                
    #print repr_has
    #print dup
    return repr_has,dup
    
def sstype_ends_coservation(gene_source,condition, filename, window, fout):
    total_junctions={}
    '''
    has_mega={'HH':{(0,0.25):0,(0.25,0.50):0,(0.50,0.75):0,(0.75,1.00):0,(1.00,1.01):0},
             'CC':{(0,0.25):0,(0.25,0.50):0,(0.50,0.75):0,(0.75,1.00):0,(1.00,1.01):0},
             'EE':{(0,0.25):0,(0.25,0.50):0,(0.50,0.75):0,(0.75,1.00):0,(1.00,1.01):0}}
    '''
    has_mega={'HH':{(0,0.1):0,(0.1,0.2):0,(0.2,0.3):0,(0.3,0.4):0,(0.4,0.5):0,(0.5,0.6):0,                    (0.6,0.7):0,(0.8,0.9):0,(0.9,1.0):0,(1.00,1.01):0},
             'CC':{(0,0.1):0,(0.1,0.2):0,(0.2,0.3):0,(0.3,0.4):0,(0.4,0.5):0,(0.5,0.6):0,\
                    (0.6,0.7):0,(0.8,0.9):0,(0.9,1.0):0,(1.00,1.01):0},
             'EE':{(0,0.1):0,(0.1,0.2):0,(0.2,0.3):0,(0.3,0.4):0,(0.4,0.5):0,(0.5,0.6):0,\
                    (0.6,0.7):0,(0.8,0.9):0,(0.9,1.0):0,(1.00,1.01):0}}
    fout.write("Helix junction freq\tGene\tExon_pair\tGene\n")
    for gene in gene_source:
        flag_p = 0
        if gene in condition:
            aaseq={}
            len_exons=len(gene_source[gene].exons)
            representative={(i,i+1):[] for i in range(0,len_exons+2)}
            representative_duplicate={(i,i+1):[] for i in range(0,len_exons+2)}
            
            for trans in gene_source[gene].transcripts:
                seq=''.join([i.seq for i in trans.exons])
                if seq not in aaseq:
                    aaseq[seq]=0
                    representative, representative_duplicate =                     conservation_junc_populator(representative,representative_duplicate, trans, trans.PI,window)
            for i in representative:
                temp_lis=[]
                if len(representative[i])>2:
                    #atleast two combinations:
                    for j in representative[i]:
                        #print j
                        seq1=j[2]
                        seq2=j[3]
                        temp_lis+=[(j[4], ss_def_customised(seq1,seq2,window))]
                if temp_lis:
                    temp_lis.sort()
                    #print temp_lis
                    parent_junction=temp_lis[-1][1]
                    if parent_junction not in total_junctions:
                        total_junctions[parent_junction]=1
                    else:
                        total_junctions[parent_junction]+=1
                    if parent_junction in ['HH','EE','CC']:
                        freq=div_fact([j[1] for j in temp_lis].count(parent_junction),len(temp_lis))
                        exon_pair_lis=[(i,j) for j in representative[i]]
                        fout.write("%s\t%s\t%s\t%s\n"%(freq,gene,exon_pair_lis,gene_source[gene].detail))
                        for vals in has_mega[parent_junction]:
                            if vals[0]<=freq<vals[1]:
                                    has_mega[parent_junction][vals]+=1
                                    break
    fout.close()
    print has_mega
    print total_junctions
    with open (filename+'exons_junctions.csv','w') as fin:
        fin.write("SSType,Total,Freq\n")
        keys=total_junctions.keys()
        keys.sort()
        total_semi=sum(total_junctions.values())
        for i in keys:
            fin.write("%s,%s,%s\n"%(i,total_junctions[i],div_fact(total_junctions[i],total_semi)))
    major_ss=has_mega.keys()
    major_ss.sort()
    for i in major_ss:
        with open (filename+'%s_conservation.csv'%i,'w') as fin:
            fin.write("Occurence,Total,Freq\n")
            mini_keys=has_mega[i].keys()
            mini_keys.sort()
            mini_sum=sum(has_mega[i].values())
            for j in mini_keys:
                fin.write("%s,%s,%s\n"%(j,has_mega[i][j],div_fact(has_mega[i][j],mini_sum)))
        
sstype_ends_coservation(has_gene, condition, results_dir_csv +'csv/ss0.6/_varying_partner_', window,fileswriter)


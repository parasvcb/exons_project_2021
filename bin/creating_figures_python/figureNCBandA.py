
import sys, os
import constructing_data.CDF as CDF
if len(sys.argv)!=2:
	print "Please type correct location of results dir having object stored"
	sys.exit()
results_dir_csv = sys.argv[1]
import re 
import numpy as np

with open (os.path.join(results_dir_csv,"matrix_NCBI.csv")) as fin:
    dat=[i for i in fin.read().split("\n") if len(i)>0 and ''.join(i[:5])!='Event']

def identity(a,b):
    tot=len(a)
    eq=0;uneq=0
    for i in range (0,len(a)):
        if a[i]==b[i]:
            eq+=1
        else:
            uneq+=1
    return round((float(eq)/tot)*100,2)

def csv_writer(lis, filename_prefix):
    # print "inside"
    template = [0, '', '', '']
    conserved_fraction = []
    posden=[]
    # first element will be total number of cases falling in this bin, secondis infininte long ss,
    # then stride and then the position of protein letters n,m,c
    has = {('-', -50): template[:], (-50, -40): template[:], (-40, -30): template[:],(-30, -20): template[:],
           (-20, -10): template[:], (-10, -5): template[:],(-5,-2):template[:], (3,5):template[:],
           (5,10):template[:],(10, 20): template[:], (20, 30): template[:], (30, 40): template[:],
           (40, 50): template[:], (50, '+'): template[:]}
    for i in lis:
        position_protein = 'N' if 0 <= i[0] <= 0.30 else 'M' if i[0] <= 0.70 else 'C'
        posden += [i[0]]
        #print i[0], position_protein
        #print i,"position"
        conserved_fraction += [identity(i[4], i[5])]
        '''
        i[1] is the length of change, i[0] is the protein position, i[2] is pssm unique, i[3] is stride 
        unique, i[4] is parentaaseq, i[5] is child aaseq
        '''
        for key in has:
            if key[0] == '-':
                if i[1] < key[1]:
                    keyres = key
                    break
            if key[1] == '+':
                #print "true"
                #print i[1], key, "i[1] > key[0]",type(i[1]), type(key[0]), i[1] > key[0], -15 > 250
                if i[1] > key[0]:
                    keyres = key
                    break
            
            if key[0] <= i[1] < key[1]:
                    keyres = key
                    break
        #print i, i[1]
        has[keyres][0] += 1
        if i[2]:
            has[keyres][1] += i[2]
        if i[3]:
            has[keyres][2] += i[3]
        has[keyres][3] += position_protein
    # print has
    with open(filename_prefix+'change_histogram.csv', 'w') as fin:
        fin.write(
            "ChangeInLength,TotalCases,Frequency,SSeqLength,Css,Hss,Ess,StrideLength,Cst,Hst,Est,Nter,Middle,Cter\n")
        totalcases = sum([i[0] for i in has.values()])
        #print lis
        #print totalcases
        keys = has.keys()
        keys.sort()
        first = keys[-1]
        del keys[-1]
        keys.insert(0, first)
        for i in keys:
            semi_total=has[i][0]
            whole_ss_seq=has[i][1]
            whole_st_seq=has[i][2]
            whole_positions=has[i][3]
            # print whole_positions
            CterFreq=round(float(whole_positions.count('C'))/len(whole_positions),2) if whole_positions else 0
            NterFreq=round(float(whole_positions.count('N'))/len(whole_positions),2) if whole_positions else 0
            MiddleFreq=round(float(whole_positions.count('M'))/len(whole_positions),2) if whole_positions else 0
            Css=round(float(whole_ss_seq.count('C'))/len(whole_ss_seq),2) if whole_ss_seq else 0
            Hss=round(float(whole_ss_seq.count('H'))/len(whole_ss_seq),2) if whole_ss_seq else 0
            Ess=round(float(whole_ss_seq.count('E'))/len(whole_ss_seq),2) if whole_ss_seq else 0
            Cst=round(float(whole_st_seq.count('C'))/len(whole_st_seq),2) if whole_st_seq else 0
            Est=round(float(whole_st_seq.count('E'))/len(whole_st_seq),2) if whole_st_seq else 0
            Hst=round(float(whole_st_seq.count('H'))/len(whole_st_seq),2) if whole_st_seq else 0
            fin.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%("_".join(map(str,i)),semi_total,            round(float(semi_total)/totalcases,2),len(whole_ss_seq),Css,Hss,Ess,len(whole_st_seq),Cst,Est,Hst,                                                                    NterFreq,MiddleFreq,CterFreq))
    with open (filename_prefix+'densityvals.csv','w') as fin:
        fin.write("position\n")
        for i in posden:
            fin.write("%s\n"%i)
    with open (filename_prefix+'densityvals_CDF.csv','w') as fin:
        data=CDF.CDF_String(posden,["Position_affected","Total","Freq"])
        fin.write("%s"%data)
    with open (filename_prefix+'conserved_fraction.csv','w') as fin:
        has_conserv={(0,25):0,(25,50):0,(50,75):0,(75,100):0,(100,101):0}
        for i in conserved_fraction:
            for j in has_conserv:
                if j[0]<=i<j[1]:
                    has_conserv[j]+=1
                    break
        fin.write("Conserved_sequence_range,Total,Frequency\n")
        keys_cons=has_conserv.keys()
        keys_cons.sort()
        total=sum(has_conserv.values())
        #print has_conserv
        for i in keys_cons:
            #print "_".join(map(str,i))
            #print has_conserv[i]
            #print round(float(has_conserv[i])/total,3)
            fin.write("%s,%s,%s\n"%("_".join(map(str,i)),has_conserv[i],round(float(has_conserv[i])/total,3)))

nlis=[];clis=[];bnlis=[];bclis=[];blis=[]
b_corr_x=[]
b_corr_y=[]
alis=[]
flis=[]
gene_type={'N':[],'C':[],'B':[],'A0':[]}
for line in dat:
    ele=line.split("\t")
    head= ele[0]
    position = ele[3]
    length=ele[4]
    if'na' not in length and 'False' not in length and not(set(map(int,length.split(","))) & set([0,1,2,-2,-1])):
        # print "->",length,map(int,length.split(","))[0]
        if re.match(r'^[-]?[\d+]\.[\d+]',position) :
            ssseq=ele[10] if ele[10]!='na' else False
            strideseq=ele[13] if ele[13]!='na' else False
            parentseq=ele[5]
            childseq=ele[6]
            #print line
            #print ele[1],ele[0], parentseq, childseq
            #print ele
            if len(parentseq)>0 and len(childseq)>0:
                #print head
                if head == 'B':
                    #print map(float,position.split(",")) 
                    value,valuen,valuec=map(int,length.split(","))
                    b_corr_x+=[valuen]
                    b_corr_y+=[valuec]
                    ssseqn,ssseqc=ele[10].split(",") if 'na' not in ele[10] else [False,False]
                    strideseqn,strideseqc=ele[13].split(",") if 'na' not in ele[13] and 'Z' not in ele[13] else [False,False]
                    parentseq=ele[5]
                    childseq=ele[6]
                    bnlis+=[(float(position),valuen,ssseqn,strideseqn,parentseq,childseq)]
                    bclis+=[(float(position),valuec,ssseqc,strideseqc,parentseq,childseq)]
                    blis+=[value]
                elif head == 'N':
                    #print "in"
                    nlis+=[(float(position),int(length),ssseq,strideseq,parentseq,childseq)]
                elif head == 'C':
                    clis+=[(float(position),int(length),ssseq,strideseq,parentseq,childseq)]
                elif head == 'A0':
                    #print "in", ele[2]
                    if ele[2].split(".")[2]=='A':
                        alis+=[float(position)]
                    else:
                        flis+=[float(position)]
                gene_type[head]+=[int(ele[1])]
            else:
                print line,ele

csv_writer(nlis,os.path.join(results_dir_csv,"Fig2_1_N_ext"))
csv_writer(clis,os.path.join(results_dir_csv,"Fig2_1_C_ext"))
csv_writer(bclis,os.path.join(results_dir_csv,"Fig2_1_BC_ext"))
csv_writer(bnlis,os.path.join(results_dir_csv,"Fig2_1_BN_ext"))

with open (os.path.join(results_dir_csv,'Fig2_1_AltA_densityvals.csv'),"w") as fin:
    fin.write("position\n")
    for i in alis:
        fin.write("%s\n"%i)
with open (os.path.join(results_dir_csv,'Fig2_1_AltA_densityvals_CDF.csv'),'w') as fin:
    data=CDF.CDF_String(alis,["Position_affected","Total","Freq"])
    fin.write("%s"%data)

with open (os.path.join(results_dir_csv,'Fig2_1_histogram_ncb.csv'),"w") as fin:
    fin.write("AS_type,Total,Freq\n")
    total=len(nlis)+len(clis)+len(bclis)
    fin.write("%s,%s,%s\n"%('N',len(nlis),round(float(len(nlis))/total,3)))
    fin.write("%s,%s,%s\n"%('C',len(clis),round(float(len(clis))/total,3)))
    fin.write("%s,%s,%s\n"%('B',len(bnlis),round(float(len(bnlis))/total,3)))
with open (os.path.join(results_dir_csv,'Fig2_1_change_in_length_all.csv'),'w') as fin:
    fin.write('Length,Type\n')
    for i in nlis:
        fin.write('%s,N\n'%(i[1]))
    for i in clis:
        fin.write('%s,C\n'%(i[1]))
    for i in blis:
        fin.write('%s,B\n'%(i))
    for i in bnlis:
        fin.write('%s,BN\n'%(i[1]))
    for i in bclis:
        fin.write('%s,BC\n'%(i[1]))


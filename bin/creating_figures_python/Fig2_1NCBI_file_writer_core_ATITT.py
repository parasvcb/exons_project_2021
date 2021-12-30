
import cPickle as pickle
import sys,os
import common.general_modules as gm
import figPanels.modules_analysis as am

import re
import numpy as np

if len(sys.argv)!=4:
	print ("Please type 1. location of object, 3. results dir having condition pickle 3, ATITrue=0, ATIFalse=1, reagrdless=-1")
	sys.exit()
prog, source_gene_object, results_dir_csv, ATIFlag = sys.argv

ATIFlag=int(ATIFlag)
fnameappend='both' if ATIFlag ==-1 else 'core' if ATIFlag==1 else 'atit' if ATIFlag==0 else 'False'

humanGeneObject=gm.readPickle(source_gene_object)
CONDITION_GENES=gm.readPickle(os.path.join(results_dir_csv,"condition_genes.pick"))

def calculatorB(parent_cood, child_cood):
    list_temp = [(child_cood, 0), (parent_cood, 1)]
    list_temp.sort()
    tfirst = list_temp[0][0][0]
    tlast = list_temp[0][0][1]
    ele = range(tfirst, (tfirst+((tlast-tfirst)+1)/3))
    prev_ele = ele[0]
    new_lis = [(list_temp[0][1], ele[0], ele[-1])]
    tfirst = list_temp[1][0][0]
    tlast = list_temp[1][0][1]
    cal_fac = int(round(float(tfirst-prev_ele)/3, 0))
    ele = range(prev_ele+cal_fac, (prev_ele+cal_fac+(tlast-tfirst)/3)+1)
    new_lis += [(list_temp[1][1], ele[0], ele[-1])]
    new_lis.sort()
    nside = new_lis[0][1]-new_lis[1][1]
    cside = new_lis[0][2]-new_lis[1][2]
    return nside, cside


def condition_checker(parent, child, specifier):
    pleft = parent.coding_span[0]
    pright = parent.coding_span[1]
    cleft = child.coding_span[0]
    cright = child.coding_span[1]
    if set(range(pleft, pright+1)) & set(range(cleft, cright+1)):
        leftdiff = abs(pleft-cleft)
        rightdiff = abs(pright-cright)
        if leftdiff == 0 and rightdiff == 0:
            if parent.seq == child.seq:
                return [False, "ATIT:Same_amino_acid"]
            else:
                return [False, "ATIT:Diff_amino_acid"]
        if specifier == "N":
            if leftdiff != 0 and rightdiff in {0, 1, 2}:
                # means left is varying and right is in consistent +-1,2 nt
                # proper N case
                return [True]
            else:
                return [False, "Variation_in_coding_coods"]
        if specifier == "C":
            if leftdiff in {0, 1, 2} and rightdiff != 0:
                return [True]
            else:
                return [False, "Variation_in_coding_coods"]
        if specifier == "B":
            if leftdiff != 0 and rightdiff != 0:
                return [True]
            elif leftdiff == 0:
                return [False, "pseudoC"]
            elif rightdiff == 0:
                return [False, "pseudoN"]
            else:
                print ("surprise!!!!!")
    else:
        return [False, "No_overlap_exists"]



def changes_specifier(N=None, C=None, parent=None, child=None):
    #print N,C
    if N != None and C != None:
        if N > 0:
            uniqueN = parent[:abs(N)]
            if C > 0:
                uniqueC = child[-1*abs(C):]
                overlapparent = parent[abs(N):]
                overlapchild = child[:-1*abs(C)]
            elif C < 0:
                uniqueC = parent[-1*abs(C):]
                overlapchild = child[:]
                overlapparent = parent[abs(N):-1*abs(C)]
            else:
                uniqueC = ''
                overlapchild = child[:]
                overlapparent = parent[abs(N):]
        elif N < 0:
            uniqueN = child[:abs(N)]
            if C > 0:
                uniqueC = child[-1*abs(C):]
                overlapparent = parent[:]
                overlapchild = child[abs(N):-1*abs(C)]
            elif C < 0:
                uniqueC = parent[-1*abs(C):]
                overlapparent = parent[:-1*abs(C)]
                overlapchild = child[abs(N):]
            else:
                uniqueC = ''
                overlapparent = parent[:]
                overlapchild = child[abs(N):]
        else:
            # N is 0
            uniqueN = ''
            if C == 0:
                uniqueC = ''
                overlapparent = parent[:]
                overlapchild = child[:]
            elif C > 0:
                uniqueC = child[-1*abs(C)]
                overlapparent = parent[:]
                overlapchild = child[:-1*abs(C)]
            else:
                uniqueC = parent[-1*abs(C)]
                overlapparent = parent[:-1*abs(C)]
                overlapchild = child[:]
        return uniqueN, uniqueC, overlapparent, overlapchild

    elif N != None and C == None:
        #print "in",N

        if N < 0:
            unique = parent[:abs(N)]
            overlapchild = child[:]
            overlapparent = parent[abs(N):]
        elif N > 0:
            #print "n",N
            unique = child[:abs(N)]
            overlapchild = child[abs(N):]
            overlapparent = parent[:]
        else:
            unique = ''
            overlapchild = child[:]
            overlapparent = parent[:]
        '''
        print parent
        print child
        print unique
        print overlapparent
        print overlapchild
        '''
    elif C != None and N == None:
        # print "why ????"
        if C < 0:
            unique = parent[-1*abs(C):]
            overlapchild = child[:]
            overlapparent = parent[:-1*abs(C)]
        elif C > 0:
            # print "everin?"
            unique = child[-1*abs(C):]
            overlapchild = child[:-1*abs(C)]
            overlapparent = parent[:]
        else:
            unique = ''
            overlapchild = child[:]
            overlapparent = parent[:]
    else:
        print ("deadly _ conduition")
        print (N, C)

    return unique, overlapparent, overlapchild



'''
This function will be doing the following 
Going to the exons in gene, checking if it is ncb or alternate exon?
once in 
search out its parent exon, if that is in the exon list,
once in2
hunto for parent and child object
get their secondary and tertiary structure stride sequenec

for case to be N
is there any intersection:
    if diff of right and left both are strict 0, and aaseq is same, then right alt initiaition site 
    an dsame or diff aaseq

    if diff of left cood varies(hasto) and diff of rigt cood also varies but in range of abs(0to2),
    then its acceptable
    right it to be pseudoB case
if not its not psudo extension case

same goes for C

for B
there should be overlap
cases can be strcit N or C and B figure them out


'''
def func1_fillerNCB(gene, has_all, PI):
    ncase = []
    ccase = []
    bcase = []
    icase = []
    for i in has_all:
        if i[0]!='0':
            #print gene, i
            parent_exon = i[2].parent
            if parent_exon:
                temp_case=i[0]
                pob = parent_exon
                cob = i[2]
                childlength = cob.length
                parentlength = pob.length
                parentaaseq = pob.seq
                childaaseq = cob.seq
                sseqp = pob.out_secondseq(PI) if pob in PI_exons else pob.out_secondseq()
                sseqc = cob.out_secondseq(has_all[i])
                #sseqc = [it for it in cob.ssseq.keys(
                #) if it is not None and it != "NULL"]
                strideseqp = pob.out_strideseqAF(PI) if pob in PI_exons else pob.out_strideseqAF()
                
                #[it for it in pob.strideseq if it !=
                #             "NULL" and it is not None and PI in pob.strideseq[it]]
                strideseqc = cob.out_strideseqAF(has_all[i])

                parentssseq = sseqp if sseqp else False
                childssseq = sseqc if sseqc else False
                #print parentssseq
                #print childssseq
                #print pob.ID, cob.ID
                #print strideseqp, strideseqc
                parentstrideseq = strideseqp if strideseqp else False
                childstrideseq = strideseqc if strideseqc else False
                len_par_aa = childlength-parentlength
                len_par_ss = len_par_aa if parentssseq and childssseq else None
                len_par_stride = len_par_aa if parentstrideseq and childstrideseq else None

                exid_pair = "%s:%s" % (cob.ID, pob.ID)
                # print exid_pair
                # print sseqp,sseqc,len_par_aa, len_par_ss, len_par_stride
                '''
                Defaults
                '''
                #rint i
                position_affected = i[3]
                overlap_aap = 'na'
                overlap_aac = 'na'
                unique_aa = 'na'
                unique_aaN='na'
                unique_aaC='na'
                overlap_ssp = 'na'
                overlap_ssc = 'na'
                overlap_stridep = 'na'
                overlap_stridec = 'na'
                unique_ss = 'na'
                unique_stride = 'na'
                unique_ssN = 'na'
                unique_ssC = 'na'
                unique_strideN = 'na'
                unique_strideC = 'na'
                nside='na'
                cside='na'
                if temp_case == "n":
                    flag_n = condition_checker(
                        parent=pob, child=cob, specifier='N')
                    # print flag_n
                    if flag_n[0]:
                        #print "isit"
                        # position_affected = i[3]
                        unique_aa,  overlap_aap, overlap_aac = changes_specifier(
                            N=len_par_aa, child=childaaseq, parent=parentaaseq)
                        # print unique_aa,  overlap_aap, overlap_aac
                        if len_par_ss is not None:
                            unique_ss, overlap_ssp, overlap_ssc = changes_specifier(
                                N=len_par_ss, child=childssseq, parent=parentssseq)
                            # print unique_ss, overlap_ssp, overlap_ssc
                        if len_par_stride is not None:
                            unique_stride, overlap_stridep, overlap_stridec = changes_specifier(
                                N=len_par_stride, child=childstrideseq, parent=parentstrideseq)

                    else:
                        #print "isit??"
                        position_affected = flag_n[1]

                    ncase += [["N", gene, exid_pair, position_affected, len_par_aa, overlap_aap, overlap_aac,
                               unique_aa, overlap_ssp, overlap_ssc, unique_ss, overlap_stridep, overlap_stridec, unique_stride, has_all[i]]]
                elif temp_case == "c":
                    flag_c = condition_checker(
                        parent=pob, child=cob, specifier='C')
                    if flag_c[0]:
                        # position_affected = 
                        unique_aa,  overlap_aap, overlap_aac = changes_specifier(
                            C=len_par_aa, child=childaaseq, parent=parentaaseq)
                        '''
                        print len_par_aa, "lenpaa"
                        print childaaseq, "caaseq"
                        print parentaaseq,"paaseq"
                        print  unique_aa, "uniqueaa"
                        print overlap_aap, "overlapaa",
                        print overlap_aac, "overlap_aac"
                        '''
                        if len_par_ss is not None:
                            unique_ss, overlap_ssp, overlap_ssc = changes_specifier(
                                C=len_par_ss, child=childssseq, parent=parentssseq)

                        if len_par_stride is not None:
                            unique_stride, overlap_stridep, overlap_stridec = changes_specifier(
                                C=len_par_stride, child=childstrideseq, parent=parentstrideseq)

                    else:
                        position_affected = flag_c[1]
                    ccase += [["C", gene, exid_pair, position_affected, len_par_aa, overlap_aap, overlap_aac,
                               unique_aa, overlap_ssp, overlap_ssc, unique_ss, overlap_stridep, overlap_stridec, unique_stride, has_all[i]]]

                elif temp_case == "b":
                    flag_b = condition_checker(
                        parent=pob, child=cob, specifier='B')
                    if flag_b[0]:

                        #position_affected = midH[pob.ID] if pob.ID in midH else 'na'
                        nside, cside = calculatorB(
                            pob.coding_span, cob.coding_span)
                        unique_aaN, unique_aaC, overlap_aap, overlap_aac = changes_specifier(C=cside, N=nside, child=childaaseq, parent=parentaaseq)
                        if len_par_ss is not None:
                            unique_ssN, unique_ssC, overlap_ssp, overlap_ssc = changes_specifier(
                                C=cside, N=nside, child=childssseq, parent=parentssseq)

                        if len_par_stride is not None:
                            unique_strideN, unique_strideC, overlap_stridep, overlap_stridec = changes_specifier(
                                C=cside, N=nside, child=childstrideseq, parent=parentstrideseq)

                    else:
                        position_affected = flag_b[1]
                    bcase += [["B", gene, exid_pair, position_affected, str(len_par_aa)+","+str(nside)+","+str(cside),
                               overlap_aap, overlap_aac, str(
                                   unique_aaN)+","+str(unique_aaC),
                               overlap_ssp, overlap_ssc, str(
                                   unique_ssN)+","+str(unique_ssC),
                               overlap_stridep, overlap_stridec, str(unique_strideN)+","+str(unique_strideC), has_all[i]]]

                    '''
                    parent_cood = pob.coding_span
                    child_cood = cob.coding_span
                    
                    overlap_aap1, overlap_aac1, unique_aa1, overlap_ssp1, overlap_ssc1, unique_ss1, overlap_stridep1, overlap_stridec1, unique_stride1 =\
                        changes_in_Nside(nside, parentaaseq, childaaseq, parentssseq,
                                         childssseq, parentstrideseq, childstrideseq)
                    unique_ss, overlap_ssp, overlap_sscoverlap_aap2, overlap_aac2, unique_aa2, overlap_ssp2, overlap_ssc2, unique_ss2, overlap_stridep2, overlap_stridec2, unique_stride2 =\
                        changes_in_Cside(cside, parentaaseq, childaaseq, parentssseq,
                                         childssseq, parentstrideseq, childstrideseq)
                    bcase += [["B", gene, exid_pair, position_affected, str(nside)+","+str(cside), str(overlap_aap1)+","+str(overlap_aap2),
                               str(overlap_aac1)+","+str(overlap_aac2), str(unique_aa1) +
                               ","+str(unique_aa2), str(overlap_ssp1) +
                               ","+str(overlap_ssp2),
                               str(overlap_ssc1)+","+str(overlap_ssc2), str(unique_ss1)+","+str(unique_ss2), str(
                                   overlap_stridep1)+","+str(overlap_stridep2), str(overlap_stridec1)+","+str(overlap_stridec2),
                               str(unique_stride1)+","+str(unique_stride2)]]
                    '''
                
                    
        else:
            cob = i[2]
            position_affected = i[3]
            icase += [["A0", gene, cob.ID, position_affected, False, False, False, cob.seq, False, False,
                               False, False, False, False, has_all[i]]]

    return ncase, ccase, bcase, icase


'''
This block is doing the following
going to gene, getting its coding eoxns in a hash with id as key and object as value

get the exons of principal isoform, and its length and record the start and end position of every exon, about 
how much fraction they are covering for the exons

1. number of total exon in princiapl soform, median mean and average, and number of alt exons in other sioforms
with same statistics

2. add the exons from other transcripts also, about their length start and end with affected protein region

3. total cases where NCB in 1 and 2, are transitioning into coding from noncoding and from noncoding to coding 
with mean and average tats per gene with range also

4. focus the data and location for coding to codng change only and the position of prtein being affected in 
exons of PI and exons of other transcripts.

5. create the core file with the following format
gene, exonsidpair, chageinlength, N,C,B or ALT, parent aaseq, overlap, chikld aaseq, overlap, overhang
and then for stride and pssm

feel ad crete histograms fro the length being changed
how many time site is divisible by 3 

create bins then 
extension happening amonsgt those bins, and frequency
underlying stride changes and ssseq achanges in those same bins
and position of protein being affected in all those bins

seems very comoprehsnive to discuss

HOMEWORK, 
include iupred and have a look at SASA and cases
merge them also


for all cases, see how many times underlying aaseq changes , (obviously in non divisible by 3 extensions)
seq remains same in divisible by 3 cases )???




thereafter the information is passed onto 
filler ncb function with arguements (geneID,starthash,endhash,exonsofGene,transcriptID)
'''

'''
oct 24 2021, 
This code will have follwing parts
A) Get the gene exons which are coding [T tag] and add them into the datex_cood 
    Posrecord: Iterate the transcripts and record which position of protein is being affected by them
    Check for exon qualifier (not retenetion and status should be not 0 but withe ro fn c or b)
    if qualified:
        if ext != '0': #['n','c','b']: 
            if cod_noncood>0:
                #exon should have T flag of not zero
                if i.parent.ID in datex_cood:
                    #factA, this exon does have its parent in gene exons
                    ultrakey=(ext,cons_alt,i,posAffected[ext])
                    #ultrakey is only attributing local properties
                    if ultrakey not in has_all:
                        has_all[ultrakey]=trans.ID
                        cod_both[ext]+=[ultrakey]
                else:
                    #fact!A, parent must be non coding
                    cod_child_noncood_parent[ext]+=[uniq_pair]
            else:
                #exon does have T flag of zero (FactB)
                if i.parent.ID in datex_cood:
                    cod_parent_noncood_child[ext]+=[uniq_pair]
                else:
                    non_cood_both[ext]+=[uniq_pair]
        else:
            #fact C of not interest but should have alternate cases without splice sites
            ultrakey=(ext,cons_alt,i,posAffected[0])
            # print ultrakey
            if ultrakey not in has_all:
                has_all[ultrakey]=trans.ID

    #give them later separate tags, if differece and its magnitude is comparatively larger or not, 
    and get to know so many varied cases

'''

n_case = []
c_case = []
b_case = []
alt_case = []
cod_both={'n':[],'c':[],'b':[]}
cod_parent_noncood_child={'n':[],'c':[],'b':[]}
cod_child_noncood_parent={'n':[],'c':[],'b':[]}
non_cood_both={'n':[],'c':[],'b':[]}
for gene in humanGeneObject:
    if gene in CONDITION_GENES:
        tExonList=am.exonScreenerBetweenTConstitutiveTransWise(humanGeneObject[gene].exons)
        # will be a list of exons of interest
        #print gene
        datex_cood={}
        datex_non_cood={}
        for exons in humanGeneObject[gene].exons:
            if re.match(r'^\w[\:\.][-]?[0-n]+\.[GFA]\.[0-n]\.[0]',exons.ID):
                if re.match(r'^T\.[0-n]',exons.ID):
                    datex_cood[exons.ID]=exons
                else:
                    datex_non_cood[exons.ID]=exons
        #above we got all the non transformed exons
       
        has_all={}
        PI_exons = [i for i in humanGeneObject[gene].PI.exons if re.match(r'\w+[\.\:][-]?[0-n]+', i.ID)]

        trans_length = sum([i.length for i in PI_exons])

        for trans in humanGeneObject[gene].transcripts:
            # print trans.ID, gene, trans.seqlen
            trans_exons = [i for i in trans.exons if re.match(r'\w+[\.\:][-]?[0-n]+', i.ID)]
            trans_length = sum([i.length for i in trans_exons])
            su = 0
            for i in trans_exons:
                # this consecutive block should have been moved below after conmdintion of not retention and other are met, however this surely should have messed with poistion ofthem affecting the transcript
                posAffected={'n':0,'c':0,'b':0}
                tempStart=round(float(su)/trans_length, 2)#start position affected
                su += i.length
                tempMiddle=round((float(su)-float(i.length)/2)/trans_length, 2) #exon length added then middle of exon length affects which region
                tempEnd=round(float(su)/trans_length, 2) #
                posAffected['n']=tempStart
                posAffected['c']=tempEnd
                posAffected['b']=tempMiddle
                
                exide=i.ID.split(".")
                numericFlag=int(i.ID.split(".")[3])
                #ATITrue=0, ATIFalse=1, reagrdless=-1
                conditionATIATT=True if ATIFlag ==-1 else numericFlag in tExonList if ATIFlag==1 else numericFlag not in tExonList if ATIFlag==0 else False
                if i.ID[0]!='R' and conditionATIATT:
                    uniq_pair=(gene,i.ID)
                    cons_alt=exide[2];ext=exide[4];cod_noncood=int(exide[1])
                    #segment A
                    if ext != '0': #['n','c','b']: 
                        #segment AB
                        if cod_noncood>0:
                            if i.parent.ID in datex_cood:
                                ultrakey=(ext,cons_alt,i,posAffected[ext])
                                #ultrakey is only attributing local properties
                                if ultrakey not in has_all:
                                    has_all[ultrakey]=trans.ID
                                    cod_both[ext]+=[ultrakey]
                            else:
                                cod_child_noncood_parent[ext]+=[uniq_pair]
                        else:
                        #segment B
                            if i.parent.ID in datex_cood:
                                cod_parent_noncood_child[ext]+=[uniq_pair]
                            else:
                                non_cood_both[ext]+=[uniq_pair]
                    else:
                        ultrakey=(ext,cons_alt,i,posAffected['n'])
                        # print ultrakey
                        if ultrakey not in has_all:
                            has_all[ultrakey]=trans.ID
                            
                
        
        a, b, c, d = func1_fillerNCB(gene, has_all, humanGeneObject[gene].PI.ID)
        n_case += a
        c_case += b
        b_case += c
        alt_case += d
        #gene ends, anothe rloop will start

def structure_val(gene):
    lis=[]
    cath=False
    for i in gene.transcripts:
        if i.cath_list:
            cath=True
        if i.coverage:
            lis+=[i.coverage]
    return max(lis) if lis else 'Cath' if cath else 'False'
            

cod_bothN={i:len(set(cod_both[i])) for i in cod_both}
cod_parent_noncood_childN={i:len(set(cod_parent_noncood_child[i])) for i in cod_parent_noncood_child}
cod_child_noncood_parentN={i:len(set(cod_child_noncood_parent[i])) for i in cod_child_noncood_parent}
non_cood_bothN={i:len(set(non_cood_both[i])) for i in non_cood_both}
with open(os.path.join(results_dir_csv,"NCB_log_"+fnameappend),"w")as fin:
    keyvals=['relevant','ATIT:Same_amino_acid','ATIT:Diff_amino_acid','Variation_in_coding_coods','pseudoC',            'pseudoN','No_overlap_exists']
    fin.write("Type\tN\tC\tB\n")
    sumtotn=cod_bothN['n']+cod_parent_noncood_childN['n']+cod_child_noncood_parentN['n']+non_cood_bothN['n']
    sumtotc=cod_bothN['c']+cod_parent_noncood_childN['c']+cod_child_noncood_parentN['c']+non_cood_bothN['c']
    sumtotb=cod_bothN['b']+cod_parent_noncood_childN['b']+cod_child_noncood_parentN['b']+non_cood_bothN['b']
    fin.write("Both_coding\t%s,%s\t%s,%s\t%s,%s\n"%(cod_bothN['n'],gm.div_fact(float(cod_bothN['n']),sumtotn),
                                                       cod_bothN['c'],gm.div_fact(cod_bothN['c'],sumtotc),
                                                       cod_bothN['b'],gm.div_fact(cod_bothN['b'],sumtotb)))
    fin.write("CParent_NCchild\t%s,%s\t%s,%s\t%s,%s\n"%(cod_parent_noncood_childN['n'],gm.div_fact(cod_parent_noncood_childN['n'],sumtotn),
                                                       cod_parent_noncood_childN['c'],gm.div_fact(cod_parent_noncood_childN['c'],sumtotc),
                                                       cod_parent_noncood_childN['b'],gm.div_fact(cod_parent_noncood_childN['b'],sumtotb)))
    fin.write("NCParent_Cchild\t%s,%s\t%s,%s\t%s,%s\n"%(cod_child_noncood_parentN['n'],gm.div_fact(cod_child_noncood_parentN['n'],sumtotn),
                                                       cod_child_noncood_parentN['c'],gm.div_fact(cod_child_noncood_parentN['c'],sumtotc),
                                                       cod_child_noncood_parentN['b'],gm.div_fact(cod_child_noncood_parentN['b'],sumtotb)))
    fin.write("Both_NonCoding\t%s,%s\t%s,%s\t%s,%s\n"%(non_cood_bothN['n'],gm.div_fact(non_cood_bothN['n'],sumtotn),
                                                       non_cood_bothN['c'],gm.div_fact(non_cood_bothN['c'],sumtotc),
                                                       non_cood_bothN['b'],gm.div_fact(non_cood_bothN['b'],sumtotb)))
        
    n_coding_both_segregation ={'ATIT:Same_amino_acid':0,'ATIT:Diff_amino_acid':0,'Variation_in_coding_coods':0,                                'pseudoC':0,'pseudoN':0,'No_overlap_exists':0, 'relevant':0}
    c_coding_both_segregation = {'ATIT:Same_amino_acid':0,'ATIT:Diff_amino_acid':0,'Variation_in_coding_coods':0,                                'pseudoC':0,'pseudoN':0,'No_overlap_exists':0, 'relevant':0}
    b_coding_both_segregation = {'ATIT:Same_amino_acid':0,'ATIT:Diff_amino_acid':0,'Variation_in_coding_coods':0,                                'pseudoC':0,'pseudoN':0,'No_overlap_exists':0, 'relevant':0}
    N_zero=0; C_zero=0; B_zero=0
    N_lis=[];B_lis=[];C_lis=[]
    N_lis_all=[];B_lis_all=[];C_lis_all=[]
    for i in n_case:
        if i[3] in n_coding_both_segregation:
            n_coding_both_segregation[i[3]]+=1
        else:
            n_coding_both_segregation['relevant']+=1
            if i[4] not in [0,1,-1,2,-2]:
                N_lis+=[(i[1],i[4],i[-1])]
            else:
                N_zero+=1
        N_lis_all+=[i[1]]
            
    for i in b_case:
        if i[3] in b_coding_both_segregation:
            b_coding_both_segregation[i[3]]+=1
        else:
            b_coding_both_segregation['relevant']+=1
            B_zero= B_zero+1 if int(i[4].split(",")[0]) in [0,1,-1,2,-2] else B_zero+0
            if int(i[4].split(",")[0]) not in [0,1,-1,2,-2]:
                B_lis+=[(i[1],i[4],i[-1])]
            else:
                B_zero+=1
        B_lis_all+=[i[1]]
    for i in c_case:
        if i[3] in c_coding_both_segregation:
            c_coding_both_segregation[i[3]]+=1
        else:
            c_coding_both_segregation['relevant']+=1
            if i[4] not in [0,1,-1,2,-2]:
                C_lis+=[(i[1],i[4],i[-1])]
            else:
                C_zero+=1
        C_lis_all+=[i[1]]
    sumN=sum(n_coding_both_segregation.values())
    sumC=sum(c_coding_both_segregation.values())
    sumB=sum(b_coding_both_segregation.values())
    fin.write("\n\nN segregation\nType\tTotal\tFreq\n")
    for keys in keyvals:
        fin.write("%s\t%s\t%s\n"%(keys,n_coding_both_segregation[keys],gm.div_fact(n_coding_both_segregation[keys],sumN)))
    fin.write("\n\nC segregation\nType\tTotal\tFreq\n")
    for keys in keyvals:
        fin.write("%s\t%s\t%s\n"%(keys,c_coding_both_segregation[keys],gm.div_fact(c_coding_both_segregation[keys],sumC)))
    fin.write("\n\nB segregation\nType\tTotal\tFreq\n")
    for keys in keyvals:
        fin.write("%s\t%s\t%s\n"%(keys,b_coding_both_segregation[keys],gm.div_fact(b_coding_both_segregation[keys],sumB)))
    fin.write("\n\nFor N from total relevant :%s cases %s(%s) have length change of -2 to 2"%    (n_coding_both_segregation['relevant'],N_zero,gm.div_fact(float(N_zero), n_coding_both_segregation['relevant'])))
    fin.write("\nFor C from total relevant :%s cases %s(%s) have length change of -2 to 2"%    (c_coding_both_segregation['relevant'],C_zero,gm.div_fact(float(C_zero), c_coding_both_segregation['relevant'])))
    fin.write("\nFor B from total relevant :%s cases %s(%s) have length change of -2 to 2"%    (b_coding_both_segregation['relevant'],B_zero,gm.div_fact(float(B_zero), b_coding_both_segregation['relevant'])))
                   
                   
    fin.write("\nAll genes in which N cases are presnet are:%s"%(len(set(N_lis_all))))
    fin.write("\nGenes in which N cases are presnet and both coding with length out of -2 to 2 arr:%s"              %(len(set(N_lis))))
    fin.write("\nAll genes in which C cases are presnet are:%s"%(len(set(C_lis_all))))
    fin.write("\nGenes in which C cases are presnet and both coding with length out of -2 to 2 arr:%s"              %(len(set(C_lis))))
    fin.write("\nAll genes in which B cases are presnet are:%s"%(len(set(B_lis_all))))
    fin.write("\nGenes in which B cases are presnet and both coding with length out of -2 to 2 arr:%s"              %(len(set(B_lis))))
    fin.write("\nAll genes in which NCB union cases are presnet are:%s"%              (len(set(N_lis_all)|set(C_lis_all)|set(B_lis_all))))
    fin.write("\nGenes in which uninon NCB cases are presnet and both coding with length out of -2 to 2 arr:%s"              %(len(set(N_lis)|set(B_lis)|set(C_lis))))
    fin.write("\nAll genes in which NCB intersection cases are presnet are:%s"%              (len(set(N_lis_all)&set(C_lis_all)&set(B_lis_all))))
    fin.write("\nGenes in which intersection NCB cases are presnet and both coding with length out of -2 to 2 arr:%s"              %(len(set(N_lis)&set(B_lis)&set(C_lis))))
    
    fin.write("\n\nGenes_list_N\nGene\tLength\tStructure\tVariants\tName\n")
    for i in set(N_lis):
        #print i
        fin.write("%s\t%s\t%s\t%s\t%s\n"%(i[0],i[1],i[2],structure_val(humanGeneObject[i[0]]),humanGeneObject[i[0]].detail))
    
    fin.write("\n\nGenes_list_C\nGene\tLength\tStructure\tVariants\tName\n")
    for i in set(C_lis):
        fin.write("%s\t%s\t%s\t%s\t%s\n"%(i[0],i[1],i[2],structure_val(humanGeneObject[i[0]]),humanGeneObject[i[0]].detail))
    
    fin.write("\n\nGenes_list_B\nGene\tLength\tStructure\tVariants\tName\n")
    for i in set(B_lis):
        fin.write("%s\t%s\t%s\t%s\t%s\n"%(i[0],i[1],i[2],structure_val(humanGeneObject[i[0]]),humanGeneObject[i[0]].detail))
    
matrix=n_case+c_case+b_case+alt_case
with open (os.path.join(results_dir_csv,"matrix_NCBI_"+fnameappend),"w") as fin:
    fin.write("Event\tGene\tEx_pair\tPosition_affected\tchangeinLength\tParent_aaseq\tChild_aaseq\t    Unique_aa_seq\tSS_parent\tSS_child\tUniqueSS\tStride_parent\tStride_child\tUniqueStride\n")
    for i in matrix:
        fin.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11],i[12],i[13]))
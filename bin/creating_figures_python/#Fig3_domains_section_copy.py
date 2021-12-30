
from posixpath import supports_unicode_filenames
import pandas as pd
import numpy as np
import re, scipy.stats, os
import cPickle as pickle
import sys
import figPanels.modules_analysis as am

print (len(sys.argv))
if 6>len(sys.argv)<8:
    print ('Please type 1. object 2. condFile 3 outputdir+append 4 pfam_or_cath, 5 ATIFlag  ATITrue=0, ATIFalse=1, regardless=-1, 6(optional) another filterChcek(pickleobject)')
    sys.exit()
if len(sys.argv)==6:
    prog,source_gene_object, filtergfname, results_dir,CONST_val, ATIFlag =sys.argv
    additionalFilter=False
else:
    prog,source_gene_object, filtergfname, results_dir,CONST_val, ATIFlag, additionalFilter=sys.argv
ATIFlag=int(ATIFlag)
fnameappend='both' if ATIFlag ==-1 else 'core' if ATIFlag==1 else 'atit' if ATIFlag==0 else 'False'

'''
additional filetr will have the follwing format
has{trans}=(gene,PI,EClis)
'''
def isoform_giver(gene_ob,CONST_val):
    '''
    Note the documenation, what was there in the file
    '''
    relev_exons=[i for i in gene_ob.PI.exons if i.length>0]
    domainsTranscript=gene_ob.PI.cath_list if CONST_val=='cath' else gene_ob.PI.pfam_list
    if domainsTranscript:
        # all the coding exons
        if CONST_val=='cath':
            relev_doms=[]
            for dom in domainsTranscript:
                flag=1
                for j in dom[2]:
                    if j<0.8:
                        #is it likely the coverage ? i dont know!!!!!!!!1
                        flag=0
                if flag:
                    relev_doms+=[dom]
        else:
            relev_doms=[i for i in gene_ob.PI.pfam_list if i[2]>=0.5]
            '''
            pfam_list looks like follwing
            [((1, 64), 'PF02295.17', 0.96, 'z-alpha', 'Domain'), ((209, 274), 'PF00035.26', 0.99, 'dsrm', 'Domain'), ((320, 385), 'PF00035.26', 0.99, 'dsrm', 'Domain'), ((432, 497), 'PF00035.26', 0.99, 'dsrm', 'Domain'), ((591, 920), 'PF02137.18', 1.0, 'A_deamin', 'Family')]
            ''' 
            relev_doms=relev_doms if relev_doms else False
        return relev_exons,relev_doms
    else:
        return False,False

def dumppickle(fname,has):
    with open (fname,'wb') as fout:
        pickle.dump(has,fout)
def loadpickle(fname):
    with open(fname,'rb') as fin:
        has=pickle.load(fin)
    return has

def div_fact(num,denom):
	try:
		return round(float(num)/denom,3)	
	except:
		return 0

#mod_out.files_writer
def update_hash(has,val,gene):
    #print has
    if val in has:
        has[val]+=[gene]
    else:
        has[val]=[gene]
    return has

def exon_screener(exList,consecCoding, ATIFlag):
    exlist2=[]
    for ex in exList:
        baseCondition=ex.length>0 # it should be coding
        #ATITrue=0, ATIFalse=1, reagrdless=-1
        #conditionATIATT=True if ATIFlag ==-1 else numericFlag in tExonList if ATIFlag==1 else numericFlag not in tExonList if ATIFlag==0 else False
        if ATIFlag==-1:
            choosenExonCondition=True
        elif ATIFlag==1:
            choosenExonCondition=int(ex.ID.split(".")[3]) in consecCoding
        else:
            #ATIFlag==0
            choosenExonCondition=int(ex.ID.split(".")[3]) not in consecCoding
        if (baseCondition and choosenExonCondition):
            exlist2+=[ex]
    return exlist2    

def domain_exon_fraction(human):
    gene_domain_annotation={}
    gene_exons_fraction={}
    #this will give me the number of const exons, number of alt exopns, 
    #how many of them code for 100% of domains, and how mnay domains are completely
    #contained within the exons
    #how mnay code for null
    # if split domains, thi inofmration will be tased outr by domain whole and the whole middle
    '''
    above two dicts will store the splitting information of domains (those are not enveloped in 1 exon), will calculate 
    [(div_fact(length_m[key[0]],sum(length_m.values())), div_fact(length_m[key[1]],sum(length_m.values())))]
    0) the fraction contribution of A exons from total exons, and 1) the fraction of G exons from total exons if key ==AG,  
    # will only store their fractions, 
    '''
    total_aa = []
    total_aa_dom = []
    #above two parameters will be used to judge the coverage amount of domains in those subsets
    #domains spans (d_aa) per single transcript of the gene will be added to total_aa_dom
    
    c = 0
    per_prot_dom=[]
    for gene in human:
        tExonList=am.exonScreenerBetweenTConstitutive(human[gene].exons)
        exlis,domlis=isoform_giver(human[gene],CONST_val)
        relev_ex=exon_screener(exlis, tExonList, ATIFlag) if domlis and exlis else False
        if domlis and relev_ex:
            #print (len(domlis))
            total_aa += [sum([ex.length for ex in relev_ex])]
            d_aa = 0
            for dom in domlis:
                domain_name=dom[1]+':'+dom[3]
                domain_namewith_cood=domain_name+'::'+str(dom[0][0])+','+str(dom[0][1])
                length={'A':0,'G':0}
                length_m={'A':0,'G':0}
                domain_constituents=''
                domain_constituents_middle=''
                span_temp=range(dom[0][0], dom[0][1]+1)
                dspan = set(span_temp)
                dspan_middle=set(span_temp[int(len(span_temp)*0.3):int(len(span_temp)*0.6)])
                d_aa += len(dspan)
                su = 0
                Flag=False
                for ex in exlis:
                    '#here the exons are individual not together'
                    espan = set(range(su, su+ex.length))
                    su += ex.length
                    interori = espan & dspan
                    intermidd = espan & dspan_middle
                    
                    dfrac = div_fact(len(interori),len(dspan))
                    efrac = div_fact(len(interori),len(espan))
                    dfrac_m=div_fact(len(intermidd),len(dspan_middle))
                    efrac_m = div_fact(len(intermidd),len(espan))
                    if (dfrac or dfrac_m) and ex in relev_ex:
                        if gene not in gene_exons_fraction:
                            gene_exons_fraction[gene]={}
                        if ex.ID not in gene_exons_fraction:
                            gene_exons_fraction[gene][ex.ID]=[]
                        
                        gene_exons_fraction[gene][ex.ID]+=[(domain_namewith_cood, dfrac, dfrac_m, efrac_m, efrac)]
                        key=ex.ID.split('.')[2] if ex.ID[0]!='R' else 'A'
                        key='G' if key == 'F' else key

                        if dfrac:
                            domain_constituents+=key
                            length[key]+=len(interori)
                        if dfrac_m:
                            Flag = True
                            domain_constituents_middle+=key
                            length_m[key]+=len(intermidd)
                #domain restarted, 
                
                #print (dspan,domain_namewith_cood,dspan_middle  )
                #print (domain_constituents, '1', domain_constituents_middle, '2')
                if Flag:
                    #means some fraction got encountered
                    keyw=list(set(domain_constituents))
                    keym=list(set(domain_constituents_middle))
                    keym.sort()
                    keyw.sort()
                    keym="".join(keym)
                    keyw="".join(keyw)
                    keyw = keyw*2 if len(keyw)==1 else keyw
                    #they will be largely, AG, GA , GG
                    #print key, length
                    
                    keym = keym*2 if len(keym)==1 else keym
                    
                    #print (keyw,keym)
                    key1contrib_W=div_fact(length[keyw[0]],sum(length.values()))
                    key2contrib_W=div_fact(length[keyw[1]],sum(length.values()))

                    key1contrib_M=div_fact(length_m[keym[0]],sum(length_m.values()))
                    key2contrib_M=div_fact(length_m[keym[1]],sum(length_m.values()))
                    if gene not in gene_domain_annotation:
                        gene_domain_annotation[gene]={}
                    if domain_namewith_cood not in gene_domain_annotation[gene]:
                        gene_domain_annotation[gene][domain_namewith_cood]={100:[],40:[]}
                    
                    sum_VW=key1contrib_W+key2contrib_W if keyw[0]!=keyw[1] else key1contrib_W
                    sum_VM=key1contrib_M+key2contrib_M if keym[0]!=keym[1] else key1contrib_M

                    if sum_VM>=0.9 and sum_VW>=0.9:                   
                        gene_domain_annotation[gene][domain_namewith_cood][100]=(keyw, key1contrib_W, key2contrib_W)
                        gene_domain_annotation[gene][domain_namewith_cood][40]=(keym, key1contrib_M, key2contrib_M)

    return gene_domain_annotation, gene_exons_fraction        




human= loadpickle(source_gene_object)
sourcedir=os.path.dirname(os.path.normpath(source_gene_object))
print (sourcedir)
filterg=loadpickle(filtergfname)

human={i:human[i] for i in human if i in filterg}
if additionalFilter:
    addfilter=loadpickle(additionalFilter)
    human={i:human[i] for i in human if i in addfilter}

dom_annot,exons_annot=domain_exon_fraction(human=human)

ex_annot_refined={}
for gene in dom_annot:
    affected_dom_list=set(list(dom_annot[gene].keys()))
    exlist=exons_annot[gene].keys()
    for ex in exlist:
        for domiter in exons_annot[gene][ex]:
            dname,dfw,dfm,efw,efm=domiter
            if dname in affected_dom_list:
                if gene not in ex_annot_refined:
                    ex_annot_refined[gene]={}
                if ex not in ex_annot_refined[gene]:
                    ex_annot_refined[gene][ex]=[]
                ex_annot_refined[gene][ex]+=[(dname,dfw,dfm,efw,efm)]

                
        

dumppickle(results_dir+'domain_annotation_%s.pickle'%fnameappend, dom_annot)
dumppickle(results_dir+'exonsWise_%s.pickle'%fnameappend, ex_annot_refined)
#python bin_analysis/domain_section.py dataset/objectsave_9606.pick results/ pfam
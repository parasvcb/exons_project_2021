
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
    relev_exons=[]
    for ex in gene_ob.PI.exons:
        key=ex.ID.split('.')[2] if ex.ID[0]!='R' else 'A'
        key='G' if key == 'F' else key
        if ex.length>0:
            relev_exons+=[[ex.length,key,ex.ID,1,0]]# -2val: include in analysis as 1 and remove in analysis as 0, -1val:0 for single exon, 1 maybe nex for merged entries
    
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

def writeRangeFile(has,fname):
    occRangeList=list(has.keys())
    occRangeList.sort()
    totalOcc=sum(has.values())
    with open (fname,'w') as fin:
        fin.write('Range\tTotal\tFreq\tCumTotal\tCumFreq\n')
        previousValue=0
        for occRange in occRangeList:
            occrangestr="_".join(map(str,list(occRange)))
            count_occrange=has[occRange]
            cumulative_count=previousValue+count_occrange
            previousValue+=count_occrange
            fin.write('%s\t%s\t%s\t%s\t%s\n'%(occrangestr,count_occrange,div_fact(count_occrange,totalOcc), cumulative_count,div_fact(cumulative_count,totalOcc)))
    return

        

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
    for ex in exList:
        #[ex.length,key,ex.ID,1,0]
        #length,consFlag,ID, inclusion Flag, merged entry or not
        baseCondition=ex[0] # it should be coding
        #ATITrue=0, ATIFalse=1, reagrdless=-1
        #conditionATIATT=True if ATIFlag ==-1 else numericFlag in tExonList if ATIFlag==1 else numericFlag not in tExonList if ATIFlag==0 else False
        if ATIFlag==-1:
            choosenExonCondition=True
        elif ATIFlag==1:
            choosenExonCondition=int(ex[2].split(".")[3]) in consecCoding
        else:
            #ATIFlag==0
            choosenExonCondition=int(ex[2].split(".")[3]) not in consecCoding
        if (baseCondition and choosenExonCondition):
            ex[3]=1
        else:
            ex[3]=0
    return exList    

def exon_screener_consecCoding(exlis, tExonList, ATIFlag, geneOB):
    consecCodingLis=geneOB.connst_togetherness_coding()
    consecCodingHas={}
    '''
    As of now there are three Choices, ATIT no G exons so a waste, other two will capture all the G exons
    '''

    # for setex in consecCodingLis:
    #     print ([i.ID for i in setex], 'setex')
    Flag=sum([len(i) for i in consecCodingLis])

    consecCodingLis=[i for i in consecCodingLis if len(i)]
    if consecCodingLis and Flag:    
        for chunk in consecCodingLis:
            chunkRep=chunk[0].ID
            length=sum([j.length for j in chunk])
            consecCodingHas[chunkRep]=[length,chunk]
        exrefined=[]
        removeLis=[]
        #print(consecCodingHas)
   
        for ex in exlis:
            #length,consFlag,ID
            if ex[2] in consecCodingHas:
                lengthT,membersT=consecCodingHas[ex[2]]
                #print (consecCodingHas[ex[2]][1])
                #[91, [<constructing_data.Classes_exons.Exon instance at 0x7fac16953e10>, <constructing_data.Classes_exons.Exon instance at 0x7fac16953b40>]]
                partnerMergedNames=[j.ID for j in membersT]
                newname=",".join(partnerMergedNames)
                exrefined+=[[lengthT,'G',newname,ex[3],1 if len(partnerMergedNames)>1 else 0]]
                removeLis+=partnerMergedNames[1:]
            else:
                exrefined+=[ex]
        #print (removeLis,'removelis')
        exrefined=[i for i in exrefined if i[2] not in removeLis]
        #print (exlis,'exlis')
        #print (exrefined,'exrefined')
    else:
        exrefined=exlis
    #this entity should remove the other exons present in the isoform.

    # for trans in geneOB.transcripts:
    #     print (trans.ID, [ex.ID for ex in trans.exons])

    # for i in consecCodingLis:
    #     print ([j.ID for j in i])
    return exrefined
    
def isolateEfracDfrac(domainContrib,exonsContrib):
    #omainContrib[ex[2]]=[domain_namewith_cood, dfrac, dfrac_m, efrac_m, efrac, espan, interori, intermidd]
    #print (1, exonsContrib)
    for ex in exonsContrib:
        if ex in domainContrib:
            interori=len(domainContrib[ex][6])
            intermid=len(domainContrib[ex][7])
            exonsContrib[ex]['W'][1]+=interori
            exonsContrib[ex]['M'][1]+=intermid
            #print (exonsContrib[ex]['W'][1],)
            exonsContrib[ex]['W'][2]=div_fact(exonsContrib[ex]['W'][1],exonsContrib[ex]['W'][0])
            exonsContrib[ex]['M'][2]=div_fact(exonsContrib[ex]['M'][1],exonsContrib[ex]['M'][0])
    domainContrib={i:[domainContrib[i][0],domainContrib[i][1],domainContrib[i][2]] for i in domainContrib}
    #print (1, exonsContrib)
    return exonsContrib,domainContrib


def domainExonIntersection(exlis,domain_namewith_cood, dspan, dspan_middle):
    lengthHas={'W':{'A':0,'G':0},'M':{'A':0,'G':0}}
    domain_constituents={'W':'','M':''}
    contained40=False
    contained100=False
    su=0
    domainContrib={}
    #print ('START')
    #print (dspan)
    #print (dspan_middle)
    for ex in exlis:
        espan = set(range(su, su+ex[0]))
        #print (espan, 'espan', ex[2])
        su += ex[0]
        interori = espan & dspan
        intermidd = espan & dspan_middle
        dfrac = div_fact(len(interori),len(dspan))
        efrac = div_fact(len(interori),len(espan))
        dfrac_m=div_fact(len(intermidd),len(dspan_middle))
        efrac_m = div_fact(len(intermidd),len(espan))
        #print (interori, intermidd, dfrac, dfrac_m)
        condition_exon_passsed_filters=ex[3]
        if (dfrac or dfrac_m) and condition_exon_passsed_filters:
            #print ('yes')
            domainContrib[ex[2]]=[domain_namewith_cood, dfrac, dfrac_m, efrac_m, efrac, espan, interori, intermidd]
            key=ex[1]
            if dfrac:
                domain_constituents['W']+=key
                lengthHas['W'][key]+=len(interori)
                if dfrac>=0.95:
                    contained100=key
            if dfrac_m:
                domain_constituents['M']+=key
                lengthHas['M'][key]+=len(intermidd)
                if dfrac_m>=0.95:
                    contained40=key

    sum_VM, sum_VW, keyw, keym, key1contrib_W, key2contrib_W, key1contrib_M, key2contrib_M = [False]*8

    if domain_constituents['W']:
        keyw=list(set(domain_constituents['W']))
        keyw.sort()
        keyw="".join(keyw)
        keyw = keyw*2 if len(keyw)==1 else keyw
        key1contrib_W=div_fact(lengthHas['W'][keyw[0]],sum(lengthHas['W'].values()))
        key2contrib_W=div_fact(lengthHas['W'][keyw[1]],sum(lengthHas['W'].values()))
        sum_VW=key1contrib_W+key2contrib_W if keyw[0]!=keyw[1] else key1contrib_W

    if domain_constituents['M']:
        keym=list(set(domain_constituents['M']))
        keym.sort()
        keym="".join(keym)
        keym = keym*2 if len(keym)==1 else keym
        key1contrib_M=div_fact(lengthHas['M'][keym[0]],sum(lengthHas['M'].values()))
        key2contrib_M=div_fact(lengthHas['M'][keym[1]],sum(lengthHas['M'].values()))
        sum_VM=key1contrib_M+key2contrib_M if keym[0]!=keym[1] else key1contrib_M

    
    return sum_VW, domainContrib, contained40, contained100, keyw, keym, key1contrib_W, key2contrib_W, key1contrib_M, key2contrib_M
    

def domain_exon_fraction(human):
    gene_exons_domains_fractions={}
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
    total_aa_dom_middle=[]
    #above two parameters will be used to judge the coverage amount of domains in those subsets
    #domains spans (d_aa) per single transcript of the gene will be added to total_aa_dom
    
    c = 0
    per_prot_dom={(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}
    per_prot_dom_middle={(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}
    changeinCONSTExonONMerger={(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}

    #testl=[131076]
    #for gene in testl:
    for gene in human:
    
        tExonList=am.exonScreenerBetweenTConstitutive(human[gene].exons)
        exlis,domlis=isoform_giver(human[gene],CONST_val)
        
        
        #print (gene, bool(domlis), bool(exlis))
        if domlis and exlis:
            flag_exons=False
            flag_exons_conseCoding=False
            c+=1
            domcovered_w=0
            domcovered_m=0
            exlis=exon_screener(exlis, tExonList, ATIFlag) if domlis and exlis else False
            #print (exlis)
            consecCoding=exon_screener_consecCoding(exlis, tExonList, ATIFlag, human[gene])
            
            #print (consecCoding)
            
            proteinLength=sum([ex[0] for ex in exlis if ex[3]])
            exonsContrib_conseCoding={ex[2]:{'W':[ex[0],0,0],'M':[ex[0],0,0]} for ex in consecCoding}
            exonsContrib={ex[2]:{'W':[ex[0],0,0],'M':[ex[0],0,0]} for ex in exlis}

            fracExchange=round(float(len([i for i in exlis if i[3]])-len([i for i in consecCoding if i[3]]))/len(exlis),2)
            changeinCONSTExonONMerger=am.has_range_adder(changeinCONSTExonONMerger,fracExchange)
            
            #print (len(domlis))
            #exonsContrib={}
            for dom in domlis:
                domain_name=dom[1]+':'+dom[3]
                domain_namewith_cood=domain_name+'::'+str(dom[0][0])+','+str(dom[0][1])
                lengthHas={'W':{'A':0,'G':0},'M':{'A':0,'G':0}}
                # domain_constituents=''
                # domain_constituents_middle=''
                span_temp=range(dom[0][0], dom[0][1]+1)
                dspan = set(span_temp)
                dspan_middle=set(span_temp[int(len(span_temp)*0.3):int(len(span_temp)*0.6)])
                # total_aa_dom+=len(dspan)
                # total_aa_dom_middle+=len(dspan_middle)

                sum_VW, domainContrib, contained40, contained100,keyw, keym, key1contrib_W, key2contrib_W, key1contrib_M, key2contrib_M=domainExonIntersection(exlis, domain_namewith_cood, dspan, dspan_middle)
                #print (domainContrib)
                domcovered_w+=sum([len(domainContrib[ex][6]) for ex in domainContrib])
                #print ('normal')

                sum_VW_conseCoding, domainContrib_conseCoding, contained40_conseCoding, contained100_conseCoding,keyw_conseCoding, keym_conseCoding, key1contrib_W_conseCoding, key2contrib_W_conseCoding, key1contrib_M_conseCoding, key2contrib_M_conseCoding=domainExonIntersection(consecCoding, domain_namewith_cood, dspan, dspan_middle)
                domcovered_m+=sum([len(domainContrib_conseCoding[ex][7]) for ex in domainContrib_conseCoding])
                #print ('consec')

                #domainContrib[ex[2]]=[domain_namewith_cood, dfrac, dfrac_m, efrac_m, efrac, espan, interori, intermidd]
                
                #print (sum_VW, sum_VW_conseCoding)
                
                if sum_VW>=0.9:
                    if gene not in gene_exons_domains_fractions:
                        gene_exons_domains_fractions[gene]={'domains':{},'exons':{},'domains_conseCoding':{},'exons_conseCoding':{}}
                    exonsContrib,domainContrib=isolateEfracDfrac(domainContrib,exonsContrib)
                    gene_exons_domains_fractions[gene]['domains'][domain_namewith_cood]=[(keyw, key1contrib_W, key2contrib_W),(keym, key1contrib_M, key2contrib_M) ,contained40, contained100 ,domainContrib]
                    flag_exons=True

                    #print (domainContrib, 'dt')
                    #print (exonsContrib, 'et')

                if sum_VW_conseCoding>=0.9:
                    exonsContrib_conseCoding,domainContrib_conseCoding=isolateEfracDfrac(domainContrib_conseCoding,exonsContrib_conseCoding)
                    gene_exons_domains_fractions[gene]['domains_conseCoding'][domain_namewith_cood]=[(keyw_conseCoding, key1contrib_W_conseCoding, key2contrib_W_conseCoding),(keym_conseCoding, key1contrib_M_conseCoding, key2contrib_M_conseCoding) ,contained40_conseCoding, contained100_conseCoding ,domainContrib_conseCoding]
                    flag_exons_conseCoding=True
        
            #print (exonsContrib)
            #print ('\n')
            #print (exonsContrib_conseCoding)
            if flag_exons:
                gene_exons_domains_fractions[gene]['exons']=exonsContrib
            if flag_exons_conseCoding:
                gene_exons_domains_fractions[gene]['exons_conseCoding']=exonsContrib_conseCoding
        
            total_aa_dom+=[domcovered_w]
            total_aa_dom_middle+=[domcovered_m]
            total_aa+=[proteinLength]
            frac_per_prot_cover=div_fact(domcovered_w,proteinLength)
            frac_per_prot_cover_m=div_fact(domcovered_w,proteinLength)
            per_prot_dom=am.has_range_adder(per_prot_dom,frac_per_prot_cover)
            per_prot_dom_middle=am.has_range_adder(per_prot_dom_middle,frac_per_prot_cover_m)

        
        #break
    print ('Total_genes=%s, included cases=%s, total_aa=%s, total_aa_dom=%s, total_aa_dom_middle=%s, fracCoveerd_wdom=%s, fraccoveerd_m_dom=%s'%(len(human),c,sum(total_aa), sum(total_aa_dom), sum(total_aa_dom_middle), div_fact(sum(total_aa_dom), sum(total_aa)), div_fact(sum(total_aa_dom_middle), sum(total_aa))))

    return gene_exons_domains_fractions, per_prot_dom, per_prot_dom_middle, changeinCONSTExonONMerger



def test(geneob):
    relev_exons=[]
    #print ([t.ID for t in geneob.transcripts])
    for t in geneob.transcripts:
        print (t.ID, [ex.ID for ex in t.exons])
    for ex in geneob.PI.exons:
        print (geneob.PI.ID,ex.ID, ex.length)
        key=ex.ID.split('.')[2] if ex.ID[0]!='R' else 'A'
        key='G' if key == 'F' else key
        #print (key)
        if ex.length>0:
            relev_exons+=[[ex.length,key,ex.ID,1,0]]# -2val: include in analysis as 1 and remove in analysis as 0, -1val:0 
    print (relev_exons)

    domainsTranscript=geneob.PI.cath_list if CONST_val=='cath' else geneob.PI.pfam_list
    print (domainsTranscript)

human= loadpickle(source_gene_object)
sourcedir=os.path.dirname(os.path.normpath(source_gene_object))
print (sourcedir)
filterg=loadpickle(filtergfname)

human={i:human[i] for i in human if i in filterg}
if additionalFilter:
    addfilter=loadpickle(additionalFilter)
    human={i:human[i] for i in human if i in addfilter}

#test(human[5328])

#sys.exit()

dom_ex_annot, w_per_prot, m_per_prot, changeinCONSTExonONMerger=domain_exon_fraction(human=human)
writeRangeFile(w_per_prot,results_dir+'F3_per_proteinDomain_covered_W_%s.csv'%fnameappend)
writeRangeFile(m_per_prot,results_dir+'F3_per_proteinDomain_covered_M_%s.csv'%fnameappend)
writeRangeFile(changeinCONSTExonONMerger,results_dir+'F3_exonsChangeFractionOnMerge_%s.csv'%fnameappend)

print (fnameappend)
dumppickle(results_dir+'domain_exon_annotation_%s.pickle'%fnameappend, dom_ex_annot)
#python bin_analysis/domain_section.py dataset/objectsave_9606.pick results/ pfam
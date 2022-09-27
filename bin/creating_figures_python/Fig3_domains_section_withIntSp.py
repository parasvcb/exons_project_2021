import pandas as pd
import numpy as np
import re, scipy.stats, os
import cPickle as pickle
import sys
import figPanels.modules_analysis as am

print (len(sys.argv))
if len(sys.argv)!=5:
    print ('Please type 1. object 2. condFile 3 outputdir+append 4 pfam_or_cath')
    sys.exit()

prog,source_gene_object, filtergfname, results_dir,CONST_val =sys.argv
additionalFilter=False

# ATIFlag=int(ATIFlag)
# fnameappend='both' if ATIFlag ==-1 else 'core' if ATIFlag==1 else 'atit' if ATIFlag==0 else 'False'

'''
additional filetr will have the follwing format
has{trans}=(gene,PI,EClis)
'''
def isoform_giver(gene_ob,CONST_val,per_dom_coverage):
    '''
    Note the documenation, what was there in the file
    '''
    relev_exons=[]
    for ex in gene_ob.PI.exons:
        key=ex.ID.split('.')[2] if ex.ID[0]!='R' else 'A'
        key='A' if key == 'F' else key
        if ex.length>0:
            relev_exons+=[[ex.length,key,ex.ID,1]]
            # -1 val: if presnet in core
    
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
            relev_doms=[]
            for i in gene_ob.PI.pfam_list:
                per_dom_coverage=am.has_range_adder(per_dom_coverage,i[2])
                if i[2]>=0.5:
                    relev_doms+=[i]
            
            '''
            pfam_list looks like follwing
            [((1, 64), 'PF02295.17', 0.96, 'z-alpha', 'Domain'), ((209, 274), 'PF00035.26', 0.99, 'dsrm', 'Domain'), ((320, 385), 'PF00035.26', 0.99, 'dsrm', 'Domain'), ((432, 497), 'PF00035.26', 0.99, 'dsrm', 'Domain'), ((591, 920), 'PF02137.18', 1.0, 'A_deamin', 'Family')]
            ''' 
            relev_doms=relev_doms if relev_doms else False
        return relev_exons,relev_doms,per_dom_coverage
    else:
        return False,False,per_dom_coverage

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

def exon_screener(exList,consecCoding):
    consecCodingrange=(consecCoding[0],consecCoding[-1]) if consecCoding else False
    for ex in exList:
        ExonCoreCondition=consecCodingrange[0]<=int(ex[2].split(".")[3])<=consecCodingrange[-1] if consecCodingrange else False
        if ExonCoreCondition:
            ex[3]=1
        else:
            ex[3]=0
    return exList    

def exon_screener_consecCoding(exlis, tExonList, geneOB):
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
                exrefined+=[[lengthT,'G',newname,ex[3]]]
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
            #coreExon=domainContrib[ex][8]
            exonsContrib[ex]['W'][2]+=interori
            exonsContrib[ex]['M'][2]+=intermid
            #print (exonsContrib[ex]['W'][1],)
            exonsContrib[ex]['W'][3]=div_fact(exonsContrib[ex]['W'][1],exonsContrib[ex]['W'][0])
            exonsContrib[ex]['M'][3]=div_fact(exonsContrib[ex]['M'][1],exonsContrib[ex]['M'][0])
    domainContrib={i:[domainContrib[i][1],domainContrib[i][2],domainContrib[i][8]] for i in domainContrib}
    #print (1, exonsContrib)
    return exonsContrib,domainContrib


def domainExonIntersection(exlis,domain_namewith_cood, dspan, dspan_middle,per_dom_dfrac_range_Full=False):
    domain={
        'considered':'',
        #in order of CORE,ATIT,BOTH
        'containedM':['','',''],
        'containedW':['','',''],
        'exContrib':{},
        #in order of CORE,ATIT,BOTH
        'splitM':[{'A':[0,''],'G':[0,'']},{'A':[0,''],'G':[0,'']},{'A':[0,''],'G':[0,'']}],
        'splitW':[{'A':[0,''],'G':[0,'']},{'A':[0,''],'G':[0,'']},{'A':[0,''],'G':[0,'']}],
        #+ in the end [combcategory[tuple(val1[0])], val1[0],val1[1], bothCateg[tuple(val2[0])],val2[0],val2[1]]
        }
  
    su=0
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
        if (dfrac or dfrac_m):
            coreExon=ex[3]
            ind=0 if coreExon else 1
            domain['considered']=True
            #ex[2] is id and ex[1] is 'A','G' flag
            domain['exContrib'][ex[2]]=[domain_namewith_cood, dfrac, dfrac_m, efrac_m, efrac, espan, interori, intermidd,coreExon]
            key=ex[1]
            if dfrac:
                if per_dom_dfrac_range_Full:
                    per_dom_dfrac_range_Full=am.has_range_adder(per_dom_dfrac_range_Full,dfrac)
                domain['splitW'][ind][key][0]+=len(interori)
                domain['splitW'][2][key][0]+=len(interori) #both should also be added to 
                proceedW=True
                if dfrac>=0.95:
                    domain['containedW'][ind]=key
                    domain['containedW'][2]=key #both should also be added to 
            if dfrac_m:
                domain['splitM'][ind][key][0]+=len(interori)
                domain['splitM'][2][key][0]+=len(interori) #both should also be added to 
                proceedM=True
                if dfrac_m>=0.95:
                    domain['containedM'][ind]=key
                    domain['containedM'][2]=key #both should also be added to
                    
    combcategory={(1,'A'):'1COREA',(1,'A',1,'G'):'1COREAG',(1,'G'):'1COREG',(2,'A'):'2ATITA',(2,'A',2,'G'):'2ATITAG',(2,'G'):'2ATITG',(1,'A',2,'A'):'1CORE_A__2ATIT_A',(1,'A',2,'G'):'1CORE_A__2ATIT_G',(1,'G',2,'G'):'1CORE_G__2ATIT_G',(1,'G',2,'A'):'1CORE_G__2ATIT_A',(1,'A',1,'G',2,'A'):'1CORE_AG__2ATIT_A',(1,'A',1,'G',2,'G'):'1CORE_AG__2ATIT_G',(1,'A',1,'G',2,'A',2,'G'):'1CORE_AG__2ATIT_AG',(1,'A',2,'A',2,'G'):'1CORE_A__2ATIT_AG',(1,'G',2,'A',2,'G'):'1CORE_G__2ATIT_AG'}
    bothCateg={(3,'A'):'3BOTHA',(3,'A',3,'G'):'3BOTHAG',(3,'G'):'3BOTHG'}
    if proceedW:
        totRes=sum([domain['splitW'][2][junc][0] for junc in domain['splitW'][2]])
        # 2 is both catgory auto symmation of the ATIT and core
        keyComb=''
        keyBoth=''
        val1=[[],[]]
        val2=[[],[]]
        for ind,atitOrCoreOrBoth in enumerate(domain['splitW']):
            for junc in atitOrCoreOrBoth:
                atitOrCoreOrBoth[junc][1]=div_fact(atitOrCoreOrBoth[junc][0],totRes)
            if ind!=2:
                for junc in atitOrCoreOrBoth:
                    if atitOrCoreOrBoth[junc][0]:
                        val1[0]+=[ind+1]
                        val1[0]+=[junc]
                        val1[1]+=[atitOrCoreOrBoth[junc][1]]
            else:
                for junc in atitOrCoreOrBoth:
                    if atitOrCoreOrBoth[junc][0]:
                        val2[0]+=[ind+1]
                        val2[0]+=[junc]
                        val2[1]+=[atitOrCoreOrBoth[junc][1]]
        domain['splitW']+=[combcategory[tuple(val1[0])], val1[0],val1[1], bothCateg[tuple(val2[0])],val2[0],val2[1]]


    if proceedM:
        totRes=sum([domain['splitM'][2][junc][0] for junc in domain['splitM'][2]])
        # 2 is both catgory auto symmation of the ATIT and core
        keyComb=''
        keyBoth=''
        val1=[[],[]]
        val2=[[],[]]
        for ind,atitOrCoreOrBoth in enumerate(domain['splitM']):
            for junc in atitOrCoreOrBoth:
                atitOrCoreOrBoth[junc][1]=div_fact(atitOrCoreOrBoth[junc][0],totRes)
            if ind!=2:
                for junc in atitOrCoreOrBoth:
                    if atitOrCoreOrBoth[junc][0]:
                        val1[0]+=[ind+1]
                        val1[0]+=[junc]
                        val1[1]+=[atitOrCoreOrBoth[junc][1]]
            else:
                for junc in atitOrCoreOrBoth:
                    if atitOrCoreOrBoth[junc][0]:
                        val2[0]+=[ind+1]
                        val2[0]+=[junc]
                        val2[1]+=[atitOrCoreOrBoth[junc][1]]
        domain['splitM']+=[combcategory[tuple(val1[0])],val1[1], bothCateg[tuple(val2[0])],val2[1]]
    
    return domain, per_dom_dfrac_range_Full
    

def formatter(gname,domname,exonsLis,codingexonslis,domcontribNormal, excontribNormal,domcontribcons, excontribcons):
    
    print ('gid:%s'%gname)
    print ('domname:%s'%domname)
    print ('exonList:')
    cumlen=0
    for i in exonsLis:
        print (i,cumlen+i[0])
        cumlen+=i[0]
    print ('ConsecCodingexonsList:')
    cumlen=0
    for i in codingexonslis:
        print (i, cumlen+i[0])
        cumlen+=i[0]
    normL=sum([i[0] for i in exonsLis])
    consL=sum([i[0] for i in codingexonslis])
    print('lengthmatch',normL==consL)
    print ('DomainPerspective:')
    exDom=domcontribNormal.keys()
    exDom=[[int(ex.split('.')[3]), ex] for ex in exDom]
    exDom.sort()
    print (exDom)
    for ex in exDom:
        print (ex[1], domcontribNormal[ex[1]])
    print ('\n')
    print ('DomainPerspective_conseCoding:')
    exDom=domcontribcons.keys()
    exDom=[[int(ex.split('.')[3]), ex] for ex in exDom]
    exDom.sort()
    for ex in exDom:
        print (ex[1],domcontribcons[ex[1]])
    print ('\n')

    
    print ('exonPerspective:')
    exDom=excontribNormal.keys()
    exDom=[[int(ex.split('.')[3]), ex] for ex in exDom]
    exDom.sort()
    for ex in exDom:
        print (ex[1],excontribNormal[ex[1]])
    print ('\n')
    print ('exonPerspective_conseCoding:')
    exDom=excontribcons.keys()
    exDom=[[int(ex.split('.')[3]), ex] for ex in exDom]
    exDom.sort()
    for ex in exDom:
        print (ex[1],excontribcons[ex[1]])
    print ('\n**')

        


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
    per_dom_coverage={(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.701,0.8):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}
    per_dom_dfrac_range_Full={(0.0,0.80):0,(0.801,0.9):0,(0.901,0.950):0,(0.951,0.990):0, (0.991,1.01):0}
    per_dom_dfrac_range_Full2={(0.0,0.80):0,(0.801,0.9):0,(0.901,0.950):0,(0.951,0.990):0, (0.991,1.01):0}
    per_prot_dom={(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.701,0.8):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}
    per_prot_dom_middle={(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.701,0.8):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}
    changeinCONSTExonONMerger={(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.701,0.8):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}

    # setdeepanshi=[7173]
    # for gene in setdeepanshi:
    for gene in human:
        tExonList=am.exonScreenerBetweenTConstitutive(human[gene].exons)
        exlis,domlis, per_dom_coverage=isoform_giver(human[gene],CONST_val,per_dom_coverage)
        #print ('tex',tExonList)
        
        #print (gene, bool(domlis), bool(exlis))
        if domlis and exlis:
            flag_exons=False
            flag_exons_conseCoding=False
            c+=1
            domcovered_w=0
            domcovered_m=0
            exlis=exon_screener(exlis, tExonList) if domlis and exlis else False
            #print ('exlis',exlis)
            consecCoding=exon_screener_consecCoding(exlis, tExonList, human[gene])
            #print ('consecCoding',consecCoding)
            
            proteinLength=sum([ex[0] for ex in exlis])
            exonsContrib_conseCoding={ex[2]:{'W':[ex[0],ex[3],0,0],'M':[ex[0],ex[3],0,0]} for ex in consecCoding}
            exonsContrib={ex[2]:{'W':[ex[0],ex[3],0,0],'M':[ex[0],ex[3],0,0]} for ex in exlis}

            fracExchange=round(float(len([i for i in exlis if i[3]])-len([i for i in consecCoding if i[3]]))/len(exlis),2)
            changeinCONSTExonONMerger=am.has_range_adder(changeinCONSTExonONMerger,fracExchange)
            
            #print (len(domlis))
            #exonsContrib={}
            for dom in domlis:
                #print ("dom",dom)
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

                domainContrib, per_dom_dfrac_range_Full=domainExonIntersection(exlis, domain_namewith_cood, dspan, dspan_middle, per_dom_dfrac_range_Full)
                # domainContrib={
                # 'considered':'',
                # #in order of ATIT,CORE,BOTH
                # 'containedM':['','',''],
                # 'containedW':['','',''],
                # 'exContrib':{}, #exon wise list of list
                #                   [domain_namewith_cood, dfrac, dfrac_m, efrac_m, efrac, espan, interori, intermidd, coreExon]
                # #in order of ATIT,CORE,BOTH
                # 'splitM':[{'A':[0,''],'G':[0,'']},{'A':[0,''],'G':[0,'']},{'A':[0,''],'G':[0,'']}],
                # 'splitW':[{'A':[0,''],'G':[0,'']},{'A':[0,''],'G':[0,'']},{'A':[0,''],'G':[0,'']}],
                # }
                domcovered_w+=sum([len(domainContrib['exContrib'][ex][6]) for ex in domainContrib['exContrib']])
                #print ('normal')

                domainContrib_conseCoding,per_dom_dfrac_range_Full2=domainExonIntersection(consecCoding, domain_namewith_cood, dspan, dspan_middle,per_dom_dfrac_range_Full2)
                domcovered_m+=sum([len(domainContrib_conseCoding['exContrib'][ex][7]) for ex in domainContrib_conseCoding['exContrib']])

                exonsContrib_conseCoding,domainContrib_conseCoding['exContrib']=isolateEfracDfrac(domainContrib_conseCoding['exContrib'],exonsContrib_conseCoding)
                exonsContrib,domainContrib['exContrib']=isolateEfracDfrac(domainContrib['exContrib'],exonsContrib)
                
                if gene not in gene_exons_domains_fractions:
                    gene_exons_domains_fractions[gene]={'domains':[],'exons':{},'domains_conseCoding':[],'exons_conseCoding':{}}
                gene_exons_domains_fractions[gene]['domains']+=[[domain_namewith_cood,domainContrib]]
                gene_exons_domains_fractions[gene]['domains_conseCoding']+=[[domain_namewith_cood,domainContrib_conseCoding]]

            #print (';')
            gene_exons_domains_fractions[gene]['exons']=exonsContrib
            gene_exons_domains_fractions[gene]['exons_conseCoding']=exonsContrib_conseCoding

            
            total_aa_dom+=[domcovered_w]
            total_aa_dom_middle+=[domcovered_m]
            total_aa+=[proteinLength]
            frac_per_prot_cover=div_fact(domcovered_w,proteinLength)
            frac_per_prot_cover_m=div_fact(domcovered_w,proteinLength)
            per_prot_dom=am.has_range_adder(per_prot_dom,frac_per_prot_cover)
            per_prot_dom_middle=am.has_range_adder(per_prot_dom_middle,frac_per_prot_cover_m)
    if 0:
            formatter(gene,domain_namewith_cood,exlis,consecCoding,domainContrib, exonsContrib,domainContrib_conseCoding, exonsContrib_conseCoding)
    print ('Total_genes=%s, included cases=%s, total_aa=%s, total_aa_dom=%s, total_aa_dom_middle=%s, fracCoveerd_wdom=%s, fraccoveerd_m_dom=%s'%(len(human),c,sum(total_aa), sum(total_aa_dom), sum(total_aa_dom_middle), div_fact(sum(total_aa_dom), sum(total_aa)), div_fact(sum(total_aa_dom_middle), sum(total_aa))))

    return gene_exons_domains_fractions, per_prot_dom, per_prot_dom_middle, changeinCONSTExonONMerger, per_dom_coverage, per_dom_dfrac_range_Full



def test(geneob):
    relev_exons=[]
    #print ([t.ID for t in geneob.transcripts])
    for t in geneob.transcripts:
        print (t.ID, [ex.ID for ex in t.exons])
    for ex in geneob.PI.exons:
        print (geneob.PI.ID,ex.ID, ex.length)
        key=ex.ID.split('.')[2] if ex.ID[0]!='R' else 'A'
        key='A' if key == 'F' else key
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
test(human[2051])

human={i:human[i] for i in human if i in filterg}
if additionalFilter:
    addfilter=loadpickle(additionalFilter)
    human={i:human[i] for i in human if i in addfilter}

fnameappend='super'
# dom_ex_annot, w_per_prot, m_per_prot, changeinCONSTExonONMerger, per_dom_coverage, per_dom_dfrac_range_Full =domain_exon_fraction(human=human)
# writeRangeFile(w_per_prot,results_dir+'F3_per_proteinDomain_covered_W_%s.csv'%fnameappend)
# writeRangeFile(m_per_prot,results_dir+'F3_per_proteinDomain_covered_M_%s.csv'%fnameappend)
# writeRangeFile(changeinCONSTExonONMerger,results_dir+'F3_exonsChangeFractionOnMerge_%s.csv'%fnameappend)
# writeRangeFile(per_dom_dfrac_range_Full, results_dir+'F3_domFractionRange80to100_%s.csv'%fnameappend)
# writeRangeFile(per_dom_coverage, results_dir+'F3_domCoveragePI_doms_%s.csv'%fnameappend)


# print (fnameappend)
# dumppickle(results_dir+'domain_exon_annotation_%s.pickle'%fnameappend, dom_ex_annot)



#python bin_analysis/domain_section.py dataset/objectsave_9606.pick results/ pfam
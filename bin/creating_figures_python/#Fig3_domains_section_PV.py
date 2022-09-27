
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
        key='A' if key == 'F' else key
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
    consecCodingrange=(consecCoding[0],consecCoding[-1]) if consecCoding else False
    for ex in exList:
        #[ex.length,key,ex.ID,1,0]
        #length,consFlag,ID, inclusion Flag, merged entry or not
        baseCondition=ex[0] # it should be coding
        #ATITrue=0, ATIFalse=1, reagrdless=-1
        #conditionATIATT=True if ATIFlag ==-1 else numericFlag in tExonList if ATIFlag==1 else numericFlag not in tExonList if ATIFlag==0 else False
        if ATIFlag==-1:
            choosenExonCondition=True
        elif ATIFlag==1:#core
            choosenExonCondition=consecCodingrange[0]<=int(ex[2].split(".")[3])<=consecCodingrange[-1] if consecCodingrange else False
        else:
            #ATIFlag==0
            choosenExonCondition=not (consecCodingrange[0]<=int(ex[2].split(".")[3])<=consecCodingrange[-1]) if consecCodingrange else True
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
    per_prot_dom={(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}
    per_prot_dom_middle={(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}
    changeinCONSTExonONMerger={(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}

    setdeepanshi=[7173]
    for gene in setdeepanshi:
    #for gene in human:
        tExonList=am.exonScreenerBetweenTConstitutive(human[gene].exons)
        exlis,domlis=isoform_giver(human[gene],CONST_val)
        #print ('tex',tExonList)
        
        #print (gene, bool(domlis), bool(exlis))
        if domlis and exlis:
            flag_exons=False
            flag_exons_conseCoding=False
            c+=1
            domcovered_w=0
            domcovered_m=0
            exlis=exon_screener(exlis, tExonList, ATIFlag) if domlis and exlis else False
            #print ('exlis',exlis)
            consecCoding=exon_screener_consecCoding(exlis, tExonList, ATIFlag, human[gene])
            #print ('consecCoding',consecCoding)
            
            proteinLength=sum([ex[0] for ex in exlis if ex[3]])
            exonsContrib_conseCoding={ex[2]:{'W':[ex[0],0,0],'M':[ex[0],0,0]} for ex in consecCoding}
            exonsContrib={ex[2]:{'W':[ex[0],0,0],'M':[ex[0],0,0]} for ex in exlis}

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
    if 1:
                formatter(gene,domain_namewith_cood,exlis,consecCoding,domainContrib, exonsContrib,domainContrib_conseCoding, exonsContrib_conseCoding)
        
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

human={i:human[i] for i in human if i in filterg}
if additionalFilter:
    addfilter=loadpickle(additionalFilter)
    human={i:human[i] for i in human if i in addfilter}

#test(human[5328])
setdeepanshi=[90113, 8193, 8202, 19, 22, 131096, 64175, 84656, 34, 41, 8239, 8242, 204851, 8458, 57409, 57410, 87, 8284, 95, 101, 114791, 131177, 57451, 109, 114803, 118, 8314, 8318, 57472, 57475, 114822, 57486, 57496, 161, 118813, 57520, 177, 57526, 284021, 204, 8226, 210, 221, 57576, 57578, 57579, 57582, 57583, 240, 8437, 57597, 57599, 57600, 8450, 57610, 729359, 8468, 285, 286, 57631, 57636, 57642, 114987, 57644, 57646, 57647, 8502, 57655, 57658, 317, 57662, 321, 8517, 8518, 333, 8526, 8539, 357, 8550, 57704, 8567, 8573, 57728, 10989, 405, 406, 409, 8603, 246175, 147872, 421, 8618, 8624, 8625, 8626, 8631, 65980, 65981, 65983, 139716, 65989, 460, 123355, 490, 492, 498, 66035, 147968, 115201, 57862, 8711, 147991, 8729, 545, 546, 547, 8751, 8764, 575, 576, 8786, 8792, 8808, 103, 8813, 23318, 352909, 8301, 23481, 64282, 8863, 672, 8867, 8880, 8887, 8894, 117, 8908, 8911, 8912, 8913, 57466, 90853, 123624, 745, 8938, 197358, 8943, 752, 762, 773, 775, 776, 778, 780, 8986, 800, 131873, 58155, 816, 818, 11059, 831, 9024, 840, 58190, 9044, 9051, 9053, 865, 9058, 115584, 898, 919, 923, 124401, 9156, 90273, 9162, 984, 987, 9181, 123872, 9185, 996, 1014, 9209, 1024, 9221, 9223, 1036, 115727, 1048, 132612, 9254, 9256, 1066, 222256, 9266, 9267, 9274, 1089, 9289, 9295, 1107, 1108, 1135, 1138, 26133, 9344, 58498, 58499, 1159, 9356, 1175, 124056, 91289, 1181, 9381, 25766, 91304, 1196, 140462, 25777, 9401, 25792, 91351, 165082, 25827, 25831, 9448, 4306, 9455, 1265, 9459, 9462, 25854, 140545, 9474, 1286, 1294, 9488, 1298, 9499, 9524, 25909, 50488, 79412, 25914, 91452, 222545, 9559, 9567, 9569, 9570, 1384, 1387, 26173, 9584, 25970, 1395, 116085, 25989, 9610, 9611, 25998, 116113, 9618, 9619, 26009, 1439, 1440, 9639, 1452, 26030, 1456, 26034, 9654, 9656, 1468, 9665, 26054, 222663, 26057, 26059, 50640, 5710, 9689, 9693, 9696, 1510, 26088, 9706, 26092, 9711, 26097, 1523, 9719, 9732, 9734, 26123, 91662, 84910, 26135, 9754, 9757, 9759, 26146, 1579, 9774, 5725, 9776, 149041, 8509, 9786, 116285, 11189, 9798, 1607, 1609, 2999, 9811, 9815, 83549, 1630, 91746, 1638, 1641, 91754, 9847, 26234, 1659, 9854, 9859, 9863, 9866, 9870, 9877, 26262, 26266, 9885, 9887, 9892, 9898, 9905, 91828, 83637, 124602, 9921, 9605, 9931, 1742, 1755, 1760, 9963, 83694, 83697, 9971, 9972, 1781, 222967, 83706, 83707, 1795, 9989, 9993, 10001, 10014, 10016, 8498, 1841, 51005, 51006, 1855, 124739, 51019, 7138, 1871, 1877, 10076, 124773, 1894, 57661, 51059, 51061, 1911, 51070, 51078, 8513, 51085, 51093, 143684, 26523, 26528, 51105, 83875, 10149, 10152, 59307, 149428, 10165, 1981, 51147, 8525, 334, 51160, 26586, 2011, 2015, 83939, 2036, 2037, 10232, 10236, 51203, 10250, 10256, 83985, 10260, 2070, 2071, 51225, 83999, 92211, 2104, 387129, 10300, 9909, 10307, 10312, 2122, 51276, 157777, 2130, 51291, 26298, 51295, 2146, 84067, 10347, 84077, 51310, 84083, 51317, 51322, 149628, 2178, 125061, 2185, 26762, 157855, 2210, 84131, 2212, 51366, 10411, 10413, 51377, 2232, 10426, 116931, 84164, 84166, 51399, 84172, 2253, 10450, 84179, 2261, 2264, 10457, 10458, 79567, 10460, 51421, 10462, 51430, 51433, 10492, 10495, 84224, 84225, 2314, 10509, 51474, 125206, 51479, 2329, 84250, 84251, 2332, 387357, 2334, 84256, 10531, 4486, 84271, 51517, 84287, 10565, 256329, 7223, 10580, 51542, 10585, 51548, 10603, 256364, 51569, 10625, 10630, 51593, 27019, 84364, 51599, 10642, 55023, 27037, 27039, 27040, 117154, 10659, 51629, 10675, 1780, 117178, 51643, 134218, 51654, 51657, 2515, 256471, 10712, 10713, 1786, 10721, 149986, 1787, 51684, 10725, 51692, 10733, 10734, 84465, 84467, 27124, 2553, 51720, 27147, 27148, 10771, 2580, 27158, 27165, 51742, 84515, 51750, 158248, 84532, 199223, 2635, 2648, 2649, 27229, 27231, 10858, 494188, 10869, 27255, 2687, 84626, 10900, 84629, 10902, 84632, 10908, 10910, 10919, 10926, 2736, 27324, 84669, 27327, 84675, 10948, 84678, 256714, 10956, 84691, 166614, 27352, 84699, 84700, 84705, 10988, 84717, 10990, 10992, 84722, 150274, 2821, 338699, 84749, 11043, 27429, 11047, 11055, 27443, 81033, 2873, 84804, 2889, 2890, 84818, 2905, 11107, 11113, 93034, 11118, 11119, 11129, 11131, 84867, 84871, 11148, 84879, 11168, 11169, 84899, 2983, 11176, 2990, 60343, 11194, 11196, 11202, 23014, 84941, 84951, 3032, 84954, 3035, 84958, 5970, 3054, 3055, 84978, 3067, 3070, 11273, 3099, 85021, 199713, 203054, 166968, 11338, 3161, 3177, 3187, 3189, 535, 93349, 126123, 60592, 126129, 60626, 3292, 388325, 23762, 167153, 60680, 3340, 153478, 200010, 200014, 142678, 56034, 126306, 3430, 3480, 93594, 339366, 85417, 93627, 85444, 85445, 85451, 85459, 201294, 3551, 126432, 3594, 3603, 3609, 3614, 3619, 3636, 3643, 126526, 3654, 3655, 642636, 3663, 3679, 3683, 3708, 54955, 3783, 3786, 3796, 3801, 151258, 118491, 3815, 143098, 10197, 93986, 255967, 10211, 100527963, 10497, 348013, 94081, 3978, 94097, 3996, 4000, 4008, 7499, 4041, 4052, 4054, 4057, 4058, 2044, 339965, 645121, 2051, 348180, 127002, 4123, 4125, 4133, 4137, 4140, 4147, 4152, 135228, 4162, 4185, 4205, 4208, 4215, 4216, 4217, 4233, 53407, 4261, 4293, 4301, 159963, 4327, 200942, 4820, 4356, 4361, 151827, 127281, 28988, 28996, 151888, 10297, 4440, 23607, 282973, 29028, 53616, 4485, 29062, 225689, 29085, 340390, 29102, 29117, 340419, 4582, 4595, 4602, 4603, 4628, 148229, 4648, 4649, 4659, 127544, 10335, 4669, 53832, 4686, 4715, 4753, 4756, 53919, 4779, 4781, 4790, 4795, 4796, 4798, 2165, 4802, 389840, 152273, 135892, 4832, 144100, 4839, 4849, 4850, 4861, 85461, 4864, 6272, 4867, 152330, 4882, 57447, 4897, 4915, 2186, 4927, 4928, 4929, 4931, 4938, 54103, 127833, 4978, 283455, 4998, 102724488, 201627, 5023, 5024, 5026, 283578, 127943, 5094, 5129, 136227, 5160, 5168, 29760, 5210, 29800, 5244, 5256, 5264, 29841, 78997, 5270, 5273, 5275, 23635, 5283, 54436, 54437, 54439, 54443, 5293, 1540, 144568, 54457, 54461, 29888, 5315, 29893, 414918, 29896, 54482, 29911, 5338, 79068, 29924, 5358, 128239, 79090, 201973, 54520, 29945, 84181, 29952, 29954, 54531, 29959, 54536, 54537, 54538, 91010, 5394, 29974, 79133, 29988, 29991, 29994, 7732, 5434, 5445, 144715, 5452, 79184, 54621, 5475, 284013, 5493, 5495, 728642, 5528, 79258, 152992, 202151, 284086, 153020, 5571, 51196, 284111, 5588, 5599, 284129, 54760, 5609, 54776, 54778, 80128, 153090, 2305, 54799, 54802, 54808, 54815, 54822, 54823, 54836, 54840, 54847, 54862, 2317, 5720, 54877, 65125, 54880, 54882, 2322, 5742, 5745, 54903, 5754, 5757, 374403, 5781, 5783, 153241, 161436, 5789, 5794, 5798, 10524, 5803, 79573, 284339, 54970, 54971, 54972, 54974, 5826, 5830, 5832, 5833, 54991, 79574, 5858, 79595, 79598, 5871, 284403, 79607, 91851, 79627, 5900, 390928, 5909, 5910, 55063, 100652824, 5913, 5915, 79649, 5926, 5934, 5937, 55605, 55107, 55108, 993, 284498, 5976, 5981, 5985, 55139, 5989, 5990, 79722, 55616, 55173, 55183, 55184, 55619, 55193, 55198, 55200, 79781, 55206, 79783, 55210, 79796, 6091, 6092, 22861, 6098, 79828, 79829, 79830, 79867, 79869, 79872, 79883, 79885, 22872, 79892, 6490, 79915, 6197, 55367, 79960, 30817, 30819, 30827, 79980, 6257, 129138, 30835, 79990, 80000, 80006, 161931, 374928, 137362, 80021, 6294, 6305, 26993, 55471, 6323, 6328, 6331, 6335, 6336, 80067, 6344, 55511, 317662, 1062, 6386, 80115, 112885, 375033, 22911, 6400, 55553, 55556, 55558, 51585, 55561, 653583, 284958, 6432, 22818, 55591, 55593, 6442, 55687, 80174, 22837, 22838, 6455, 6464, 317761, 55691, 22852, 6472, 80204, 80205, 22864, 55634, 80216, 22873, 22874, 22880, 6497, 80227, 80232, 55657, 55660, 22893, 55666, 22899, 63893, 6528, 88455, 55689, 80267, 22925, 22926, 55695, 22931, 55700, 22933, 55704, 22937, 440730, 22941, 6559, 6560, 55714, 55718, 55726, 55728, 55729, 63925, 55738, 55740, 203197, 55743, 22976, 6596, 6597, 6598, 6601, 55755, 80332, 10658, 22990, 22998, 154075, 55777, 63976, 6640, 23025, 23028, 23031, 23032, 23033, 23037, 55806, 23039, 23040, 55809, 55812, 55814, 55815, 23048, 23054, 6672, 203286, 23065, 23066, 162333, 23072, 6689, 6692, 399909, 80305, 23082, 55860, 64062, 64066, 64067, 23113, 6733, 23120, 23122, 23125, 23126, 23130, 64093, 64094, 23136, 23140, 55914, 23155, 6773, 64118, 55930, 23164, 23175, 23176, 23181, 285331, 23189, 64151, 22980, 6810, 6813, 23205, 64174, 23215, 6833, 203447, 23240, 23241, 6861, 64210, 23252, 23254, 23255, 6872, 23258, 6875, 23262, 23266, 23271, 23272, 23276, 64241, 137970, 51667, 23284, 23286, 23294, 6504, 23301, 23303, 23305, 23306, 23307, 6925, 6928, 6929, 6934, 6938, 219931, 23325, 146206, 23328, 121643, 23350, 23355, 23358, 23367, 23371, 64333, 23380, 23383, 56155, 23389, 7006, 220001, 7010, 56163, 23397, 23400, 23404, 7027, 7029, 23418, 220032, 100129669, 7046, 7059, 285596, 7073, 55793, 23464, 64427, 7088, 7089, 7090, 56243, 55796, 23492, 23500, 23505, 23518, 375775, 56288, 23522, 7139, 7145, 56301, 55805, 23542, 23543, 23544, 64506, 23549, 23550, 23552, 23555, 23041, 654346, 171019, 1198, 154661, 252983, 158219, 23620, 23623, 23624, 23626, 7248, 7249, 203859, 130132, 64599, 23647, 23660, 146547, 7297, 728194, 7318, 80169, 56478, 29894, 7337, 64682, 162989, 10781, 79048, 23731, 7378, 10787, 23774, 7392, 23787, 7404, 79058, 7408, 64757, 64760, 64770, 7433, 64784, 146713, 7455, 146722, 7481, 7490, 5345, 56649, 253260, 64849, 64854, 64857, 7536, 146802, 343413, 7555, 64900, 8086, 138639, 286097, 121274, 64919, 79087, 64926, 155061, 286148, 130507, 7629, 146894, 54521, 94120, 65010, 286204, 146956, 10840, 7702, 65056, 7716, 23133, 56884, 7746, 7748, 56904, 24140, 81490, 56916, 56917, 7772, 65117, 7781, 56945, 56946, 56947, 56948, 56957, 56965, 81545, 56970, 23149, 56980, 56987, 7849, 81578, 7862, 57019, 7869, 89790, 646851, 6774, 81607, 253650, 23162, 401124, 6780, 7917, 65267, 65268, 89848, 8148, 92799, 122622, 7936, 57096, 253738, 7994, 57154, 89941, 57178, 139105, 57187, 401258, 57205, 57211, 8061, 81794, 57223, 253832, 55959, 57231, 221074, 221078, 81831, 23457, 286676]

#test(human[setdeepanshi[0]])

#sys.exit()

dom_ex_annot, w_per_prot, m_per_prot, changeinCONSTExonONMerger=domain_exon_fraction(human=human)
# writeRangeFile(w_per_prot,results_dir+'F3_per_proteinDomain_covered_W_%s.csv'%fnameappend)
# writeRangeFile(m_per_prot,results_dir+'F3_per_proteinDomain_covered_M_%s.csv'%fnameappend)
# writeRangeFile(changeinCONSTExonONMerger,results_dir+'F3_exonsChangeFractionOnMerge_%s.csv'%fnameappend)

# print (fnameappend)
# dumppickle(results_dir+'domain_exon_annotation_%s.pickle'%fnameappend, dom_ex_annot)



#python bin_analysis/domain_section.py dataset/objectsave_9606.pick results/ pfam
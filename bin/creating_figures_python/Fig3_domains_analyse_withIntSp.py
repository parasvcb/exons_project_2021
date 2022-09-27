import sys,os
import pandas as pd
import cPickle as pickle
from collections import Counter
if len(sys.argv)!=3:
    print ('Please enter correct args 1. geneExondata 2. outputdir')
    sys.exit()
prog,geneEx,outdir=sys.argv

fappend=geneEx.split('_')[-1].split('.')[0]
outdir=outdir+fappend+'_'
def div_fact(num,denom):
	try:
		return round(float(num)/denom,3)	
	except:
		return 0

def has_range_adder(has,val):
    val=round(val,2)
    for rangeV in has:
        if rangeV[0]<=val<=rangeV[1]:
            has[rangeV]+=1
            break
    return has

def has_range_adder_val(has,val, ex):
    val=round(val,2)
    for rangeV in has:
        if rangeV[0]<=val<=rangeV[1]:
            has[rangeV]+=ex
            break
    return has

'''
gene_exons_domains_fractions[gene]={'domains':[],'exons':{},'domains_conseCoding':[],'exons_conseCoding':{}}
    gene_exons_domains_fractions[gene]['domains']=[[domainsName,domainContrib]]
    gene_exons_domains_fractions[gene]['domains_conseCoding']=....
    
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

    domainContrib['exContrib']=[dnamewithcood,dfrac,dfracm, exonCore]

    gene_exons_domains_fractions[gene]['exons']=exonsContrib
    gene_exons_domains_fractions[gene]['exons_conseCoding']=exonsContrib_conseCoding
    
    # exonsContrib[ex]['W'][2]+=interori
    # exonsContrib[ex]['M'][2]+=intermid
    # exonsContrib[ex]['W'][3]=round(float(exonsContrib[ex]['W'][1]/exonsContrib[ex]['W'][0]),2)
    # 0th ele is G/A key, 1st is exonCore, second element is domain span and third is complete fractio coding for such regions

L0: domain presence conundrum, [only domaincontrib fraction]
L1: How many domains undergo split and how many are contained within the exons (X% of D undergo split and Y% is contained within)
L2: For the domains undergoing split (X%), what is the fraction of the per exon contribution [barplot of 0-0.10,,,,, 0.9-1.0, lower panel their exon types]
L3: How many split instances (Gx% of X) are split within the GG exons, and how many are split within other exon types (AGx% and AAx%)
 For the AGx%, draw fraction contained within Gs (cumulative frequency distribution)
Analyse with ATI/ATT and core, repeat with middle 40% of region also

'''
def addvaltoHaslis(has,key,val):
    if key not in has:
        has[key]=[]
    if val not in has[key]:
        has[key]+=[val]
    return has


def summarise(ob):
    domKeys=['domains','domains_conseCoding']
    for keys in domKeys:
        print ('key:',keys)
        for domains in ob[keys]:
            print ('\n',domains,ob[keys][domains][:3])
            for exons in ob[keys][domains][4]:
                print (exons,ob[keys][domains][4][exons])
        print ('\n')
def l1(has, fout, category,WorM):
    #category={'domains_conseCoding':'ConsecCoding', 'domains':'normal'}
    #one of those but in list
    #WorM is Whole or middle
    completeTag='containedW' if WorM=='W' else 'containedM'
    splitTag='splitW' if WorM=='W' else 'splitM'
    total={}
    contained={'2ATIT':{'A':{},'G':{}},'3BOTH':{'A':{},'G':{}},'1CORE':{'A':{},'G':{}}}
    
    split={
        '2ATIT':{'A':{},'G':{},'AG':{}},
        '3BOTH':{'A':{},'G':{},'AG':{}},
        '1CORE':{'A':{},'G':{},'AG':{}},
        '1CORE+2ATIT':{'1CORE_A__2ATIT_A':{},'1CORE_A__2ATIT_G':{},'1CORE_G__2ATIT_G':{},'1CORE_G__2ATIT_A':{},'1CORE_AG__2ATIT_A':{},'1CORE_AG__2ATIT_G':{},'1CORE_AG__2ATIT_AG':{},'1CORE_A__2ATIT_AG':{},'1CORE_G__2ATIT_AG':{}}
        }
    
    for gene in has:
      #if gene==405:
            for domain in has[gene][category[0]]:
                domname, domainContrib=domain
                # domainContrib={
                # 'considered':'',
                # #in order of CORE,ATIT,BOTH
                # 'containedM':['','',''], || 'containedW':['','',''],
                # 'exContrib':{}, #exon wise list of list
                #                   [domain_namewith_cood, dfrac, dfrac_m, efrac_m, efrac, espan, interori, intermidd, coreExon]
                # 'splitM':[{'A':[0,''],'G':[0,'']},{'A':[0,''],'G':[0,'']},{'A':[0,''],'G':[0,'']}], || 'splitW' # }
                # + in the end [categoryTupleATITCORE, theActualTuple, fractionsThere, categoryTupleBOTH, theActualTuple, fractionsThere]
          
                localCateg=category[1]
                # will give value normal/consecCoding
                
                total=addvaltoHaslis(total,gene,domname)
                comp=[i for i in domainContrib[completeTag] if i]
                #print (gene,domname)
                #print (domainContrib)
                if comp:
                    key='1CORE' if domainContrib[completeTag][0] else '2ATIT'
                    junc=domainContrib[completeTag][0] if domainContrib[completeTag][0] else domainContrib[completeTag][1]
                    #print ('Yes',gene,domname, domainContrib[completeTag],comp,key)
                    contained[key][junc]=addvaltoHaslis(contained[key][junc],gene,domname)
                    contained['3BOTH'][domainContrib[completeTag][2]]=addvaltoHaslis(contained['3BOTH'][domainContrib[completeTag][2]],gene,domname)
                else:
                    CORE,ATIT,BOTH,COREATIT_cat,COREATIT_realTup,COREATIT_frac,BOTH_cat,BOTH_realTup,BOTH_frac=domainContrib[splitTag]
                    #splitW': [{'A': [0, 0.0], 'G': [91, 0.711]}, {'A': [37, 0.289], 'G': [0, 0.0]}, {'A': [37, 0.289], 'G': [91, 0.711]}, '1CORE_G__2ATIT_A', [1, 'G', 2, 'A'], [0.711, 0.289], '3BOTHAG', [3, 'A', 3, 'G'], [0.289, 0.711]]
                    #print(BOTH_cat)
                    key=BOTH_cat[:5]
                    junc=BOTH_cat[5:]
                    #print (key,junc,'split')
                    split[key][junc]=addvaltoHaslis(split[key][junc],gene,domname)
                    if key=='2ATIT':
                        print (domainContrib)
                        sys.exit()

                    if len(COREATIT_cat) in [6,7]:
                        key=COREATIT_cat[:5]
                        junc=COREATIT_cat[5:]
                        split[key][junc]=addvaltoHaslis(split[key][junc],gene,domname)
                    else:
                        key='1CORE+2ATIT'
                        junc=COREATIT_cat
                        split[key][junc]=addvaltoHaslis(split[key][junc],gene,domname)

                    '''
                    '1COREA','1COREAG','1COREG','2ATITA','2ATITAG','2ATITG','1CORE_A__2ATIT_A','1CORE_A__2ATIT_G','1CORE_G__2ATIT_G','1CORE_G__2ATIT_A','1CORE_AG__2ATIT_A','1CORE_AG__2ATIT_G','1CORE_AG__2ATIT_AG','1CORE_A__2ATIT_AG','1CORE_G__2ATIT_AG'
                    '''
    #print (contained['A'])
    
    #print (split['1CORE']['G'])
    fout.write('Category: %s, domains:%s\n'%(category,WorM))
    
    fout.write('Broadcat\texScope\tTotalDomain_Core\tUniqueDom_Core\tGenes_Core\tTotalDomain_ATIT\tUniqueDom_ATIT\tGenes_ATIT\tTotalDomain_Both\tUniqueDom_Both\tGenes_Both\tTotalDomain_ATIT+Core\tUniqueDom_ATIT+Core\tGenes_ATIT+Core\n')
    genes=total.keys()
    domains=[i for j in total.values() for i in j]
    print (domains[:5])
    uniqueDom=set(domains)
    fout.write('Total\t.\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\t.\t.\t.\n'%(len(domains), len(uniqueDom), len(genes)))
    
    vals3=['1CORE','2ATIT','3BOTH']
    for exJunc in ['A','G']:
        string=''
        for ind,ttype in enumerate(vals3):   
            genes=contained[ttype][exJunc].keys()
            #print (len(genes))
            domains=[i for j in contained[ttype][exJunc].values() for i in j]
            uniqueDom=set(domains)
            string+='\t%s\t%s\t%s'%(len(domains),len(uniqueDom),len(genes))
            #print (string,exJunc,ttype)
        fout.write('Contained\t%s%s\t.\t.\t.\n'%(exJunc,string))
    
    for exJunc in ['A','G','AG']:
        string=''
        for ind,ttype in enumerate(vals3):    
            genes=split[ttype][exJunc].keys()
            domains=[i for j in split[ttype][exJunc].values() for i in j]
            uniqueDom=set(domains)
            string+='\t%s\t%s\t%s'%(len(domains),len(uniqueDom),len(genes))
        fout.write('Split\t%s%s\t.\t.\t.\n'%(exJunc,string))
    
    

    for exJunc in ['1CORE_A__2ATIT_A','1CORE_A__2ATIT_G','1CORE_G__2ATIT_G','1CORE_G__2ATIT_A','1CORE_AG__2ATIT_A','1CORE_AG__2ATIT_G','1CORE_AG__2ATIT_AG','1CORE_A__2ATIT_AG','1CORE_G__2ATIT_AG']:
        string=''
        for ind,ttype in enumerate(['1CORE+2ATIT']):
            genes=split[ttype][exJunc].keys()
            domains=[i for j in split[ttype][exJunc].values() for i in j]
            uniqueDom=set(domains)
            string+='\t%s\t%s\t%s'%(len(domains),len(uniqueDom),len(genes))
        fout.write('Split\t%s\t.\t.\t.\t.\t.\t.\t.\t.\t%s\n'%(exJunc,string))

    fout.write('\n\n')
    return fout



def writetoFileRangesandStrings(hasVals,hasStrings,fname):
    #print (hasStrings)
    occRangeList=list(hasVals.keys())
    occRangeList.sort()
    totalOcc=sum(hasVals.values())
    
    with open (fname,'w') as fin:
        fin.write('Range\tTotal\tFreq\tTag\n')
        #previousValue=0
        
        for occRange in occRangeList:
            occrangestr="_".join(map(str,list(occRange)))
            count_occrange=hasVals[occRange]
            #cumulative_count=previousValue+count_occrange
            #previousValue+=count_occrange
            
            totalStr=hasStrings[occRange]
            lenStr=len(totalStr)
            GCount=totalStr.count('G')
            GFreq=div_fact(GCount,lenStr)
            fin.write('%s\t%s\t%s\t%s\n'%(occrangestr,count_occrange,div_fact(count_occrange,totalOcc),'exonsContribution'))
            fin.write('%s\t%s\t%s\t%s\n'%(occrangestr,GCount,GFreq,'GfracColumn'))
    return

def writetoFileRanges(has,fname):
    occRangeList=list(has.keys())
    occRangeList.sort()
    totalOcc=sum(has.values())
    
    with open (fname,'w') as fin:
        fin.write('Range\tTotal\tFreq\tTag\n')
        #previousValue=0
        
        for occRange in occRangeList:
            occrangestr="_".join(map(str,list(occRange)))
            count_occrange=has[occRange]
            #cumulative_count=previousValue+count_occrange
            #previousValue+=count_occrange
            
            fin.write('%s\t%s\t%s\t%s\n'%(occrangestr,count_occrange,div_fact(count_occrange,totalOcc),'exonsContribution'))
    return

def l2(has,fnames):
    '''
    L2: For the domains undergoing split (X%), what is the fraction of the per exon contribution [barplot of 0-0.10,,,,, 0.9-1.0, lower panel their exon types]
    L3: How many split instances (Gx% of X) are split within the GG exons, and how many are split within other exon types (AGx% and AAx%)
    For the AGx%, draw fraction contained within Gs (cumulative frequency distribution)
    Analyse with ATI/ATT and core, repeat with middle 40% of region also
    '''
    category={'domains_conseCoding':'ConsecCoding', 'domains':'normal'}
    

    domainPerspective={'normal':{'W':{(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0},'M':{(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}},'ConsecCoding':{'W':{(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0},'M':{(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}}}

    domainPerspective_EX={'normal':{'W':{(0.0,0.01):'',(0.01,0.1):'',(0.101,0.2):'',(0.201,0.3):'',(0.301,0.4):'',(0.401,0.5):'',(0.501,0.6):'',(0.601,0.7):'',(0.801,0.9):'',(0.901,0.990):'', (0.991,1.01):''},'M':{(0.0,0.01):'',(0.01,0.1):'',(0.101,0.2):'',(0.201,0.3):'',(0.301,0.4):'',(0.401,0.5):'',(0.501,0.6):'',(0.601,0.7):'',(0.801,0.9):'',(0.901,0.990):'', (0.991,1.01):''}},'ConsecCoding':{'W':{(0.0,0.01):'',(0.01,0.1):'',(0.101,0.2):'',(0.201,0.3):'',(0.301,0.4):'',(0.401,0.5):'',(0.501,0.6):'',(0.601,0.7):'',(0.801,0.9):'',(0.901,0.990):'', (0.991,1.01):''},'M':{(0.0,0.01):'',(0.01,0.1):'',(0.101,0.2):'',(0.201,0.3):'',(0.301,0.4):'',(0.401,0.5):'',(0.501,0.6):'',(0.601,0.7):'',(0.801,0.9):'',(0.901,0.990):'', (0.991,1.01):''}}}


    exonsPerspective={'normal':{'W':{(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0},'M':{(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}},'ConsecCoding':{'W':{(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0},'M':{(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}}}

    exonsPerspective_EX={'normal':{'W':{(0.0,0.01):'',(0.01,0.1):'',(0.101,0.2):'',(0.201,0.3):'',(0.301,0.4):'',(0.401,0.5):'',(0.501,0.6):'',(0.601,0.7):'',(0.801,0.9):'',(0.901,0.990):'', (0.991,1.01):''},'M':{(0.0,0.01):'',(0.01,0.1):'',(0.101,0.2):'',(0.201,0.3):'',(0.301,0.4):'',(0.401,0.5):'',(0.501,0.6):'',(0.601,0.7):'',(0.801,0.9):'',(0.901,0.990):'', (0.991,1.01):''}},'ConsecCoding':{'W':{(0.0,0.01):'',(0.01,0.1):'',(0.101,0.2):'',(0.201,0.3):'',(0.301,0.4):'',(0.401,0.5):'',(0.501,0.6):'',(0.601,0.7):'',(0.801,0.9):'',(0.901,0.990):'', (0.991,1.01):''},'M':{(0.0,0.01):'',(0.01,0.1):'',(0.101,0.2):'',(0.201,0.3):'',(0.301,0.4):'',(0.401,0.5):'',(0.501,0.6):'',(0.601,0.7):'',(0.801,0.9):'',(0.901,0.990):'', (0.991,1.01):''}}}
    
    #this is domain perspective
    
    for gene in has:
        #print (has[gene].keys())
        #sys.exit()
        for exonsCateg in category:
            localCateg=category[exonsCateg]
            for domain in has[gene][exonsCateg]:
                W_contri,M_contri,contained40, contained100 ,domainContrib=has[gene][exonsCateg][domain]
                #[(keyw, key1contrib_W, key2contrib_W),(keym, key1contrib_M, key2contrib_M) ,contained40, contained100 ,domainContrib]
                    
                for ex in domainContrib:
                    keyID=ex.split('.')[2] if ex[0]!='R' else 'A'
                    keyID='G' if keyID == 'F' else keyID
                    dname,dfrac,dfracm=domainContrib[ex]
                    if not contained100:
                        domainPerspective[localCateg]['W']=has_range_adder(domainPerspective[localCateg]['W'],dfrac)
                        domainPerspective_EX[localCateg]['W']=has_range_adder_val(domainPerspective_EX[localCateg]['W'],dfrac, keyID)
                    if not contained40:
                        domainPerspective[localCateg]['M']=has_range_adder(domainPerspective[localCateg]['M'],dfracm)
                        domainPerspective_EX[localCateg]['M']=has_range_adder_val(domainPerspective_EX[localCateg]['M'],dfracm, keyID)
    

    category2={'exons_conseCoding':'ConsecCoding', 'exons':'normal'}
    # #exonsWise
    # gene_exons_domains_fractions[gene]['exons']=exonsContrib
    # gene_exons_domains_fractions[gene]['exons_conseCoding']=exonsContrib_conseCoding
    # exonsContrib[ex]['W'][1]+=interori
    # exonsContrib[ex]['M'][1]+=intermid
    # exonsContrib[ex]['W'][2]=round(float(exonsContrib[ex]['W'][1]/exonsContrib[ex]['W'][0]),2)
    
    #print (exonsPerspective)
    for gene in has:
        # print (gene)
        # print (has[gene].keys())
        # print (has[gene])
        # sys.exit()
        for exonsCateg in category2:
            localCateg=category2[exonsCateg]
            #print (localCateg,exonsCateg,1)
            exonHhas=has[gene][exonsCateg]
            #print (gene,exonHhas)
            for ex in exonHhas:
                keyID=ex.split('.')[2] if ex[0]!='R' else 'A'
                keyID='G' if keyID == 'F' else keyID
                for localGlobalScope in exonHhas[ex]:
                    lengthex,span,frac=exonHhas[ex][localGlobalScope]
                    #print (localCateg,exonsCateg,frac)
                    # if localCateg=='normal':
                    #     print ('yes')
                    #     print (frac)
                    # a condition cab be imposed heer to check if underlying domain is fully conatined or not, lets not worry about it as of now
                    exonsPerspective[localCateg][localGlobalScope]=has_range_adder(exonsPerspective[localCateg][localGlobalScope],frac)
                    exonsPerspective_EX[localCateg][localGlobalScope]=has_range_adder_val(exonsPerspective_EX[localCateg][localGlobalScope],frac, keyID)
    #print (exonsPerspective)

    for exonsCateg in domainPerspective:
        for domaincateg in domainPerspective[exonsCateg]:
            writetoFileRangesandStrings(domainPerspective[exonsCateg][domaincateg], domainPerspective_EX[exonsCateg][domaincateg] ,fnames+'domainsPerspective_%s_%s.csv'%(exonsCateg,domaincateg))
    
    for exonsCateg in exonsPerspective_EX:
        for domaincateg in exonsPerspective_EX[exonsCateg]:
            writetoFileRangesandStrings(exonsPerspective[exonsCateg][domaincateg], exonsPerspective_EX[exonsCateg][domaincateg] ,fnames+'exonsPerspective_%s_%s.csv'%(exonsCateg,domaincateg))

   


def l3(has,fnames, foutLog):
    category={'domains_conseCoding':'ConsecCoding', 'domains':'normal'}

    genFrac={'normal':{'W':{'AG':0,'GG':0,'AA':0},'M':{'AG':0,'GG':0,'AA':0}},'ConsecCoding':{'W':{'AG':0,'GG':0,'AA':0},'M':{'AG':0,'GG':0,'AA':0}}}


    AGFrac={'normal':{'W':{(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0},'M':{(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}},'ConsecCoding':{'W':{(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0},'M':{(0.0,0.01):0,(0.01,0.1):0,(0.101,0.2):0,(0.201,0.3):0,(0.301,0.4):0,(0.401,0.5):0,(0.501,0.6):0,(0.601,0.7):0,(0.801,0.9):0,(0.901,0.990):0, (0.991,1.01):0}}}


    for gene in has:
        for exonsCateg in category:
            localCateg=category[exonsCateg]
            for domain in has[gene][exonsCateg]:
                W_contri,M_contri,contained40, contained100 ,domainContrib=has[gene][exonsCateg][domain]
                #[(keyw, key1contrib_W, key2contrib_W),(keym, key1contrib_M, key2contrib_M) ,contained40, contained100 ,domainContrib]
                keym, key1contrib_M, key2contrib_M=M_contri
                keyw, key1contrib_W, key2contrib_W=W_contri
                if not contained100 and keyw:
                    genFrac[localCateg]['W'][keyw]+=1
                    if keyw=='AG':
                        #AGFrac[localCateg]['W']+=1
                        AGFrac[localCateg]['W']=has_range_adder(AGFrac[localCateg]['W'],key2contrib_W)
                
                if not contained40 and keym:
                    genFrac[localCateg]['M'][keym]+=1
                    if keym=='AG':
                        #AGFrac[localCateg]['M']+=1
                        #print (keym, key1contrib_M, key2contrib_M)
                        AGFrac[localCateg]['M']=has_range_adder(AGFrac[localCateg]['M'],key2contrib_W)
                

    junc=['AA','AG','GG']
    foutLog.write('\n**section about split cases and totality of them\n')
    foutLog.write('localGlobalScope\tdomIntegrity\texonJunction\tcount\tfreq\n')
    
    for localGlobalScope in genFrac:
        #foutLog.write('\nfor the %s exonmerger\n'%(localGlobalScope))
        for domIntegrity in genFrac[localGlobalScope]:
            totalDom=sum(genFrac[localGlobalScope][domIntegrity].values())
            #foutLog.write('\tthe %s domainType, total split cases are %s\n'%(domIntegrity,totalDom))
            for exJunc in junc:
                count=genFrac[localGlobalScope][domIntegrity][exJunc]
                #foutLog.write('\t\tfor junctype %s total are %s (%s) \n'%(exJunc,count, div_fact(count,totalDom)))
                foutLog.write('%s\t%s\t%s\t%s\t%s\n'%(localGlobalScope,domIntegrity,exJunc,count,div_fact(count,totalDom)))
    
    
    
    for exonsCateg in AGFrac:
        for domaincateg in AGFrac[exonsCateg]:
            writetoFileRanges(AGFrac[exonsCateg][domaincateg],fnames+'_%s_%s.csv'%(exonsCateg,domaincateg))

def readpickle(fname):
    with open (fname,'rb') as fin:
        dat=pickle.load(fin)
    return dat


geneEx=readpickle(geneEx)

logfile=outdir+'domainsSectionStatsLog'
logf=open (logfile,'w')


print (geneEx[2051])
#l1(geneEx,logf, category=['domains','normal'], WorM='W')
#l1(geneEx,logf, category=['domains_conseCoding','ConsecCoding'], WorM='W')


# l2(geneEx,outdir+'perFractionDomains_')
# l3(geneEx,outdir+'FractionAG_', logf )
logf.close()

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
gene_exons_domains_fractions[gene]={'domains':{},'exons':{},'domains_conseCoding':{},'exons_conseCoding':{}}
    gene_exons_domains_fractions[gene]['domains'][domain_namewith_cood]=[(keyw, key1contrib_W, key2contrib_W),(keym, key1contrib_M, key2contrib_M) ,contained40, contained100 ,domainContrib]
    gene_exons_domains_fractions[gene]['domains_conseCoding']=....
    #domainContrib[ex]=[dnamewithcood,dfrac,dfracm]


    gene_exons_domains_fractions[gene]['exons']=exonsContrib
    gene_exons_domains_fractions[gene]['exons_conseCoding']=exonsContrib_conseCoding
    
    # exonsContrib[ex]['W'][1]+=interori
    # exonsContrib[ex]['M'][1]+=intermid
    # exonsContrib[ex]['W'][2]=round(float(exonsContrib[ex]['W'][1]/exonsContrib[ex]['W'][0]),2)
    # first element is domain span and secod is complete fractio coding for such regions


gene_exons_domains_fractions[gene][domain_namewith_cood]={100:(keyw, key1contrib_W, key2contrib_W),40:(keym, key1contrib_M, key2contrib_M),'containedM':contained40,'containedW':contained100 ,'contrib':exContrib}

exContrib[ex.ID]+=[(domain_namewith_cood, dfrac, dfrac_m, efrac_m, efrac)]

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

def l1(has, fout):
    total={'normal':{},'ConsecCoding':{}}
    contained={'normal':{'W':{},'M':{}},'ConsecCoding':{'W':{},'M':{}}}
    split={'normal':{'W':{},'M':{}},'ConsecCoding':{'W':{},'M':{}}}
    category={'domains_conseCoding':'ConsecCoding', 'domains':'normal'}

    for gene in has:
        for exonsCateg in category:
            for domain in has[gene][exonsCateg]:
                W_contri,M_contri,contained40, contained100 ,domainContrib=has[gene][exonsCateg][domain]
                #[(keyw, key1contrib_W, key2contrib_W),(keym, key1contrib_M, key2contrib_M) ,contained40, contained100 ,domainContrib]

                localCateg=category[exonsCateg]
                total[localCateg]=addvaltoHaslis(total[localCateg],gene,domain)
                if contained100:
                    contained[localCateg]['W']=addvaltoHaslis(contained[localCateg]['W'],gene,domain)
                else:
                    split[localCateg]['W']=addvaltoHaslis(split[localCateg]['W'],gene,domain)
                if contained40:
                    contained[localCateg]['M']=addvaltoHaslis(contained[localCateg]['M'],gene,domain)
                else:
                    split[localCateg]['M']=addvaltoHaslis(split[localCateg]['M'],gene,domain)
    
    localCateg=category.values()
    fout.write('*Total domain stats and contained cases summary\n')
    for cat in total:
        genes=total[cat].keys()
        domains=[i for j in total[cat].values() for i in j]
        uniqueDom=set(domains)
        fout.write('\tTotal Domains for categ %s is %s (%s unique), gene count is %s\n  '%(cat,len(domains), len(uniqueDom), len(genes)))
    
    fout.write('\n* For contained cases')
    for cat in contained:
        for localGlobalScope in contained[cat]:
            genes=contained[cat][localGlobalScope].keys()
            domains=[i for j in contained[cat][localGlobalScope].values() for i in j]
            uniqueDom=set(domains)
            fout.write('\tContained Domains for categ %s and scope %s is %s (%s unique), gene count is %s\n  '%(cat,localGlobalScope,len(domains), len(uniqueDom), len(genes)))

    fout.write('\n* For split cases')
    for cat in split:
        for localGlobalScope in split[cat]:
            genes=split[cat][localGlobalScope].keys()
            domains=[i for j in split[cat][localGlobalScope].values() for i in j]
            uniqueDom=set(domains)
            fout.write('\tsplit Domains for categ %s and scope %s is %s (%s unique), gene count is %s\n  '%(cat,localGlobalScope,len(domains), len(uniqueDom), len(genes)))
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
    foutLog.write('\n**section about split cases and totality of them ')
    for h in genFrac:
        foutLog.write('\nfor the %s exonmerger\n'%(h))
        for WM in genFrac[h]:
            totalDom=sum(genFrac[h][WM].values())
            foutLog.write('\tthe %s domainType, total split cases are %s\n'%(WM,totalDom))
            for ju in junc:
                count=genFrac[h][WM][ju]
                foutLog.write('\t\tfor junctype %s total are %s (%s) \n'%(ju,count, div_fact(count,totalDom)))
    
    
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
l1(geneEx,logf)


l2(geneEx,outdir+'perFractionDomains_')
l3(geneEx,outdir+'FractionAG_', logf )
logf.close()
#python bin_analysis/analyse_domains.py results/exonsWise.pickle results/domain_annotation.pickle results/
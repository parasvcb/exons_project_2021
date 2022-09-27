import sys,os
import pandas as pd
import cPickle as pickle
import constructing_data.Classes_exons as Classes_exons
sys.modules['Classes_exons'] = Classes_exons
from collections import Counter
if len(sys.argv)!=4:
    print ('Please enter correct args 1. geneExondata 2 humangeneob 3. outputdir')
    sys.exit()
prog,geneEx,human,outdir=sys.argv

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
def addvaltoHaslis(has,key,val,exon=False):
    if exon:
        val=val+'?'+exon
    if key not in has:
        has[key]=[]
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
def l1(has, fout, human, category,WorM):
    #category={'domains_conseCoding':'ConsecCoding', 'domains':'normal'}
    #one of those but in list
    #WorM is Whole or middle
    completeTag='containedW' if WorM=='W' else 'containedM'
    splitTag='splitW' if WorM=='W' else 'splitM'
    domindex=1 if WorM=='W' else 2
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
                exContrib=domainContrib['exContrib']
                if comp:
                   
                    exid=''
                    for exon in exContrib:
                        if exContrib[exon][domindex]>=0.95:
                            exid=exon
                            break

                    key='1CORE' if domainContrib[completeTag][0] else '2ATIT'
                    junc=domainContrib[completeTag][0] if domainContrib[completeTag][0] else domainContrib[completeTag][1]
                    #print ('Yes',gene,domname, domainContrib[completeTag],comp,key)
                    contained[key][junc]=addvaltoHaslis(contained[key][junc],gene,domname,exon=exid)

                    contained['3BOTH'][domainContrib[completeTag][2]]=addvaltoHaslis(contained['3BOTH'][domainContrib[completeTag][2]],gene,domname,exon=exid)
                else:
                    CORE,ATIT,BOTH,COREATIT_cat,COREATIT_realTup,COREATIT_frac,BOTH_cat,BOTH_realTup,BOTH_frac=domainContrib[splitTag]
                    #splitW': [{'A': [0, 0.0], 'G': [91, 0.711]}, 
                    #          {'A': [37, 0.289], 'G': [0, 0.0]}, 
                    #          {'A': [37, 0.289], 'G': [91, 0.711]}, 
                    #           '1CORE_G__2ATIT_A', [1, 'G', 2, 'A'], [0.711, 0.289], 
                    #           '3BOTHAG', [3, 'A', 3, 'G'], [0.289, 0.711]]
                    #print(BOTH_cat)
                    key=BOTH_cat[:5]
                    junc=BOTH_cat[5:]
                    #print (key,junc,'split')
                    exid=''
                    for exon in exContrib:
                        if exContrib[exon][domindex]:
                            exid+=exon+':'+str(exContrib[exon][domindex])+', '
                    split[key][junc]=addvaltoHaslis(split[key][junc],gene,domname,exon=exid)
                    
                    if len(COREATIT_cat) in [6,7]:
                        key=COREATIT_cat[:5]
                        junc=COREATIT_cat[5:]
                        exid+='-'+("_".join(map(str,COREATIT_frac)))
                        split[key][junc]=addvaltoHaslis(split[key][junc],gene,domname,exon=exid)
                        
                    else:
                        key='1CORE+2ATIT'
                        junc=COREATIT_cat
                        exid+='-'+("_".join(map(str,COREATIT_frac)))
                        split[key][junc]=addvaltoHaslis(split[key][junc],gene,domname,exon=exid)

                    '''
                    '1COREA','1COREAG','1COREG','2ATITA','2ATITAG','2ATITG','1CORE_A__2ATIT_A','1CORE_A__2ATIT_G','1CORE_G__2ATIT_G','1CORE_G__2ATIT_A','1CORE_AG__2ATIT_A','1CORE_AG__2ATIT_G','1CORE_AG__2ATIT_AG','1CORE_A__2ATIT_AG','1CORE_G__2ATIT_AG'
                    '''
    #print (contained['A'])
    
    #print (split['1CORE']['G'])
    #fout.write('Category: %s, domains:%s\n'%(category,WorM))
    
    #fout.write('Broadcat\texScope\tTotalDomain_Core\tUniqueDom_Core\tGenes_Core\tTotalDomain_ATIT\tUniqueDom_ATIT\tGenes_ATIT\tTotalDomain_Both\tUniqueDom_Both\tGenes_Both\tTotalDomain_ATIT+Core\tUniqueDom_ATIT+Core\tGenes_ATIT+Core\n')
    genes=total.keys()
    domains=[i for j in total.values() for i in j]
    print (domains[:5])
    uniqueDom=set(domains)
    #fout.write('Total\t.\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\t.\t.\t.\n'%(len(domains), len(uniqueDom), len(genes)))
    
    string='DomainID\tDomainName\tDomLength\tCoordinates\tGeneID\tGeneName\tPILength\tExonType\texonID\tContained\tCOREorATIT\tFracContrib\n'
    vals3=['1CORE','2ATIT']
    #print(contained['1CORE']['A'])
    for ind,ttype in enumerate(vals3):
        for exJunc in ['A','G']:   
            genes=contained[ttype][exJunc].keys()
            for gene in genes:
                for domain in contained[ttype][exJunc][gene]:
                    #print (domain)
                    dom,exon=domain.split('?')
                    dom,coordinates=dom.split('::')
                    domID,domName=dom.split(':')
                    st,end=coordinates.split(',')
                    domLength=(int(end)-int(st))+1
                    pilen=human[gene].PI.seqlen
                    string+='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tcontained\t%s\tNULL\n'%(domID,domName,domLength,coordinates,gene,human[gene].detail,pilen,exJunc,exon,ttype)
    '''
    split={
        '2ATIT':{'A':{},'G':{},'AG':{}},
        '3BOTH':{'A':{},'G':{},'AG':{}},
        '1CORE':{'A':{},'G':{},'AG':{}},
        '1CORE+2ATIT':{'1CORE_A__2ATIT_A':{},'1CORE_A__2ATIT_G':{},'1CORE_G__2ATIT_G':{},'1CORE_G__2ATIT_A':{},'1CORE_AG__2ATIT_A':{},'1CORE_AG__2ATIT_G':{},'1CORE_AG__2ATIT_AG':{},'1CORE_A__2ATIT_AG':{},'1CORE_G__2ATIT_AG':{}}
        }
    '''
    vals3=['2ATIT','1CORE']
    for ind,ttype in enumerate(vals3):
        for exJunc in ['A','G','AG']:   
            genes=split[ttype][exJunc].keys()
            for gene in genes:
                for domain in split[ttype][exJunc][gene]:
                    #print (domain)
                    dom,exon=domain.split('?')
                    exon,fraction=exon.split('-')
                    dom,coordinates=dom.split('::')
                    domID,domName=dom.split(':')
                    st,end=coordinates.split(',')
                    domLength=(int(end)-int(st))+1
                    pilen=human[gene].PI.seqlen
                    string+='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tsplit\t%s\t%s\n'%(domID,domName,domLength,coordinates,gene,human[gene].detail,pilen,exJunc,exon,ttype,fraction)
    
    exJuncLis=split['1CORE+2ATIT'].keys()
    
    for exJunc in exJuncLis:  
            genes=split['1CORE+2ATIT'][exJunc].keys()
            for gene in genes:
                for domain in split['1CORE+2ATIT'][exJunc][gene]:
                    #print (domain)
                    dom,exon=domain.split('?')
                    exon,fraction=exon.split('-')
                    dom,coordinates=dom.split('::')
                    domID,domName=dom.split(':')
                    st,end=coordinates.split(',')
                    domLength=(int(end)-int(st))+1
                    pilen=human[gene].PI.seqlen
                    string+='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tsplit\t%s\t%s\n'%(domID,domName,domLength,coordinates,gene,human[gene].detail,pilen,exJunc,exon,'1CORE+2ATIT',fraction)
    
    return string
                    
   


def readpickle(fname):
    with open (fname,'rb') as fin:
        dat=pickle.load(fin)
    return dat


geneEx=readpickle(geneEx)
human=readpickle(human)

logfile=outdir+'domainsSectionStatsLog'
logf=open (logfile,'w')


print (geneEx[2051])
str1=l1(geneEx,logf,human, category=['domains','normal'], WorM='W')
str2=l1(geneEx,logf,human, category=['domains_conseCoding','ConsecCoding'], WorM='W')

def getgenes(str1,human):

    # DomainID	DomainName	DomLength	Coordinates	GeneID	GeneName	PILength	ExonType	exonID	Contained	COREorATIT	FracContrib
    # PF07645.15	EGF_CA	43	796,838	7173	TPO thyroid peroxidase	933	A	T.1.A.14.0.0	contained	1CORE	NULL
    # PF00096.26	zf-C2H2	25	266,290	57862	ZNF410 zinc finger protein 410	516	A	T.1.F.8.0.0	contained	1CORE	NULL
    # PF00
    dat=[i for i in str1.split('\n')[1:] if len(i)>10]
    genes=[]
    for i in dat:
        #print (i.split('\t'))
        DomainID,DomainName,DomLength,Coordinates,GeneID,GeneName,PILength,ExonType,exonID,Contained,COREorATIT,FracContrib=i.split('\t')
        genes+=[(GeneID,GeneName)]
    genes=set(genes)
    
    strnew='Gene\tGeneName\tIsoformCount\tPI\tpilen\tPIexons\tGeneExons\n'
    for gene in genes:
        su=0
        exs=''
        for ex in human[int(gene[0])].PI.exons:
            exs+='%s (%s,%s) : '%(ex.ID,su,ex.length)
            su+=ex.length
        exsGene=''
        for ex in human[int(gene[0])].exons:
            exsGene+='%s (%s,%s) : '%(ex.ID,su,ex.length)
            su+=ex.length
        
        strnew+='%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(gene[0],gene[1],len(human[int(gene[0])].transcripts),human[int(gene[0])].PI.ID,human[int(gene[0])].PI.seqlen,exs,exsGene)
    return strnew

def getsets(str1):

    # DomainID	DomainName	DomLength	Coordinates	GeneID	GeneName	PILength	ExonType	exonID	Contained	COREorATIT	FracContrib
    # PF07645.15	EGF_CA	43	796,838	7173	TPO thyroid peroxidase	933	A	T.1.A.14.0.0	contained	1CORE	NULL
    # PF00096.26	zf-C2H2	25	266,290	57862	ZNF410 zinc finger protein 410	516	A	T.1.F.8.0.0	contained	1CORE	NULL
    # PF00
    '''
    Contained are the cases when its only contained, and not split, 
    '''
    dat=[i for i in str1.split('\n')[1:] if len(i)>10]
    hasCount={'contained':{'A':{},'G':{}},'split':{'G':{},'A':{},'AG':{}}}
    hasGenes={'contained':{'A':{},'G':{}},'split':{'G':{},'A':{},'AG':{}}}

    for i in dat:
        #print (i.split('\t'))
        DomainID,DomainName,DomLength,Coordinates,GeneID,GeneName,PILength,ExonType,exonID,Contained,COREorATIT,FracContrib=i.split('\t')
        keyexon=ExonType if len(ExonType)<=2 else ''
        Gene=(GeneID,GeneName,PILength,DomLength,Coordinates)
        DomainNameID=(DomainName,DomainID)
        if not keyexon:
            junc=''
            
            one,two=ExonType.split('__')
            ele1=one.split('_')[1]
            ele2=two.split('_')[1]
            junc+=ele1+ele2
            junc =set(junc)
            if len(junc)=='2':
                keyexon='AG'
            else:
                keyexon=list(junc)[0]
        if DomainNameID not in hasCount[Contained][keyexon]:
            hasCount[Contained][keyexon][DomainNameID]=0
        hasCount[Contained][keyexon][DomainNameID]+=1
        if DomainNameID not in hasGenes[Contained][keyexon]:
            hasGenes[Contained][keyexon][DomainNameID]={}
        if Gene not in hasGenes[Contained][keyexon][DomainNameID]:
            hasGenes[Contained][keyexon][DomainNameID][Gene]=0
        hasGenes[Contained][keyexon][DomainNameID][Gene]+=1

    containedDomain={domain:0 for ex in hasCount['contained'] for domain in hasCount['contained'][ex]}
    SplitDomain={domain:0 for ex in hasCount['split'] for domain in hasCount['split'][ex]}
    print (len(containedDomain),'containedDomain')
    print (len(SplitDomain),'splitDomain')
    pureContain=set(containedDomain.keys())-set(SplitDomain.keys())
    pureSplit=set(SplitDomain.keys())-set(containedDomain.keys())
    print (len(pureContain),'purecontain')
    print (len(pureSplit),'pureSplit')
    
    hasPure={'contained':pureContain,'split':pureSplit}

    strnew='domainname\tDomainID\tLength\tCoordinates\tExonType\tContained\tGenesCount\ttotalDomainCount\tcountinGene\tpilength\tGene\n'
    for cat in hasCount:
        for ex in hasCount[cat]:
            for domain in hasCount[cat][ex]:
                if domain in hasPure[cat]:
                    genes=len(hasGenes[cat][ex][domain].keys())
                    domainCount=hasCount[cat][ex][domain]
                    geneIds={}
                    for gene in hasGenes[cat][ex][domain]:
                        idg,name,pilength,domlength,coods=gene
                        if idg not in geneIds:
                            geneIds[idg]=[(name,pilength),[(domlength,coods)]]
                        else:
                            geneIds[idg][1]+=[(domlength,coods)]
                    for geneIter in geneIds:
                        gname,pilength=geneIds[geneIter][0]
                        for instances in geneIds[geneIter][1]:
                            length,coods=instances
                            strnew+='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(domain[0],domain[1],length,coods,ex,cat,len(geneIds),domainCount,len(geneIds[geneIter][1]),pilength,gname+':'+str(geneIter))
    return strnew

with open(outdir+'domainsWrapped_normal.csv','w')as fout:
    fout.write('%s\n'%str1)

strnew=getgenes(str1,human)
with open (outdir+'domainsWrapped_normal_genes.csv','w')as fout:
    fout.write('%s\n'%strnew)

with open(outdir+'domainsWrapped_conseCoding.csv','w')as fout:
    fout.write('%s\n'%str2)

string=getsets(str1)
with open (outdir+'domainsWrapped_normal_pureSplitContained.csv','w')as fout:
    fout.write('%s\n'%string)

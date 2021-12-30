import sys,os
import pandas as pd
import cPickle as pickle
from collections import Counter
if len(sys.argv)!=5:
    print ('Please enter correct args 1. geneExondata 2. geneDomain data 4. genePDBobject 3. outputdir')
    sys.exit()
prog,geneEx,geneDom,genePDBobject,outdir=sys.argv

fappend=geneEx.split('_')[-1].split('.')[0]
def div_fact(num,denom):
	try:
		return round(float(num)/denom,3)	
	except:
		return 0

'''
Have the geneEx in format of  has[gene][exID]=[(dname,dfrac,dfrac_m,efrac_m,efrac),]
and geneDom in has[gene][domname]=
                                  {'whole':
                                            (exontype[AG/AA/GG], firstlet frac, second let frac),
                                    'core':
                                            (exontype[AG/AA/GG], firstlet frac, second let frac),
                                  }
1st level:
read the second and startcovering the genes, if that is either of GG # but this wont give if individual exon wrapped it or multiple such cases
so read the geneEx data, and for every exon check
             if dfrac is 1 or dfrac_m is 1, add them separtely,
                keep a track of the genes and domain names in has
                
             else these events will be split events
                for split events, keep a track of the gene and domname, and then query later the anothe has for their fraction,

'''


def readpickle(fname):
    with open (fname,'rb') as fin:
        dat=pickle.load(fin)
    return dat

def addvaltoHaslis(has,key,val):
    if key not in has:
        has[key]=[]
    if val not in has[key]:
        has[key]+=[val]
    return has

def writeLevelled(has,fname):
    with open (fname,'w') as fout:
        fout.write('Gene\tDomain\tIntegrity\tsplitorNot\n')
        for splitorNot in has:
            for domIntegrity in has[splitorNot]:
                for gene in has[splitorNot][domIntegrity]:
                    for doms in has[splitorNot][domIntegrity][gene]:
                        fout.write('%s\t%s\t%s\t%s\n'%(gene,doms,domIntegrity,splitorNot))

def toDataframe(hasSplitCases,fname):
    exJunc='AG'
    for Integrity in hasSplitCases:
        f_name=fname+str(Integrity)+'.csv'
        if exJunc in hasSplitCases[Integrity]:
            with open (f_name,'w') as fout:
                fout.write('Gene,dname,A,G,pdbInfo\n')
                print (hasSplitCases[Integrity].keys())
                for i in hasSplitCases[Integrity][exJunc]:
                    gene,dname,a,g=i
                    name=dname.split('::')[0]
                    pdbinfo=genePDB[int(gene)] if int(gene) in genePDB else False
                    fout.write('%s,%s,%s,%s,%s\n'%(gene,name,a,g,pdbinfo))

def writeDFfractions(has,fname):
    total=sum([has[i]['totalDom'] for i in has])
    for i in has:
        has[i]['domFrac']=div_fact(has[i]['totalDom'],total)
    df= pd.DataFrame.from_dict(has,orient='index')
    df.to_csv(fname, sep='\t')
    with open (fname) as fin:
        dat=fin.read().split('\n')
        dat[0]="\t".join(['category']+dat[0].split())
    with open (fname,'w') as fin:
        for i in dat:
            fin.write('%s\n'%(i))
    
def level1(domhas,exhas, outname):
    '''
    it will parse the exhas having the general gene exon wise data, and based on the dfrac and dfracm values it will segregate them to the contained and split cases
    in has named has
    every individaul domain per gene if got split was analyzed for its fraction and that will be written below to a file, quadrant will be analysed,
    now two ultimate returns are needed, 
    domains that got split
        40 
        100
    domains that got contained 
        40
        100
    but with only the names 
    '''
    has={'contained':{
        100:{},40:{}},
        'split':{
        100:{},40:{}}}
    for gene in exhas:
        for exid in exhas[gene]:
            for splittingevents in exhas[gene][exid]:
                dname,dfrac,dfrac_m,efrac_m,efrac=splittingevents
                if dfrac==1:
                    has['contained'][100]=addvaltoHaslis(has['contained'][100],gene,dname)
                else:
                    has['split'][100]=addvaltoHaslis(has['split'][100],gene,dname)
                if dfrac_m==1:
                    has['contained'][40]=addvaltoHaslis(has['contained'][40],gene,dname)
                else:
                    has['split'][40]=addvaltoHaslis(has['split'][40],gene,dname)
    '''
    Transforming the default gene wise exon fractions/per domain to above format
    keys in fasion of contained/split category > Integrity Check of 100/40, then add gene as key and in end unique domains
    always contained and domains that got split are modified.
    '''
    #writeLevelled(has,outname)
    
    has_return_split={}
    has_return_enveloped={}
    for splitype in has:
        has_choosen=has_return_split if splitype=='split' else has_return_enveloped        
        for Integrity in has['split']:
            if Integrity not in has_choosen:
                has_choosen[Integrity]=set()
            for gene in has['split'][Integrity]:
                for dnames in has['split'][Integrity][gene]:
                    dname_simple=dnames.split('::')[0]
                    has_choosen[Integrity]|=set([dname_simple])
    '''
    Using the organised has dict above to specifically tag the doomains names that will undergo split or containment, this will be returned and wont of our interest but to marking domain sets undegroing always split or contained.
    '''

    has_split_analysis={}
    for Integrity in has['split']:
        #starting from integrity of 100/40 as splittype was mentioned apriori
        if Integrity not in has_split_analysis:
            has_split_analysis[Integrity]={}
        for gene in has['split'][Integrity]:
            for dnames in has['split'][Integrity][gene]:
                dname_simple=dnames.split('::')[0]
                exonJuncType,frac1,frac2=domhas[gene][dnames][Integrity]
                if exonJuncType not in has_split_analysis[Integrity]:
                    has_split_analysis[Integrity][exonJuncType]=[]
                #print (exonJuncType)
                has_split_analysis[Integrity][exonJuncType]+=[(gene,dname_simple,frac1,frac2)]

    toDataframe(has_split_analysis,outname+'dataframeAGjunction')
    '''
    This is the analysis of subset undergoing split and their fraction codings, specifically the interest would be get the histogram of such events and then to later demarcate the exact contribution per exon or gene
    -> has to write to a file, and histogram is pending
    '''
    hasrefined={}
    for splittingevents in has:
        #contained or split
        if splittingevents not in hasrefined:
            hasrefined[splittingevents]={}
        for Integrity in has[splittingevents]:
            #40, 100
            if Integrity not in hasrefined[splittingevents]:
                hasrefined[splittingevents][Integrity]=[]
            domlis=[]
            for gene in has[splittingevents][Integrity]:
                for domains in has[splittingevents][Integrity][gene]:
                    domlis+=[domains.split('::')[0]]
            hasrefined[splittingevents][Integrity]=domlis
    '''
    has_refined was created to get a concise summary for writeup only, but for the fractions, sets should not be used as  they may remove redundancy for multiple domaisn in single protein and everywhere
    '''
    dataframe_100={'contained':{'domFrac':0,'totalDom':0,'uniqueDom':0,'geneCount':0},
                    'AG':{'domFrac':0,'totalDom':0,'uniqueDom':0,'geneCount':0},
                    'GG':{'domFrac':0,'totalDom':0,'uniqueDom':0,'geneCount':0},
                    'AA':{'domFrac':0,'totalDom':0,'uniqueDom':0,'geneCount':0},
                    }
    dataframe_40={'contained':{'domFrac':0,'totalDom':0,'uniqueDom':0,'geneCount':0},
                    'AG':{'domFrac':0,'totalDom':0,'uniqueDom':0,'geneCount':0},
                    'GG':{'domFrac':0,'totalDom':0,'uniqueDom':0,'geneCount':0},
                    'AA':{'domFrac':0,'totalDom':0,'uniqueDom':0,'geneCount':0},
                    }
    string=''
    for i in hasrefined:
        #contained,split
        string+='For %s cases\n'%i
        for j in hasrefined[i]:
            #40,100
            df = dataframe_100 if j==100 else dataframe_40
            string+='\tif %spcnt domain was considered\n'%j
            #print (len(has[i][j]))
            genes=len(has[i][j])
            domains=sum([len(has[i][j][k]) for k in has[i][j]])
            string+='\t\t %s genes and %s domains participated\n'%(genes,domains)
            #print (hasrefined[i][j])
            consolidate=Counter(hasrefined[i][j])
            consolidate=[[consolidate[k],k] for k in consolidate]
            consolidate.sort(reverse=True)
            #print (consolidate)
            unique=len(consolidate)
            if i=='contained':    
                df[i]['totalDom']=domains
                df[i]['uniqueDom']=unique
                df[i]['geneCount']=genes
                
            string+='\t\t %s of those domains were unique\n' %(unique)
            string+='\t\t top 30 of them were\n%s'%('\t\t\t\n'.join([str(k[0])+','+k[1] for k in consolidate [:30] ]) )
            string +='\n ---------\n\n'
        string +=' ********\n'

   

    string+='->For split cases in general\n'
    for Integrity in has_split_analysis:
        #100,40
        df = dataframe_100 if Integrity==100 else dataframe_40
        string+='\tif %s domains were considered\n'%Integrity
        for exjunc in has_split_analysis[Integrity]:
            #print (exjunc)
            string+='\t\tFor ex_junc type %s\n'%exjunc
            geneslis=[i[0] for i in has_split_analysis[Integrity][exjunc]]
            dnames=[i[1] for i in has_split_analysis[Integrity][exjunc]]
            consolidate=Counter(dnames)
            consolidate=[[consolidate[k],k] for k in consolidate]
            consolidate.sort(reverse=True)
            #print (df)
            df[exjunc]['totalDom']=len(dnames)
            df[exjunc]['uniqueDom']=len(set(dnames))
            df[exjunc]['geneCount']=len(set(geneslis))
            string+='\t\t\t%s genes involved\n'%(len(set(geneslis)))
            string+='\t\t\t%s domains involved (%s unique)\n'%(len(dnames),len(set(dnames)))
            string+='\t\t\t Top 30 of them were %s\n'%("\n".join([str(k[0])+','+k[1] for k in consolidate [:30] ]))
            string +='\n ---------\n\n'
        string +=' ********\n'
    #has_split_analysis[Integrity]+=[(gene,dname_simple,frac1,frac2)]
    with open (outname+'.summ','w') as fout:
        fout.write('%s'%(string))

    writeDFfractions(dataframe_100, outname+'100fractionDomains.csv')
    writeDFfractions(dataframe_40, outname+'40fractionDomains.csv')
    
    return has_return_enveloped,has_return_split

geneEx=readpickle(geneEx)
geneDom=readpickle(geneDom)
genePDB=readpickle(genePDBobject)


def updateGeneEx(geneEx,gene,exID,val):
    if gene not in geneEx:
        geneEx[gene]={}
    if exid not in geneEx[gene]:
        geneEx[gene][exid]=[]
    geneEx[gene][exid]+=[val]
    return geneEx


h_ret_env,h_ret_spl=level1(geneDom,geneEx,os.path.join(outdir,'Fig3D_%s_level1_broadview'%fappend))
geneExneverSplit={}
geneExalwaysSplit={}

#1{
for gene in geneEx:
    for exid in geneEx[gene]:
        for splittingevents in geneEx[gene][exid]:
            dname,dfrac,dfrac_m,efrac_m,efrac=splittingevents
            dname_simple=dname.split('::')[0]
                    
            if dfrac==1:
                #complete domain (envelope)
                if dname_simple not in h_ret_spl[100]:
                    geneExneverSplit=updateGeneEx(geneExneverSplit,gene,exid,splittingevents)

            else:
                #this event splits
                if dname_simple not in h_ret_env[100]:
                    geneExalwaysSplit=updateGeneEx(geneExalwaysSplit,gene,exid,splittingevents)

            if dfrac_m==1:
                #complete domain (envelope)
                if dname_simple not in h_ret_spl[40]:
                    geneExneverSplit=updateGeneEx(geneExneverSplit,gene,exid,splittingevents)
            else:
                #this event splits
                if dname_simple not in h_ret_env[40]:
                    geneExalwaysSplit=updateGeneEx(geneExalwaysSplit,gene,exid,splittingevents)
#}1 This block deals and only stores the domain names which undergoes split and other events

h_ret_env,h_ret_spl=level1(geneDom,geneExneverSplit,os.path.join(outdir,'Fig3D_%s_level2_notsplitting'%fappend))
h_ret_env,h_ret_spl=level1(geneDom,geneExalwaysSplit,os.path.join(outdir,'Fig3D_%s_level2_alwayssplitting'%fappend))


#python bin_analysis/analyse_domains.py results/exonsWise.pickle results/domain_annotation.pickle results/
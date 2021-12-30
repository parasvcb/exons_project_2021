import cPickle as pickle
import common.general_modules as gm
import figPanels.modules_analysis as am
import sys,os
if len(sys.argv)!=4:
        print ("Please give 1 object 2 conditionFile 3 results dir 4 ATIFlag  ATITrue=0, ATIFalse=1, reagrdless=-1")
        sys.exit()
prog,source_gene_object,gene_conFile,results_dir_csv = sys.argv


fileswriter=open(results_dir_csv+"ss.log","w")
window=3

humanGeneObject=gm.readPickle(source_gene_object)

CONDITION_GENES=gm.readPickle(gene_conFile)



exonsATI={(0,0.10):0,(0.101,0.20):0,(0.201,0.30):0,(0.301,0.40):0,(0.401,0.50):0,(0.501,0.60):0,(0.601,0.70):0,(0.701,0.80):0,(0.801,0.90):0,(0.901,1.01):0}
protATI={(0,0.10):0,(0.101,0.20):0,(0.201,0.30):0,(0.301,0.40):0,(0.401,0.50):0,(0.501,0.60):0,(0.601,0.70):0,(0.701,0.80):0,(0.801,0.90):0,(0.901,1.01):0}

for gene in CONDITION_GENES:
    PI=[i for i in humanGeneObject[gene].transcripts if i.PI][0]
    if PI:
        tExonList=am.exonScreenerBetweenTConstitutive(humanGeneObject[gene].exons)
        PI_representative=[]
        exonsTotal=0
        exonsCore=0
        protTotal=0
        protCore=0
        for ex in PI.exons:
            numericFlag=int(ex.ID.split(".")[3])
            if numericFlag in tExonList:
                exonsCore+=1
                protCore+=ex.length
            exonsTotal+=1
            protTotal+=ex.length
        exCoreFrac=round(float(exonsCore)/exonsTotal,2)
        protCoreFrac=round(float(protCore)/protTotal,2)
        flage=True
        for rangeVal in exonsATI:
            if rangeVal[0]<=exCoreFrac<=rangeVal[1]:
                flage=False
                exonsATI[rangeVal]+=1
        if flage:
            print (exCoreFrac)
        flagp=True
        for rangeVal in protATI:
            if rangeVal[0]<=protCoreFrac<=rangeVal[1]:
                flagp=False
                protATI[rangeVal]+=1
        if flagp:
            print (protCoreFrac)
list_WEF=list(protATI.keys())
list_WEF.sort()
f1=open (os.path.join(results_dir_csv+'ATICoreExonsProteinFraction.csv'),'w')
f1.write('FractionCovered\tTotal\tFreq\tCumFreq\tTag\n')

countvar_exons=sum(protATI.values())
Tag='ProteinFrac'
print (Tag,countvar_exons)
cumFreq=0
c=0
for WEFrange in list_WEF:
    count=protATI[WEFrange]
    c+=count
    freq=gm.div_fact(count,countvar_exons)
    cumFreq+=freq
    f1.write('%s\t%s\t%s\t%s\t%s\n'%("_".join(map(str,list(WEFrange))),count,freq,cumFreq, Tag))
print (c)
cumFreq=0
countvar_exons=sum(exonsATI.values())
print (Tag,countvar_exons)
c=0
Tag='ExonsFrac'
for WEFrange in list_WEF:
    count=exonsATI[WEFrange]
    c+=count
    freq=gm.div_fact(count,countvar_exons)
    cumFreq+=freq
    f1.write('%s\t%s\t%s\t%s\t%s\n'%("_".join(map(str,list(WEFrange))),count,freq,cumFreq, Tag))
print (c)
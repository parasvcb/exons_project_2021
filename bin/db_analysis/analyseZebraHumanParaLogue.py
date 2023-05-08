import re
with open ("one_to_many_zebrafish_ortholog_data") as fin:
    dat=[i for i in fin.read().split('\n')[1:] if len(i)>0]

# 28990:4,['NP_001275879.1', 'NP_054784.2', 'XP_024309254.1', 'XP_016861755.1']	100003699:1,['XP_001343186.3']	100006826:1,['NP_001092900.1']	
# 23440:1,['NP_115485.1']	30316:2,['XP_005165562.1', 'NP_571175.1']	560759:2,['XP_005170306.1', 'NP_001122175.1']
dictData = {}

lis=[]
total=0
morethan1 =0
for i in dat:
    total+=1
    j=re.sub('\[.+?\]', '', i)
    j=re.sub('\,', ' ', j)
    j=j.split()
    humanId,humanTranscripts= j[0].split(':')
    totalFishTrans=0
    if len(j[1:])>1:
        morethan1+=1
        for m in j[1:]:
            FishId,FishTranscripts= m.split(':')
            totalFishTrans+=int(FishTranscripts)
        diff = int(humanTranscripts) -  totalFishTrans
        frac = round(float(diff)/float(humanTranscripts),3)
        lis+=[(humanId, humanTranscripts, totalFishTrans, diff, frac)]

print (total,morethan1)
with open ("tableViolin.tsv",'w') as fout:
    fout.write('GeneId\thumanIsf\tFishIsf\tdiff\tfrac\n')
    for i in lis:
        fout.write('%s\t%s\t%s\t%s\t%s\n'%(i[0],i[1],i[2],i[3],i[4]))

    
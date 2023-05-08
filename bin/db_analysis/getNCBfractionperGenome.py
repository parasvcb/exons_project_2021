import sys
import common.general_modules as cm
import figPanels.modules_analysis as ca
import common.general_modules as cm

if len(sys.argv)!=3:
        print ("Please give 1 object 2 condFile")
        sys.exit()
prog, source_gene_object, conditionFile = sys.argv

humanGeneObject=cm.readPickle(source_gene_object)
condition = cm.readPickle(conditionFile)

n,c,b = [0,0,0]
nT,cT,bT = [0,0,0]
nG,cG,bG = [[],[],[]]

# with open ("fishGeneCondition2",'w') as fout:
#      for i in condition:
#           fout.write('%s\n'%(i))

totalExons=0
totalCexons=0
totalCGexons=0
totalCAexons=0
totalCFexons=0
for gene in humanGeneObject:
    if gene in condition: # and gene == 11259:
        exonMatrix = ca.positionalExonMatrix_forExCharacterization(humanGeneObject[gene])
        #{1: ['U', 'F', [2, 0, 2], 0, 1, 0, 0], 2: ['U', 'A', [0, 0, 0], 0, 1, 0, 0], 3: ['T', 'A', [0, 0, 0], 0, 0, 0, 0], 4: ['T', 'F', [0, 1, 0], 0, 0, 0, 0]}
        
        # exonMatrix[placeholder] = ['','',[0,0,0], 0,0,0, 0, [0,0,0]] #0, 1,    2   , 3,4,5, 6, 7
        # 0th ele, UTMRD tag, 
        # 1st ele, AGF tag, 
        # 2nd records children order of ncb, and 
        # 3rd occurrences of the aa change on this position, (should not be evaluated if values of the ncb counterparts is true
        # 4th coding intron retention events starting from this
        # 5th non coding intron retention events starting from this
        # 6th whether aa has been ever removed from this sequence (yes consider all the ncb variation also but not retention)
        # 7th whether ncb variations ever had change in aa case, like FMR1 gene does have, in that case it will be binary and value 1 means atleast 1 sduch case exists, 
        totalExons += len(exonMatrix)
        for pos in exonMatrix:
            utmrd, agf, ncb, t1, t2, t3, t4, t5 = exonMatrix[pos]
            if utmrd =='T':
                totalCexons+=1
                if agf=='A':
                    totalCAexons+=1
                elif agf=="F":
                    totalCFexons+=1
                else:
                    totalCGexons+=1

                if agf in ["A","F"]:
                    if ncb[0]>0:
                        n+=1
                        nT+=ncb[0]
                        nG+=[gene]
                    if ncb[1]>0:
                        c+=1
                        cT+=ncb[1]
                        cG+=[gene]
                    if ncb[2]>0:
                        b+=1
                        bT+=ncb[2]
                        bG+=[gene]
                
nG=set(nG)
cG=set(cG)
bG=set(bG)
print ('genesTotalandCondition: %s\t%s'%(len(humanGeneObject), len(condition)))
print ('TotalExons:%s|totalCodingExons:%s|totalCodingGexons:%s|totalCodingAexons:%s|totalCodingFexons:%s' %(totalExons, totalCexons, totalCGexons, totalCAexons, totalCFexons))
print ('PosAffected\t%s\t%s\t%s'%(n,c,b))
print ('TotPolymInPos\t%s\t%s\t%s'%(nT,cT,bT))
print ('GenesAffected\t%s\t%s\t%s'%(len(nG),len(cG),len(bG)))
print ('GenesAffected\t%s\t%s\t%s'%(cm.div_fact(len(nG), len(condition)),cm.div_fact(len(cG), len(condition)),cm.div_fact(len(bG), len(condition))))
print ('TotalGenesAffected\t%s' %(len(nG|cG|bG)))
print ('TotalGenesAffected\t%s' %(cm.div_fact(len(nG|cG|bG), len(condition))))


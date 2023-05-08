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

fname=open("MasterFile_%s.tsv"%org,'w')
fname.write('Gene\tName\tExons\tTranscript\tPILength\tUA\tUAss\tUG\tUF\tDA\tDAch\tDGch\tDFch\tDAss\tDG\tDF\tTA\tTAch\tTAss\tTG\tTGch\tTF\tTFch\n')
for gene in humanGeneObject:
    if gene in condition:
        exonsName = [i.ID for i in humanGeneObject[gene].exons]
        transInterest = [(i.PI, i) for i in humanGeneObject[gene].transcripts]
        transInterest.sort()
        length = sum([i.length for i in transInterest[0][1].exons])


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

        Strict=[i for i in exonMatrix if exonMatrix[i][0]==utmrdtag]
    MajorlyConstitutive=[i for i in exonMatrix if exonMatrix[i][0]==utmrdtag and exonMatrix[i][1]=='F']   
    Alternate=[i for i in exonMatrix if exonMatrix[i][0]==utmrdtag and exonMatrix[i][1]=='A' and sum(exonMatrix[i][2]) ==0]
    AlternateWithSS=[i for i in exonMatrix if exonMatrix[i][0]==utmrdtag and exonMatrix[i][1]=='A' and sum(exonMatrix[i][2])>0]  
    Constitutive=[i for i in exonMatrix if exonMatrix[i][0]==utmrdtag and exonMatrix[i][1]=='G']

    for tagListUtmrdConAlt in [("D","G"), ("D","A"), ("D","F"), ("T","G"), ("T","A"), ("T","F")]:
        utmrdtag, constAltTag = tagListUtmrdConAlt
        ConstitutiveWaachange = [i for i in exonMatrix if exonMatrix[i][0]==utmrdtag and exonMatrix[i][1]==constAltTag and (exonMatrix[i][3]>0 or sum(exonMatrix[i][7])>0)]

        if "G" not in ",".join(exonsName) and "F" not in ",".join(exonsName):
             fnameAlt.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(org, gene, length, len(transInterest), len(exonsName),humanGeneObject[gene].detail))
        else:
             fnameConst.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(org, gene, length, len(transInterest), len(exonsName),humanGeneObject[gene].detail))
fnameConst.close()
fnameAlt.close()
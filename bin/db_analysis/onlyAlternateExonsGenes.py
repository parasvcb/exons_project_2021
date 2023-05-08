import sys
import common.general_modules as cm

if len(sys.argv)!=4:
        print ("Please give 1 object 2 condFile 3 orgTag")
        sys.exit()
prog, source_gene_object, conditionFile, org = sys.argv

humanGeneObject=cm.readPickle(source_gene_object)
condition = cm.readPickle(conditionFile)


fnameAlt=open("AlternateOnly_%s.tsv"%org,'w')
fnameConst=open("constandFOnly_%s.tsv"%org,'w')
fnameAlt.write('org\tGene\tPILength\tIsoforms\texonsCount\tgeneName\n')
fnameConst.write('org\tGene\tPILength\tIsoforms\texonsCount\tgeneName\n')
for gene in humanGeneObject:
    if gene in condition:
        exonsName = [i.ID for i in humanGeneObject[gene].exons]
        transInterest = [(i.PI, i) for i in humanGeneObject[gene].transcripts]
        # print (transInterest)
        # sys.exit()
        transInterest.sort()
        length = sum([i.length for i in transInterest[0][1].exons])
        if "G" not in ",".join(exonsName) and "F" not in ",".join(exonsName):
             fnameAlt.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(org, gene, length, len(transInterest), len(exonsName),humanGeneObject[gene].detail))
        else:
             fnameConst.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(org, gene, length, len(transInterest), len(exonsName),humanGeneObject[gene].detail))
fnameConst.close()
fnameAlt.close()

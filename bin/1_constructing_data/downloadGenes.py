import constructing_data.retriever as ret
import sys,os

if len (sys.argv)!=4:
    print ('Please enter correct args 1. flatfile having genes infor, where to store genes and, wher to store the proteins extracted')
prog,flatfile,outdirGenes,outdirProtein=sys.argv
ret.entrez_retriever(lis=[],nature="gene", directory=outdirGenes,fname=flatfile)
proteins=ret.id_extractorGT(gene_add=outdirGenes, nature='protein')

prevProteinSaveDir=os.path.dirname(os.path.normpath(outdirProtein))

print ("done genes")
with open (os.path.join(prevProteinSaveDir,"proteinsinGenome"),'w') as fin:
    for i in proteins:
        fin.write('%s\n'%i)

ret.entrez_retriever(lis=proteins,nature="protein", directory=outdirProtein)

print ('Done proteins')
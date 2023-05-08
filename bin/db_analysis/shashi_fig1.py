'''
doc:
23-01-2023
its important to figure out and plot the UTR variations laong in exondata and transdarta and make panel figure differently
'''
import common.general_modules as cm
import seaborn as sns
import scipy, os
import numpy as np
import pandas as pd
import figPanels.modules_analysis as ca

def variations(geneob):
    def numberUniqueExonsCoding(geneOrTrancript_ob):
        exonMatrix = ca.positionalExonMatrix_forExCharacterization(geneob)
        # { 1: ['U', 'F', [2, 0, 2], 0, 1, 0, 0], 
        #   2: ['U', 'A', [0, 0, 0], 0, 1, 0, 0], 
        #   3: ['T', 'A', [0, 0, 0], 0, 0, 0, 0], 
        #   4: ['T', 'F', [0, 1, 0], 0, 0, 0, 0]
        # }
        # 0th ele, UTMRD tag, 
        # 1st ele, AGF tag, 
        # 2nd records children order of ncb, and 
        # 3rd occurrences of the aa change on this position, (should not be evaluated if values of the ncb counterparts is true
        # 4th coding intron retention events starting from this
        # 5th non coding intron retention events starting from this
        # 6th whether aa has been ever removed from this sequenc
        pos_Coding = [i for i in exonMatrix if exonMatrix[i][0] in ['T','D']]
        return (len(pos_Coding), len(exonMatrix))
    UniqueCodingExonCountGene,UniqueAllExonCountGene = numberUniqueExonsCoding(geneob)
    dictIsfSeqAsKeyTransObAsValue = {}
    for transcripts in geneob.transcripts:
        seqTranscriptaa=''.join(i.seq for i in transcripts.exons)
        if seqTranscriptaa not in dictIsfSeqAsKeyTransObAsValue:
            dictIsfSeqAsKeyTransObAsValue[seqTranscriptaa]=[]
        dictIsfSeqAsKeyTransObAsValue[seqTranscriptaa]+=[transcripts]
    return ([UniqueCodingExonCountGene,UniqueAllExonCountGene, geneob.detail,[len(dictIsfSeqAsKeyTransObAsValue[i]) for i in dictIsfSeqAsKeyTransObAsValue],len(dictIsfSeqAsKeyTransObAsValue)])

def transcripts_type(has,organism,fin):
    for gene in has:
        info=variations(has[gene])
        UniqueCodingExonCountGene, UniqueAllExonCountGene, name, lisAllvar, uniqueProt = info
        Allvars = sum(lisAllvar)
        fin.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(gene,Allvars,uniqueProt,UniqueAllExonCountGene,UniqueCodingExonCountGene,organism))


obLis = {
    '/home/paras/exonsdrive/project/protein_splicing/projectDir/outdir/latest_5_org_2022_27_nov/objectsToDB/objectsave_6239__noStructure.pick':'worm',
    '/home/paras/exonsdrive/project/protein_splicing/projectDir/outdir/latest_5_org_2022_27_nov/objectsToDB/objectsave_7227__noStructure.pick':'fly',
    '/home/paras/exonsdrive/project/protein_splicing/projectDir/outdir/latest_5_org_2022_27_nov/objectsToDB/objectsave_7955__noStructure.pick':'fish',
    '/home/paras/exonsdrive/project/protein_splicing/projectDir/outdir/latest_5_org_2022_27_nov/objectsToDB/objectsave_10090__noStructure.pick':'mouse',
    '/home/paras/exonsdrive/project/protein_splicing/projectDir/outdir/latest_5_org_2022_27_nov/objectsToDB/objectsave_9606__noStructure.pick':'human',
}

with open ('geneExonTransUpdated.tab','w') as fin:
    fin.write("Gene\tTotalVariants\tUniqueProteinVariants\tTotalExons\tUniqueCodingExons\tOrganism\n")
    for i in obLis:
        org = obLis[i]
        humanGeneObject=cm.readPickle(i)
        transcripts_type(humanGeneObject,org,fin)
        
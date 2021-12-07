import cPickle as pickle
import common.general_modules as gm
import figPanels.modules_analysis as am
import sys
if len(sys.argv)!=4:
        print ("Please give 1 object 2 conditionFile 3 results dir")
        sys.exit()
prog,source_gene_object,gene_conFile,results_dir_csv = sys.argv

window=3

humanGeneObject=gm.readPickle(source_gene_object)
CONDITION_GENES=gm.readPickle(gene_conFile)

desiredlistofgenes={i for i in humanGeneObject if i in CONDITION_GENES}

for gene in desiredlistofgenes:
        PI=[i for i in gene_source[gene].transcripts if i.PI][0]
        if PI:
            pi_seq=''.join([i.seq for i in PI.exons])
            ss_seq=''
            for i in PI.exons:
                if i.length>0:
                    ssres=i.out_secondseq(PI.ID)
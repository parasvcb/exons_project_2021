import cPickle as pickle
import os
import re
from progress.bar import Bar
from shutil import rmtree
import subprocess
import sys

import constructing_data.retriever as retriever
import constructing_data.prediction_assignment_ss as prediction_assignment_ss
import constructing_data.mappings_refseqProtein as mappings_refseqProtein
import constructing_data.pfam_files as pfam_files
import constructing_data.cath_files_multi as cath_files_multi
import constructing_data.ncbia_raw_coods as ncbia_raw_coods
import constructing_data.structure_refseq as structure_refseq
import constructing_data.multifasta_to_fasta as multifasta_to_fasta
from constructing_data.Classes_exons import Gene as Gene
from constructing_data.Classes_exons import Exon as Exon
from constructing_data.Classes_exons import Transcript as Transcript


source_data = os.path.join(parent_dir,"source_data")
pfam_dir = os.path.join(source_data,"pfam")
ssp_raw = os.path.join(source_data,"SSP")
pid_add = os.path.join(source_data,'Refseq_protein')
pid_raw_multifasta = os.path.join(source_data,'%s_ref.faa' % organism)
ens_raw_multifasta = os.path.join(source_data,'%s_ens.faa.gz' % organism)
swiss_raw_multifasta = os.path.join(source_data,'%s_swiss.faa.gz' % organism)
ens_pid = os.path.join(source_data,'ensemble_data')
swiss_pid = os.path.join(source_data,'swissprot_data')
gene_add = os.path.join(source_data,'gene_tables')

derived_data = os.path.join(outdir1,"derived_data/")
pickle_add = os.path.join(derived_data,"pickles")
results_dir = os.path.join(derived_data,"results")
save_structure_dir = derived_data + 'structure_data/PDB_structure_derived_data/'
save_structure_dirAF = derived_data + 'structure_data/PDB_structure_derived_dataAF/'
#-> added this time
stride_reformat_dir = derived_data + "stride_pssm_form/"
stride_exons_dir = derived_data + "exons_wise/stride_exons/"
stride_reformat_dirAF = derived_data + "stride_pssm_formAF/"
stride_exons_dirAF = derived_data + "exons_wise/stride_exonsAF/"

pfam_der_writer = derived_data + "domains/pfam/"
cath_der_writer = derived_data + "domains/cath/"
ss_exons_dir = derived_data + "exons_wise/ss_exons/"
ss_exons_dir06 = derived_data + "exons_wise/ss_exons_0.6/"
aaseq_exons_dir = derived_data + "exons_wise/aaseq_exons/"
pfam_der_writer = derived_data + "domains/pfam/"
cath_der_writer = derived_data + "domains/cath/"
faaDir = os.path.join(derived_data,'structure_data/fasta_file_pdb')
lisDir = os.path.join(derived_data,'structure_data/res_no_pdb')
faaDirAF = os.path.join(derived_data,'structure_data/fasta_file_pdbAF')
lisDirAF = os.path.join(derived_data,'structure_data/res_no_pdbAF')

'''
Source data needs to be set
organism (9606)/
    source_data/
        pfam/
        SSP/
        Refseq_protein/
        ensemble_data/ #ensembleProteins(single)
        swissprot_data/ 
        gene_tables/ #geneTables
'''
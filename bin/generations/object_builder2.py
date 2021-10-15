import cPickle as pickle
import os
import re
from progress.bar import Bar
from shutil import rmtree
import subprocess
import sys

import modules.retriever as retriever
import modules.prediction_assignment_ss as prediction_assignment_ss
import modules.mappings_refseqProtein as mappings_refseqProtein
import modules.pfam_files as pfam_files
import modules.cath_files_multi as cath_files_multi
import modules.ncbia_raw_coods as ncbia_raw_coods
import modules.structure_refseq as structure_refseq
import modules.multifasta_to_fasta as multifasta_to_fasta
from Classes_exons import Gene as Gene
from Classes_exons import Exon as Exon
from Classes_exons import Transcript as Transcript

if len(sys.argv) != 2:
    print "enter the corret arguements, object_builder.py followed by aruguments.txt\n, please modify arguements.txt file"
    with open("arguments.txt", "w") as fin:
        fin.write("""OrganismFolder:\tOrganismID:\tPDB_source_dir:\tCommon_files:\t
        \nCores:\tstartafresh:1or0\tanalyse:yes\tanalysis_scripts_dir:dir
        # Dont_add_space_after_column
        # freshStart values: 0 means to start from last run state with those variables pickled
        # 1 means to delete the pickles and fresh start
        # 2 means to delet whole derived data and results directory and start from only source files
        """)
    sys.exit()

def part(a, b):
    print ("->Part:%s, Section:%s" % (a, b))

def makedir(dirname):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

def file_pickle_dumper(has, filename):
    with open(filename, "w") as fin:
        pickle.dump(has, fin, protocol=pickle.HIGHEST_PROTOCOL)

def file_loader(filename):
    with open(filename) as fin:
        has = pickle.load(fin)
    return has


'''
Create a organism txid folder
in this keep following in following folders and just pass the addrress of that folder
results and rest files will be created automatically
ask for organism taxid and address of that folder folder having mutual files
create files recursively
Download ensemble data as well as source data
'''


def results_dir_analysis(makeOrDestroy):
    if makeOrDestroy == 'Destroy':
        if os.path.isdir(results_dir):
            rmtree(results_dir)
            results_dir_analysis('make')
    else:
        makedir(results_dir+"csv/domains/pfam")
        makedir(results_dir+"csv/domains/cath")
        makedir(results_dir+'csv/general/RAW')
        makedir(results_dir+"csv/ncbi")
        makedir(results_dir+'csv/ss')
        makedir(results_dir+'csv/ss0.6')
        makedir(results_dir+'csv/stride')
        makedir(results_dir+"plots/pfam")
        makedir(results_dir+"plots/cath")
        makedir(results_dir+'plots/general')
        makedir(results_dir+"plots/ncbi")
        makedir(results_dir+'plots/ss')
        makedir(results_dir+'plots/ss0.6')
        makedir(results_dir+'plots/stride')

with open(sys.argv[1]) as fin:
    raw_arg_fileLines = fin.read().split("\n")
    print (raw_arg_fileLines)
    line_1 = raw_arg_fileLines[0].split()
    parent_dir = line_1[0].split(":")[1]
    print (parent_dir)
    #sys.exit()
    organism = line_1[1].split(":")[1]
    struct_repo_dir = line_1[2].split(":")[1]
    struct_repo_dirAF = os.path.join(os.path.abspath(os.path.dirname(os.path.normpath(struct_repo_dir))),'Alphafold')
    common_data = line_1[3].split(":")[1]
    outdir1=line_1[4].split(":")[1]
    print raw_arg_fileLines[1].split()[0].split(":")
    cores = int(raw_arg_fileLines[1].split()[0].split(":")[1])
    freshStart = int(raw_arg_fileLines[1].split()[1].split(":")[1])
    #analysis_needed = True if raw_arg_fileLines[1].split()[2].split(":")[1] in ['Yes', 'Y', 'yes', 'y'] else False
    analysis_needed = False
    analysis_dir = raw_arg_fileLines[1].split()[3].split(":")[1]
    if not os.path.isdir(analysis_dir):
        print "analysis_dir doesnt exists"
        sys.exit()
    organism_fil = parent_dir.split("/")[-2]

int_patt = re.compile(r'^\d+$')
build_structure = 1

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
makedir(lisDir)
makedir(faaDir)
makedir(lisDirAF)
makedir(faaDirAF)

makedir(pickle_add)
makedir(results_dir)
makedir(save_structure_dir)
makedir(save_structure_dirAF)
makedir(stride_reformat_dir)
makedir(stride_reformat_dirAF)

makedir(pfam_der_writer)
makedir(cath_der_writer)
makedir(ss_exons_dir)
makedir(ss_exons_dir06)
makedir(stride_exons_dir)
makedir(stride_exons_dirAF)

makedir(aaseq_exons_dir)
makedir(pfam_der_writer)
makedir(cath_der_writer)
makedir(results_dir)
results_dir_analysis('make')

if freshStart:
    def remove_file(directory):
        if os.path.isdir(directory):
            rmtree(directory)
    if freshStart == 2:
        remove_file(derived_data)
        remove_file(swiss_pid)
        makedir(swiss_pid)
        remove_file(ens_pid)
        makedir(ens_pid)
        remove_file(pid_add)
        makedir(pid_add)
        makedir(pickle_add)
        results_dir_analysis('Destroy')
        print "Data_cleaned"
        multifasta_to_fasta.refseqP_write_seq(pid_raw_multifasta, pid_add)
        multifasta_to_fasta.swissprot_write_seq(swiss_raw_multifasta, swiss_pid)
        multifasta_to_fasta.ensemble_write_seq(ens_raw_multifasta, ens_pid)
    elif freshStart == 1:
        remove_file(pickle_add)
        makedir(pickle_add)
    else:
        print "wrong fresh start argument, please enter 1,2 or 0"
        sys.exit()

#source fils below
cath_domall_file = os.path.join(common_data,"cath-domain-description-file.txt")
sift_db_file = os.path.join(common_data,"pdb_chain_uniprot.tsv.gz")
#-> should i change it too ? (sift_db_file)
gene2ens_file = os.path.join(common_data,"gene2ensembl.gz")
gene2go_file = os.path.join(common_data,"gene2go.gz")
tabdelim_alias = os.path.join(common_data,"aliases/%s.txt" % organism)
cx = os.path.join(common_data,"uniprot_cx_%s.tab.gz" % organism)
#-> going to change the above, to customized file

# softwares
stride_location = '/home/paras/bin/./stride'
align_location = '/home/paras/bin/./align'

print ('Done')
#sys.exit()
def retConFor(func, file_int, args):
    # retConFor=returner_control_forwarder
    # args here is the list of parameters to be passed aheaddto the functions_flow
    has = file_loader(file_int) if os.path.isfile(file_int) else func(*args)
    if not os.path.isfile(file_int):
        file_pickle_dumper(has, file_int)
    return has
    # wow i used to write good codes too ( :), sep 13 2021) 

def retConForDual(func, file_int1, file_int2, args):
    # retConFor=returner_control_forwarder
    # args here is the list of parameters to be passed aheaddto the functions_flow
    if os.path.isfile(file_int1) and os.path.isfile(file_int2):
        has1 = file_loader(file_int1)
        has2 = file_loader(file_int2)
        return has1, has2
    else:
        has1, has2 = func(*args)
        file_pickle_dumper(has1, file_int1)
        file_pickle_dumper(has2, file_int2)
        return has1, has2

def functions_imp():
    # os.remove(pickle_add + "has_gene_var_%s.pick" % organism)
    has_gene_var = retConFor(retriever.hash_gene_var_list, os.path.join(pickle_add, "has_gene_var_%s.pick" % organism), [gene_add, pid_add])
    part(1, "has_gene_var_reltionship")
    has_NCBI_Ensemble_gene = retConFor(mappings_refseqProtein.ensemble_ncbi_genes, os.path.join(pickle_add,"has_NCBI_Ensemble_gene_%s.pick" % organism), [gene2ens_file, organism])
    part(2, "mapping refseq to swissprot and ensemble")
    has_refseq_ens, has_refseq_swiss = retConForDual(mappings_refseqProtein.refseqEnsemble, os.path.join(pickle_add,"has_refseq_ens_%s.pick" % organism), os.path.join(pickle_add,"has_refseq_swiss_%s.pick" % organism),[pid_add, ens_pid, swiss_pid])
    part(3, "gene description")
    has_gene_description = retConFor(mappings_refseqProtein.gene_description,
                                     os.path.join(pickle_add,"has_gene_decription_short_%s.pick" % organism), [gene_add])
    part(4, "Gene ontology")
    has_gene_go, has_gene_components = retConForDual(mappings_refseqProtein.go_terms_and_cellular_location, os.path.join(pickle_add,"has_gene_gonumbers_%s.pick" % organism), os.path.join(pickle_add,"has_gene_components_%s.pick" % organism),[gene2go_file, organism])
    part(5, "storing exonic coordinates from files")
    has_gene_cood = retConFor(retriever.hash_gene_coordinates, os.path.join(pickle_add,"has_gene_cood_%s.pick" % organism),[gene_add, pid_add])
    part(6, "pfam section")
    pfam_fil_has = retConFor(pfam_files.pfam_runner, os.path.join(pickle_add,"pfam_final_has_%s.pick" %organism), [pfam_der_writer,common_data, pickle_add, pfam_dir, has_gene_var])
    
    part(7, "seq to exons")
    prediction_assignment_ss.aaseq(has_gene_cood, pid_add, aaseq_exons_dir)
    has_pi = retConFor(prediction_assignment_ss.principal_isoform, os.path.join(pickle_add,"has_pi_%s.pick" % organism), [has_gene_var, aaseq_exons_dir])
    makedir(parent_dir+"unassigned_ss/")
    part(8, "ss to exons and assignment")
    #-> uncommneted them below
    '''
    prediction_assignment_ss.assigner(ssp_raw, pid_add, pid_add, ssp_raw, organism, os.path.join(parent_dir,"unassigned_ss"))
    prediction_assignment_ss.ss_to_exons(has_gene_cood, ssp_raw, ss_exons_dir, 0)
    prediction_assignment_ss.ss_to_exons(has_gene_cood, ssp_raw, ss_exons_dir06, 0.6)
    '''
    part(9, "cath_domains")
    cath_fil_has = retConFor(cath_files_multi.cath_runner,os.path.join(pickle_add,"cath_final_has_%s.pick" % organism), [pickle_add, cath_domall_file, cath_der_writer, sift_db_file, swiss_pid, struct_repo_dir, has_refseq_swiss, has_gene_var,common_data, align_location, cores])
    part(10, "structures")
    #-> keep this the same, domains should be assigned based on only the coordinates of original pdb files (do they have to)
    return has_gene_var, has_pi, has_NCBI_Ensemble_gene, has_refseq_ens, has_refseq_swiss, has_gene_description, has_gene_cood, has_gene_go, has_gene_components, pfam_fil_has, cath_fil_has

def functions_flow(has_gene_var, has_gene_cood):
    # has_gene_cood = retConFor(retriever.hash_gene_coordinates,pickle_add+"has_gene_cood_%s.pick" % organism,[gene_add, pid_add])
    # prediction_assignment_ss.assigner(ssp_raw, pid_add, pid_add, ssp_raw)
    # prediction_assignment_ss.ss_to_exons(has_gene_cood, ssp_raw, ss_exons_dir)
    if not os.path.isfile(os.path.join(pickle_add,"has_var_pdbinfo_%s.pick" % organism)):
        structure_refseq.horse(tabdelim_alias, cx, has_gene_var, faaDir, lisDir,struct_repo_dir, save_structure_dir, pid_add, align_location, cores)
        #->this is important, however this only creates files having the name of the vanriant, so create one more file here
        prediction_assignment_ss.prerequisite_stride(stride_location, has_gene_var, save_structure_dir, stride_reformat_dir, pid_add)
        # -> this needs to be tweaked and 
        prediction_assignment_ss.stride_assigner( has_gene_cood, stride_reformat_dir, stride_exons_dir)
        
        # -> (latest) first tsructures were stored and then their stride files like NP_000025.1.ss_3 and then their variant info like in exons in last
    else:
        print (os.path.join(pickle_add,"has_var_pdbinfo_%s.pick" % organism), " this file exists, so not reforming structures, \n delete this file for refreshing strcutres")
    has_var_pdbinfo = retConFor(mappings_refseqProtein.proteins_info_var, os.path.join(pickle_add,"has_var_pdbinfo_%s.pick" % organism),[save_structure_dir])

    if not os.path.isfile(os.path.join(pickle_add,"has_var_pdbinfoAF_%s.pick" % organism)):
    #if 1:
        pass
        # structure_refseq.AF_generator(has_gene_var,has_refseq_swiss, faaDirAF, lisDirAF,struct_repo_dirAF, save_structure_dirAF, pid_add, align_location, cores)
        # #->this is important, however this only creates files having the name of the vanriant, so create one more file here
        # prediction_assignment_ss.prerequisite_stride(stride_location, has_gene_var, save_structure_dirAF, stride_reformat_dirAF, pid_add)
        # #-> this needs to be tweaked and 
        # prediction_assignment_ss.stride_assigner( has_gene_cood, stride_reformat_dirAF, stride_exons_dirAF)
        
        # -> (latest) first tsructures were stored and then their stride files like NP_000025.1.ss_3 and then their variant info like in exons in last
    else:

        print (os.path.join(pickle_add,"has_var_pdbinfoAF_%s.pick" % organism), " this file exists, so not reforming structures, \n delete this file for refreshing strcutres")
    has_var_pdbinfoAF = retConFor(mappings_refseqProtein.proteins_info_var, os.path.join(pickle_add,"has_var_pdbinfoAF_%s.pick" % organism),[save_structure_dirAF])
    
    return has_var_pdbinfo, has_var_pdbinfoAF


def stride_sspred(var, template, ss_ex_dir, stride_ex_dir, AF_exon_dir):
    template.sort()
    # template has coding exons
    # print var, template, len(template)
    # print "insdies"
    ss = {}
    st = {}
    stAF = {}
    # print "hello"
    # try:
    ssFile=os.path.join(ss_ex_dir,"%s" % var)
    stFile=os.path.join(stride_ex_dir,"%s" % var)
    stAFfile=os.path.join(AF_exon_dir,"%s" % var)
    if os.path.isfile(ssFile):

        with open(ssFile) as fin:
            dat = {i[1]: i[0] for i in pickle.load(fin)}
            # print "dat", len(dat), dat, " insitridesspred"
        for i in template:
            ss[i] = dat[i]

    else:
        for i in range(0, len(template)):
            ss[template[i]] = "NULL"
    # print "ss", len(ss), ss
    if os.path.isfile(stFile):
        with open(stFile) as fin:
            dat = {i[1]: i[0] for i in pickle.load(fin)}
        for i in template:
            st[i] = dat[i]
    else:
        for i in range(0, len(template)):
            st[template[i]] = "NULL"

    if os.path.isfile(stAFfile):
        with open(stAFfile) as fin:
            dat = {i[1]: i[0] for i in pickle.load(fin)}
        for i in template:
            stAF[i] = dat[i]
    else:
        for i in range(0, len(template)):
            stAF[template[i]] = "NULL"
    return ss, st, stAF


def redef(ssp, sst, sstAF, exids):
    # print "exids", exids
    ssp_2 = {}
    sst_2 = {}
    sstAF_2 = {}
    exid_2 = {}
    for i in exids:
        if i not in exid_2:
            exid_2[i] = exids[i]
    # KEY IS SPAN AND VALUE IS ID
    for i in ssp:
        ssp_2[exid_2[i]] = ssp[i]
    for i in sst:
        sst_2[exid_2[i]] = sst[i]
    for i in sstAF:
        sstAF_2[exid_2[i]] = sstAF[i]
    
    # print ssp_2, sst_2
    return ssp_2, sst_2, sstAF_2


# functions_flow()
has_gene_var, has_pi, has_NCBI_Ensemble_gene, has_refseq_ens, has_refseq_swiss, \
    has_gene_description, has_gene_cood, has_gene_go, has_gene_components, pfam_fil_has, cath_fil_has = functions_imp()
# sys.exit()
has_var_pdbinfo, has_var_pdbinfoAF  = functions_flow(has_gene_var, has_gene_cood)

print ('done')
#sys.exit()

def refineFcategory(has):
    for gene in has:
        exhas={}
        Tcount=len(has[gene].transcripts)
        for transcripts in has[gene].transcripts:
            for exons in transcripts.exons:
                exID=exons.ID
                if exID[0]!='R':
                    ele=exID.split('.')
                    seqNumber=ele[3]
                    flag=ele[2]
                    if flag!='G':    
                        if seqNumber not in exhas:
                            exhas[seqNumber]=[]
                        exhas[seqNumber]+=[transcripts.ID]
        changeNeeded={i:0 for i in exhas if len(set(exhas[i]))==Tcount}
        if changeNeeded:
            for transcripts in has[gene].transcripts:
                for exons in transcripts.exons:
                    exID=exons.ID
                    if exID[0]!='R':
                        ele=exID.split('.')
                        if ele[3] in changeNeeded and ele[2]!='G':
                            ele[2]='F'
                            exons.ID=".".join(ele)
            changedHas={}
            for exons in has[gene].exons:
                exID=exons.ID
                if exID[0]!='R':
                    ele=exID.split('.')
                    if ele[3] in changeNeeded and ele[2]!='G':
                        ele[2]='F'
                        changedHas[exons.ID]=".".join(ele)
                        exons.ID=".".join(ele)
            for exons in has[gene].exons:
                exID=exons.ID
                if exID[0]=='R':
                    for changeDFrom in changedHas:
                        if changeDFrom in exID:
                            exons.ID=re.sub(r'%s'%(exID),'%s'%(changedHas[changeDFrom]))            
    return has


def gene_object_creator(gene_ob_name, ss_location):
    # if not os.path.isfile(gene_ob_name):
    if 1:
        has_objects = {}
        c = 0
        # print "before main program"
        bar = Bar('Processing genes', max=len(has_gene_var))
        for gene in has_gene_var:
            try:
                bar.next()
                c += 1
                if gene:
                    # 37, 103, 27, 199
                    # exids = ncbia_events.ncbia(has_gene_var[gene], has_pi[gene], aaseq_exons_dir)
                    exids = ncbia_raw_coods.ncbia(
                        has_gene_var[gene], has_pi[gene], aaseq_exons_dir)
                    #continue
                    #print (1)
                    #sys.exit()
                    has_exons = {}
                    oGens = has_NCBI_Ensemble_gene[gene] if gene in has_NCBI_Ensemble_gene else False
                    oLoc = (has_gene_components[gene],
                            True) if gene in has_gene_components else (False, False)
                    # ((0-1),Bool) if gene in hash, else, (False, False) 1 is membranous
                    oGo = has_gene_go[gene] if gene in has_gene_go else False
                    oDes = has_gene_description[gene] if gene in has_gene_description else False
                    geneObject = Gene(organism, gene, oGens, oLoc, oDes, oGo)
                    translisGene = []
                    #print "22"
                    '''
                    # {((1812, 2108), 'AAPPPPVLMHHGESSQVLHPGNKVTLTCVAPLSGVDFQLRRGEKELLVPRSSTSPDRIFFHLNAVALGDGGHYTCRYRLHDNQNGWSGDSAPVELILSD'): ['500', 1, ((1813, 2109),)]
                    '''
                    # print has_gene_var[gene]
                    for var in has_gene_var[gene]:
                        # if var == "NP_116573.1":
                        # print var
                        with open(aaseq_exons_dir+"%s" % var) as fin:
                            dat = pickle.load(fin)
                        # print "dat", dat
                        # print var
                        ssp, sst, sstAF = stride_sspred(
                            var, [i[1] for i in dat if i[3] == "C"], ss_location, stride_exons_dir, stride_exons_dirAF)
                        #print (23)
                        # -> needs a change, include sstAF
                        '''
                        sending here vague id's are givfing vage rsuklts
                        send only hash of raw_coods and correspodning genomic id's
                        '''
                        upd_exids = {
                            mexo[1]: exids[(mexo[1], mexo[0], mexo[3])][0] for mexo in dat}
                        # print "updexids", upd_exids
                        ssp, sst, sstAF = redef(ssp, sst, sstAF, upd_exids)
                        # print "redef", ssp
                        #print (24)
                        '''
                            # create and complete the transcript object here only, and add it to the list one by one
                        '''
                        # print "::33::"
                        opistat = True if has_pi[gene] == var else False
                        olentrans = len("".join([i[0]
                                                 for i in dat if i[3] == 'C']))
                        oens = has_refseq_ens[var] if var in has_refseq_ens else False
                        oswiss = has_refseq_swiss[var] if var in has_refseq_swiss else False
                        #print ":2"
                        transObject = Transcript(
                            var, opistat, olentrans, oens, oswiss)
                        # print "1"
                        # oStructure, oStruclen = has_var_pdbinfo[var][0], has_var_pdbinfo[var][1] if var in has_var_pdbinfo else False, False
                        # print has_var_pdbinfo[var][0]
                        oStructure = has_var_pdbinfo[var][0] if var in has_var_pdbinfo else False
                        oStruclen = has_var_pdbinfo[var][1] if var in has_var_pdbinfo else False
                        oStructureAF = has_var_pdbinfoAF[var][0] if var in has_var_pdbinfoAF else False
                        oStruclenAF = has_var_pdbinfoAF[var][1] if var in has_var_pdbinfoAF else False
                        # print oStructure, oStruclen
                        #print "2"
                        if oStructure:
                            transObject.pdb_info(oStructure, oStruclen)
                        if oStructureAF:
                            #print (oStructureAF)
                            transObject.pdb_infoAF(oStructureAF, oStruclenAF)
                        transexlis = []
                        # print dat, "dat"
                        # print ssp
                        for miniexons in dat:
                            # print miniexons, "me"
                            non_coding_span = miniexons[1]
                            coding_span = miniexons[2]
                            aaseq = miniexons[0]
                            exStatus = miniexons[3]
                            tupsql = (non_coding_span, aaseq, exStatus)
                            exon_id = exids[tupsql][0]
                            
                            if exon_id not in has_exons:
                                has_exons[exon_id] = Exon(
                                    id=exon_id, aaseq=aaseq, coding_span=coding_span, nc_span=non_coding_span)
                            
                            #print ('e1')
                            has_exons[exon_id].second_seq(
                                ssp[exon_id], var) if exon_id in ssp else has_exons[exon_id].second_seq(None, var)
                            #print ('e2')
                            has_exons[exon_id].stride_seq(
                                sst[exon_id], var) if exon_id in sst else has_exons[exon_id].stride_seq(None, var)
                            #print ('e3')
                            #print (sstAF)
                            has_exons[exon_id].stride_seqAF(
                                sstAF[exon_id], var) if exon_id in sstAF else has_exons[exon_id].stride_seqAF(None, var)
                            #print('**Outaa Here')
                            transexlis += [has_exons[exon_id]]
                        transObject.exons(transexlis)
                        transObject.intron_status()
                        #print ('e4')
                        if gene in pfam_fil_has:
                            if var in pfam_fil_has[gene]:
                                transObject.pfam_update(
                                    pfam_fil_has[gene][var])
                        if gene in cath_fil_has:
                            if var in cath_fil_has[gene]:
                                transObject.cath_update(
                                    cath_fil_has[gene][var])
                        translisGene += [transObject]
                    geneObject.transcripts(translisGene)
                    geneObject.exons(has_exons.values())
                    geneObject.constitutive_alternate_with_freq()
                    
                    has_objects[gene] = geneObject
                    

            except Exception as E:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                print "ERROR: main()", gene, E, exc_tb.tb_lineno
                break
            
        print ("in main program")
        
        has_objects=refineFcategory(has_objects)
        with open(gene_ob_name, "wb") as fin:
           pickle.dump(has_objects, fin, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        print "gene object %s is already exisiting" % gene_ob_name


# gene_object_creator(gene_ob_name=results_dir+"objectsave_%s.pick" %
#                     organism, ss_location=ss_exons_dir)
# part(12, "Normal_processing_done")
gene_object_creator(gene_ob_name=os.path.join(results_dir,"objectsave_%s_0.6.pick" %
                   organism), ss_location=ss_exons_dir06)
part(13, "helix_confidence_processing_done")

if analysis_needed and 0:
    results_dir_analysis('make')
    '''
    subprocess.check_output(['python', analysis_dir+'general_stats_p1.py', results_dir, organism])
    subprocess.check_output(['python', analysis_dir+'NCBI_file_writer_core.py', results_dir, organism])
    subprocess.check_output(['python', analysis_dir+'figureNCBandA.py', results_dir, organism])
    subprocess.check_output(['python', analysis_dir+'figure_ss_0.6.py', results_dir, organism])
    subprocess.check_output(['python', analysis_dir+'figure_ss.py', results_dir, organism])
    subprocess.check_output(['python', analysis_dir+'figure_stride.py', results_dir, organism])
    subprocess.check_output(['python', analysis_dir+'domain_section.py', results_dir, organism, 'pfam'])
    subprocess.check_output(['python', analysis_dir+'domain_section.py', results_dir, organism, 'cath'])

    os.system(analysis_dir+'general_stats_p1.py %s %s' % (results_dir, organism))
    os.system(analysis_dir+'NCBI_file_writer_core.py %s %s' % (results_dir, organism))
    os.system(analysis_dir+'figureNCBandA.py %s %s' % (results_dir, organism))
    os.system(analysis_dir+'figure_ss_0.6.py %s %s' % (results_dir, organism))
    os.system(analysis_dir+'figure_ss_new.py %s %s' % (results_dir, organism))
    os.system(analysis_dir+'figure_stride.py %s %s' % (results_dir, organism))
    os.system(analysis_dir+'domain_section.py %s %s cath' % (results_dir, organism))
    os.system(analysis_dir+'domain_section.py %s %s pfam' % (results_dir, organism))

    os.system('general_stats_p1.py %s %s' % (results_dir, organism))
    os.system('NCBI_file_writer_core.py %s %s' % (results_dir, organism))
    os.system('figureNCBandA.py %s %s' % (results_dir, organism))
    os.system('figure_ss_0.6.py %s %s' % (results_dir, organism))
    os.system('figure_ss_new.py %s %s' % (results_dir, organism))
    os.system('figure_stride.py %s %s' % (results_dir, organism))
    os.system('domain_section.py %s %s cath' % (results_dir, organism))
    os.system('domain_section.py %s %s pfam' % (results_dir, organism))
    '''
    with open("%s/temp_bash_%s.sh" % (analysis_dir, organism_fil), "w") as fin:
        fin.write("#!/bin/bash")
        fin.write('''
        python %s/general_stats_p1.py %s %s
        python %s/NCBI_file_writer_core.py %s %s
        python %s/figureNCBandA.py %s %s
        python %s/figure_ss_0.6.py %s %s
        python %s/figure_ss_new.py %s %s
        python %s/figure_stride.py %s %s
        python %s/domain_section.py %s %s pfam
        python %s/domain_section.py %s %s cath
        ''' % (analysis_dir, results_dir, organism,
               analysis_dir, results_dir, organism,
               analysis_dir, results_dir, organism,
               analysis_dir, results_dir, organism,
               analysis_dir, results_dir, organism,
               analysis_dir, results_dir, organism,
               analysis_dir, results_dir, organism,
               analysis_dir, results_dir, organism))
    subprocess.check_output(
        ['chmod', 'a+x', '%s/temp_bash_%s.sh' % (analysis_dir, organism_fil)])
    out = subprocess.check_output(
        ['%s/./temp_bash_%s.sh' % (analysis_dir, organism_fil)])
    # os.remove('temp_bash.sh')
    with open("%s/temp_plots_%s.sh" % (analysis_dir, organism_fil), "w") as fin:
        fin.write("#!/bin/bash\n")
        fin.write("""
        Rscript /home/paras/project/protein_splicing/scripts/plots/pssm.R %s/csv/stride/ %s/plots/stride/
        Rscript /home/paras/project/protein_splicing/scripts/plots/N.R %s/csv/ncbi/ %s/plots/ncbi/
        Rscript /home/paras/project/protein_splicing/scripts/plots/B.R %s/csv/ncbi/ %s/plots/ncbi/
        Rscript /home/paras/project/protein_splicing/scripts/plots/C.R %s/csv/ncbi/ %s/plots/ncbi/
        Rscript /home/paras/project/protein_splicing/scripts/plots/density_plot.R %s/csv/ncbi/ %s/plots/ncbi/
        Rscript /home/paras/project/protein_splicing/scripts/plots/density_plot_CDF.R %s/csv/ncbi/ %s/plots/ncbi/
        Rscript /home/paras/project/protein_splicing/scripts/plots/domains_new.R %s/csv/domains/pfam/ %s/plots/pfam/
        Rscript /home/paras/project/protein_splicing/scripts/plots/domains_new.R %s/csv/domains/cath/ %s/plots/cath/
        Rscript /home/paras/project/protein_splicing/scripts/plots/general_fig1.R %s/csv/general/ %s/plots/general/
        Rscript /home/paras/project/protein_splicing/scripts/plots/pssm.R %s/csv/ss/ %s/plots/ss/
        Rscript /home/paras/project/protein_splicing/scripts/plots/pssm.R %s/csv/ss0.6/ %s/plots/ss0.6/
        """ % (results_dir, results_dir,
               results_dir, results_dir,
               results_dir, results_dir,
               results_dir, results_dir,
               results_dir, results_dir,
               results_dir, results_dir,
               results_dir, results_dir,
               results_dir, results_dir,
               results_dir, results_dir,
               results_dir, results_dir,
               results_dir, results_dir,
               ))
    subprocess.check_output(
        ['chmod', 'a+x', '%s/temp_plots_%s.sh' % (analysis_dir, organism_fil)])
    # out = subprocess.check_output(['./temp_plots.sh'], shell=True)
part(14, "congratulations, analysis done too")
# print "output_log:%s" % out

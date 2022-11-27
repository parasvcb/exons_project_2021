import cPickle as pickle
import os
import re
from progress.bar import Bar
from shutil import rmtree
import sys
import multiprocessing
import subprocess
import constructing_data.retriever as retriever
import constructing_data.prediction_assignment_ss as prediction_assignment_ss
import constructing_data.mappings_refseqProtein as mappings_refseqProtein
import constructing_data.pfam_files as pfam_files
import constructing_data.cath_files_multi as cath_files_multi
import constructing_data.ncbia_raw_coods_error_free as ncbia_raw_coods
import constructing_data.structure_refseq as structure_refseq
import constructing_data.multifasta_to_fasta as multifasta_to_fasta
from constructing_data.Classes_exons import Transcript as Transcript
from constructing_data.Classes_exons import Exon as Exon
from constructing_data.Classes_exons import Gene as Gene
import object_module as om
'''
Very big file to document completely,
half of its steps are done in prep.bash and other halves are peding, I thiught i didnt knew i had modified them in this file
This file can be segmented in
    0. parsing args
    1. fileStructure generation, reseting it
    2. setting core folders and file names
    3. making sure all of them exixting and have contents (pfam, ssp, disoredr and strcutre regions again)
    4. general mapping, very much NCBI intensive, will see if similar supplementary files can be generated in the Ensembl Build also
    5. making structure build unofficial
    6. improving args inout format
    7. different module and efficient analysis_component


# ensembl filew styructre inside willbe liottle complictaed or should that be done differently ?
yes make it diffrent 
# file naming for object should also be differet and name of file and direrctories should also be different


'''

currCmdString = """
parent_dir:/home/paras/project/protein_splicing/10090/ #OrganismFolder
organism:10090 #txid
module: 0 # 0 for nom, 1, for nom + ss + dom+ disorder, 2 prev + strcutre data aspects 
struct_repo_dir:/media/CSB/pdbs_splicing/new_pdb/ #PDB_source_dir
common_data:/home/paras/project/protein_splicing/common_files/
freshStart:0
structureandCath:-1 #want to analyse structre, 1 yes, -1 no
analysis_needed:yes
analysis_dir:/home/paras/project/protein_splicing/scripts/analysis/

#cores used would be cpu count -4 , change in main program
"""

# ############################################################# <Part 0> ##############################################################################

if len(sys.argv) != 2:
    print ("enter the corret arguements, object_builder.py followed by aruguments.txt\n, please modify arguements.txt file")
    with open("arguments.txt", "w") as fin:
        fin.write(currCmdString)
    sys.exit()

params = {}
with open(sys.argv[1]) as fin:
    raw_arg_fileLines = [i for i in fin.read().split(
        "\n") if len(i) > 3 and i[0] != '#']

    for i in raw_arg_fileLines:
        key, value = i.split(':')
        key = key.strip()
        value = value.split("#")[0].strip()
        params[key] = value
    
    freshStart = int(params['freshStart'])

    parent_dir = params['parent_dir']
    outdir1 = parent_dir
    organism = params['organism']
    common_data = params['common_data']
    analysis_needed = int(params['analysis_needed'])
    analysis_dir = params['analysis_dir']

    sslevel = float(params['sslevel']) if params['sslevel'] else False
    
    struct_repo_dir = params['struct_repo_dir']
    struct_repo_dirAF = os.path.join(os.path.abspath(
        os.path.dirname(os.path.normpath(struct_repo_dir))), 'Alphafold')
    structureandCath = int(params['structureandCath'])
    # if above is true, then only handle them otherwise skip for now
    
    cores = multiprocessing.cpu_count() - 4
"""
This is the module, there should be two/three aspects of this, 
source, is that NCBI or ensembl ?
    Basic functionality (nomewnclaure only, )
    DB functionality (nomenclature, ss, disorder, pfam)
    advaned functionality (above and cath, pdb structure, alphfaold structre and plDDT significance)
"""

# ############################################################# <Part 0/> ##############################################################################

# ############################################################# <Part 1> ##############################################################################
# setting names of the directorues and eseting
source_data, pfam_dir, ssp_raw, pid_add, pid_raw_multifasta, ens_raw_multifasta, swiss_raw_multifasta, ens_pid, swiss_pid, gene_add, derived_data, disp_raw, pickle_add, results_dir, save_structure_dir, save_structure_dirAF, stride_reformat_dir, stride_exons_dir, stride_reformat_dirAF, stride_exons_dirAF, pfam_der_writer, cath_der_writer, ss_exons_dir, ss_exons_dir06, dis_exons_dir, aaseq_exons_dir, pfam_der_writer, cath_der_writer, faaDir, lisDir, faaDirAF, lisDirAF = om.setnames(parent_dir, analysis_needed)
#-> this should also be updated to has dir level


if freshStart:
    def remove_file(directory):
        if os.path.isdir(directory):
            om.rmtree(directory)
    if freshStart == 2:
        remove_file(derived_data)
        remove_file(swiss_pid)
        om.makedir(swiss_pid)
        remove_file(ens_pid)
        om.makedir(ens_pid)
        remove_file(pid_add)
        om.makedir(pid_add)
        om.makedir(pickle_add)
        results_dir_analysis('Destroy')
        print ("Data_cleaned")
        multifasta_to_fasta.refseqP_write_seq(pid_raw_multifasta, pid_add)
        multifasta_to_fasta.swissprot_write_seq(
            swiss_raw_multifasta, swiss_pid)
        multifasta_to_fasta.ensemble_write_seq(ens_raw_multifasta, ens_pid)
    elif freshStart == 1:
        remove_file(pickle_add)
        makedir(pickle_add)
    else:
        print ("wrong fresh start argument, please enter 1,2 or 0")
        sys.exit()

results_dir_analysis('make')
# -> apparently there is reduinda cuy in code aboea nd in prep.bash segment, i think both needs to underho some change, keep it preent here and emov from the prep.bash
# ############################################################# <Part 1/> ##############################################################################

# ############################################################# <Part 2> ##############################################################################
# SOURCE FILES, DIRECTORIES AND PROGRAM REFERENCE


cath_domall_file, sift_db_file, gene2ens_file, gene2go_file, tabdelim_alias, cx, stride_location, align_location = om.setsourcefnames (common_data)
# -> as i will be shifting sooner to mapping on the laphafold or colved collabfold entries, there is not point acallingthe akign program, nor their is any need to update with teh gene2ens files

# ############################################################# <Part 2/> ##############################################################################


# ############################################################# <Part 3> ##############################################################################
# CHECKING IF THE MAJOR FILES AND ORTGANIZAION IS PRESENT OR NOT,
# A) CORE FILES (GENE TABLES, PROTEINS, ENSEMBL AND SWISSPROT PROTEOMES, GENE2ENSEMBL AND MORE FILES
# B) PREDICTIONS DONE (DOMAINS (PFAM), SSP, DISORDER
# C) STRUCTURE CROSS MAPPING
# APPROACH IS PENDING, SHOULD REFER TO PROGRAMS USED FOR SUCH PROGRAMS IN DIFFERNT SYSTEMS AND THEIR DOCUMENTATION ALSO
# ############################################################# <Part 3/> ##############################################################################

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


def functions_imp(structureandCath):
    '''
    ->
    ENSEMBL: 
    has_gene_var is easy,
    -has_NCBI_Ensemble_gene, do i need it, havent used anywhere so far but could be used
    *has_refseq_ens, *has_refseq_swiss, call the ensembl mapper file,
    has_gene_description can also be done
    has_gene_go, has_gene_components can also be done, match how they have and ensembl has listed the values 
    has_gene_cood can also be done
    pfam section is different
    ss section is different
    
    -> module division could be crerated easily here
    '''
    # os.remove(pickle_add + "has_gene_var_%s.pick" % organism)
    has_gene_var = retConFor(retriever.hash_gene_var_list, os.path.join(
        pickle_add, "has_gene_var_%s.pick" % organism), [gene_add, pid_add])
    part(1, "has_gene_var_reltionship")
    has_NCBI_Ensemble_gene = retConFor(mappings_refseqProtein.ensemble_ncbi_genes, os.path.join(
        pickle_add, "has_NCBI_Ensemble_gene_%s.pick" % organism), [gene2ens_file, organism])
    part(2, "mapping refseq to swissprot and ensemble")
    # above can be skipped
    has_refseq_ens, has_refseq_swiss = retConForDual(mappings_refseqProtein.refseqEnsemble, os.path.join(
        pickle_add, "has_refseq_ens_%s.pick" % organism), os.path.join(pickle_add, "has_refseq_swiss_%s.pick" % organism), [pid_add, ens_pid, swiss_pid])
    # with bigger pool, more memory will be required, not a very efficinet way,
    part(3, "gene description")
    has_gene_description = retConFor(mappings_refseqProtein.gene_description,
                                     os.path.join(pickle_add, "has_gene_decription_short_%s.pick" % organism), [gene_add])
    part(4, "Gene ontology")
    has_gene_go, has_gene_components = retConForDual(mappings_refseqProtein.go_terms_and_cellular_location, os.path.join(
        pickle_add, "has_gene_gonumbers_%s.pick" % organism), os.path.join(pickle_add, "has_gene_components_%s.pick" % organism), [gene2go_file, organism])
    # its a hash of the go numbers and has having values, where value close to one is membranous and opposite, checka ny alanysis section to see how that was done, (likely in gi1 generat stats)

    part(5, "storing exonic coordinates from files")
    has_gene_cood = retConFor(retriever.hash_gene_coordinates, os.path.join(
        pickle_add, "has_gene_cood_%s.pick" % organism), [gene_add, pid_add])

    part(6, "pfam section")

    pfam_fil_has = retConFor(pfam_files.pfam_runner, os.path.join(pickle_add, "pfam_final_has_%s.pick" % organism), [pfam_der_writer, common_data, pickle_add, pfam_dir, has_gene_var])
    # pfam_fil_has = []
    part(7, "seq to exons")
    prediction_assignment_ss.aaseq(has_gene_cood, pid_add, aaseq_exons_dir)
    has_pi = retConFor(prediction_assignment_ss.principal_isoform, os.path.join(
        pickle_add, "has_pi_%s.pick" % organism), [has_gene_var, aaseq_exons_dir])
    makedir(parent_dir + "unassigned_ss/")
    part(8, "ss to exons and assignment")
    # -> uncommneted them below

    prediction_assignment_ss.assigner(
        ssp_raw, pid_add, pid_add, ssp_raw, os.path.join(parent_dir, "unassigned_ss"))

    # -> fucntion equivalent arg names, (done_dir, fastadir_done, pool_of_left, dir_to_write, where_to_write_rep_of_unassigned) screens the ssp_raw for already done prediction, and stores aaseq iun memory (from fastadirDone), screens pool_ofleft(pid_add: refseq dir) and check aaseq if predicted, if doesnt then ask for the prediction and write to unassigned_ss dir

    prediction_assignment_ss.ss_to_exons(
        has_gene_cood, ssp_raw, ss_exons_dir, None)
    # -> write sspred to exons format (matrix) and pickles them, threshold used is for confidence
    prediction_assignment_ss.ss_to_exons(
        has_gene_cood, ssp_raw, ss_exons_dir06, 0.6)
    prediction_assignment_ss.disorder_to_exons(
        has_gene_cood, disp_raw, dis_exons_dir)
    # don't worry about the its being done twice, onl;y on eth eidrectory wikl be used below

    part(9, "cath_domains")
    if structureandCath:
        cath_fil_has = retConFor(cath_files_multi.cath_runner, os.path.join(pickle_add, "cath_final_has_%s.pick" % organism), [
            pickle_add, cath_domall_file, cath_der_writer, sift_db_file, swiss_pid, struct_repo_dir, has_refseq_swiss, has_gene_var, common_data, align_location, cores])
    else:
        cath_fil_has = []
    # -> keep this the same, domains should be assigned based on only the coordinates of original pdb files (do they have to)
    return has_gene_var, has_pi, has_NCBI_Ensemble_gene, has_refseq_ens, has_refseq_swiss, has_gene_description, has_gene_cood, has_gene_go, has_gene_components, pfam_fil_has, cath_fil_has


def functions_flow(has_gene_var, has_gene_cood, structureandCath):
    if structureandCath:
        if not os.path.isfile(os.path.join(pickle_add, "has_var_pdbinfo_%s.pick" % organism)):
            structure_refseq.horse(tabdelim_alias, cx, has_gene_var, faaDir, lisDir,
                                   struct_repo_dir, save_structure_dir, pid_add, align_location, cores)
            # ->this is important, however this only creates files having the name of the vanriant, so create one more file here
            prediction_assignment_ss.prerequisite_stride(
                stride_location, has_gene_var, save_structure_dir, stride_reformat_dir, pid_add)
            # -> this needs to be tweaked and
            prediction_assignment_ss.stride_assigner(
                has_gene_cood, stride_reformat_dir, stride_exons_dir)

            # -> (latest) first tsructures were stored and then their stride files like NP_000025.1.ss_3 and then their variant info like in exons in last
        else:
            print (os.path.join(pickle_add, "has_var_pdbinfo_%s.pick" % organism),
                   " this file exists, so not reforming structures, \n delete this file for refreshing strcutres")
        has_var_pdbinfo = retConFor(mappings_refseqProtein.proteins_info_var, os.path.join(
            pickle_add, "has_var_pdbinfo_%s.pick" % organism), [save_structure_dir])

        if not os.path.isfile(os.path.join(pickle_add, "has_var_pdbinfoAF_%s.pick" % organism)):
            # if 1:
            pass
            # structure_refseq.AF_generator(has_gene_var,has_refseq_swiss, faaDirAF, lisDirAF,struct_repo_dirAF, save_structure_dirAF, pid_add, align_location, cores)
            # #->this is important, however this only creates files having the name of the vanriant, so create one more file here
            # prediction_assignment_ss.prerequisite_stride(stride_location, has_gene_var, save_structure_dirAF, stride_reformat_dirAF, pid_add)
            # #-> this needs to be tweaked and
            # prediction_assignment_ss.stride_assigner( has_gene_cood, stride_reformat_dirAF, stride_exons_dirAF)

            # -> (latest) first tsructures were stored and then their stride files like NP_000025.1.ss_3 and then their variant info like in exons in last
        else:

            print (os.path.join(pickle_add, "has_var_pdbinfoAF_%s.pick" % organism),
                   " this file exists, so not reforming structures, \n delete this file for refreshing strcutres")
        has_var_pdbinfoAF = retConFor(mappings_refseqProtein.proteins_info_var, os.path.join(
            pickle_add, "has_var_pdbinfoAF_%s.pick" % organism), [save_structure_dirAF])
    else:
        has_var_pdbinfo = []
        has_var_pdbinfoAF = []
    return has_var_pdbinfo, has_var_pdbinfoAF




# functions_flow()
has_gene_var, has_pi, has_NCBI_Ensemble_gene, has_refseq_ens, has_refseq_swiss, \
    has_gene_description, has_gene_cood, has_gene_go, has_gene_components, pfam_fil_has, cath_fil_has = functions_imp(
        structureandCath)
# sys.exit()
has_var_pdbinfo, has_var_pdbinfoAF = functions_flow(
    has_gene_var, has_gene_cood, structureandCath)

# structureandCath if not 1, then var will be empty/false cath_fil_has, has_var_pdbinfo, has_var_pdbinfoAF (last two will be t7urn false, no problem at all), same for the former
print ('done')

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
                # if gene == 319701:
                if gene:
                    # for var in has_gene_cood[gene]:
                    #     for exon in has_gene_cood[gene][var]:
                    #         print (var, exon)
                    # # 37, 103, 27, 199
                    # exids = ncbia_events.ncbia(has_gene_var[gene], has_pi[gene], aaseq_exons_dir)
                    exids = ncbia_raw_coods.ncbia(
                        has_gene_var[gene], has_pi[gene], aaseq_exons_dir)
                    has_exons = {}
                    oGens = has_NCBI_Ensemble_gene[gene] if gene in has_NCBI_Ensemble_gene else False
                    oLoc = (has_gene_components[gene],
                            True) if gene in has_gene_components else (False, False)
                    # ((0-1),Bool) if gene in hash, else, (False, False) 1 is membranous
                    oGo = has_gene_go[gene] if gene in has_gene_go else False
                    oDes = has_gene_description[gene] if gene in has_gene_description else False
                    geneObject = Gene(organism, gene, oGens, oLoc, oDes, oGo)
                    translisGene = []
                    # print "22"
                    '''
                    # {((1812, 2108), 'AAPPPPVLMHHGESSQVLHPGNKVTLTCVAPLSGVDFQLRRGEKELLVPRSSTSPDRIFFHLNAVALGDGGHYTCRYRLHDNQNGWSGDSAPVELILSD'): ['500', 1, ((1813, 2109),)]
                    '''
                    # print has_gene_var[gene]
                    for var in has_gene_var[gene]:
                        # if var == "NP_116573.1":
                        with open(aaseq_exons_dir + "%s" % var) as fin:
                            dat = pickle.load(fin)
                        ssp, diso, sst, sstAF = stride_sspred_disorder(
                            var, [i[1] for i in dat if i[3] == "C"], ss_location, dis_exons_dir, stride_exons_dir, stride_exons_dirAF)
                        '''
                        sending here vague id's are givfing vage rsuklts
                        send only hash of raw_coods and correspodning genomic id's
                        '''
                        upd_exids = {
                            mexo[1]: exids[(mexo[1], mexo[0], mexo[3])][0] for mexo in dat}
                        ssp, diso, sst, sstAF = redef(
                            ssp, diso, sst, sstAF, upd_exids)
                        '''
                            # create and complete the transcript object here only, and add it to the list one by one
                        '''
                        opistat = True if has_pi[gene] == var else False
                        olentrans = len("".join([i[0]
                                                 for i in dat if i[3] == 'C']))
                        oens = has_refseq_ens[var] if var in has_refseq_ens else False
                        oswiss = has_refseq_swiss[var] if var in has_refseq_swiss else False
                        transObject = Transcript(
                            var, opistat, olentrans, oens, oswiss)
                        oStructure = has_var_pdbinfo[var][0] if var in has_var_pdbinfo else False
                        oStruclen = has_var_pdbinfo[var][1] if var in has_var_pdbinfo else False
                        oStructureAF = has_var_pdbinfoAF[var][0] if var in has_var_pdbinfoAF else False
                        oStruclenAF = has_var_pdbinfoAF[var][1] if var in has_var_pdbinfoAF else False
                        if oStructure:
                            transObject.pdb_info(oStructure, oStruclen)
                        if oStructureAF:
                            transObject.pdb_infoAF(oStructureAF, oStruclenAF)
                        transexlis = []
                        for miniexons in dat:
                            
                            non_coding_span = miniexons[1]
                            coding_span = miniexons[2]
                            aaseq = miniexons[0]
                            exStatus = miniexons[3]
                            tupsql = (non_coding_span, aaseq, exStatus)
                            exon_id = exids[tupsql][0]

                            if exon_id not in has_exons:
                                has_exons[exon_id] = Exon(
                                    id=exon_id, aaseq=aaseq, coding_span=coding_span, nc_span=non_coding_span)

                            has_exons[exon_id].second_seq(
                                ssp[exon_id], var) if exon_id in ssp else has_exons[exon_id].second_seq(None, var)
                            has_exons[exon_id].disord_seq(
                                diso[exon_id], var) if exon_id in diso else has_exons[exon_id].disord_seq(None, var)
                            has_exons[exon_id].stride_seq(
                                sst[exon_id], var) if exon_id in sst else has_exons[exon_id].stride_seq(None, var)
                            has_exons[exon_id].stride_seqAF(
                                sstAF[exon_id], var) if exon_id in sstAF else has_exons[exon_id].stride_seqAF(None, var)

                            transexlis += [has_exons[exon_id]]
                            # print (has_exons)
                        transObject.exons(transexlis)
                        transObject.intron_status()
                        # print ('e4')
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
                print ("ERROR: main()", gene, E, exc_tb.tb_lineno)
                break

        print ("in main program")
        with open(gene_ob_name, "wb") as fin:
            pickle.dump(has_objects, fin, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        print ("gene object %s is already exisiting") % gene_ob_name


# gene_object_creator(gene_ob_name=results_dir+"objectsave_%s.pick" %
#                     organism, ss_location=ss_exons_dir)
# part(12, "Normal_processing_done")

strss = '0.6' if sslevel == 0.6 else ''
ss_location = ss_exons_dir06 if sslevel == 0.6 else ss_exons_dir
strstr = 'withStructure' if structureandCath else 'noStructure'
print (strss, ss_location, strstr, organism)
gene_ob_name = os.path.join(
    results_dir, "objectsave_%s_%s_%s.pick" % (organism, strss, strstr))

gene_object_creator(gene_ob_name=gene_ob_name, ss_location=ss_location)
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
    with open("%s/temp_bash_%s.sh" % (analysis_dir, organism), "w") as fin:
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
        ['chmod', 'a+x', '%s/temp_bash_%s.sh' % (analysis_dir, organism)])
    out = subprocess.check_output(
        ['%s/./temp_bash_%s.sh' % (analysis_dir, organism)])
    # os.remove('temp_bash.sh')
    with open("%s/temp_plots_%s.sh" % (analysis_dir, organism), "w") as fin:
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
        ['chmod', 'a+x', '%s/temp_plots_%s.sh' % (analysis_dir, organism)])
    # out = subprocess.check_output(['./temp_plots.sh'], shell=True)
part(14, "congratulations, analysis done too")
# print "output_log:%s" % out

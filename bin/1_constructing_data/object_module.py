
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

# ############################################################# <Part 1> ##############################################################################
# DIRECTORY STRUCTURE GENERATION AND RESETTING IT


def results_dir_analysis(makeOrDestroy):
    if makeOrDestroy == 'Destroy':
        if os.path.isdir(results_dir):
            rmtree(results_dir)
            results_dir_analysis('make')
    else:
        makedir(results_dir + "csv/domains/pfam")
        makedir(results_dir + "csv/domains/cath")
        makedir(results_dir + 'csv/general/RAW')
        makedir(results_dir + "csv/ncbi")
        makedir(results_dir + 'csv/ss')
        makedir(results_dir + 'csv/ss0.6')
        makedir(results_dir + 'csv/stride')
        makedir(results_dir + "plots/pfam")
        makedir(results_dir + "plots/cath")
        makedir(results_dir + 'plots/general')
        makedir(results_dir + "plots/ncbi")
        makedir(results_dir + 'plots/ss')
        makedir(results_dir + 'plots/ss0.6')
        makedir(results_dir + 'plots/stride')

int_patt = re.compile(r'^\d+$')

def setnames(parent_dir, analysis_needed):
    source_data = os.path.join(parent_dir, "source_data")
    pfam_dir = os.path.join(source_data, "pfam")
    ssp_raw = os.path.join(source_data, "SSP")
    pid_add = os.path.join(source_data, 'Refseq_protein')
    pid_raw_multifasta = os.path.join(source_data, '%s_ref.faa' % organism)
    ens_raw_multifasta = os.path.join(source_data, '%s_ens.faa.gz' % organism)
    swiss_raw_multifasta = os.path.join(source_data, '%s_swiss.faa.gz' % organism)
    ens_pid = os.path.join(source_data, 'ensemble_data')
    swiss_pid = os.path.join(source_data, 'swissprot_data')
    gene_add = os.path.join(source_data, 'gene_tables')

    derived_data = os.path.join(outdir1, "derived_data/")
    disp_raw = os.path.join(derived_data, "disorder_data/")
    pickle_add = os.path.join(derived_data, "pickles")
    results_dir = os.path.join(derived_data, "results")
    save_structure_dir = os.path.join(
        derived_data + 'structure_data/PDB_structure_derived_data/')
    save_structure_dirAF = os.path.join(
        derived_data + 'structure_data/PDB_structure_derived_dataAF/')

    # -> added this time
    stride_reformat_dir = os.path.join(derived_data, "stride_pssm_form/")
    stride_exons_dir = os.path.join(derived_data, "exons_wise/stride_exons/")
    stride_reformat_dirAF = os.path.join(derived_data, "stride_pssm_formAF/")
    stride_exons_dirAF = os.path.join(derived_data, "exons_wise/stride_exonsAF/")

    pfam_der_writer = os.path.join(derived_data, "domains/pfam/")
    cath_der_writer = os.path.join(derived_data, "domains/cath/")
    ss_exons_dir = os.path.join(derived_data, "exons_wise/ss_exons/")
    ss_exons_dir06 = os.path.join(derived_data, "exons_wise/ss_exons_0.6/")

    dis_exons_dir = os.path.join(derived_data, "exons_wise/dis_exons/")
    aaseq_exons_dir = os.path.join(derived_data, "exons_wise/aaseq_exons/")
    pfam_der_writer = os.path.join(derived_data, "domains/pfam/")
    cath_der_writer = os.path.join(derived_data, "domains/cath/")
    faaDir = os.path.join(derived_data, 'structure_data/fasta_file_pdb')
    lisDir = os.path.join(derived_data, 'structure_data/res_no_pdb')
    faaDirAF = os.path.join(derived_data, 'structure_data/fasta_file_pdbAF')
    lisDirAF = os.path.join(derived_data, 'structure_data/res_no_pdbAF')

    makedir(lisDir)
    makedir(faaDir)
    makedir(lisDirAF)
    makedir(faaDirAF)
    if analysis_needed:
        makedir(analysis_dir)

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
    makedir(dis_exons_dir)
    makedir(stride_exons_dir)
    makedir(stride_exons_dirAF)

    makedir(aaseq_exons_dir)
    makedir(pfam_der_writer)
    makedir(cath_der_writer)
    makedir(results_dir)
    
    names = [source_data, pfam_dir, ssp_raw, pid_add, pid_raw_multifasta, ens_raw_multifasta,
    swiss_raw_multifasta, ens_pid, swiss_pid, gene_add, derived_data ,
    disp_raw, pickle_add, results_dir, save_structure_dir,
    save_structure_dirAF,
    stride_reformat_dir, stride_exons_dir, stride_reformat_dirAF, stride_exons_dirAF,
    pfam_der_writer, cath_der_writer, ss_exons_dir, ss_exons_dir06,
    dis_exons_dir, aaseq_exons_dir, pfam_der_writer, cath_der_writer,
    faaDir, lisDir, faaDirAF, lisDirAF
    ]
    return names 


results_dir_analysis('make')


# ############################################################# <Part 2> ##############################################################################
# SOURCE FILES, DIRECTORIES AND PROGRAM REFERENCE

def setsourcefnames():
    cath_domall_file = os.path.join(common_data, "cath-domain-description-file.txt")
    sift_db_file = os.path.join(common_data, "pdb_chain_uniprot.tsv.gz")
    # -> should i change it too ? (sift_db_file)
    gene2ens_file = os.path.join(common_data, "gene2ensembl.gz")
    gene2go_file = os.path.join(common_data, "gene2go.gz")
    tabdelim_alias = os.path.join(common_data, "aliases/%s.txt" % organism)
    cx = os.path.join(common_data, "uniprot_cx_%s.tab.gz" % organism)
    # -> going to change the above, to customized file
    # softwares
    stride_location = '/home/paras/bin/./stride'
    align_location = '/home/paras/bin/./align'

    return cath_domall_file, sift_db_file, gene2ens_file, gene2go_file, tabdelim_alias, cx, stride_location, align_location

# ############################################################# <Part 2/> ##############################################################################

def stride_sspred_disorder(var, template, ss_ex_dir, dis_ex_dir, stride_ex_dir, AF_exon_dir):
    template.sort()
    # template has coding exons
    ss = {}
    st = {}
    stAF = {}
    dis = {}
    ssFile = os.path.join(ss_ex_dir, "%s" % var)
    disfile = os.path.join(dis_ex_dir, "%s" % var)
    stFile = os.path.join(stride_ex_dir, "%s" % var)
    stAFfile = os.path.join(AF_exon_dir, "%s" % var)

    def mini_wroking_Code(filename):
        ob = {}
        if os.path.isfile(filename):
            with open(filename) as fin:
                dat = {i[1]: i[0] for i in pickle.load(fin)}
            for i in template:
                ob[i] = dat[i]

        else:
            for i in range(0, len(template)):
                ob[template[i]] = "NULL"
        return ob

    ss = mini_wroking_Code(ssFile)
    dis = mini_wroking_Code(disfile)
    st = mini_wroking_Code(stFile)
    stAF = mini_wroking_Code(stAFfile)
    return ss, dis, st, stAF


def redef(ssp, dis, sst, sstAF, exids):
    # print "exids", exids
    ssp_2 = {}
    dis_2 = {}
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
    for i in dis:
        dis_2[exid_2[i]] = dis[i]
    return ssp_2, dis_2, sst_2, sstAF_2

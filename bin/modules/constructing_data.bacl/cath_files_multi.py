import ProDomStruc_multiprocess as pmb
import os
import cPickle as pickle
import re
import numpy as np



def readPickle(fname):
    with open (fname) as fin:
        has=pickle.load(fin)
    return has

def writePickle (fname,has):
    with open(fname, "w") as fin:
        pickle.dump(has, fin, protocol=pickle.HIGHEST_PROTOCOL)


# download_pdb_flag = 1


def feeder(gene, var, sp, has, coods, cov, topol, domid):
    if gene not in has:
        has[gene] = {}
        has[gene][var] = [(domid, coods, cov, topol)]
    else:
        if var not in has[gene]:
            has[gene][var] = [(domid, coods, cov, topol)]
        else:
            has[gene][var] += [(domid, coods, cov, topol)]
    return has


def sorter(has):
    child_has = {}
    for i in has:
        for j in has[i]:
            if i not in child_has:
                child_has[i] = {}
            child_has[i][j] = list(set(has[i][j]))
    return child_has


def preparation(pickle_add, domall_file, domall_seqres_format, sift, cases_refined, pdbdir, spdir, humanspshash, has_gene_var,
                download_pdb_flag, cases_before_align, cases_after_align, error_logs, cores, common_data, align_add):
    print ('in cath preparation')
    if not os.path.isfile(os.path.join(pickle_add,"cath_domall_coods.pick")) or not os.path.isfile(os.path.join(pickle_add,"cath_domall_topology.pick")):
        domall_hash, topology = pmb.domall2seqres(
            domall=domall_file, filsave=domall_seqres_format)
        # will write to file in this format
        '''
        #domoainID      PDBCood
        4uk9F02 122,172
        4uk9F03 175,245
        4uk9F01 2,117
        4o8rA00 5,224
        5k1zA00 2,230
        2vhaB01 3,111:208,279
        '''
        writePickle(os.path.join(pickle_add,"cath_domall_coods.pick"),domall_hash)
        writePickle(os.path.join(pickle_add,"cath_domall_topology.pick"),topology)

    else:
        domall_hash=readPickle(os.path.join(pickle_add,"cath_domall_coods.pick"))
        topology=readPickle(os.path.join(pickle_add,"cath_domall_topology.pick"))

    print (os.path.isfile(cases_after_align), cases_after_align)
    if not os.path.isfile(cases_after_align):
        print 'True'
        pmb.cath(domall_hash, sift, pdbdir, spdir, humanspshash,
                 download_pdb_flag, cases_before_align, cases_after_align, error_logs, cores, common_data, align_add)
    '''
    format of below will be like this
    "#Uniprot\tDomID\tCoordinates\tIdentity\tUnip_Coods
    '''

    with open(cases_after_align) as fin:
        dat = [i.split("\t") for i in fin.read().split(
            "\n") if len(i) > 10 and i[0] != '#']

    humansps = humanspshash
    if not os.path.isfile(cases_refined):
        has_merged = {}
        for i in dat:
            if i[0] not in has_merged:
                has_merged[i[0]] = [i[1:]]
            else:
                has_merged[i[0]] += [i[1:]]
        # print has_merged.keys()[:5]
        cases_with_unip = 0
        cases_single = 0
        cases_morethan_one = 0
        gene_rif = {}
        with open(cases_refined, "w") as fin:
            fin.write(
                "#Gene\tVar\tSP\tDomID\tIdentity\tPDBCoods\tUnip_Coods\tPseudo_cov\tTopology\n")
            for i in has_gene_var:
                for j in has_gene_var[i]:
                    if j in humansps:
                        temp_a = humansps[j].split(",")
                        len_fac = len(temp_a)
                        if (len_fac) == 1:
                            cases_single += 1
                        if (len_fac) > 1:
                            cases_morethan_one += 1
                        for k in temp_a:
                            k = k.split(".")[0]
                            if k in has_merged:
                                gene_rif[i] = 0
                                cases_with_unip += 1
                                for l in has_merged[k]:
                                    fin.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                                        i, j, k, l[0], l[2], l[1], l[3], l[4], topology[l[0]]))
        print len(gene_rif), "no of genes"
        print cases_with_unip, cases_single, cases_morethan_one
        print len(has_merged)
    if not os.path.isfile(pickle_add,"cath_all_domains_hash.pick"):
        with open(cases_refined) as fin:
            dat = [i.split("\t") for i in fin.read().split(
                "\n") if len(i) > 10 and i[0] != '#']
        has_one_dom_per_gene_all = {}
        '''
        has[gene]=[(uniprot,ncbi,[[intcoods]],topology)]
        '''
        # Gene\tVar\tSP\tLength\tPI\tDomID\tIdentity\tPDBCoods\tUnip_Coods\tTopology\n"
        for i in dat:
            gene = int(i[0])
            uniprot = i[2]
            domid = i[3]
            ncbi = i[1]
            coods = i[6]
            cov = i[7]
            topol = i[8]
            has_one_dom_per_gene_all = feeder(
                gene, ncbi, uniprot, has_one_dom_per_gene_all, coods, cov, topol, domid)

        # print has_one_dom_per_gene_continuous[16]
        has_gene_all_dom = sorter(has_one_dom_per_gene_all)

        with open(os.path.join(pickle_add,"cath_all_domains_hash.pick"), "w") as fin:
            pickle.dump(has_gene_all_dom, fin,
                        protocol=pickle.HIGHEST_PROTOCOL)
        return has_gene_all_dom
    else:
        with open(os.path.join(pickle_add,"cath_all_domains_hash.pick")) as fin:
            has_gene_all_dom = pickle.load(fin)
        return has_gene_all_dom


def overlapper(has):
    child_has = {}
    for gene in has:
        for var in has[gene]:
            longest = (0, {-1}, 0)
            doms = []
            # setlists = []
            for iter1 in has[gene][var]:
                if 'undef' not in iter1[1]:
                    span1 = iter1[1]
                    val = set()
                    for span in span1.split(":"):
                        ele = map(int, (span.split(",")))
                        val |= set(range(ele[0], ele[1]+1))
                    doms += [(len(val), val, iter1)]
                    if len(val) > longest[0]:
                        longest = (len(val), val, iter1)
            # got the longest domain
            if doms:
                setlists = [longest]
                # print setlists
                # print doms
                # print i
                for val in doms:
                    counter = 0
                    for dom in setlists:
                        intersec = dom[1] & val[1]
                        if intersec:
                            # print intersec, "intersec"
                            # print dom, len(dom)
                            if float(len(intersec))/dom[0] < 0.25:
                                counter += 1
                            else:
                                # if domain is larger ?
                                if val[0] > dom[0]:
                                    setlists = [val if ix ==
                                                dom else ix for ix in setlists]
                        else:
                            counter += 1
                    if counter == len(setlists):
                        setlists += [val]
                if (0, {-1}, 0) in setlists:
                    setlists.remove((0, {-1}, 0))
                if gene not in child_has:
                    child_has[gene] = {}
                if var not in child_has[gene]:
                    # print setlists
                    child_has[gene][var] = [(tuple(map(int, re.findall(
                        r'[-]?\d+', i[2][1]))), i[2][0], tuple(map(float, re.findall(r'[-]?\d+\.\d+', i[2][2]))), i[2][3]) for i in setlists]
                else:
                    child_has[gene][var] += [(tuple(map(int, re.findall(
                        r'[-]?\d+', i[2][1]))), i[2][0], tuple(map(float, re.findall(r'[-]?\d+\.\d+', i[2][2]))), i[2][3]) for i in setlists]
    return child_has


def np_to_csv(arr, filename, cols):
    values = arr[0]
    bins = arr[1]
    su = sum(values)
    with open(filename+".csv", "w") as fin:
        fin.write("%s\n" % (",".join(cols)))
        for i in range(0, len(bins)-1):
            fin.write("%s,%s,%s\n" %
                      ("-".join(map(str, [round(bins[i], 2), round(bins[i+1], 2)])), values[i], round(float(values[i])/su, 3)))


def summary(all_doms, total_has_merged, filename):
    with open(filename, "w") as fin:
        fin.write(
            "Gene\tVariant\tSP\tPI\tDomain_count\tDomain_length_all\tProtein_length\tFold\n")
        for i in all_doms:
            gene = i
            variant = total_has_merged[i][3]
            sp = total_has_merged[i][2]
            pi = total_has_merged[i][0]
            domain_count = len(all_doms[i])
            summary_length = set()
            for j in all_doms[i]:
                summary_length |= j
            domain_length = len(summary_length)
            Protein_length = total_has_merged[i][1]
            fold = ":".join(total_has_merged[i][5])
            fin.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                      (gene, variant, sp, pi, domain_count, domain_length, Protein_length, fold))


# res_dir_cath = "/home/paras/mysql/project1/jupyter/results/"


def domdistribution(arr, ranges):
    np_arr = np.histogram(arr, bins=ranges)
    np_to_csv(np_arr, res_dir_cath+"csv/cath/Domains_cath_distribution",
              cols=["No.OfDomains", "Total", "Freq"])


def cath_runner(pickle_add, domall_file, derived_data_add, sift, sp_dir, pdbdir, humanspshash, has_gene_var,common_data, align_add, cores):
    '''
    pickle_add, domall_file, domall_seqres_format, sift, cases_refined, pdbdir, spdir, humanspshash, has_gene_var,
                    download_pdb_flag, cases_before_align, cases_after_align, error_logs
    '''
    print ('main cath')

    download_pdb_flag = 0
    domall_seqres_format = os.path.join(derived_data_add,"cath_domain_all_human_readable.txt")

    cases_before_align = os.path.join(derived_data_add,"structural_domains_uniprot_predicted_multi.txt")
    cases_after_align = os.path.join(derived_data_add,"structural_domains_uniprot_assigned_multi.txt")
    error_logs = os.path.join(derived_data_add,"cath_structural_domains_multicases_error.txt")
    cases_refined = os.path.join(derived_data_add,"cath_structural_domains_uniprot_gene_var_CX.txt")
    # $results_dir = "/home/paras/mysql/project1/res/"

    total_has_merged = preparation(pickle_add, domall_file, domall_seqres_format, sift, cases_refined, pdbdir, sp_dir, humanspshash, has_gene_var, download_pdb_flag, cases_before_align, cases_after_align, error_logs, cores, common_data, align_add)
    all_doms = overlapper(total_has_merged)
    return all_doms
    # t_doms_all = [len(all_doms[i]) for i in all_doms]
    '''
    domdistribution(t_doms_all, range(min(t_doms_all), max(t_doms_all)+1, 1))
    summary(all_doms, total_has_merged, filename=results_dir +
            "csv/cath/cath_all_domain_summary_statistics.log")
    summary(cont_doms, cont_has, filename=results_dir +
            "csv/cath/cath_continuous_domain_summary_statistics.log")
    summary(discont_doms, discont_has, filename=results_dir +
            "csv/cath/cath_discontinuous_domain_summary_statistics.log")
    # print sum([len(cont_doms[i]) for i in cont_doms])
    # print sum([len(discont_doms[i]) for i in discont_doms])

    # res_dir_cath = "/home/paras/mysql/project1/jupyter/results/"
    ProDomStruc_multiprocess.discontinuous_domains(
        3, condition_genes, discont_doms, discont_has, source_gene_object, res_dir_cath + "csv/cath/")
    ProDomStruc_multiprocess.fraction_domain_ind(
        3, condition_genes, cont_doms, cont_has, source_gene_object, res_dir_cath + "csv/cath/")
    ProDomStruc_multiprocess.fraction_domain_pair(
        3, condition_genes, cont_doms, cont_has, source_gene_object, res_dir_cath + "csv/cath/")
    '''

import os
import cPickle as pickle
import gzip
import re
import numpy as np
import ProDomStruc_multiprocess

def preparation(src_dir, common_data, pickle_add, has_gene, pfam_dir):
    if not os.path.isfile(src_dir+"pfam_processed_core.flat.gz") or not os.path.isfile(os.path.join(pickle_add,"pfam_domains_hash.pick")):
        if not os.path.isfile(os.path.join(pickle_add,"pfam_length_hmm.pick")):
            #print "1inside"
            has_pfam_hmm = {}
            with open(os.path.join(common_data,"Pfam-A.hmm")) as fin:
                dat = fin.read()
                dat2 = dat.split("//")
                for i in dat2:
                    if len(i) > 10:
                        acc = re.search(r'^ACC.*$', i, re.MULTILINE)
                        leng = re.search(r'^LENG.*$', i, re.MULTILINE)
                        if acc and leng:
                            # print acc.group(0), leng.group(0)
                            has_pfam_hmm[acc.group(0).split()[1]] = int(
                                leng.group(0).split()[1])
                        else:
                            print i[:250]
            print len(has_pfam_hmm)
            with open(os.path.join(pickle_add,"pfam_length_hmm.pick"), 'w') as fin:
                pickle.dump(has_pfam_hmm, fin,
                            protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(os.path.join(pickle_add,"pfam_length_hmm.pick")) as fin:
                has_pfam_hmm = pickle.load(fin)
        with gzip.open(src_dir+"pfam_processed_core.flat.gz", "wb") as fout1:
            fout1.write(
                "Gene\tNCBIVar\tPfamID\tCoverage\tstart\tend\thmmacc\ttype\tE-val\tClan\n")
            pfam_all = {}
            # print has_gene, "has_gene"
            for gene in has_gene:
                # print gene, type(gene)
                for var in has_gene[gene]:
                    # print pfam_dir, var

                    with open(os.path.join(pfam_dir,"%s" % var)) as fin:
                        dat = fin.read()
                    res_temp = re.findall(
                        r'^[XNY]P.*$', dat, flags=re.MULTILINE)
                    if res_temp:
                        for j in res_temp:
                            ele = j.split()
                            if float(ele[-3]) <= 0.01:
                                start = ele[3]
                                end = ele[4]
                                lenpfam = (int(end)-int(start))+1
                                pfam_id = ele[5]
                                cov = round(float(lenpfam) /
                                            has_pfam_hmm[pfam_id], 2)
                                hmmtype = ele[6]
                                dtype = ele[7]
                                evalue = ele[-3]
                                clan = ele[-1]
                                coods = (int(start), int(end))
                                fout1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                                            (gene, var, pfam_id, cov, start, end, hmmtype, dtype, evalue, clan))
                                if gene not in pfam_all:
                                    pfam_all[gene] = {}
                                if var not in pfam_all[gene]:
                                    pfam_all[gene][var] = []
                                pfam_all[gene][var] += [(pfam_id, coods, cov, hmmtype, dtype)]
                                

        with open(os.path.join(pickle_add,"pfam_domains_hash.pick"), "w") as fin:
            pickle.dump(pfam_all, fin, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(os.path.join(pickle_add,"pfam_domains_hash.pick")) as fin:
            pfam_all = pickle.load(fin)
        # print pfam_all[19]
    return pfam_all


def overlapper(has):
    child_has = {}
    for gene in has:
        if gene:
            for var in has[gene]:
                # if [1 for i in has[gene][var] if (i[1][1]-i[1][0]+1) > 39:
                if 1:
                    #print gene, var
                    #print has[gene][var]
                    longest = (0, 0, 0)
                    #these three values will correspond to length of domain, domainspanset, and domain itself
                    domsori = []
                    setlists = []
                    '''
                    iterate the domains, for every domain add to the list of original domains as [lengthofdomain, its spanrange in set and, domain], 
                    and save the longest sepeartelt
                    '''
                    for iter1 in has[gene][var]:
                        span1 = iter1[1]
                        val = set(range(span1[0], span1[1]+1))
                        domsori += [(len(val), val, iter1)]
                    domsori.sort()
                    longest=tuple(domsori[-1])
                    doms = domsori[:]
                    
                    '''
                    for every domain in the doms
                        compare it with setlist (which initially  contains only the l;argest domain)
                            if there is no span intersection and if it is and is less than 0.25% then count this domain in counter
                        check if couter is exactly eual to the total number of the largest domains (setlist), 
                            add this to the domains list
                    '''
                    setlists = [longest]
                    for val in doms:
                        counter = 0
                        for dom in setlists:
                            # print dom
                            # print val
                            intersec = dom[1] & val[1]
                            if intersec:
                                if float(len(intersec))/dom[0] < 0.25:
                                    counter += 1
                            else:
                                counter += 1
                        if counter == len(setlists):
                            setlists += [val]

                    if gene not in child_has:
                        child_has[gene] = {}
                    # pfam_id, coods, hmmtype, dtype)]
                    # child_has[gene][var] = [{'pfam_id': i[0], 'coods':i[1], 'hmmtype': i[2], 'dtype':i[3]} for i in setlists]
                    # print setlists

                    if var not in child_has[gene]:
                        #print setlists
                        child_has[gene][var] = [
                            (i[2][1], i[2][0], i[2][2], i[2][3], i[2][4]) for i in setlists]
                        #easch set list comtains domain iteration, represented by i above, its second elemnt contains complete information in forms of [(pfam_id, coods, cov, hmmtype, dtype)]
                    else:
                        child_has[gene][var] += [
                            (i[2][1], i[2][0], i[2][2], i[2][3], i[2][4]) for i in setlists]
                # break
    return child_has


def summary(all_doms, total_has_merged, filename):
    with open(filename, "w") as fin:
        fin.write(
            "Gene\tVariant\tDomain_count\tDomain_length_all\tProtein_length\tFold\tType\n")
        for i in all_doms:
            gene = i
            variant = total_has_merged[i][0]
            domain_count = len(all_doms[i])
            summary_length = set()
            for j in all_doms[i]:
                summary_length |= j
            domain_length = len(summary_length)
            Protein_length = total_has_merged[i][1]
            fold = ":".join(total_has_merged[i][3])
            dtype = ":".join(total_has_merged[i][4])
            # pfam_all[gene] = [pi, length, [coods], [hmmtype], [dtype]]
            fin.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                      (gene, variant, domain_count, domain_length, Protein_length, fold, dtype))


def domdistribution(arr, ranges, filename):
    np_arr = np.histogram(arr, bins=ranges)
    np_to_csv(np_arr, filename, cols=["No.OfDomains", "Total", "Freq"])


def np_to_csv(arr, filename, cols):
    values = arr[0]
    bins = arr[1]
    su = sum(values)
    with open(filename, "w") as fin:
        fin.write("%s\n" % (",".join(cols)))
        for i in range(0, len(bins)-1):
            fin.write("%s,%s,%s\n" %
                      ("-".join(map(str, [round(bins[i], 2), round(bins[i+1], 2)])), values[i], round(float(values[i])/su, 3)))


def pfam_runner(src_dir,common_data, pickle_add, pfam_dir, has_gene_var):

    # print "gained_control"
    # print src_dir, pickle_add, pfam_dir

    pfam_comprehensive = preparation(
        src_dir, common_data, pickle_add, has_gene_var, pfam_dir)

    # print "out"
    pfam_doms_refined = overlapper(pfam_comprehensive)

    return pfam_doms_refined
    '''
    pfam_doms_all = [len(pfam_doms_refined[i]) for i in pfam_doms_refined]
    domdistribution(pfam_doms_all, [1, 2, 3, 4, 5, 6, 10, 15, 25, 50, 100, 150, 500], filename=results_dir +
                    "csv/pfam/pfam_domain_distribution.csv")
    summary(pfam_doms_refined, pfam_comprehensive, filename=results_dir +
            "csv/pfam/pfam_all_domain_summary_statistics.log")

    ProDomStruc_multiprocess.fraction_domain_ind(
        0, condition_genes, pfam_doms_refined, pfam_comprehensive, source_gene_object, results_dir + "csv/pfam/")
    ProDomStruc_multiprocess.fraction_domain_pair(
        0, condition_genes, pfam_doms_refined, pfam_comprehensive, source_gene_object, results_dir + "csv/pfam/")

    print pfam_comprehensive[26]
    '''

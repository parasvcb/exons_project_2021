import os
import multiprocessing
from Bio import Entrez
import datetime
import re
import gzip
from progress.bar import Bar
from Bio.PDB import PDBList
print ("---> In Retriever Module")


def description():
    print ('''
    ::: This Module has follwoing functions
    1. entrez_retriever(lis, nature, directory)
        Retrieves proteins and genetables from NCBI and stores them on the directory listed
        Uses API keys and nature type should be  ["protein", "mRNA", "gene"]
    2. id_extractorGT(gene_add, nature):
        extracta protein mRNA ids's from gee tables and return them as list
    3. def hash_gene_coordinates(gene_add, pid_add, flatfile=None):
        return hash having gene,var and then value as tuple (exon seq, nc cood, coding cood and coding staus)
    4. def hash_gene_var_list(gene_add, pid_add, flatfile=None):
        return has having gene and list of varinats in side it
    5. def mmcif_pdb_retriever(lis,directory)
    ''')


def hidden_retriever1(lis, dbtype, address, has):
    if dbtype == "gene":
        fmt = "gene_table"
        ext = ".txt"
    else:
        fmt = "fasta"
        ext = ".faa"
    for j in lis:
        if not(os.path.isfile(os.path.join(address,"%s" % j+ext))):
            try:
                Entrez.email = "parasvcb@gmail.com"
                handle = Entrez.efetch(db=dbtype, id=str(
                    j), api_key="d13f309e2f3022c4a996b71b453258233907", rettype=fmt, retmode="text")
                b = handle.read()
                b = b.encode("ascii")
                with open(address+"%s" % j+ext, "w") as fout:
                    fout.write("%s" % b)
            except Exception as E:
                has[j] = E


def hidden_retriever2(lis, dbtype, address, has):
    if dbtype == "gene":
        fmt = "gene_table"
        ext = ".txt"
    else:
        fmt = "fasta"
        ext = ".faa"

    for j in lis:
        if not(os.path.isfile(os.path.join(address,"%s" % j+ext))):
            try:
                Entrez.email = "paras.verma@live.com"
                handle = Entrez.efetch(db=dbtype, id=str(
                    j), api_key="bce08d2b16ea6ff6cf78dcbb16b9ef683a09", rettype=fmt, retmode="text")
                b = handle.read()
                b = b.encode("ascii")
                with open(address+"%s" % (j)+ext, "w") as fout:
                    fout.write("%s" % b)
            except Exception as E:
                has[j] = E


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def entrez_retriever(lis, nature, directory, fname=False):
    if fname:
        with open (fname) as fin:
            lis=[i for i in fin.read().split('\n')]
    if nature not in ["protein", "mRNA", "gene"]:
        return "Exiting: please add second arguemnet as either 'gene','protein' or 'mRNA' only"
    manager = multiprocessing.Manager()
    has = manager.dict()
    cpus = 10
    sublists = list(chunks(lis, cpus))
    processes = []
    for i in range(0, cpus):
        if i < 5:
            # send processes to one api_key
            p = multiprocessing.Process(target=hidden_retriever1, args=(
                sublists[i], nature, directory, has,))
            processes.append(p)
            p.start()
        else:
            p = multiprocessing.Process(target=hidden_retriever2, args=(
                sublists[i], nature, directory, has,))
            processes.append(p)
            p.start()
    for proc in processes:
        proc.join()
    # retrying missed
    has = manager.dict()
    p = multiprocessing.Process(target=hidden_retriever1, args=(lis, nature, directory, has,))
    p.start()
    p.join()
    print ("data_successfully_retrived\n******************")
    if has:
        addr = "/".join(directory.split("/")[:-2])+"/"
        print ("untreieved logs are present in %s" % addr)
        dat_t = str(datetime.datetime.now())
        dat_t = re.sub(r'\:\d+\.\d+$', '', dat_t)
        dat_t = re.sub(r':', '_', dat_t)
        filname = "unretrieved_log_%s_%s.log" % (nature, dat_t)
        with open(addr+"%s" % filname, "w") as fout1:
            for i in has:
                fout1.write("%s:\t%s\n" % i, has[i])

def id_extractorGT(gene_add, nature):
    if nature not in ["protein", "mRNA"]:
        return "Exiting: please add second arguemnet as either 'protein' or 'mRNA' only"
    bar = Bar('Processing gene tables:', max=len(os.listdir(gene_add)))
    pattern_core = re.compile(
        r'[m]?RNA .*\s+.*[XN]P_[\d]+[.][\d]+[\s]?[\(]?[\w]*.[\w]*[\)]?[,]? \b[\d]+\s+coding\s+exon[s]?,.*\s+.*\s+.*\s+[-]+[\d\n\-\s]*')
    # this pattern will search out the transcripts that will hacve portein counterpart, excluding ncRNA and lncRNA
    if nature == "protein":
        mol_pat = re.compile("[WNXY]P_\d+.\d+")
    else:
        mol_pat = re.compile("[WNXY]M_\d+.\d+")
    err = ""
    ids = []
    for f in os.listdir(gene_add):
        try:
            with open(gene_add+"%s" % f, "r") as fin:
                dat2 = fin.read()
            gene_id = f.split(".")[0]
            matches = pattern_core.findall(dat2, re.IGNORECASE)
            # below code wont go for ncRNA genes or genes have only ncRNA
            for i in range(0, len(matches)):
                prots = mol_pat.search(matches[i])
                ids += [prots.group(0)]
        except Exception as e:
            err += "%s\t%s\n" % (gene_id, e)
        bar.next()
    bar.finish()
    if err:
        addr = "/".join(gene_add.split("/")[:-2])+"/"
        print ("Error retrieving logs are present in %s" % addr)
        dat_t = str(datetime.datetime.now())
        dat_t = re.sub(r'\:\d+\.\d+$', '', dat_t)
        dat_t = re.sub(r':', '_', dat_t)
        filname = "error_patt_%s_%s.log" % (nature, dat_t)
        with open(addr+"%s" % filname, "w") as fout1:
            fout1.write("%s" % err)
    return ids

def hash_gene_coordinates(gene_add, pid_add, flatfile=None):
    # print "hello"
    has = {}
    pid_dir_hash = {i: 0 for i in os.listdir(pid_add)}
    bar = Bar('Processing coordinates from gene tables:', max=len(os.listdir(gene_add)))
    pid = [i.strip() for i in os.listdir(pid_add)]
    if len(pid) == 0:
        return "Exiting: empty PID list"

    pattern_core = re.compile(
        r'[m]?RNA .*\s+.*[XN]P_[\d]+[.][\d]+[\s]?[\(]?[\w]*.[\w]*[\)]?[,]? \b[\d]+\s+coding\s+exon[s]?,.*\s+.*\s+.*\s+[-]+[\d\n\-\s]*')
    pattern_exon = re.compile(r'\d+\-\d+.*\n')
    prot = re.compile("[NX]P_\d+.\d+")
    err = ""
    # print len(os.listdir(gene_add))
    for f in os.listdir(gene_add):
        try:
            with open(os.path.join(gene_add,"%s" % f)) as fin:
                dat2 = fin.read()
            gene_id = f.split(".")[0]
            gene_id = int(gene_id)
            matches = pattern_core.findall(dat2, re.IGNORECASE)
            has[gene_id] = {}
            # this is creating the gene as hash
            # below code wont go for ncRNA genes or genes have only ncRNA
            for i in range(0, len(matches)):
                ex_mat = pattern_exon.findall(matches[i], re.IGNORECASE)
                prots = prot.search(matches[i])
                # prots = prots.group(0)
                # print prots
                if prots.group(0) in pid_dir_hash:
                    # no need to extract whats not avialabe as sequence
                    has[gene_id][prots.group(0)] = []  # has to hash on basis of pid
                    for j in range(0, len(ex_mat)):
                        spl_exo = ex_mat[j].split()
                        fil_len = len(spl_exo)
                        if fil_len in [3, 4]:
                            # this is a exon UTR
                            tuple_nc = tuple(map(int, spl_exo[1].split("-")))
                            # in this case, when length is no there,
                            # 1st refersto gene location
                            has[gene_id][prots.group(0)] += [[j, tuple_nc, (0, 0), "N"]]
                        elif fil_len in [6, 7]:
                            # this is a coding exon
                            # gene information is avaialbel 2nd exon
                            tuple_nc = tuple(map(int, spl_exo[2].split("-")))
                            tuple_c = tuple(map(int, spl_exo[3].split("-")))
                            has[gene_id][prots.group(0)] += [[j, tuple_nc, tuple_c, "C"]]
                            # now pushing value [(nccood),(ccood),"N"or"C"]

        except Exception as e:
            err += "%s\t%s\n" % (f, e)
            print (e, f)
        bar.next()
    bar.finish()
    if err:
        addr = "/".join(gene_add.split("/")[:-2])+"/"
        print ("Error extracting coods logs are present in %s" % addr)
        dat_t = str(datetime.datetime.now())
        dat_t = re.sub(r'\:\d+\.\d+$', '', dat_t)
        dat_t = re.sub(r':', '_', dat_t)
        filname = "error_cood_extraction_%s.log" % (dat_t)
        with open(os.path.join(addr,"%s" % filname), "w") as fout1:
            fout1.write("%s" % err)
    if flatfile:
        with gzip.open(flatfile+".gz"+"w") as fin:
            fin.write("Gene\tProtein\tEx_no\t\Non_coding_cood\tCoding_cood\tCoding_status\n")
            for gene in has:
                for protein in has[gene]:
                    cood_tup = has[gene][protein]
                    fin.write("%s\t%s\t%s\t%s\t%s\t%s\n" %
                              (gene, protein, cood_tup[0], cood_tup[1], cood_tup[2], cood_tup[3]))
    return has


def hash_gene_var_list(gene_add, pid_add, flatfile=None):
    #print "pid_add", pid_add
    has = {}
    pid_dir_hash = {i: 0 for i in os.listdir(pid_add)}
    bar = Bar('Processing gene_var_lits from gene tables:', max=len(os.listdir(gene_add)))
    pid = [i.strip() for i in os.listdir(pid_add)]
    if len(pid) == 0:
        return "Exiting: empty PID list"

    pattern_core = re.compile(
        r'[m]?RNA .*\s+.*[XN]P_[\d]+[.][\d]+[\s]?[\(]?[\w]*.[\w]*[\)]?[,]? \b[\d]+\s+coding\s+exon[s]?,.*\s+.*\s+.*\s+[-]+[\d\n\-\s]*')
    pattern_exon = re.compile(r'\d+\-\d+.*\n')
    prot = re.compile("[NX]P_\d+.\d+")
    err = ""
    for f in os.listdir(gene_add):
        if True:
            try:
                with open(os.path.join(gene_add,"%s" % f)) as fin:
                    dat2 = fin.read()
                # print dat2
                gene_id = int(f.split(".")[0])
                matches = pattern_core.findall(dat2, re.IGNORECASE)
                has[gene_id] = []
                # this is creating the gene as hash
                # below code wont go for ncRNA genes or genes have only ncRNA
                for i in range(0, len(matches)):
                    ex_mat = pattern_exon.findall(matches[i], re.IGNORECASE)
                    prots = prot.search(matches[i])
                    #print prots
                    if prots.group(0) in pid_dir_hash:
                        # no need to extract whats not avialabe as sequence
                        has[gene_id] += [prots.group(0)]
                # print prots.group(0)
                # print matches
                # break
            except Exception as e:
                err += "%s\t%s\n" % (f, e)
                print (e)
            bar.next()
    bar.finish()
    if err:
        addr = "/".join(gene_add.split("/")[:-2])+"/"
        print ("Error extracting gene_var_rltnship logs are present in %s" % addr)
        dat_t = str(datetime.datetime.now())
        dat_t = re.sub(r'\:\d+\.\d+$', '', dat_t)
        dat_t = re.sub(r':', '_', dat_t)
        filname = "error_gene_var_extraction_%s.log" % (dat_t)
        with open(addr+"%s" % filname, "w") as fout1:
            fout1.write("%s" % err)
    if flatfile:
        with gzip.open(flatfile+".gz"+"w") as fin:
            fin.write("Gene\tProtein\tEx_no\t\Non_coding_cood\tCoding_cood\tCoding_status\n")
            for gene in has:
                for protein in has[gene]:
                    l=has[gene][protein] #imaginary
                    fin.write("%s\t%s\n" % (gene, protein, l[0], l[1], l[2], l[3]))
    # print has['103']
    has_ref = {i: has[i] for i in has if has[i]}
    return has_ref


def mmcif_pdb_retriever(lis, directory):
    pdbl = PDBList()
    c = 0
    for i in lis:
        print (c, len(lis))
        c += 1
        if not (os.path.isfile(directory+"%s.cif" % i) or os.path.isfile(directory+"pdb%s.ent" % i)):
            try:
                pdbl.retrieve_pdb_file(i, pdir=directory, file_format="mmCif")
            except Exception as E:
                print ("%s: getting file through pdbformat for %s" % (E, i))
                pdbl.retrieve_pdb_file(i, pdir=directory, file_format="pdb")

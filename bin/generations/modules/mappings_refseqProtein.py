import gzip
import re
import os
from progress.bar import Bar
from Bio.PDB import PDBParser, Selection
parser = PDBParser()
print '''
---> In mappings_module
'''


def description():
    print '''
    This module has follwing functions and utility
    1. ensemble_ncbi_genes (flat_file, txid, flatadd=None):
      flat file from ncbi, can gzipped or uncompressed and txid for which relatiosnhip needs to be extracted
    2. def refseqEnsemble(ref, ens, swiss, flataddens=None, flataddswiss=None):
       returns two mappings, ncbi to ensemble and ncbi to swiisprot
       will take 3 inputs major, refseq repostry, sensemble repositury and swissprot repostry
    3. def gene_description(gene_table_dir):
        from gene table reads first line as decsription nand return a hash, with key as gene id and value as gene_description
    4. def go_terms_and_cellular_location(go_term_file, txid):
        reads the go term file and extracts info for requested txid,
        returns two hashes,
        first, having gene id as key and then corresponding list of go numbers
        second, geneid as key and then value, where 1 corresponds membranous protein, 0 as globular and rest in between both
    5. def proteins_info_var(protein_structure_dir):
        returns the hash with key as variant NP and then list having two values, first filename and then length of sequence

    Future improvements :
    will add multiprocessing capabilities and split func 2
    will add gene description mapping, [ADDED]
    go terms and cellular location of genes [ADDED]'''


def ensemble_ncbi_genes(flat_file, txid, flatadd=None):
    '''
    #tax_id	GeneID	Ensembl_gene_identifier	RNA_nucleotide_accession.version	Ensembl_rna_identifier	protein_accession.version	Ensembl_protein_identifier
    7227	30970	FBgn0040373	NM_001297803.1	FBtr0340207	NP_001284732.1	FBpp0309182
    '''
    print (flat_file)
    if re.search(r"\.gz$", flat_file):
        with gzip.open(flat_file, "rb") as fin:
            dat3 = [i for i in fin.read().split("\n") if re.search(r'^%s' % txid, i)]
    else:
        with open(flat_file, "r") as fin:
            dat3 = [i for i in fin.read().split("\n") if re.search(r'^%s' % txid, i)]
    gen = {}
    # print len(dat3), "dat3"
    for i in dat3:
        a = i.split("\t")
        gen[int(a[1])] = a[2]
    if flatadd:
        lis = map(int, gen.keys())
        lis.sort()
        with gzip.open(flatadd+".gz", "wb") as fin:
            fin.write("EntrezID\tEnsembleID\n")
            for i in lis:
                fin.write("%s\t%s\n" % (i, gen[str(i)]))
    return gen


def refseqEnsemble(ref, ens, swiss, flataddens=None, flataddswiss=None):
    ncbi_has = {}
    ens_has = {}
    swiss_has = {}
    ncbi_ens = {}
    ncbi_swiss = {}
    ncbi = os.listdir(ref)
    for i in ncbi:
        with open(os.path.join(ref,"%s" % i)) as fin:
            seq = "".join(fin.read().split("\n")[1:]).strip()
        if seq not in ncbi_has:
            ncbi_has[seq] = [i]
        else:
            ncbi_has[seq] += [i]
    # ncbi_hashes_in_memory
    if ens:
        for i in os.listdir(ens):
            with open(os.path.join(ens,"%s") % (i)) as fin:
                seq = "".join(fin.read().split("\n")[1:]).strip()
            if seq not in ens_has:
                ens_has[seq] = [i]
            else:
                ens_has[seq] += [i]
        for key in ncbi_has:
            if key in ens_has:
                for vals1 in ncbi_has[key]:
                    ncbi_ens[vals1] = ",".join(ens_has[key])
        if flataddens:
            lis = map(int, ncbi_ens.keys())
            lis.sort()
            with gzip.open(flataddens+".gz", "wb") as fin:
                fin.write("Refseq\tEnsembleID\n")
                for i in lis:
                    fin.write("%s\t%s\n" % (i, ncbi_ens[str(i)]))

    if swiss:
        for i in os.listdir(swiss):
            with open(os.path.join(swiss,"%s") % (i)) as fin:
                seq = "".join(fin.read().split("\n")[1:]).strip()
            if seq not in swiss_has:
                swiss_has[seq] = [i]
            else:
                swiss_has[seq] += [i]

        for key in ncbi_has:
            if key in swiss_has:
                for vals1 in ncbi_has[key]:
                    ncbi_swiss[vals1] = ",".join(swiss_has[key])
        if flataddswiss:
            lis = map(int, ncbi_swiss.keys())
            lis.sort()
            with gzip.open(flataddswiss+".gz", "wb") as fin:
                fin.write("Refseq\tSP\n")
                for i in lis:
                    fin.write("%s\t%s\n" % (i, ncbi_swiss[str(i)]))
    return ncbi_ens, ncbi_swiss


def gene_description(gene_table_dir):
    has = {}
    for i in os.listdir(gene_table_dir):
        with open(os.path.join(gene_table_dir,"%s") % i) as fin:
            dat = fin.read().split("\n")[0].split("[")[0]
        has[int(i.split(".")[0])] = dat
    return has


def go_terms_and_cellular_location(go_term_file, txid):
    if re.search(r"\.gz$", go_term_file):
        with gzip.open(go_term_file, "rb") as fin:
            dat3 = [i for i in fin.read().split("\n") if re.search(r'^%s' % txid, i)]
    else:
        with open(go_term_file, "r") as fin:
            dat3 = [i for i in fin.read().split("\n") if re.search(r'^%s' % txid, i)]
    gen_go = {}
    gen_comp = {}
    for i in dat3:
        ele = i.split("\t")
        gene = ele[1]
        gene = int(gene)
        if gene not in gen_go:
            gen_go[gene] = [ele[2]]
        else:
            gen_go[gene] += [ele[2]]
        if ele[7] == "Component":
            if gene not in gen_comp:
                gen_comp[gene] = [ele[5]]
            else:
                gen_comp[gene] += [ele[5]]
    gen_comp_revise = {}
    for i in gen_comp:
        valc = gen_comp[i]
        ctemp = 0
        for j in valc:
            if re.search(r'membran', j, re.IGNORECASE):
                ctemp += 1
        gen_comp_revise[i] = round(float(ctemp)/len(valc), 3)
    return gen_go, gen_comp_revise


def proteins_info_var(protein_structure_dir):
    bar2 = Bar('Parsing_structure for length:', max=len(os.listdir(protein_structure_dir)))
    has = {}
    for i in os.listdir(protein_structure_dir):
        stru = parser.get_structure("", os.path.join(protein_structure_dir,"%s") % i)
        res_list = Selection.unfold_entities(stru, 'R')
        var = "_".join(i.split("_")[:2])
        has[var] = [i, len(res_list)]
        bar2.next()
    bar2.finish()
    return has
    '''
    will return has, var as key, file as val [0 and length as val[1]
    '''

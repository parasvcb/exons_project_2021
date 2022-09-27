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

def updateHas(has,seq,fidentifier):
    if seq not in has:
        has[seq]=[]
    has[seq]+=[fidentifier]
    return has

def mappingstwoHasandFlatFile(hasp,hasc,flatf):
    haspc={}
    for key in hasp:
        if key in hasc:
            for vals1 in hasp[key]:
                haspc[vals1] = ",".join(hasc[key])
    if flatf:
        lis = map(int, haspc.keys())
        lis.sort()
        with gzip.open(flatf+".gz", "wb") as fin:
            fin.write("Refseq\tSP\n")
            for i in lis:
                fin.write("%s\t%s\n" % (i, haspc[str(i)]))
    return haspc

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
        ncbi_has=updateHas(ncbi_has,seq,i)
    # ncbi_hashes_in_memory
    if ens:
        for i in os.listdir(ens):
            with open(os.path.join(ens,"%s") % (i)) as fin:
                seq = "".join(fin.read().split("\n")[1:]).strip()
            ens_has=updateHas(ens_has,seq,i)
        ncbi_ens=mappingstwoHasandFlatFile(ncbi_has,ens_has,flataddens)
        
    if swiss:
        for i in os.listdir(swiss):
            with open(os.path.join(swiss,"%s") % (i)) as fin:
                seq = "".join(fin.read().split("\n")[1:]).strip()
            swiss_has=updateHas(swiss_has,seq,i)
        ncbi_swiss=mappingstwoHasandFlatFile(ncbi_has,swiss_has,flataddens)
    return ncbi_ens, ncbi_swiss

def mapping_alreadyDone_vs_Pending(unassignedDir,alreadydoneDir,copyDonetoUnassigned,pendingToRun,createTheMapFile):
    unassignedSeq=os.listdir(unassignedDir)
    alreadydoneSeq=[i for i in os.listdir(alreadydoneDir) if i[-7:]=='.dat.ss']
    print (alreadydoneSeq[:5])
    print (len(alreadydoneSeq))
    unassignedSeqHas={}
    alreadydoneSeqHas={}
    pendingHas={(0,300):['less300'],(301,600):['300to600'],(601,1500):['600to1500'],(1501,3000):['1500to3000']}
    for i in unassignedSeq:
        with open(os.path.join(unassignedDir,i)) as fin:
            seq = "".join(fin.read().split("\n")[1:]).strip()
        unassignedSeqHas=updateHas(unassignedSeqHas,seq,i)
    for i in alreadydoneSeq:
        #NP_000005.2.dat.ss
        #print (os.path.join(alreadydoneDir,i))
        with open(os.path.join(alreadydoneDir,i)) as fin:
            seq = "".join([j.split()[1] for j in fin.read().split("\n")[1:] if len(j)]).strip()
        alreadydoneSeqHas=updateHas(alreadydoneSeqHas,seq,i)
    #copyTheDone and storing pending
    print (len(alreadydoneSeqHas))
    mapfout= open(createTheMapFile,'w')
    megafile=[]
    alDoneTot=0
    alDonenr=0
    runTot=0
    runnr=0
    for i in unassignedSeqHas:
        if i in alreadydoneSeqHas:
            #map found
            alDonenr+=1
            alDoneTot+=len(unassignedSeqHas[i])
            mapFile=alreadydoneSeqHas[i][0]
            with open (os.path.join(alreadydoneDir,mapFile)) as fin:
                dat=fin.read()
            filesToCopy=unassignedSeqHas[i]    
            for iter in filesToCopy:
                with open (os.path.join(copyDonetoUnassigned,iter+'.dat.ss'),'w') as fout:
                    fout.write('%s'%dat)
        else:
            length=len(i)
            flag=False
            for j in pendingHas:
               if j[0]<=length<=j[1]:
                   pendingHas[j]+=[unassignedSeqHas[i][0]]
                   flag=True
                   break
            if flag:
                runnr+=1
                runTot+=len(unassignedSeqHas[i])
                string=",".join(unassignedSeqHas[i][1:])
                mapfout.write("%s\t%s\t%s"%(unassignedSeqHas[i][0],pendingHas[j],string))
            else:
                megafile+=[[unassignedSeqHas[i]]]
    mapfout.close()
    #makingDirs

    for i in pendingHas:
        dirname=pendingHas[i][0]
        if len(pendingHas[i])>1:
            os.makedirs(os.path.join(pendingToRun,dirname))

    for i in pendingHas:
        for files in pendingHas[i][1:]:
            with open (os.path.join(unassignedDir,files)) as fin:
                dat=fin.read()
            with open (os.path.join(pendingToRun,pendingHas[i][0],files),'w') as fin:
                fin.write('%s'%dat)

    print ('Total Were %s entries (%s non red)'%(len(unassignedSeq),len(unassignedSeqHas)))
    print ('>3000 Were %s entries (%s non red)'%(len([i for j in megafile for i in j]),len(megafile)))
    print ('PreAss Were %s entries (%s non red)'%(alDoneTot,alDonenr))
    print ('RunSet Were %s entries (%s non red)'%(runTot,runnr))

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

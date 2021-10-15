import os
import sys
import re
from math import ceil
from Bio.PDB import *
import cPickle as pickle
import subprocess
import multiprocessing
from Bio import SeqIO
import Bio.SeqUtils as three2one
import common_modules as cm
pdb_io = PDBIO()
parserpdb = PDBParser(QUIET=True)
parsermmcif = MMCIFParser(QUIET=True)
import tqdm
import gc
pat_colon = re.compile(r':+')
word_count = re.compile(r'\w+')
pat = re.compile(r'\d+\t\w+\.\w+\t\w+\t\d+\t\d+\t\w+\t\w{4}\n?')
# this pattern gets the first block of output from the align program
se_pat = re.compile(r'Sequence identity.*')
length_se2_pat = re.compile(r'Length of sequence 2:.*')
core_three_matcher = re.compile(r'[\-A-Z]+\n[\s\:]+\n[\-A-Z]+')


def sing_l(three):
    try:
        form = three[0].upper()+three[1:].lower()
        return three2one.IUPACData.protein_letters_3to1[form]
    except:
        return "X"


def pdb_writer(has_pseq_pdb, needapdbfileopenener, variant, struct_repo_dir, save_structure_dir, AF=False):
    #print ('Ever here ?', needapdbfileopenener)
    needs_reform = 0
    mod_id = 0

    try:
        new = []
        old = {}
        i = needapdbfileopenener
        if AF:
            #means alphafold strcutre
            a = i.split("_") # AF-P04217-F1-model_v1_A (file name with chain)
            pdb = "_".join(a[:-1])
            chain_d = a[-1]
        else:
            a = i.split("_")
            pdb = a[0]
            chain_d = a[1]
        #print (a,pdb,chain_d)
        if len(chain_d) > 1:
            needs_reform = 1
        structure = returnme_structure(pdb, struct_repo_dir)
        #print (structure)
    except Exception as e:
        #fout.write("in pdb writer,structre editing part, variant:%s %s %s\n" % (variant, gene, e))
        print ("in pdb writer,structre editing part, variant:%s %s %s\n" % (variant, gene, e))
    #print ('out', needs_reform)

    try:
        # pt=structure[mod_id][chain_d]
        c = -100000
        #print ('0.1')
        #print (structure, dir(structure))
        ch = structure[mod_id][chain_d]
        #print ('0.2')
        if needs_reform:
            try:
                #print (1)
                temp_mod = structure[mod_id]
                temp_mod.detach_child("A")
                #print ('okay')
            except:
                #print ('yes ERR')
                # means that chain doesnt exist
                pass
            chain_d = "A"
            ch.id = "A"
        for residue in ch:
            res_id = list(residue.id)
            res_id[1] = res_id[1]+10000
            residue.id = tuple(res_id)

        for residue in ch:
            if is_aa(residue, standard=True):
                if residue.id[2].isspace():
                    res_id = list(residue.id)
                    if (res_id[1]-10000) in has_pseq_pdb:
                        res_id[1] = has_pseq_pdb[res_id[1]-10000]
                        new += [res_id[1]]
                        residue.id = tuple(res_id)
                    else:
                        res_id[1] = c
                        residue.id = tuple(res_id)
                        c += 1
                else:
                    # means insertion code
                    res_id = list(residue.id)
                    if (str(res_id[1]-10000)+"_"+res_id[2]) in has_pseq_pdb:
                        res_id[1] = has_pseq_pdb[(
                            str(res_id[1]-10000) + "_"+res_id[2])]
                        new += [res_id[1]]
                        residue.id = tuple(res_id)
                    else:
                        res_id[1] = c
                        residue.id = tuple(res_id)
                        c += 1
            #print (residue,'residue')
    except Exception as e:
        print ("ipdbw:6 err::%s %s %s\n" % (variant, gene, e))
    #print ('progreess ??')
    try:
        class ResSelect(Select):

            def accept_residue(self, res):
                # print res.get_parent().parent.id
                # print dir(res)
                # rint "resid:",res.id,"resparentid:",res.parent.id,
                # "res_spuer_parent",res.get_parent().parent.id
                if res.id[1] in new and res.parent.id == chain_d and res.get_parent().parent.id == mod_id:
                    return True
                else:
                    return False

        pdb_io.set_structure(structure)
        #print('fnamePDB',os.path.join(save_structure_dir,"%s.pdb" %(variant+"_" + needapdbfileopenener)))
        pdb_io.save(os.path.join(save_structure_dir,"%s.pdb" %(variant+"_" + needapdbfileopenener)), ResSelect())
        # print save_structure_dir + \    "%s.pdb" % (variant+"_" + needapdbfileopenener)
        # print "saved seemslike"
    except Exception as e:
        print ("inpdw:8err:%s %s %s\n" % (variant, gene, e))
        raise


def gapped_regions_from_spans(spans_t, upper_t, lower_t):
    try:

        gapped_resgions_temp = []
        for i in range(0, len(spans_t)-1):
            gapped_resgions_temp += [(spans_t[i][1]+1,
                                      spans_t[i+1][0]-1)]
        # adding a front and end gapped sequence based on \
        # spans of pdb seq in below script
        # print "1P"
        before_span_lower_word_count = list(
            word_count.            finditer(upper_t[0:spans_t[0][0]]))
        after_span_lower_word_count = list(
            word_count.            finditer(upper_t[spans_t[-1][1]:len(lower_t)]))
        # print "2"
        if before_span_lower_word_count:
            for temp_addition in before_span_lower_word_count:
                gapped_resgions_temp += [temp_addition.span()]
        # print "@"
        if after_span_lower_word_count:
            add_factor = spans_t[-1][1]
            for temp_addition in after_span_lower_word_count:
                temp_val_in_tup = temp_addition.span()
                temp_tup = (temp_val_in_tup[0]+add_factor,
                            temp_val_in_tup[1]+add_factor)
                gapped_resgions_temp += [temp_tup]
        #print gapped_resgions_temp, "gapped_resgions_temp"
        return gapped_resgions_temp
    except Exception as e:
        print "in gapp()2:", e
        raise


def cond1_2(gapped_regions, upper, lower):
    try:
        proceed_t = 1
        for i in gapped_regions:
            if upper[i[0]:i[1]+1].count("-") > 0 or lower[i[0]:i[1]+1].count("-") != i[1]-i[0]+1:
                proceed_t = 0

        return proceed_t
    except Exception as e:
        raise
        return 0


def cond3_4(file_to_be_opened_temp, gapped_reg_temp1, lower, lisDir):
    '''
    This fucntion is called to perform the following functions
    1st:
    to load residue numbers from pickle in to a list, if that contains
    str in it modify it to consective numbers list else dont tdo anythng and retrun both
    2: modify the gapped regions to maintain the numbering of lower pdb sequnec
    equivalent to that of structure that is remove gaps, nad equate it equivalent
    to that pdb seq list
    3: a condition that checks if gaps in pdb are there , they arer exactly the missing
    residues numbering wise,stop the process if not
    '''
    try:
        proceed_flag = 1
        with open(os.path.join(lisDir,"%s.pck") % file_to_be_opened_temp, "rb") as fin:
            lis_pdb_res_map = pickle.load(fin)
        # print lis_pdb_res_map
        check_temp = not all(isinstance(xt, int) for xt in lis_pdb_res_map)
        # print check_temp
        if check_temp:
            # this will replace residue from insertion to a number
            # print "5"
            lis_var_mod = lis_pdb_res_map[:]
            for ind, val in enumerate(lis_var_mod):
                if type(val) == str:
                    lis_var_mod[ind] = lis_var_mod[ind-1]+1
                    for i in range(ind+1, len(lis_var_mod)):
                        if type(lis_var_mod[i]) == str:
                            break
                        else:
                            lis_var_mod[i] = lis_var_mod[i]+1
        else:
            lis_var_mod = lis_pdb_res_map[:]
    except Exception as e:
        print "cond3_4 p1"
        raise
    '''
    2 is done below
    '''
    try:
        gap_reg_mod = []
        for i in gapped_reg_temp1:
            count_dash_low = lower[0:i[0]].count("-")
            i2 = (i[0]-count_dash_low, i[1]-count_dash_low)
            gap_reg_mod += [i2]

        # print "6"
        # print len(lis_var_mod),"len)"
        # print gap_reg_mod,"gap_reg_mod"
        # count_dash_low=lower[0:spans[0][0]].count("-")
        # print count_dash_low

        # checking 1 as sub part
    except Exception as e:
        print "cond3_4 p2"
        raise
    '''
    3rd is a condition that
    '''
    try:
        for i2 in gap_reg_mod:
            if i2[0] != 0 and i2[1] != len(lower):
                # interested in between gaps of pdb, which can be disordered
                diff = i2[1]-i2[0]+1
                # print diff
                sub_l = list(
                    range(lis_var_mod[i2[0]-1], lis_var_mod[i2[0]-1]+diff+1))[1:]
                if set(sub_l) & set(lis_var_mod):  # means no disordered, they are present
                    proceed_flag = 0  # dont proceed ahead
        # print "7"

        return proceed_flag, gap_reg_mod, lis_pdb_res_map
    except exception as e:
        print "3_4 p3"
        raise


def pdb_setter(check, spanr, upper, lower, file_to_be_opened_wo_ext, lisDir):
    #print (check, file_to_be_opened_wo_ext, 'PDB_SETTER')
    if check:  # means only one match
        try:
            '''
                no condition ever evlautaed for this block , whoch makes it a part of pain
            '''
            # if 1:
            #print "isis"
            if not (re.match(r'[A-Z]+', lower[0:spanr[0]]) or upper[0:spanr[0]].count('-') > 0 or upper[spanr[1]:len(upper)].count('-') > 0):
                    # protseq_subm=upper[spanr[0]:spanr[1]]
                    # pdbseq_subm=lower[spanr[0]:spanr[1]]
                #print ('if in')
                count_dash_upp = upper[0:spanr[0]].count("-")
                count_dash_low = lower[0:spanr[0]].count("-")
                protseq_ind = (spanr[0]-count_dash_upp,
                               spanr[1]-count_dash_upp)
                # above has a tuple of length 2 having exact indices match of
                # prot seq of matched substring
                pdbseq_ind = (spanr[0]-count_dash_low, spanr[1]-count_dash_low)
                # but actual pdb mapping of resdue id's in pickle dumed file,
                # lets open it
                # print "protseq_ind pdb_seq_ind", protseq_ind,pdbseq_ind
                with open(os.path.join(lisDir,"%s.pck" % file_to_be_opened_wo_ext), "rb") as fin:
                    lis_pdb_res_map = pickle.load(fin)
                pdbseq_ind_1 = lis_pdb_res_map[pdbseq_ind[0]:pdbseq_ind[1]]
                protseq_ind_1 = list(range(protseq_ind[0], protseq_ind[1]))
                has_pseq_pdb = {}
                for t_i in range(0, len(protseq_ind_1)):
                    # pdb residue no as keys
                    has_pseq_pdb[pdbseq_ind_1[t_i]] = protseq_ind_1[t_i]
                # and prot seq NP resid as value
                # the pdb can be writen gracefully with keys as resideu number in sleection for modifictaion
                # but have to pay attetion to the inserted residue and its keys
                # print "**",has_pseq_pdb
                #print ('ifout')
                return has_pseq_pdb
            else:
                return False
            #print ('done >')
        except Exception as e:
            print (E, 'pdbsetter')
            
            fout.write("in pdb_setter first if, variant:%s %s %s\n" %
                       (variant, gene, e))
        #print ('setter 1')
    else:
        try:
            # 1.if they are in chunks then the difference of chunks
            # must be gaps in pdb,
            # 2.there should not be any gaps in the primary sequnece
            # for the difference string span,
            # it should contain alphabets
            # 3.if satisfied then they must correspond to the
            # disoredered regions, else kill them with grace
            # notice spanr is a iterator now
            proceed = 1
            spans = []
            for i in spanr:
                spans += [i.span()]
            # matching spans are stored in spans

            gapped_regions = gapped_regions_from_spans(spans, upper, lower)
            # this will give list of tuples having gapped regions
            # checking codition 2 and 1 in a function consectuvelily
            cond1_2_flag = cond1_2(gapped_regions, upper, lower)
            # print "2?"
            if cond1_2_flag:
                # going in if 1 and 2 is true
                # 3.first chevk if the residue list has str variable
                # in form of insertion coded residue it,
                # if yes, mutate it to another variable coounting for
                # the continuity
                # 4. and create modified gapped regions by substarcting nay
                # gaps before the very fisrt matches
                # cheking 1
                # print spans

                # now gaps regions are modified after substarcting gaps before alignment,
                # and but still not equal to pdb mapping
                # print "4"
                proceed_flagm, gap_reg_mod, lis_pdb_res_map = cond3_4(
                    file_to_be_opened_wo_ext, gapped_regions, lower, lisDir)

                if proceed_flagm:  # means all the gaps were disordered residues

                    count_dash_upp = upper[0:spans[0][0]].count("-")
                    # count_dash_low=lower[0:spans[0][0]].count("-")
                    # print count_dash_low
                    has_pseq_pdb = {}
                    temp_sub_fact = 0
                    upd_span_prot = []
                    for temp_i2 in range(0, len(spans)):
                        if temp_i2 == 0:
                            temp_sub_fact += lower[0:spans[temp_i2]
                                                   [0]].count("-")
                            upd_span_prot += [(spans[temp_i2][0]-temp_sub_fact,
                                               spans[temp_i2][1]-temp_sub_fact)]
                        else:
                            temp_sub_fact += lower[spans[temp_i2-1][1]:
                                                   spans[temp_i2][0]].count("-")
                            upd_span_prot += [(spans[temp_i2][0]-temp_sub_fact,
                                               spans[temp_i2][1]-temp_sub_fact)]
                    partial_pdbseq_list = []
                    # print "?"
                    for temp_i3 in upd_span_prot:
                        partial_pdbseq_list += list(
                            range(temp_i3[0], temp_i3[1]))
                    partial_NP_lis = []
                    # print "??"
                    for temp_i3 in spans:
                        partial_NP_lis += list(range(temp_i3[0]-count_dash_upp,
                                                     temp_i3[1]-count_dash_upp))
                    # print "???"
                    # print partial_NP_lis
                    # print partial_pdbseq_list
                    # print len(partial_NP_lis),len(partial_pdbseq_list)
                    for t_i in range(0, len(partial_NP_lis)):
                        # now there are gaps also in the sequence how that
                        # will be accounted ?
                        has_pseq_pdb[lis_pdb_res_map[partial_pdbseq_list[t_i]]
                                     ] = partial_NP_lis[t_i]
                    # print "**in pdb setter else",file_to_be_opened_wo_ext
                    return has_pseq_pdb
        except Exception as E:
            print ("in pdb_setter file last part_:%s %s\n" %
                   (file_to_be_opened_wo_ext, E))
            raise
    #print ('outofsetter')       

def pre_worker(tabdelim_alias, cx, has_gene_var):
    has_alias_to_id = {}
    with open(tabdelim_alias) as fin:
        '''
        tax_id  Org_name        GeneID  CurrentID       Status  Symbol  Aliases description     other_designations      map_location    chromosome      genomic_nucleotide_accession.version      start_position_on_the_genomic_accession end_position_on_the_genomic_accession   orientation     exon_count      OMIM
        9606    Homo sapiens    7157    0       live    TP53    BCC7, BMFS5, LFS1, P53, TRP53   tumor protein p53       cellular tumor antigen p53|antigen NY-CO-13|mutant tumor protein 53|p53 tumor suppressor|phosphoprotein p53|transformation-related protein 53|tumor protein 53|tumor supressor p53        17p13.1 17      NC_000017.11    7668402 7687550 minus   12191170
        9606    Homo sapiens    1956    0       live    EGFR    ERBB, ERBB1, HER1, NISBD2, PIG61, mENA  epidermal growth factor receptor        epidermal growth factor receptor|avian erythroblastic leukemia viral (v-erb-b) oncogene homolog|cell growth inhibiting protein 40|cell proliferation-inducing protein 61|epidermal growth factor receptor tyrosine kinase domain|erb-b2 receptor tyrosine kinase 1|proto-oncogene c-ErbB-1|receptor tyrosine-protein kinase erbB-1  7p11.2  7       NC_000007.14    55019021        55208080        plus    31131550
        '''
        dat = [i for i in fin.read().split('\n') if len(i) > 10]
    columns = dat[0].split('\t')
    # print columns
    sym_ind = columns.index('Symbol')
    alias_ind = columns.index('Aliases')
    gene_ind = columns.index('GeneID')
    for entry in dat[1:]:
        ele = entry.split('\t')
        valEsymInd = ele[sym_ind].strip()
        if valEsymInd not in has_alias_to_id:
            has_alias_to_id[valEsymInd] = []
        has_alias_to_id[valEsymInd] += [int(ele[gene_ind])]
        for alias in ele[alias_ind].split(','):
            indent = alias.strip()
            if indent not in has_alias_to_id:
                has_alias_to_id[indent] = []
            has_alias_to_id[indent] += [int(ele[gene_ind])]
    
    dat=cm.readFileSepslashN(cx,5)
    '''
    will read and return the file
    '''
    columns = dat[0].split("\t")
    np_ind = columns.index('Cross-reference (RefSeq)')
    pdb_ind = columns.index('Cross-reference (PDB)')
    entrez_ind = columns.index('Cross-reference (GeneID)')
    gene_name = columns.index('Gene names')
    has_id = {}
    has_npxp = {}
    cv = 0
    for entry in dat[1:]:
        ele = entry.split('\t')
        cv += 1
        if len(ele[pdb_ind]) > 0:
            # print ele[0], ele[pdb_ind]
            pdbx = [i.lower() for i in ele[pdb_ind].split(';') if len(i) > 0]
            # print pdbx
            # means pdb record exits
            if len(ele[np_ind]) > 0:
                # mean refseq npxp record exists:
                # print '1'
                for refseqId in re.findall(r'[NXY]P\_\w+\.\w+', ele[np_ind]):
                    if refseqId not in has_npxp:
                        has_npxp[refseqId] = set()
                    has_npxp[refseqId].update(pdbx)
                    # print '1!'
            if len(ele[entrez_ind]) > 0:
                for j in ele[entrez_ind].strip().split(';'):
                    if len(j) > 0:
                        # print '2'
                        #j = j.strip()
                        j = int(j)
                        if j not in has_id:
                            has_id[j] = set()
                        has_id[j].update(pdbx)
                        # print '2!'
            if len(ele[gene_name]) > 0:
                # print '3'
                for j in ele[gene_name].split():
                    j = j.strip()
                    if j in has_alias_to_id:
                        for k in has_alias_to_id[j]:
                            if k not in has_id:
                                has_id[k] = set()
                            has_id[k].update(pdbx)
                # print '3!'
    # print 'out'
    # print has_
    has_to_be_returned = {}
    # print 'his started'
    for gene in has_gene_var:
        if gene in has_id:
            if gene not in has_to_be_returned:
                has_to_be_returned[gene] = set()
            has_to_be_returned[gene].update(has_id[gene])
        for var in has_gene_var[gene]:
            if var in has_npxp:
                if gene not in has_to_be_returned:
                    has_to_be_returned[gene] = set()
                has_to_be_returned[gene].update(has_npxp[var])
    # for gene in has_to_be_returned:
        # has_to_be_returned[gene] = list(set(has_to_be_returned[gene]))
    has_to_be_returned = {gene: [i.lower() for i in has_to_be_returned[gene]] for gene in has_to_be_returned}
    # print has_to_be_returned.values()[:5]
    unique_pdbs = [j for i in has_to_be_returned.values() for j in i]
    # print 'has returned', unique_pdbs[:5]
    return has_to_be_returned, unique_pdbs

    '''
    This function will take in the refernce uniprot file with pdb CX
    then assign pdb directly to refseqs based on any Np mathcH if any
    thereafter assign pdb based on gene id and then based on gene alias
    so has hasingene alias to geneid match, is also needed
    return 1:
    gene and strcutres, no variant information in between
    '''


def chunks(l, n):
    div_fac = int(ceil(float(len(l))/n))
    # print div_fac
    """Yield successive n-sized chunks from l."""
    for i in range(0, n+1):
        yield l[i*div_fac: (i+1)*div_fac]


def returnme_structure(pdbid, repo_dir):
    #print (pdbid,repo_dir)
    id1 = pdbid+'.pdb'
    id2 = pdbid+'.cif'
    id3 = "pdb%s.ent" % (pdbid)
    filenames = [id1, id2, id3]
    for i in filenames:
        fname=os.path.join(repo_dir,'%s'%i)
        #print (fname,os.path.isfile(fname))
        if os.path.isfile(fname):
            if '.cif' in i:
                structure = parsermmcif.get_structure(pdbid, fname)
                return structure
            else:
                structure = parserpdb.get_structure(pdbid, fname)
                return structure
    return False


def list_pickler(pdbid, faaDir, lisDir, rep_dir):
    #print (pdbid, faaDir, lisDir, rep_dir)
    structure_ob = returnme_structure(pdbid, rep_dir)
    if structure_ob:
        #print(pdbid)
        for model in structure_ob:
            # print (model.id)
            for ch in model:
                # print (ch.id)
                chainPdbFile=os.path.join(faaDir,"%s_%s.faa" %(pdbid, ch.id))
                # print (chainPdbFile)
                # print ('why')
                if not os.path.isfile(chainPdbFile):
                    atoms = ch.get_atoms()
                    lisatom = set([ida.get_name() for ida in atoms if is_aa(
                        ida.get_parent()) and ida.get_parent().id[0].isspace()])
                    # must have been used to check if structure has good resolution
                    resnames = []
                    resids = []
                    if len(lisatom) > 4:
                        for res in ch:
                            if is_aa(res, standard=True):
                                if res.id[2].isspace():
                                    resnames += [sing_l(res.resname)]
                                    resids += [res.id[1]]
                                else:
                                    resnames += [sing_l(res.resname)]
                                    resids += [str(res.id[1])+res.id[2]]
                    # print (len(resnames))
                    if len(resnames) > 39:
                        # write them in a pdb fasta format and in pickle dumped list
                        try:
                            with open(chainPdbFile, "w") as fin:
                                fin.write(">%s_%s\n%s" % (pdbid, ch.id, cm.fastareturn(''.join(resnames))))
                            
                            with open(os.path.join(lisDir,"%s_%s.pck") % (pdbid, ch.id), "wb") as fin:
                                pickle.dump(
                                    resids, fin, protocol=pickle.HIGHEST_PROTOCOL)
                        except Exception as E:
                            print "Err: lis pickle, pdb id: %s, pdb chain: %s, error: %s" % (pdbid, ch.id, E)
            break
        
def pendingFilesPDB(inpDirFaa,sourceList,pdbRep, pdb=True):
    #print ('pdb',pdb)
    # here to read the file and gets in chain IDS as list
    #sourceList has name of PDBFiles, inpDirFaa contaiosn name of files +'_' chainId, pdbRep
    has_done_pdb_faa={}
    for doneFaa in os.listdir(inpDirFaa):
        if pdb:
            pdbid,chain=doneFaa.split('.')[0].split('_')
        else:
            chain=doneFaa.split('.')[0].split('_')[-1]
            pdbid="_".join(doneFaa.split('.')[0].split('_')[:-1])
        
        if pdbid not in has_done_pdb_faa:
            has_done_pdb_faa[pdbid]=set()
        has_done_pdb_faa[pdbid]|=set([chain])
    #print (sourceList[:5], has_done_pdbid_faa.keys()[:5])
    return [i for i in has_done_pdb_faa if i not in sourceList]





def fasta_dis_list(pdb_lis_organism, faaDir, lisDir, rep_dir, cores):
    def update(a):
        pbar.update()
    cpu = cores
    # precheck if the lists are already there
    # assuming name without extension represents pdb file name
    # pdb_lis_organism = [i for i in pdb_lis_organism if i == '4uxr']
    if pdb_lis_organism:
        sublis = list(cm.chunks_based_on_element_size(pdb_lis_organism, 1000))
        count=0
        for work in sublis:
            #print ('count pool pdbwriter: %s/%s'%(count,len(sublis)))
            pbar = tqdm.tqdm(total=len(work))
            pool = multiprocessing.Pool()
            for pdbFile in work:
                #print (pdbFile, faaDir, lisDir, rep_dir)
                pool.apply_async(list_pickler, args=(pdbFile, faaDir, lisDir, rep_dir,), callback=update)
            pool.close()
            pool.join()
            gc.collect()
            count+=1
        print "lis_writing_part_is_done"
    else:
        print "fasta of pdbs were written already"


def structure_miner(gene, faaDir, lisDir, struct_repo_dir, genes_with_structure,
                    has_gene_var, pdb_to_fasta_link, file_in_pdb,
                    file_in_cif, PID_add, save_structure_dir, align_program, AFFlag =False):
    #print (gene, faaDir, lisDir, struct_repo_dir, PID_add, save_structure_dir, align_program)
    '''
    lis, faaDir, lisDir, struct_repo_dir,
                                  genes_with_structure, has_gene_var,
                                  pdb_to_fasta_link, file_in_pdb,
                                  file_in_cif,pid_dir, save_structure_dir
    '''
    for variant in has_gene_var[gene]:
        if 1:
            lis_to_pro = []
            #print gene, variant, genes_with_structure[gene]
            for pdbs in genes_with_structure[gene]:
                if pdbs in pdb_to_fasta_link:
                    lis_to_pro.extend(pdb_to_fasta_link[pdbs])
                    # this will have the list of fasta file sto be compared
            lis_to_pro = set(lis_to_pro)
            #print (lis_to_pro)
            # we have now the set of files to align with refseq sequence
            # print lis_to_pro, "lis_to_pro"
            lis_sco = []
            # print lis_to_pro,"lis_to_pro"
            for two_var_file in lis_to_pro:
                if 0:
                    print two_var_file, variant
                    print PID_add+"%s" % variant
                    print faaDir+"%s" % two_var_file
                try:
                    result_seq = subprocess.check_output(
                        [align_program, os.path.join(PID_add,"%s" % variant),
                            os.path.join(faaDir,"%s") % two_var_file, '0'])
                except:
                    result_seq = 0
                
                if bool(result_seq):
                    # print 'in'
                    #print ('1,',result_seq)
                    # sequence identity pattern
                    mat = se_pat.findall(result_seq)
                    identity = float(mat[0].split(
                        ":")[-1].split("(")[0].strip())  # score identity
                    len_mat_se2 = length_se2_pat.findall(result_seq)
                    len_seq_2 = int(len_mat_se2[0].split(":")[
                                    1].split("-")[0])
                    core_aligned_reg = result_seq.split("\n\n")[
                        1].split("\n")
                    length_factor = len(core_aligned_reg)/4
                    tamp_var = 0
                    upp1t = "".join(core_aligned_reg[0:length_factor])
                    midd1t = "".join(
                        core_aligned_reg[length_factor:length_factor*2])
                    low1t = "".join(
                        core_aligned_reg[length_factor*2:length_factor*3])
                    upp = re.sub(r"\s+", "", upp1t)
                    low = re.sub(r"\s+", "", low1t)

                # print "^^^^^^^^^^",variant, three_var_file, identity
                    # print identity, len_seq_2
                    if identity >= 0.95:
                        # print "***hurray*****",variant, three_var_file, identity,len_seq_2
                        lis_sco += [[len_seq_2, identity,
                                        two_var_file, [upp1t, midd1t, low1t]]]
                    # this has lengthofprotein,%seq_identity,name_offile,
                # chain.faa,and matched strings upper,middle,lower
                # in three different array elements
            #print ('liscso',lis_sco)
            if bool(lis_sco):
                # print "yes lis_sco"
                # fout.write("in pdb writer variant:%s %s %s\n"%(variant,gene,e))
                lis_sco.sort(reverse=True)

                # print lis_sco,"********lis_sco"
                # print lis_sco[0]

                for pdb_ch_hits in lis_sco:
                    temp_file_name_3 = pdb_ch_hits[2].split(".")[0]
                    # temp_matches = re.finditer(r'\w+', low1t)
                    temp_matches = pat_colon.finditer(
                        pdb_ch_hits[3][1])  # one is middle
                    # where [pdbexists
                    # print temp_middle,"temp_middle"
                    # print temp_upper,"temp_upper"
                    # print temp_lower,"temp_lower"
                    # print temp_file_name_3, variant, "temp_file_name_3"

                    temp_mat_list = list(temp_matches)
                    # print len(temp_mat_list)
                    if len(temp_mat_list) == 1:
                        spanrange = temp_mat_list[0].span()
                        # print "lessthan1_spanrange", spanrange
                        returning_has = pdb_setter(
                            1, spanrange, pdb_ch_hits[3][0], pdb_ch_hits[3][2], temp_file_name_3, lisDir)
                        # last three are upper and lower and filename
                    else:
                        # print "!lessthan1_spanrange", temp_mat_list
                        # print "yes interesting******"
                        returning_has = pdb_setter(0, temp_mat_list,
                                                    pdb_ch_hits[3][0], pdb_ch_hits[3][2], temp_file_name_3, lisDir)
                        # print temp_upper
                        # print temp_middle
                        # print temp_lower
                    # print "returning has", returning_has
                    if returning_has:
                        #print 'yes',returning_has
                        try:
                            pdb_writer(returning_has, temp_file_name_3,
                                        variant, struct_repo_dir, save_structure_dir, AF=AFFlag )
                            
                        except Exception as E:
                            raise
                            print gene, E

                        # print "pdb_ch_hits"

                        break


def horse(tabdelim_alias, cx, has_gene_var, faaDir, lisDir,
          struct_repo_dir, save_structure_dir, pid_dir, align_location, cores):
    def update(a):
        pbar.update()
    
    '''
    this will have athe list of genes,. whose structure
    can vbe solved or resolved
    this list will be dive in the number of cores
    that system can handle -5 in my case
    rest of the fucntions will be callable from this very isinstance
    other fyunmctions will be, aligner, condition verifer, and then setter

    '''
    genes_with_structure, pdb_lis_organism = pre_worker(
        tabdelim_alias, cx, has_gene_var)
    # print pdb_lis_organism[:5]
    # print genes_with_structure
    '''
    create a list and fast asequence file fro only one model that is 0th
    iterate the list so that all needed fasta and
    residue numbers will be stored locally
    '''
    #-> make above convenient for the variat or the gene may still work
    pdb_lis_organism = pendingFilesPDB(faaDir,pdb_lis_organism, struct_repo_dir)    
    print (len(pdb_lis_organism))
    #sys.exit()
    
    fasta_dis_list(pdb_lis_organism, faaDir, lisDir, struct_repo_dir, cores)
    
    pdb_to_fasta_link = {}
    for fil in os.listdir(faaDir):
        pdb_fil_name = fil[0:4]
        if pdb_fil_name not in pdb_to_fasta_link:
            pdb_to_fasta_link[pdb_fil_name] = []
        pdb_to_fasta_link[pdb_fil_name] += [fil]
        #->making link towards the pdb_to_its chains fasta seqeunces
    
    file_in_pdb = {}
    file_in_cif = {}

    for fil_temp in os.listdir(struct_repo_dir):
        if ".pdb" in fil_temp:
            file_in_pdb[fil_temp.split(".")[0].lower()] = 0
        
    #cpu = multiprocessing.cpu_count()-34
    cpu = cores
    # verify once, if the transcript strcutre is available then dont count that transcipt to be proceed ahead,
    # else
    sublis = list(cm.chunks_based_on_element_size(list(genes_with_structure.keys()), 1000))
    count=0
    print ('going to model now')
    for lis in sublis:
        print ('count pool : %s/%s'%(count,len(sublis)))
        pbar = tqdm.tqdm(total=len(lis))
        pool = multiprocessing.Pool()
        for gene in lis:
            pool.apply_async(structure_miner, args=(gene, faaDir, lisDir, struct_repo_dir, genes_with_structure, has_gene_var, pdb_to_fasta_link, file_in_pdb, file_in_cif, pid_dir, save_structure_dir, align_location,), callback=update)
        pool.close()
        pool.join()
        gc.collect()
        count+=1
        #sys.exit()
    print "donw ith care"


def AF_generator(has_gene_var, has_refseq_swiss, faaDirAF, lisDirAF,
          struct_repo_dirAF, save_structure_dirAF, pid_dir, align_location, cores):
    #sys.exit()
    def update(a):
        pbar.update()
    swissPdb_models={}
    for AFpdbFile in os.listdir(struct_repo_dirAF):
        #AF-X6R8D5-F1-model_v1.pdb.gz
        unip=AFpdbFile.split('-')[1]
        #print (AFpdbFile,unip)

        if unip not in swissPdb_models:
            swissPdb_models[unip]=[]
        swissPdb_models[unip]+=[AFpdbFile.split('.')[0]]
        #break
    genes_with_structure={gene:[] for gene in has_gene_var}
    for gene in has_gene_var:
        for var in has_gene_var[gene]:
            if var in has_refseq_swiss:
                swissGenome=has_refseq_swiss[var]
                #print (swissGenome)
                #sys.exit()
                for indEntry in swissGenome.split(','):
                    if indEntry in swissPdb_models:
                        genes_with_structure[gene]+=swissPdb_models[indEntry]
    genes_with_structure={gene:genes_with_structure[gene] for gene in genes_with_structure if len(genes_with_structure[gene])>0}
    print (genes_with_structure[genes_with_structure.keys()[0]])
    print (len(genes_with_structure), 'genes have associtaed tsructres')
    print (genes_with_structure.values()[:5])
    pdb_lis_organism=list(set([j for i in genes_with_structure.values() for j in i]))
    pdb_lis_organism = pendingFilesPDB(faaDirAF,pdb_lis_organism, struct_repo_dirAF,pdb=False)
    print (len(pdb_lis_organism),len(genes_with_structure))
    #sys.exit()
    #print (ab)
    fasta_dis_list(pdb_lis_organism, faaDirAF, lisDirAF, struct_repo_dirAF, cores)
    #capitalized
    
    pdb_to_fasta_link = {}
    for fil in os.listdir(faaDirAF):
        pdb_fil_name = fil[:-6] # removing chain name also
        if pdb_fil_name not in pdb_to_fasta_link:
            pdb_to_fasta_link[pdb_fil_name] = []
        pdb_to_fasta_link[pdb_fil_name] += [fil]
        #->making link towards the pdb_to_its chains fasta seqeunces
    
    file_in_pdb = {}
    file_in_cif = {}

    for fil_temp in os.listdir(struct_repo_dirAF):
        if ".pdb" in fil_temp:
            file_in_pdb[fil_temp.split(".")[0]] = 0
        elif ".ent" in fil_temp:
            file_in_pdb[fil_temp.split(".")[0][-4:]] = 0
        else:
            file_in_cif[fil_temp.split(".")[0]] = 0
    #cpu = multiprocessing.cpu_count()-34
    cpu = cores
    # verify once, if the transcript strcutre is available then dont count that sublis = list(cm.chunks_based_on_element_size(list(genes_with_structure.keys()), 1000))
    count=0
    sublis = list(cm.chunks_based_on_element_size(list(genes_with_structure.keys()), 1000))
    print ('going to model now')
    for lis in sublis:
        print ('count pool : %s/%s'%(count,len(sublis)))
        pbar = tqdm.tqdm(total=len(lis))
        pool = multiprocessing.Pool()
        for gene in lis:
            pool.apply_async(structure_miner, args=(gene, faaDirAF, lisDirAF, struct_repo_dirAF, genes_with_structure, has_gene_var, pdb_to_fasta_link, file_in_pdb, file_in_cif, pid_dir, save_structure_dirAF, align_location, True,), callback=update)
        pool.close()
        pool.join()
        gc.collect()
        count+=1
    print "donw ith care"


'''
tabdelim_alias = '/home/paras/project/protein_splicing/raw_files/alias/9606.txt'
cx = '/home/paras/project/protein_splicing/raw_files/uniprot-organism__Homo+sapiens+(Human)+[9606]_.tab.gz'
with open("/home/paras/mysql/project1/src/cpickles/has_gene_var_9606.pick") as fin:
    has_gene_var = pickle.load(fin)
faaDir = '/home/paras/project/protein_splicing/9606/derived_data/structure_data/fasta_file_pdb/'
lisDir = '/home/paras/project/protein_splicing/9606/derived_data/structure_data/res_no_pdb/'
struct_repo_dir = '/media/CSB/pdbs_splicing/new_pdb/'
save_structure_dir = '/home/paras/project/protein_splicing/9606/derived_data/structure_data/PDB_structure_derived_data/'
pid_dir = '/home/paras/mysql/Refseq_protein/9606/'
align_location = '/home/paras/mysql/project1/bin/temp/./align'
horse(tabdelim_alias, cx, has_gene_var, faaDir, lisDir,
      struct_repo_dir, save_structure_dir, pid_dir, align_location)
'''

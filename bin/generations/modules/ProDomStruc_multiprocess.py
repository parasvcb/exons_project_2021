import cPickle as pickle
import gzip
import multiprocessing
import os
import re
import retriever
from Bio.PDB import *
import subprocess
import numpy as np
import Bio.SeqUtils as three2one
from Bio import SeqIO
from progress.bar import Bar
import more_itertools as mit
from math import ceil
import warnings
import sys
import tqdm
import gc
import common_modules as cm
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
'''
pfam is a sequence repository
return total residues/ residues covered in doamins, with e val criteria and without that criterai, with overlap conditions

cath domains file
sift map
pdb map (supposed to be generated)

flow will be like this, from cath domains all file, get the domain id and domain seqres chopping file
get to know associated uniprot of the file from shift map , if both are from humans,
check sequence identity and match and flagt it as correct and get updated uniprot coordinates

once done, see , if this is continuous or discontiuous, and the information to object
'''


def domall2seqres(domall, filsave):
    with open(domall) as fin:
        dat = [i for i in fin.read().split("//") if len(i) > 10]
    has1 = {}
    has_top = {}
    # print dat[0]
    # print dat[-1]
    for i in dat:
        span1 = []
        span2 = []
        try:
            domid = re.search(r'^DOMAIN\s+.*$', i,
                              re.MULTILINE).group().split()[1]
        except:
            print i, "i::"
        allmat = re.finditer(
            r'^SRANGE\s+START=\-?\d+\w?\s+STOP=\-?\d+\w?$', i, re.MULTILINE)
        alltop = re.search(r'^TOPOL\s+.*$', i, re.MULTILINE)
        if alltop:
            # print alltop.group()
            acttop = " ".join(alltop.group().split()[1:])
            has_top[domid] = acttop
        for j in allmat:
            testres = j.group()
            testres = re.sub(r'SRANGE\s+', '', testres)
            testres = re.sub('START=', '', testres)
            testres = re.sub('STOP=', '', testres)
            testres = re.sub(r'\s+', ',', testres)
            span1 += [testres]
            # print testres
            # span2 += [testres.split(",")]
        has1[domid] = ":".join(span1)
        
    with open(filsave, "w") as fin:
        fin.write("#domoainID\tPDBCood\n")
        for i in has1:
            fin.write("%s\t%s\n" % (i, has1[i]))
    return has1, has_top


def sing_l(three):
    try:
        form = three.capitalize()
        return three2one.IUPACData.protein_letters_3to1[form]
    except Exception as E:
        print E
        return "X"


def pdb_setter(result_seq, lis_pdb_res_map):
    #print (result_seq)
    #print (lis_pdb_res_map)
    core_aligned_reg = result_seq.split("\n\n")[1].split("\n")
    length_factor = len(core_aligned_reg)/4
    upp1t = "".join(core_aligned_reg[0:length_factor])
    # midd1t = "".join(core_aligned_reg[length_factor:length_factor*2])
    low1t = "".join(core_aligned_reg[length_factor*2:length_factor*3])
    '''
    temp_matches = re.finditer(r'\w+', low1t)  # matching where pdb exists
    spans = [i.span() for i in temp_matches]
    # print spans, "spans"
    has_pseq_pdb = {}
    for spanr in spans:
        print spanr
        count_dash_upp_before = upp1t[0:spanr[0]].count("-")
        count_dash_upp_after = upp1t[0:spanr[1]].count("-")
        count_dash_low = low1t[0:spanr[0]].count("-")
        print count_dash_upp_before, 'count_dash_uppbefore', count_dash_upp_after, "count_dash_upp_after", count_dash_low, 'count_dash_low'
        protseq_ind = (spanr[0]-count_dash_upp_before,
                       spanr[1]-count_dash_upp_after)
        print protseq_ind, 'protseq_ind'
        pdbseq_ind = (spanr[0]-count_dash_low, spanr[1]-count_dash_low)
        print pdbseq_ind, 'pdbseq_ind'
        pdbseq_ind_1 = lis_pdb_res_map[pdbseq_ind[0]:pdbseq_ind[1]]
        protseq_ind_1 = list(range(protseq_ind[0], protseq_ind[1]))
        print pdbseq_ind_1
        for t_i in range(0, len(protseq_ind_1)):
            # pdb residue no as keys
            print t_i, "t_i", 'pdbseq_ind_1[t_i]', pdbseq_ind_1[t_i], 'protseq_ind_1[t_i]', protseq_ind_1[t_i]

            has_pseq_pdb[pdbseq_ind_1[t_i]] = protseq_ind_1[t_i]
            # break
    return has_pseq_pdb, (upp1t, low1t)
    '''
    has_raw = {}
    for i in range(0, len(low1t)):
        if low1t[i] != '-' and upp1t[i] != '-':
            newli = i-low1t[0:i].count('-')
            newui = i-upp1t[0:i].count('-')
            has_raw[lis_pdb_res_map[newli]] = newui
    return has_raw, (upp1t, low1t)



def aligner(iteration, pdbdir, spdir, temp_add, align_add):
    out=[]
    '''
    p = multiprocessing.Process(target=aligner, args=(
        sublis[i], i, new_matrix_queue, pdbdir, spdir,))
    '''
    c = 0

    try:
        uniprot = iteration[0]
        pdb = iteration[1][: 4]
        chain = iteration[1][4]
        if os.path.isfile(pdbdir+"%s.cif" % pdb):
            pdbname = "%s.cif" % pdb
            parser = MMCIFParser()
        else:
            pdbname = "%s.pdb" % pdb
            parser = PDBParser()
        seq = os.path.join(spdir,uniprot)
        pdbloc = os.path.join(pdbdir,pdbname)
        coods = iteration[3]

        se_pat = re.compile(r'Sequence identity.*')
        structure = parser.get_structure("", pdbloc)
        reslist = structure[0][chain]
        # print len(reslist), "reslist"
        resiterator = []
        resnames = []
        resids = []
        for res in reslist:
            if is_aa(res, standard=True):
                bitresid = [sing_l(res.resname), res.id[1]] if res.id[2].isspace() else [
                    sing_l(res.resname), str(res.id[1])+res.id[2]]
                bitresnam = sing_l(res.resname)
                bitresid1 = str(res.id[1]) if res.id[2].isspace() else str(
                    res.id[1])+res.id[2]
                resnames += [bitresnam]
                resids += [bitresid1]
                resiterator += [bitresid]
                # print bitresid1, bitresnam, "bit"
        del(structure)
        if bool(resiterator):
            flagT=0
            
            fname2=os.path.join(temp_add,"test_%s.strucfaa" %(pdb+chain))
            try:    
                if not os.path.isfile(fname2):
                    with open(fname2, "w") as fin:
                        tempstring=''.join(resnames)
                        stringTowrite = cm.fastareturn(tempstring)
                        fin.write(">pdb_%s\n%s" % (chain, stringTowrite))

            except Exception as E:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                print("oops, there's error in writer: %s, pdb:%s, chain:%s, %s\n" % (
                    E, pdbloc, chain, iteration), exc_tb.tb_lineno)
            result_text = subprocess.check_output(
                [align_add, seq, fname2, "0"])
            out+=['In iteration', iteration, result_text]
            mat = se_pat.findall(result_text)  # sequence identity pattern
            identity = float(mat[0].split(
                ":")[-1].split("(")[0].strip())  # score identity

            if identity > 0.949:
                residue_uniprot_map, alignment = pdb_setter(
                    result_text, resids)
                #print ('out')
                map_2 = [i for i in resids if i in residue_uniprot_map]
                boundary = (map_2[0], map_2[-1])
                updated_coods = []
                cov = []
                dom_len_pdb = []
                coodsnew = coods.split(":")
                for coordsingle in coodsnew:
                    left, right = coordsingle.split(",")
                    if left in residue_uniprot_map and right in residue_uniprot_map:
                        left_cood = str(residue_uniprot_map[left])
                        right_cood = str(residue_uniprot_map[right])
                        updated_coods += [left_cood+"," + right_cood]
                    else:
                        log_statement = ''
                        if left not in residue_uniprot_map:
                            leftind = 0
                            left1 = left
                            while True:
                                left1 = resids[resids.index(left1)+1]
                                leftind += 1
                                if left1 in residue_uniprot_map:
                                    left_cood = str(
                                        residue_uniprot_map[left1])
                                    log_statement += 'LeftInd:+%s' % leftind
                                    break
                                if left1 == boundary[1] or left1 == right:
                                    left_cood = 'undef'
                                    break
                        else:
                            left_cood = str(residue_uniprot_map[left])
                            log_statement += 'LeftInd:in_has'
                        if right not in residue_uniprot_map:
                            rightind = 0
                            right1 = right
                            while True:
                                #right1 = str(int(re.match(r'[-]?\d+', right1).group(0))-1)
                                right1 = resids[resids.index(right1)-1]
                                rightind += 1
                                if right1 in residue_uniprot_map:
                                    right_cood = str(
                                        residue_uniprot_map[right1])
                                    log_statement += '\nRightInd:-%s' % rightind
                                    break
                                if right1 == boundary[0] or right1 == left:
                                    right_cood = 'undef'
                                    log_statement += '\nRightInd:_undef'
                                    break
                        else:
                            right_cood = str(residue_uniprot_map[right])
                            log_statement += '\nRightInd:in_has'
                        updated_coods += [left_cood+"," + right_cood]
                        log_statement += "\nAlignment:\n%s" % (
                            "\n".join(alignment))
                        log_statement += '\nPDB_sequence:\n'
                        log_statement += ",".join([iind+":"+jind
                                                    for iind, jind in zip(resids, resnames)])
                        log_statement += "\nold_coods:%s, new_coods:%s\n" % (
                            coordsingle, left_cood+"," + right_cood)
                        if log_statement:
                            out+=['returning,logstatemnet']
                            #print (out)
                            return ["Log", log_statement, uniprot, iteration[2], iteration[3], 'undef', 'undef']

                    '''
                    coverage will be domain length i npdb, and in that domain length aand exclueisobvely in that
                    residue rage, how many are covered in uniport seuqnce
                    '''
                    flag_proceed = 0
                    try:
                        residue_dom_range = resids[resids.index(
                            left):resids.index(right)+1]
                        covered_residues = sum(
                            [1 for i in residue_dom_range if i in residue_uniprot_map])
                        cov += [str(round(float(covered_residues) /
                                            (len(residue_dom_range)), 2))]
                        dom_len_pdb += [(str(len(residue_dom_range)))]
                        flag_proceed = 1
                    except:
                        cov += ['undef']

                if 'undef' not in cov:
                    out+=['returning,283']
                    #print (out)
                    return [uniprot, iteration[2], iteration[3], round(
                        identity, 3), ":".join(updated_coods), ":".join(cov), ":".join(dom_len_pdb)]
                else:
                    out+=['returning,283']
                    #print (out)
                    return ["Error", "coverage_error",
                                uniprot, iteration[2], iteration[3], 'undef', 'undef']
            else:
                out+=['None']
                #print (out)
                return ['None']
    except Exception as EE:
        print ('returning, exception', iteration)
        exc_type, exc_obj, exc_tb = sys.exc_info()
        print (EE, iteration, exc_tb.tb_lineno, "super *****************")
        return ["Error", EE, uniprot, iteration[2],
                    iteration[3], 'undef', 'undef']
    print ('Not returning')


'''
aligner([['P56817', '3ohfA', '', '14,146']], 0, [], '/home/paras/mysql/src/cath/cif_files/',
        '/home/paras/project/protein_splicing/9606/source_data/swissprot_data/')

'P56817', '3ohfA', '3ohfA01', '14,146', True
# P61626	1lmtA00	1,130	undef
# P62837	4ddgA03	1088,1229
# P34998	4k5yA01	115,1002:224,368	0.952	114,396:252,396	0.61:1.0	407:141
#['P56817', '3msjB', '3msjB01', '14,146', True]
'''


def cath(domall2seqres, sift, pdbdir, spdir, humanspshash, download_pdb_flag, temporary_iteration_saver, permanent_iteration_saver, error_log, cores, common_data, align_add):
    print ('In prodomstruc')
    print (len(domall2seqres), domall2seqres.keys()[:5], domall2seqres[domall2seqres.keys()[0]])
    '''
    will take 2 files, seqres chopping and uniprot mapp with a directory of both uniport and assocuated pdb_ssp_human_prot

    ProDomStruc.cath(domall_hash, sift, pdbdir, spdir, humansps,
                     download_pdb_flag, cases_before_align, cases_after_align, error_logs)

    '''
    def update(a):
        pbar.update()
        if a[0]!='None':
            result_list.append(a)

    scratchDir=os.path.join(common_data,'scratch')
    if not os.path.isdir(scratchDir):
        os.makedirs(scratchDir)
    

    # humansps will be hash with list of sps
    dat = humanspshash
    sps = dat.values() # string and comma sperated
    has_temp = {}
    for i in sps:
        temp_a = i.split(",")
        for j in temp_a:
            has_temp[j.split(".")[0]] = 0
    sps = has_temp # now a dict with each key as individual entity
    # print sps.keys()[: 10]
    # sps is a hash now with sps of humans stored in it
    '''
    PDB	CHAIN	SP_PRIMARY	RES_BEG	RES_END	PDB_BEG	PDB_END	SP_BEG	SP_END
    101m	A	P02185	1	154	0	153	1	154
    102l	A	P00720	1	40	1	40	1	40
    '''
    # from this sift file extract the uniprots as keys and pdb+chain as values
    # intersect it with sps and proceed ahead for those to match and then transfer the coods
    if re.search(r"\.gz$", sift):
        with gzip.open(sift, "rb") as fin:
            dat3 = [i for i in fin.read().split(
                "\n") if re.search(r'^\w{4}\s+', i)]
    else:
        with open(sift, "r") as fin:
            dat3 = [i for i in fin.read().split(
                "\n") if re.search(r'^\w{4}\s+', i)]
    print 'len', len(dat3)
    cases_interest = {}
    for i in dat3:
        ele = i.split("\t")
        if ele[2] in sps:
            if ele[2] not in cases_interest:
                cases_interest[ele[2]] = []
            cases_interest[ele[2]] += [ele[0]+ele[1]]
    print len(cases_interest), "casint"
    #cases intesrt contains uniprot as key and corrwsponding pdb chain + id
    '''
    cases_interest should now proceed ahead to align,
    BUT have a look at cases of domall file also
    Parse domall_file
    '''
    pdb_s = []
    for i in os.listdir(pdbdir):
        tpdbs = i.split(".")
        if len(tpdbs[0]) == 4:
            pdb_s += [i[:4].lower()]
        else:
            pdb_s += [i[3:7].lower()]

    total_pdb_set = [i[:4].lower() for j in cases_interest.values() for i in j]
    total_domain_set = [i[:4].lower() for i in domall2seqres.keys()]
    total_domain_set_with_pdbs = set(total_pdb_set) & set(total_domain_set)
    print (len(total_pdb_set),total_pdb_set[:5])
    print (len(total_domain_set),total_domain_set[:5])
    
    print len(total_domain_set_with_pdbs), "total_domain_set_with_pdbs"
    pending_pdbs = list(set(total_pdb_set)-set(pdb_s))
    print len(pending_pdbs), "pending"
    with open(os.path.join(common_data,"__pdb_pending.lis"), "a+") as fin:
        for i in pending_pdbs:
            fin.write("%s\n" % i)

    print "retrieving them and exiting"
    # retriever.mmcif_pdb_retriever(pending_pdbs, pdbdir)
    print "done_rretrieving_now_exiiting"
    # sys.exit()
    matrix = []
    domall2seqres_refines = {}
    # domall2seqres refines hash will have key as pbdchain and value as list of all domains inside it
    for i in domall2seqres:
        keyi = i[:5]
        if keyi not in domall2seqres_refines:
            domall2seqres_refines[keyi] = []
        domall2seqres_refines[keyi] += [(i, domall2seqres[i])]
    # pdbchain as key values as domain id and its coordinates
    for uniprot in cases_interest:
        #uniprot is uniport and value is the list of pdb and chain assocuiated wit it 
        for listofPdbChain in cases_interest[uniprot]:
            if listofPdbChain in domall2seqres_refines:
                for tupleOfcomplDomIdandSpan in domall2seqres_refines[listofPdbChain]:
                    continuous = True if ':' not in tupleOfcomplDomIdandSpan[1] else False
                    matrix += [[uniprot, listofPdbChain, tupleOfcomplDomIdandSpan[0],
                                tupleOfcomplDomIdandSpan[1], continuous]]
    matrix.sort()

    with open(temporary_iteration_saver, "w") as fin:
        fin.write("#Uniprot\tPDBChain\tDomID\tCoordinates\tContinuous\n")
        for i in matrix:
            fin.write("%s\t%s\t%s\t%s\t%s\n" % (i[0], i[1], i[2], i[3], i[4]))

    # now you have to align
    # now the task is to send the file to prgranm which will create a fasta seq for the file and save a luist in memory fro resids positions
    # for case_i in cases_interest:
    # cpu = multiprocessing.cpu_count()-1
    cpu = cores
    processes = []
    result_list=[]
    finerr = open(error_log, "w")
    findone = open(permanent_iteration_saver, "w")
    findone.write(
        "#Uniprot\tDomID\tCoordinates\tIdentity\tUnip_Coods\tPseudo_cov\tDom_len\n")
    count=0
    #dividing list into chunks_based_on_cores of 1000 and then processing in pool
    sublis = list(cm.chunks_based_on_element_size(matrix, 1000))
    for mat in sublis:
        print ('count pool : %s/%s'%(count,len(sublis)))
        pbar = tqdm.tqdm(total=len(mat))
        pool = multiprocessing.Pool()
        for i in mat:
            pool.apply_async(aligner, args=(i, pdbdir, spdir, scratchDir, align_add, ), callback=update)
        pool.close()
        pool.join()
        gc.collect()
        count+=1
    print (len(result_list))
    #aligner(['O14770', '5eg0A', '5eg0A00', '279,334', True],pdbdir, spdir, scratchDir, align_add)
    result_list.sort()
    #print (result_list)
    bar_prog = Bar('writing :: ', max=len(result_list))
    
    for i in result_list:
        bar_prog.next()
        print ('writing,',i)    
        if i[0] == "Error":
            finerr.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                         (i[0], i[1], i[2], i[3], i[4], i[5], i[6]))
        elif i[0] == "Log":
            finerr.write("\n******%s\t%s\t%s\t%s\t%s\n%s******\n" %
                         (i[0], i[2], i[3], i[4], i[5], i[1]))
        else:
            # [uniprot, ele[1], ele[2], ele[3], ele[4],round(identity, 3), ":".join(updated_coods)])
            findone.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                          (i[0], i[1], i[2], i[3], i[4], i[5], i[6]))
    for i in os.listdir(scratchDir):
        os.remove(os.path.join(scratchDir,i))
    print ('Done')
    bar_prog.finish()

    findone.close()
    finerr.close()
    

def mini_csv(has, csv):
    # print has
    print csv
    # print has.values()
    total = len("".join(has.values()).strip())
    # print total
    keys = has.keys()
    keys.sort()
    with open(csv, "w") as fin:
        fin.write("Fraction,Total,Freq,G,A\n")
        for i in keys:
            cCount = round(float(has[i].count('G'))/total, 3)
            aCount = round(float(has[i].count('A'))/total, 3)
            tot_iter = len(has[i])
            tot_iter_freq = round(float(tot_iter)/total, 3)
            ikEy = str(i)
            ikEy = re.sub(r'\(', '', ikEy)
            ikEy = re.sub(r'\)', '', ikEy)
            ikEy = re.sub(r'\s+', '', ikEy)
            imod = "-".join(ikEy.split(","))
            fin.write("%s,%s,%s,%s,%s\n" %
                      (imod, tot_iter, tot_iter_freq, cCount, aCount))


def np_to_csv(arr, filename, cols):
    values = arr[0]
    bins = arr[1]
    su = sum(values)
    with open(filename, "w") as fin:
        fin.write("%s\n" % (",".join(cols)))
        for i in range(0, len(bins)-1):
            fin.write("%s,%s,%s\n" %
                      ("-".join(map(str, [round(bins[i], 2), round(bins[i+1], 2)])), values[i], round(float(values[i])/su, 3)))


def fraction_domain_ind(varcol, genes_cond, has_dom, has_merged, source_gene_object, fileprefix):
    # fileprefic will have follwing chatracteristyics
    # cath, pfam setup ? continuous discontinuos set up
    #
    prefix_func = "Individual_"
    with open(source_gene_object) as fin:
        has_gene = pickle.load(fin)
    with open(genes_cond) as fin:
        condition = pickle.load(fin)
    exon_counter = []
    total_alt_exons = 0
    total_const_exons = 0
    ex_dom_vals = []
    exon_dom_participation = {"G": [], "A": []}
    finfile = open(fileprefix+prefix_func +
                   "Exon_comprehensive_overlapp_domains.tabdelim", "w")
    finfile.write(
        'Gene\tVar\tLenVAr\tExonCountinVar\tExid\tExlen\tDomInVar\tDomsAffected\tDomsAffctedFrac\tWEF\n')
    fraction_largest = {(0, 0.25): '', (0.25, 0.50): '',
                        (0.50, 0.75): '', (0.75, 1.00): ''}
    exon_fraction_all = {(0, 0.25): '', (0.25, 0.50): '',
                         (0.50, 0.75): '', (0.75, 1.00): ''}
    c = 0
    for gene in has_dom:
        if gene in condition:
            c += 1
            domains = has_dom[gene]
            mat_domain = [[] for i in domains]
            isoform_object = ''
            flag_trans = 0
            for transcripts in has_gene[gene].transcripts:
                if transcripts.ID == has_merged[gene][varcol]:
                    isoform_object = transcripts
                    flag_trans = 1
                    break
            if flag_trans:
                # print "in"
                ex_scann = 0
                leniso = sum([i.length for i in isoform_object.exons])
                for ex in isoform_object.exons:
                    if i.length > 0:
                        ex_cases = []
                        lenexo = ex.length
                        ex_const = i.cons_flag
                        if ex_const == "G":
                            total_const_exons += 1
                        else:
                            total_alt_exons += 1
                        keyparticipation = ex_const
                        ex_range = set(range(ex_scann, ex_scann+lenexo))
                        ex_scann += lenexo
                        for ind, val in enumerate(domains):
                            intersec = ex_range & val
                            if intersec:
                                fracval = round(
                                    float(len(intersec))/len(val), 2)
                                fracval_ex = round(
                                    float(len(intersec))/len(ex_range), 2)
                                mat_domain[ind] += [[len(intersec), lenexo,
                                                     fracval, ex_const]]
                                ex_cases += [str(fracval)]
                                ex_dom_vals += [(fracval_ex,
                                                 fracval, ex_const)]
                                exon_dom_participation[keyparticipation] += [
                                    str(gene)+ex.ID]
                            if ex_cases:
                                finfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene, isoform_object.ID, leniso, len(isoform_object.exons),
                                ex.ID, lenexo, len(domains), len(ex_cases), ",".join(ex_cases), ex.WEF))

                for i in mat_domain:
                    if i:
                        i.sort()
                        valTobeBinnedSpecific = i[-1][2]
                        valAlphabetSpecific = i[-1][-1]
                        exon_counter += [len(i)]
                        for key in fraction_largest:
                            if key[0] < valTobeBinnedSpecific <= key[1]:
                                fraction_largest[key] += valAlphabetSpecific
                                break
                        for j in i:
                            valTobeBinned = j[2]
                            valAlphabet = j[-1]
                            for key in exon_fraction_all:
                                if key[0] < valTobeBinned <= key[1]:
                                    exon_fraction_all[key] += valAlphabet
                                    break
    mini_csv(exon_fraction_all, fileprefix +
             prefix_func+"all_exons_contribution.csv")
    mini_csv(fraction_largest, fileprefix+prefix_func +
             "largest_exon_contribution.csv")
    finfile.close()
    with open(fileprefix+prefix_func+"summary.txt", "w") as fin:
        fin.write("Const_cases\nTotal_exons:%s\nDomain_overlaps:%s\n\nAlt_Cases\nTotal_exons:%s\nDomain_overlaps:%s" % (
            total_const_exons, len(set(exon_dom_participation['G'])), total_alt_exons,  len(set(exon_dom_participation['A']))))
    with open(fileprefix+prefix_func+"ex_dom_quadrants.csv", "w") as fin:
        fin.write("Exon_fraction,Domain_fraction,Exon_type\n")
        for i in ex_dom_vals:
            fin.write("%s,%s,%s\n" % (i[0], i[1], i[2]))

    np_arr1 = np.histogram(exon_counter, bins=[
        1, 2, 3, 4, 5, 6, 10, 15, 20, 25, 30, 35, 50, 101])
    np_to_csv(np_arr1, fileprefix+prefix_func+"exons_numbers_coding.csv",
              cols=["NumberOfExons", "Total", "Freq"])


def fraction_domain_pair(varcol, genes_cond, has_dom, has_merged, source_gene_object, fileprefix):
    with open(source_gene_object) as fin:
        has_gene = pickle.load(fin)
    with open(genes_cond) as fin:
        condition = pickle.load(fin)
    prefix_func = "Paired_"
    exon_counter = []
    ex_dom_vals = []
    total_alt_exons = 0
    total_const_exons = 0
    exon_dom_participation = {"G": [], "A": []}
    fraction_largest = {(0, 0.25): '', (0.25, 0.50): '',
                        (0.50, 0.75): '', (0.75, 1.00): ''}
    exon_fraction_all = {(0, 0.25): '', (0.25, 0.50): '',
                         (0.50, 0.75): '', (0.75, 1.00): ''}
    finfile = open(fileprefix+prefix_func +
                   "Exon_comprehensive_overlapp_domains.tabdelim", "w")
    finfile.write(
        'Gene\tVar\tLenVAr\tExonCountinVar(Merged)\tExid\tExlen\tDomInVar\tDomsAffected\tDomsAffctedFrac\tWEF\n')

    c = 0
    for gene in has_dom:
        # if gene == 9:

        c += 1
        domains = has_dom[gene]
        mat_domain = [[] for i in domains]
        isoform_object = ''
        flag_trans = 0
        for transcripts in has_gene[gene].transcripts:
            if transcripts.ID == has_merged[gene][varcol]:
                isoform_object = transcripts
                flag_trans = 1
                break

        if gene in condition and flag_trans:
            # print gene, c
            # print isoform_object.ID
            const_exons = has_gene[gene].connst_togetherness_coding()
            # print "const_exons"
            const_exons_ref = []
            for i in const_exons:
                templis = []
                for j in i:
                    templis += [j]
                const_exons_ref += [templis]
            const_exons = const_exons_ref[:]
            list_exons = isoform_object.exons
            # print list_exons
            indexons = {i: [i.length, i.ID, i.cons_flag] for i in list_exons}
            # has_ref = {i[0]: [i, sum([indexons[j] for j in i])] for i in const_exons}
            # print indexons,"indexons"
            has_ref = {}
            for i in const_exons:
                if i:
                    has_ref[i[0]] = [i[1:]]
                    lencons = 0
                    consid = ''
                    for j in i:
                        lencons += indexons[j][0]
                        consid += ":%s" % indexons[j][1]
                    has_ref[i[0]] += [lencons, consid]
            # print has_ref
            merged_list = []
            for i in list_exons:
                if i not in has_ref:
                    merged_list += [indexons[i]]
                else:
                    tempMerge = [has_ref[i][1], has_ref[i][2], "G"]
                    merged_list += [tempMerge]
                    for j in has_ref[i][0]:
                        list_exons.remove(j)
            # print merged_list
            # print has_ref
            # print list_exons
            # print merged_list
            leniso = sum([i[0] for i in merged_list])
            ex_scann = 0
            for ex in merged_list:
                lenexo = ex[0]
                ex_cases = []
                ex_const = ex[2]
                keyparticipation = ex_const
                if ex_const == "G":
                    total_const_exons += 1
                else:
                    total_alt_exons += 1
                ex_range = set(range(ex_scann, ex_scann+lenexo))
                ex_scann += lenexo
                for ind, val in enumerate(domains):
                    intersec = ex_range & val
                    if intersec:
                        fracval_ex = round(
                            float(len(intersec))/len(ex_range), 2)
                        fracval = round(float(len(intersec))/len(val), 2)
                        mat_domain[ind] += [[len(intersec), lenexo,
                                             fracval, ex_const]]
                        ex_cases += [str(fracval)]
                        ex_dom_vals += [(fracval_ex,
                                         fracval, ex_const)]
                        exon_dom_participation[keyparticipation] += [
                            str(gene)+ex[1]]
                if ex_cases:
                    finfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene, isoform_object.ID, leniso, len(merged_list), ex[1], lenexo, len(domains), len(ex_cases), ",".join(ex_cases), ex[2]))
            for i in mat_domain:
                if i:
                    i.sort()
                    valTobeBinnedSpecific = i[-1][2]
                    valAlphabetSpecific = i[-1][-1]
                    exon_counter += [len(i)]
                    for key in fraction_largest:
                        if key[0] < valTobeBinnedSpecific <= key[1]:
                            fraction_largest[key] += valAlphabetSpecific
                            break
                    for j in i:
                        valTobeBinned = j[2]
                        valAlphabet = j[-1]
                        for key in exon_fraction_all:
                            if key[0] < valTobeBinned <= key[1]:
                                exon_fraction_all[key] += valAlphabet
                                break
    mini_csv(exon_fraction_all, fileprefix +
             prefix_func+"all_exons_contribution.csv")
    mini_csv(fraction_largest, fileprefix+prefix_func +
             "largest_exon_contribution.csv")
    finfile.close()
    with open(fileprefix+prefix_func+"summary.txt", "w") as fin:
        fin.write("Const_cases\nTotal_exons:%s\nDomain_overlaps:%s\n\nAlt_Cases\nTotal_exons:%s\nDomain_overlaps:%s" % (
            total_const_exons, len(set(exon_dom_participation['G'])), total_alt_exons,  len(set(exon_dom_participation['A']))))
    with open(fileprefix+prefix_func+"ex_dom_quadrants.csv", "w") as fin:
        fin.write("Exon_fraction,Domain_fraction,Exon_type\n")
        for i in ex_dom_vals:
            fin.write("%s,%s,%s\n" % (i[0], i[1], i[2]))
    np_arr1 = np.histogram(exon_counter, bins=[
        1, 2, 3, 4, 5, 6, 10, 15, 20, 25, 30, 35, 50, 101])
    np_to_csv(np_arr1, fileprefix+prefix_func+"exons_numbers_coding.csv",
              cols=["NumberOfExons", "Total", "Freq"])


def discontinuous_domains(varcol, genes_cond, has_dom, has_merged, source_gene_object, fileprefix):
    # fileprefic will have follwing chatracteristyics
    # cath, pfam setup ? continuous discontinuos set up
    #
    prefix_func = "Discontinuous_"
    with open(source_gene_object) as fin:
        has_gene = pickle.load(fin)
    with open(genes_cond) as fin:
        condition = pickle.load(fin)
    total_alt_exons = 0
    total_const_exons = 0
    exon_dom_participation_ind = {"Long": '', "Short": ''}
    fraction_diff_ind = {(-1, 0.25): ['', []], (0.25, 0.50): ['', []],
                         (0.50, 0.75): ['', []], (0.75, 1.00): ['', []]}
    exon_dom_participation_const = {"Long": '', "Short": ''}
    fraction_diff_const = {(-1, 0.25): ['', []], (0.25, 0.50): ['', []],
                           (0.50, 0.75): ['', []], (0.75, 1.00): ['', []]}

    for gene in has_dom:
        if gene in condition:
            domains = has_dom[gene]
            isoform_object = ''
            flag_trans = 0
            for transcripts in has_gene[gene].transcripts:
                if transcripts.ID == has_merged[gene][varcol]:
                    isoform_object = transcripts
                    flag_trans = 1
                    break
            ex_scann = 0
            # leniso = sum([i.length for i in isoform_object.exons])
            if flag_trans:
                for ex in isoform_object.exons:
                    lenexo = ex.length
                    ex_const = ex.cons_flag
                    if ex_const == "G":
                        total_const_exons += 1
                    else:
                        total_alt_exons += 1
                    keyparticipation = ex_const
                    ex_range = set(range(ex_scann, ex_scann+lenexo))
                    ex_scann += lenexo
                    for ind, val in enumerate(domains):
                        parts_dom = [set(group)
                                     for group in mit.consecutive_groups(list(val))]
                        dis_parts = [
                            [round(float(len(i & val))/len(val), 2), i] for i in parts_dom]
                        dis_parts.sort()
                        # print dis_parts, "disparts"
                        for partval in dis_parts:
                            intersec = ex_range & partval[1]
                            if intersec:
                                fracval = round(
                                    float(len(intersec))/len(partval[1]), 2)
                                partval.insert(1, fracval)
                            else:
                                partval.insert(1, 0)
                        if dis_parts[-1][0] >= 0.60 and dis_parts[-1][1] >= 0.70:
                            if len(dis_parts) > 1:
                                valTobeBinnedSpecific = round(
                                    dis_parts[-1][1]-dis_parts[-2][1], 2)
                                valAlphabetSpecific = keyparticipation
                                for key in fraction_diff_ind:
                                    if key[0] < valTobeBinnedSpecific <= key[1]:
                                        fraction_diff_ind[key][0] += valAlphabetSpecific
                                        fraction_diff_ind[key][1] += [
                                            str(gene)+"_"+ex.ID]
                                        break
                            else:
                                print gene, "gene"
                        if dis_parts[-1][0] >= 0.60:
                            for inddom, valdom in enumerate(dis_parts[::-1]):
                                if valdom[1]:

                                    if inddom == 0:
                                        exon_dom_participation_ind['Long'] += keyparticipation
                                    else:
                                        exon_dom_participation_ind['Short'] += keyparticipation

    mini_csv({i: fraction_diff_ind[i][0] for i in fraction_diff_ind},
             fileprefix+prefix_func+"Individual_longer_fraction_difference.csv")
    mini_csv(exon_dom_participation_ind, fileprefix +
             prefix_func+"Individual_overlaptype.csv")
    with open(fileprefix+prefix_func+"Individual_longer_fraction_difference_summary_genes.csv", "w") as fin:
        for i in fraction_diff_ind:
            fin.write("%s\t%s\n" %
                      (i, ":".join(list(set(fraction_diff_ind[i][1])))))
    print "parta_done"

    for gene in has_dom:
        if gene in condition:
            domains = has_dom[gene]
            flag_trans = 0
            for transcripts in has_gene[gene].transcripts:
                if transcripts.ID == has_merged[gene][3]:
                    isoform_object = transcripts
                    flag_trans = 1
                    break

            if flag_trans:
                # print isoform_object.ID
                const_exons = has_gene[gene].connst_togetherness_coding()
                list_exons = isoform_object.exons
                # print list_exons
                indexons = {i: [i.length, i.ID, i.cons_flag]
                            for i in list_exons}
                # has_ref = {i[0]: [i, sum([indexons[j] for j in i])] for i in const_exons}
                has_ref = {}
                for i in const_exons:
                    if len(i) > 0:
                        has_ref[i[0]] = [i[1:]]
                        lencons = 0
                        consid = ''
                        for j in i:
                            lencons += indexons[j][0]
                            consid += ":%s" % indexons[j][1]
                        has_ref[i[0]] += [lencons, consid]
                # print has_ref
                merged_list = []
                for i in list_exons:
                    if i not in has_ref:
                        merged_list += [indexons[i]]
                    else:
                        tempMerge = [has_ref[i][1], has_ref[i][2], 'G']
                        merged_list += [tempMerge]
                        for j in has_ref[i][0]:
                            list_exons.remove(j)
                # print merged_list, "merged_list"
                ex_scann = 0
                for ex in merged_list:
                    # print ex, 'ex'
                    lenexo = ex[0]
                    exid = ex[1]
                    ex_const = ex[2]
                    keyparticipation = ex_const
                    if ex_const == "G":
                        total_const_exons += 1
                    else:
                        total_alt_exons += 1
                    ex_range = set(range(ex_scann, ex_scann+lenexo))
                    ex_scann += lenexo
                    for ind, val in enumerate(domains):
                        parts_dom = [set(group)
                                     for group in mit.consecutive_groups(list(val))]
                        dis_parts = [
                            [round(float(len(i & val))/len(val), 2), i] for i in parts_dom]
                        dis_parts.sort()
                        for partval in dis_parts:
                            intersec = ex_range & partval[1]
                            if intersec:
                                fracval = round(
                                    float(len(intersec))/len(partval[1]), 2)
                                partval.insert(1, fracval)
                            else:
                                partval.insert(1, 0)
                        if dis_parts[-1][0] >= 0.60 and dis_parts[-1][1] >= 0.70:
                            if len(dis_parts) > 1:
                                valTobeBinnedSpecific = round(
                                    dis_parts[-1][1]-dis_parts[-2][1], 2)
                                valAlphabetSpecific = keyparticipation
                                for key in fraction_diff_const:
                                    if key[0] < valTobeBinnedSpecific <= key[1]:
                                        fraction_diff_const[key][0] += valAlphabetSpecific
                                        fraction_diff_const[key][1] += [
                                            str(gene)+"_"+exid]
                                        break
                            else:
                                print gene, "part2"
                        if dis_parts[-1][0] >= 0.60:
                            for inddom, valdom in enumerate(dis_parts[::-1]):
                                if valdom[1]:
                                    if inddom == 0:
                                        exon_dom_participation_const['Long'] += keyparticipation
                                    else:
                                        exon_dom_participation_const['Short'] += keyparticipation

    mini_csv({i: fraction_diff_const[i][0] for i in fraction_diff_const},
             fileprefix+prefix_func+"Paired_longer_fraction_difference.csv")
    mini_csv(exon_dom_participation_const, fileprefix +
             prefix_func+"Paired_overlaptype.csv")
    with open(fileprefix+prefix_func+"Paired_longer_fraction_difference_summary_genes.csv", "w") as fin:
        for i in fraction_diff_const:
            fin.write("%s\t%s\n" %
                      (i, ":".join(list(set(fraction_diff_ind[i][1])))))

from shutil import rmtree
from progress.bar import Bar
import cPickle as pickle
import os
import re
import Bio
import subprocess
import Bio.SeqUtils
import sys
print "--> In Prediction_Assignment_module_successfully"


def description():
    print '''
    ::: This module has following functions
    1. assigner(done_dir, fastadir_done, pool_of_left, dir_to_write)
       # assigns prediction to unpredicted sequences based on already predcted pool
    2. def aaseq(has_gene_cood, aaseq_fasta_dir, aaseq_dir):
       # write aaseq of sequences as exons
    3. def ss_to_exons(has_gene_cood, ssrawdir, ssfastadir, ssExonsWriteDir):
       # write sstype of proteins as sequence shape
    4. principal_isoform(has_gene_var, aaseq_dump):
       will give transcript as largest number of exons as primary, preferably NP
    5. def prerequisite_stride(has_gene_var, pdb_dir, stride_reformat_dir, pid_dir):
        # will take pdb_dir as input directory and astride reformat dir as parent directory
        # to store files in pssp format
    6. def stride_assigner(has_gene_cood, stride_dir_ssformat, ssExonsWriteDir)
        # has_gene_cood as input and then will write files present in stride_dir_ssformat
        # to pickle in ssExonsWriteDir
    5. future functions
       stride ss writer and disorder writer as exons
    '''


def hidden_fasta_seq(add):
    with open(add) as fin:
        dat = fin.read()
        dat = re.sub(r"\r", "", dat)
        if ">" == dat[0]:
            key = "".join(dat.split("\n")[1:]).strip()
        else:
            key = "".join(dat.split("\n")).strip()
    return key


def hidden_copy_file(src, dest):
    with open(src) as fin:
        dat = fin.read()
    with open(dest, "w") as fin:
        fin.write(dat)


def assigner(done_dir, fastadir_done, pool_of_left, dir_to_write, org_id, where_to_write_rep_of_unassigned):
    '''
    read PID of done, strip .dat.ss from end
    for the read ID's read the fasta sequence and save it in the ncbi_hashes_in_memory
    has as seq value as PID
    read the parent pool,
        id shjouldnt repeat
        read the fasta, ad chcek if it has been Predicted
        then write and count,
        else store, yet to be retrieved
    '''

    unknown = {}
    copy_done = 0
    done_dir_id = {re.sub(r'.dat.ss', '', i): 0 for i in os.listdir(
        done_dir) if re.search(r'\.dat.ss$', i)}
    # print done_dir_id.keys()[:5]
    bar1 = Bar('Adding sequences to memory:', max=len(done_dir_id))
    done_hash = {}
    for i in done_dir_id:
        bar1.next()
        seq = os.path.join(fastadir_done,"%s" % i)
        try:
            key = hidden_fasta_seq(seq)
        except:
            continue
        # print i, key
        done_hash[key] = i
    bar1.finish()

    # print "isit"
    # print done_hash.keys()[:5]
    bar2 = Bar('Predicting sequences for left proteins:',
               max=len(os.listdir(pool_of_left)))
    for i in os.listdir(pool_of_left):
        bar2.next()
        # pool of left is a dir yet to be Predicted
        if i not in done_dir_id:
            seq = os.path.join(pool_of_left,"%s" % i)
            key = hidden_fasta_seq(seq)
            if key not in done_hash:
                unknown[key] = i
            else:
                copy_done += 1
                hidden_copy_file(os.path.join(done_dir,"%s.dat.ss" %
                                 done_hash[key]), os.path.join(dir_to_write,"%s.dat.ss" % i))
    bar2.finish()
    temsum = 0
    for temi in unknown:
        if len(temi) > 3000:
            temsum += 1
    print "-> Newly predicted sequences are: %s" % copy_done
    print "-> Unassigned sequences are: %s, (%s are greter than 3000)" % (len(unknown), temsum)
    if os.path.isdir(where_to_write_rep_of_unassigned):
        rmtree(where_to_write_rep_of_unassigned)
        os.makedirs(where_to_write_rep_of_unassigned)
    for i in unknown:
        if len(i) <= 3000:
            hidden_copy_file(os.path.join(pool_of_left,unknown[i]), where_to_write_rep_of_unassigned+"%s" % unknown[i])
    print "-> written pending sequences in %s" % where_to_write_rep_of_unassigned
    return unknown


def hidden_ret_array_fromss(ss_file, path_ssp_dir_raw, conf=None):
    arr = []
    try:
        with open(os.path.join(path_ssp_dir_raw,"%s" % ss_file)) as fin:

            data = fin.read().split("\n")
            del data[0]
            if not len(data[-1].split()) == 5:
                del data[-1]
            for i in data:
                a = i.split()

                if conf is not None:
                    if a[2] == 'H':
                        a[2] = 'H' if float(a[4]) >= conf else 'C'
                    arr += [a[2]]
                else:
                    arr += [a[2]]
                arr += [a[2]]
        return arr

    except Exception as E:
        print E
        return 0

def hidden_ret_array_fromdis(dis_file, path_dis_dir_raw):
    arr = []
    try:
        with open(os.path.join(path_dis_dir_raw,"%s" % dis_file)) as fin:

            data = [i for i in fin.read().split("\n") if len(i)>0 and i[0]!='#']
            for i in data:
                a = i.split()
                if round(float(a[2]),4) > 0.5000:
                    arr += "D"
                else:
                    arr += "S"
        return arr

    except Exception as E:
        print E
        return 0

def hidden_align(lis1, pid1, path_pid_raw):
    try:
        mat_c = []
        seqfile = os.path.join(path_pid_raw,"%s" % pid1)
        fas1 = hidden_fasta_seq(seqfile)
        a = 0
        b = 0
        for i in lis1:
            if i[2] == "C":
                c = (abs(i[1][1]-i[1][0])+1)/3
                # i[0] is nc cood, i[1] is coding cood
                b += c
                # mat_c+=[[fas1[a:b],tuple(i[0]),tuple(i[1]),i[2]]]

                mat_c += [[fas1[a:b], tuple(i[0]), tuple(i[1]), i[2]]]
                # mat_c+=[[fas1[a:b],i[0],tuple(i[1]),i[2]]]
                a = b
                # now mat[1][0] is aa seq,
                # mat[1][1] is tuple of coods,
                # mat[1][1][0] is tuple of nccoods's, mat[1][1][1]  is tuple of coding cood's
            else:
                # mat_c+=[["",tuple(i[0]),tuple(i[1]),i[2]]]
                mat_c += [["", tuple(i[0]), tuple(i[1]), i[2]]]
                # mat_c+=[["",i[0],tuple(i[1]),i[2]]]
        return mat_c
    except Exception as E:
        print E, "**->Hidden_align"
        return 0


def hidden_ss_aligner(ss_arr, lis1):
    mat_c = []
    fas1 = "".join(ss_arr)
    a = 0
    b = 0
    for i in lis1:
        if i[2] == "C":
            c = (abs(i[1][1]-i[1][0])+1)/3
            b += c
            mat_c += [[fas1[a:b], tuple(i[0]), tuple(i[1]), i[2]]]
            # mat_c+=[[fas1[a:b],i[0],tuple(i[1]),i[2]]]
            a = b
        else:
            mat_c += [["", tuple(i[0]), tuple(i[1]), i[2]]]
            # because tuple in nc will have only one exouic span not duplicates
            # mat_c+=[["",i[0],tuple(i[1]),i[2]]]
    return mat_c


def hidden_mod3(lis1):
    lis = lis1[:]
    # print lis1
    for i in range(0, len(lis)-1):
        if lis[i][2] == "C":
            # means coding exon
            mod = ((lis[i][1][1]-lis[i][1][0])+1) % 3
            # print lis[i][1],mod
            if mod != 0:
                if mod == 1:
                    lis[i][1][1] -= 1
                    # print "mod1 change sub ", lis[i][1]
                    for j in range(i+1, len(lis)):
                        # searching immediate coding exon
                        if lis[j][2] == "C":  # and lis[j][0]!=[0,0]:
                            lis[j][1][0] -= 1
                            # print "added ",lis[j][1]
                            # donate the boundary to next coding exon
                            if lis[i][1][1]-lis[i][1][0] == -1:
                                # chcek for condition if sum is not -1 after donation in perset exon.
                                # if it is, pass the nc rna cood to very next nnon ncRNA and make thi 0,0
                                lis[i][2] = "M"
                            break
                        break

                else:
                    lis[i][1][1] += 1
                    # print "mod 2change sub ", lis[i][1]
                    for j in range(i+1, len(lis)):
                        if lis[j][2] == "C":  # and lis[j][0]!=[0,0]:
                            lis[j][1][0] += 1
                            if lis[j][1][1]-lis[j][1][0] == -1:
                                # print ";"
                                lis[j][2] = "M"
                                break
                            break

    # after processing now reove last three nuclotides from last coding eoxn in lis
    for i in range(len(lis)-1, -1, -1):
        # reverse reading a list
        if lis[i][2] == "C":
            # print "lis[i]",lis[i]
            lis[i][1][1] -= 3
            # print "lis[i]",lis[i]
            if lis[i][1][1]-lis[i][1][0] == -1:
                # print "?"
                lis[i][2] = "M"
                break
            break

    lis_upd = lis[:]

    for i in lis_upd:
        i[0].sort()

    return lis_upd
    # hoping it to be correct


def aaseq(has_gene_cood, aaseq_fasta_dir, aaseq_dir):
    bar = Bar('Processing aaseq into exons, gene_wise:',
              max=len(has_gene_cood))
    aaseq_dir_hash = {i: 0 for i in os.listdir(aaseq_dir)}
    var_inside = 0
    total_var = 0
    for gene in has_gene_cood:
        bar.next()
        try:
            pid = has_gene_cood[gene].keys()
            # get pid variants of gene ids
            for var in pid:
                total_var += 1
                if var not in aaseq_dir_hash:
                    # if gene == 100534673 and var == "XP_009302357.2":
                    exon1temp = []
                    k1_keys = has_gene_cood[gene][var]
                    # print k1_keys
                    # k1keys will have exon seq no,, its coods, then coding sattasitus
                    k1_keys.sort()
                    # k1_keys will have [[1,(10,20)(18,20),"N"] and so on in increasng order]
                    for coods in k1_keys:
                        # template=has[f_f][pid[k1]][k1_keys[k_k1]]
                        exon1temp += [[list(coods[1]),
                                       list(coods[2]), coods[3]]]
                    exon1mod = hidden_mod3(exon1temp)
                    # print exon1temp, "temp"
                    # print exon1mod, "mod"
                    mat_c = hidden_align(exon1mod, var, aaseq_fasta_dir)
                    # print "\n%s\n" % mat_c
                    if mat_c:
                        # print mat_c
                        with open(os.path.join(aaseq_dir,"%s" % var), "w") as outter:
                            pickle.dump(mat_c, outter,
                                        protocol=pickle.HIGHEST_PROTOCOL)
                        # break
                    else:
                        print "err", gene, var
                else:
                    var_inside += 1

        except Exception as E:
            print "in seq_to_exons()", E
        # break
    bar.finish()
    print "Out of total %s variants, %s were already written, use freshStart mode 2 for rewrite all" % (
        total_var, var_inside)


def ss_to_exons(has_gene_cood, ssrawdir, ssExonsWriteDir, threshold):
    # count_left1 = []
    # count_left2 = 0
    # print ssrawdir
    # print ssExonsWriteDir
    # return
    bar = Bar('Processing ssseq into exons, gene_wise:',
              max=len(has_gene_cood))
    ssseq_dir_hash = {i: 0 for i in os.listdir(ssExonsWriteDir)}
    print len(ssseq_dir_hash), "files in pickle"
    print len([i for i in (os.listdir(ssrawdir))
               if ".dat.ss" in i]), "raw _set"
    for gene in has_gene_cood:
        bar.next()
        try:
            pid = has_gene_cood[gene].keys()
            # get pid variants of gene ids
            for var in pid:
                if var not in ssseq_dir_hash:
                    # count_left1 += [var]
                    ss_file = var+".dat.ss"
                    if os.path.isfile(os.path.join(ssrawdir,"%s" % ss_file)):
                        # print "in"
                        # count_left2 += 1
                        exon1temp = []
                        k1_keys = has_gene_cood[gene][var]
                        # k1keys will have exon seq no,, its coods, then coding sattasitus
                        k1_keys.sort()
                        # k1_keys will have [[1,(10,20)(18,20),"N"] and so on in increasng order]
                        for coods in k1_keys:
                            # template=has[f_f][pid[k1]][k1_keys[k_k1]]
                            exon1temp += [[list(coods[1]),
                                           list(coods[2]), coods[3]]]
                        exon1mod = hidden_mod3(exon1temp)
                        array_ss_type = hidden_ret_array_fromss(
                            ss_file, ssrawdir, threshold)
                        mat_ss = hidden_ss_aligner(array_ss_type, exon1mod)
                        if mat_ss:
                            # print "in2"
                            with open(os.path.join(ssExonsWriteDir,"%s" % var), "w") as outter:
                                pickle.dump(mat_ss, outter,
                                            protocol=pickle.HIGHEST_PROTOCOL)

        except Exception as e:
            print e, gene, "sseq_writer"
        # break
    bar.finish()

def disorder_to_exons(has_gene_cood, disrawdir, disExonsWriteDir):
    bar = Bar('Processing disorder into exons, gene_wise:',
              max=len(has_gene_cood))
    disseq_dir_hash = {i: 0 for i in os.listdir(disExonsWriteDir)}
    print len(disseq_dir_hash), "files in pickle"
    print len(os.listdir(disrawdir)), "raw _set"
    for gene in has_gene_cood:
        bar.next()
        try:
            pid = has_gene_cood[gene].keys()
            # get pid variants of gene ids
            for var in pid:
                if var not in disseq_dir_hash:
                    # count_left1 += [var]
                    dis_file = var
                    if os.path.isfile(os.path.join(disrawdir,"%s" % dis_file)):
                        exon1temp = []
                        k1_keys = has_gene_cood[gene][var]
                        # k1keys will have exon seq no,, its coods, then coding sattasitus
                        k1_keys.sort()
                        # k1_keys will have [[1,(10,20)(18,20),"N"] and so on in increasng order]
                        for coods in k1_keys:
                            # template=has[f_f][pid[k1]][k1_keys[k_k1]]
                            exon1temp += [[list(coods[1]),
                                           list(coods[2]), coods[3]]]
                        exon1mod = hidden_mod3(exon1temp)
                        array_dis_type = hidden_ret_array_fromdis(
                            dis_file, disrawdir)

                        mat_dis = hidden_ss_aligner(array_dis_type, exon1mod)
                        if mat_dis:
                            # print "in2"
                            with open(os.path.join(disExonsWriteDir,"%s" % var), "w") as outter:
                                pickle.dump(mat_dis, outter,
                                            protocol=pickle.HIGHEST_PROTOCOL)

        except Exception as e:
            print e, gene, "disseq_writer"
        # break
    bar.finish()


def principal_isoform(has_gene_var, aaseq_dump):
    pi = {}
    for gene in has_gene_var:
        if 1:
            # print gene
            alli = []
            for var in has_gene_var[gene]:
                # print gene, var
                with open(os.path.join(aaseq_dump,"%s" % var), "rb") as fin:
                    # print var
                    # print pickle.load(fin)
                    dat = []
                    for i in pickle.load(fin):
                        if i[3] == 'C':
                            dat += [i[0]]
                    # dat = [i[0] for i in pickle.load(fin) if i[3] == "C"]
                length = 0
                # print "here", gene
                # print dat
                for i in dat:
                    length += len(i)
                if "NP" in var:
                    alli += [[1, len(dat), length, var]]
                else:
                    alli += [[0, len(dat), length, var]]

            alli.sort()
            # print alli
            pi[gene] = alli[-1][-1]
    # sys.exit()
    return pi


def hidden_modify_output_stride_pss_form(res, pid, pid_add, stride_dir):
    p_seq = hidden_fasta_seq(os.path.join(pid_add,pid))
    store_seq_original = {}
    core = res.split(
        "REM  |---Residue---|    |--Structure--|   |-Phi-|   |-Psi-|  |-Area-|      ~~~~")[-1]
    for lin in core.split("\n"):
        vals = lin.split()
        if len(vals) > 8:
            aa3 = vals[1]
            converteraa_val = Bio.SeqUtils.IUPACData.protein_letters_3to1[aa3[0].upper(
            )+aa3[1].lower()+aa3[2].lower()]
            store_seq_original[int(vals[3])] = (converteraa_val, vals[5])
    hidden_matcherandedit(store_seq_original, p_seq, pid, stride_dir)


def hidden_matcherandedit(sto_seq, prot_seq, pid, stride_dir):
    try:
        ss_trans = {"H": "H", "G": "H", "I": "H", "B": "E",
                    "T": "C", "S": "C", "C": "C", "E": "E", "b": "E"}
        with open(os.path.join(stride_dir,"%s.ss_3" % pid), "w") as fout_3:
            for i in range(0, len(prot_seq)):
                if i in sto_seq:
                    fout_3.write("%5s %s %s 0 0 0\n" %
                                 (i, sto_seq[i][0], ss_trans[sto_seq[i][1]]))
                else:
                    fout_3.write("%5s %s Z 0 0 0\n" % (i, prot_seq[i]))
    except Exception as EM:
        print "  in matchandedit funct: %s\t%s" % (pid, EM)


def prerequisite_stride(stride_location, has_gene_var, save_structure_dir, stride_reformat_dir, pid_dir):
    '''
    create a one more function, that will run striode on the created pdb's and add output to the
    requested files to further proceeed, in both orgonal and ss fromat os psspred
    '''
    # save_structure_dir must have file formatted in manner "var+chain+pdb.pdb"
    print (os.listdir(save_structure_dir)[:5])
    file_list_string_hash = {
        "_".join(i.split("_")[:2]): i for i in os.listdir(save_structure_dir)}
    for gene in has_gene_var:
        flag=0
        for var in has_gene_var[gene]:
            if var in file_list_string_hash:
                # means pdb for that file exists
                # check iof stride output for that exists
                # else run stride on them
                if not os.path.isfile(os.path.join(stride_reformat_dir,var+".ss")):
                    try:
                        #print ([stride_location, os.path.join(save_structure_dir,file_list_string_hash[var])])
                        #print (gene,var,)
                        result_stride = subprocess.check_output(
                            [stride_location, os.path.join(save_structure_dir,file_list_string_hash[var])])
                        hidden_modify_output_stride_pss_form(
                            result_stride, var, pid_dir, stride_reformat_dir)
                    except Exception as E:
                        print "Raising Exception******* in prediction_assignment_ss module, () is prerequisite_stride,\n"
                        print "during running stride: %s" % E,
                        print gene, var

        


def stride_assigner(has_gene_cood, stride_dir_ssformat, ssExonsWriteDir):
    '''
    has to first check if stride for all thepdb's are made and keopt
    '''
    bar = Bar('Processing sstride into exons, gene_wise:',
              max=len(has_gene_cood))
    stride_seq_dir_hash = {i: 0 for i in os.listdir(ssExonsWriteDir)}
    for gene in has_gene_cood:
        bar.next()
        try:
            for var in has_gene_cood[gene]:
                if os.path.isfile(os.path.join(stride_dir_ssformat,"%s.ss_3" % var)) and var not in stride_seq_dir_hash:
                    exon1temp = []
                    k1_keys = has_gene_cood[gene][var]
                    k1_keys.sort()
                    for coods in k1_keys:
                        exon1temp += [[list(coods[1]),
                                       list(coods[2]), coods[3]]]
                    exon1mod = hidden_mod3(exon1temp)
                    ss_file1 = var+".ss_3"
                    array_ss_type1 = hidden_ret_array_fromss(
                        ss_file1, stride_dir_ssformat)
                    mat_ss1 = hidden_ss_aligner(array_ss_type1, exon1mod)
                    if mat_ss1:
                        with open(os.path.join(ssExonsWriteDir,"%s" % var), "w") as outter:
                            pickle.dump(mat_ss1, outter,
                                        protocol=pickle.HIGHEST_PROTOCOL)
        except Exception as EE:
            print "error in strideassigner(), %s" % EE

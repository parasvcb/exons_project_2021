import re
import more_itertools as mit
int_patt = re.compile(r'^\d+$')


class Gene:
    #print ('gene imported')

    def __init__(self, txid, entrezid, ensembleid, localization, description, go):
        self.txid = txid
        self.ID = entrezid
        self.ensembleId = ensembleid
        self.localization = localization
        # 1 is membrabous protein, 0 in globular protein
        self.detail = description
        self.go = go

    def transcripts(self, translist):
        self.transcripts = translist
        for i in translist:
            if i.PI:
                self.PI = i
                break

    def exons(self, exonlist):
        self.exons = exonlist
        # print ("ROUTINE classes exons")
        # print [i.ID for i in exonlist]
        id_ex_has = {
            ".".join(i.ID.split(".")[3:]): i for i in exonlist if i.ID[0] != 'R'}
        # print id_ex_has
        for i in exonlist:
            # print i.ID, "bef"
            if re.match(r'^[TUMD]\.[+-]?[0-n]+\.[AFG]\.[0-n]+\.[ncb]', i.ID):
                # print i.ID, "aaft"
                parent_id = ".".join(
                    re.sub(r'[ncb]\.[0-n]$', '0.0', i.ID).split(".")[3:])
                # print i.ID, parent_id, "in"
                if parent_id in id_ex_has:
                    i.parent = id_ex_has[parent_id]
                else:
                    # print parent_id, i.ID, "not"
                    i.parent = False
                    print "ERR*******, self.ID is:%s, parent err" % self.ID
            else:
                i.parent = False
        # print "doone"

    def description(self):
        print ("GeneName: %s, Id: %s, Txid: %s, Ensemble: %s, Localization: %s, Go are: %s") % (
            self.detail, self.ID, self.txid, self.ensembleId, self.localization, self.go)
        print ("It has follwing transcripts: (Total: %s)") % (
            len(self.transcripts))
        print ([i.description() for i in self.transcripts])
        print ([i.description() for i in self.exons])

    def connst_togetherness_coding(self):
        transcripts_exons_matrix = [
            j.exons for j in self.transcripts]
        # going for a transcript if it doesnt ave rrtained intron in it
        if len(transcripts_exons_matrix) <= 1:
            return False
        mat = []
        cons_pairs_total = []
        # print (len(transcripts_exons_matrix),'lenmat')
        # print transcripts_exons_matrix
        for var in transcripts_exons_matrix:
            cons_exon_id = []
            cons_exon_no = []
            for ind, val in enumerate(var):
                # this wil iterate the exons of very first transcript
                if re.match(r'^T\.[1]\.G', val.ID):
                    if val.length > 0:
                        cons_exon_no += [ind]
                    else:
                        cons_exon_no += [-10]

                else:
                    cons_exon_no += [-10]
                cons_exon_id += [val]

            pairs = [list(group)
                     for group in mit.consecutive_groups(cons_exon_no)]

            pairs_ref = [i for i in pairs if sum(i) >= 0]

            pairs_id = []
            for i in pairs_ref:
                temp = []
                for j in i:
                    temp += [cons_exon_id[j]]
                pairs_id += [temp]
            mat += [pairs_id]
        '''
        for i in mat:
            for j in i:
                print "p", [k.ID for k in j]
            print "\n\n"
        '''
        if mat:
            if len([1 for i in mat if i == mat[0]]) == len(mat):
                # all elements if matrix are same
                return mat[0]
            else:
                # print "true"
                trans = mat[0]
                '''
                for i in trans:
                    print "pp", [j.ID for j in i]
                #print "trans", [j.ID for j in trans]
                '''
                template = [set(tempi) for ti in mat[1:] for tempi in ti]
                # this is not sorted
                '''
                print len(template), "template"
                for j in template:
                    print "j", j
                    print "i", [i.ID for i in j]
                '''
                for pairs in trans:
                    subs = 0
                    sube = len(pairs)
                    Flag = True
                    while (Flag):
                        if len([1 for i2 in template if set(pairs[subs:sube]) <= i2]) == len(mat[1:]):
                            Flag = False
                            # print "1"
                            cons_pairs_total += [pairs[subs:sube]]
                        elif len([1 for i2 in template if set(pairs[subs+1:sube]) <= i2]) == len(mat[1:]):
                            Flag = False
                            # print "2"
                            cons_pairs_total += [pairs[subs+1:sube]]
                        elif len([1 for i2 in template if set(pairs[subs:sube-1]) <= i2]) == len(mat[1:]):
                            Flag = False
                            # print "3"
                            cons_pairs_total += [pairs[subs:sube-1]]
                        subs += 1
                        sube -= 1
                        if not pairs[subs:sube]:
                            break
                return cons_pairs_total
        else:
            return False

    def constitutive_alternate_with_freq(self):
        mat = []
        for i in self.transcripts:
            mat += [i.exons]
        '''
        b = set(mat[0])
        for i in range(1, len(mat)):
            b = b & set(mat[i])
        for i in b:
            i.Exon_characteristics(1)
        '''
        for i in self.exons:
            count = len([1 for c in mat if i in c])
            i.Exon_characteristics(round(float(count)/len(mat), 3))
    '''
    to this add symbol, go term, molecular location
    function, intron retention,
    const exons finder,
    alt exons frequency
    exons merger if any const exons needed
    '''


class Transcript():
    '''swissprot, ensemble, pdb filename , seqlen, pdb seqlen, retained intron, PI, FK'''
    TranscriptCount = 0

    def __init__(self, id, pi, seqlen, ensemble, swissprot):
        self.ID = id
        self.PI = pi
        self.seqlen = seqlen
        self.swiss_pid = swissprot
        self.ensemble_pid = ensemble
        self.transcriptkey = Transcript.TranscriptCount+1
        Transcript.TranscriptCount += 1
        self.Intron_status = False
        self.structure = False
        self.structure_len = False
        self.structure_file = False
        self.structureAF = False
        self.structure_lenAF = False
        self.structure_fileAF = False
        self.coverage = False
        self.pfam_list = False
        self.cath_list = False
        #

    def pfam_update(self, lis):
        self.pfam_list = lis
        self.pfam_list.sort()
    # lis[0] if coodspan, lis[1] is pfamid, lis[2] i scoveragelis2 is hmmtype, and then dtype

    def cath_update(self, lis):
        self.cath_list = lis
        self.cath_list.sort()

    def exons(self, exonlist):
        self.exons = exonlist

    def intron_status(self):
        self.Intron_status = True if [
            i for i in self.exons if i.IntronFlag] else False

    def description(self):
        print "---->In transcript Object"
        print "---->ID:%s, PI: %s, Length: %s, SP: %s, ENS: %s, IntronR: %s, StructureFile: %s" % (
            self.ID, self.PI, self.seqlen, self.swiss_pid, self.ensemble_pid, self.Intron_status, self.structure_file)
        print "**********> exons of this transcript are follows:"
        print[i.ID for i in self.exons]
        print "**********> ssseq is as follows %s" % ([i[1]
                                                       for i in Transcript.personal_ss_list(self)])
        print "**********> strideseq is as follows %s" % ([i[1]
                                                           for i in Transcript.personal_stride_list(self)])

    def personal_ss_list(self):
        lis = []
        for i in self.exons:
            if i.coding_span:
                for key in i.ssseq:
                    if self.ID in i.ssseq[key]:
                        lis += [[i, key]]
                        break
        return lis

    def personal_stride_list(self):
        lis = []
        for i in self.exons:
            if i.coding_span:
                for key in i.strideseq:
                    if self.ID in i.strideseq[key]:
                        lis += [[i, key]]
                        break
        return lis

    def pdb_info(self, pdb, pdblen):
        self.structure = True
        self.structure_file = pdb
        self.structure_len = pdblen
        ele = self.structure_file.split("_")
        # print ele, pdb
        self.pdb = ele[2]
        self.chain = ele[3].split('.')[0]
        self.coverage = round(float(pdblen)/self.seqlen, 2)
    
    def pdb_infoAF(self, pdb, pdblen):
        self.structureAF = True
        self.structure_fileAF = pdb
        self.structure_lenAF = pdblen
        ele = self.structure_fileAF.split(".")[0].split('_')
        # print ele, pdb
        self.pdb = "_".join(ele[:-1])
        self.chain=ele[-1]
        self.coverage = round(float(pdblen)/self.seqlen, 2)

    '''
    to this add transcript structure information, pdb file if any, length, chain id
    method, if any of of the exon has retained intron it will be called the same
    '''


class Exon():

    def __init__(self, id, aaseq, coding_span, nc_span):
        self.ID = id
        self.const = True if self.ID.split(".")[2] == 'G' else False
        self.seq = aaseq
        self.coding_span = coding_span
        self.nc_span = nc_span
        self.length = len(aaseq)
        self.ssseq = {}
        self.disorderseq={}
        # self.ssseq[None] = []
        self.strideseq = {}
        self.strideseqAF= {}
        # self.strideseq[None] = []
        self.IntronFlag = True if self.ID[0] == 'R' else False
        ele = id.split(".")
        if ele[0] == 'R':
            self.cons_flag = 'A'
        elif ele[2] == 'G':
            self.cons_flag = 'G'
        elif ele[2] == 'F':
            self.cons_flag = 'F'
        else:
            self.cons_flag = 'A'

    def out_secondseq(self, transID=None):
        '''
        will need, transcript ID not object
        '''
        if self.ssseq:
            var_lis = []
            key = {}
            for sseq in self.ssseq:
                if sseq is not None and sseq != "NULL":
                    for var in self.ssseq[sseq]:
                        var_lis += [var]
                        key[var] = sseq
            if transID is None:
                if len(var_lis) > 0:
                    var_lis.sort()
                    return key[var_lis[0]]
                else:
                    return False
            else:
                return key[transID] if transID in key else False
        else:
            return False
    def out_disorderseq(self, transID=None):
        '''
        will need, transcript ID not object
        '''
        if self.disorderseq:
            var_lis = []
            key = {}
            for diso in self.disorderseq:
                if diso is not None and diso != "NULL":
                    for var in self.disorderseq[diso]:
                        var_lis += [var]
                        key[var] = diso
            if transID is None:
                if len(var_lis) > 0:
                    var_lis.sort()
                    return key[var_lis[0]]
                else:
                    return False
            else:
                return key[transID] if transID in key else False
        else:
            return False


    def out_strideseq(self, transID=None):
        '''
        will need, transcript ID not object
        '''

        if self.strideseq:
            var_lis = []
            key = {}
            for strideseq in self.strideseq:
                if strideseq is not None and strideseq != "NULL":
                    for var in self.strideseq[strideseq]:
                        var_lis += [var]
                        key[var] = strideseq
            if transID is None:
                if len(var_lis) > 1:
                    var_lis.sort()
                    return key[var_lis[0]]
            else:
                return key[transID] if transID in key else False
        else:
            return False
    
    def out_strideseqAF(self, transID=None):
        '''
        will need, transcript ID not object
        '''

        if self.strideseqAF:
            var_lis = []
            key = {}
            for strideseqAF in self.strideseqAF:
                if strideseqAF is not None and strideseqAF != "NULL":
                    for var in self.strideseqAF[strideseqAF]:
                        var_lis += [var]
                        key[var] = strideseqAF
            if transID is None:
                if len(var_lis) > 1:
                    var_lis.sort()
                    return key[var_lis[0]]
            else:
                return key[transID] if transID in key else False
        else:
            return False

    def second_seq(self, sstring, var):
        if sstring not in self.ssseq:
            self.ssseq[sstring] = [var]
        else:
            self.ssseq[sstring] += [var]

    def stride_seq(self, sstring, var):
        if sstring not in self.strideseq:
            self.strideseq[sstring] = [var]
        else:
            self.strideseq[sstring] += [var]
    
   

    def stride_seqAF(self, sstring, var):
        if sstring not in self.strideseqAF:
            self.strideseqAF[sstring] = [var]
        else:
            self.strideseqAF[sstring] += [var]
    def disord_seq(self, distring, var):
        if distring not in self.disorderseq:
            self.disorderseq[distring] = [var]
        else:
            self.disorderseq[distring] += [var]

    def description(self):
        print "+++++++Exon class"
        print self.ID, self.coding_span, self.nc_span, self.length, self.seq, len(
            self.ssseq), self.ssseq.keys(), len(self.strideseq), self.strideseq.keys(), self.const, self.IntronFlag

    def Exon_characteristics(self, WEF):
        self.WEF = WEF
    '''
    add coding type, length of sequence
    strie seqyeuinces and disorder sequence addition are notthe subject of this revision
    '''

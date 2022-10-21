import os
import cPickle as pickle
import MySQLdb
import re
organismdict = {
    '9606': "Homo sapiens", '7227': "Drosophila melanogaster", '10090': "Mus musculus", '7955': "Danio rerio", '6239': "Caenorhabditis elegans"}
alpha = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
         'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
colorcode = {
    0: "aqua",
    1: "lawngreen",
    2: "black",
    3: "blueviolet",
    4: "brown",
    5: "hotpink",
    6: "cadetblue",
    7: "crimson",
    8: "darkcyan",
    9: "darkgreen",
    10: "darkkhaki",
    11: "darkolivegreen",
    12: "gold",
    13: "slategrey",
    14: "burlywood",
    15: "indianred",
    16: "bisque",
    23: "lightgreen",
    17: "mediumorchid",
    18: "mediumseagreen",
    19: "slateblue",
    20: "goldenrod",
    21: "springgreen",
    22: "tan",
    24: "teal",
    25: "salmon",
    26: "sandybrown",
    27: "turquoise",
    28: "yellowgreen",
    29: "tomato",
    30: "thistle",
    31: "turquoise",
    32: "wheat",
    33: "blanchedalmond"
}

colcode = []
for i in range(0, 25):
    colcode += [(alpha[i], colorcode[i])]


# In[ ]:

'''
have a look at changes made again to calss exon genes
run the programs with new configuration saved

rename previous database to nextrapv1_1

rename the updated object names to iterate
'''

# In[2]:


def geneObcreator(genob, genepk):
    name = genob.detail
    name = re.sub(r'\'', '\'\'', name)
    entrezid = genob.ID
    organism = organismdict[genob.txid]
    txid = genob.txid
    sql_iter = "INSERT INTO exonapp_gene (id, name, entrezid, organism, txid)    VALUES ('%d', '%s', '%d', '%s', '%s')" % (
        genepk, name, entrezid, organism, str(txid))
    return sql_iter


def exonObcreator(genefk, exon, exonpk, trans_fkhash, exons_domainseq, exonssspk, exonsdispk, exonsdompk):
    sqliterations = []
    wef = exon.WEF
    exId = exon.ID
    length = exon.length
    codst = exon.coding_span[0] if exon.coding_span else 0
    codend = exon.coding_span[1] if exon.coding_span else 0
    sstype = len(exon.ssseq) if exon.ssseq else 0
    rawst = exon.nc_span[0] if exon.nc_span else 0
    rawend = exon.nc_span[1] if exon.nc_span else 0
    aaseq = exon.seq
    parent = exon.parent.ID if exon.parent else "False"

    sqliterations += ["INSERT INTO exonapp_exongenes (parent,id, wef, exId, aaseq, length, codst, codend, rawst, rawend, sstype, gene_id) \
                    VALUES ('%s','%d', '%s', '%s', '%s', '%d', '%d', '%d', '%d', '%d', '%s', '%d')" %
                      (parent, exonpk, str(wef), exId, aaseq, length, codst, codend, rawst,
                       rawend, sstype, genefk)]
    # print "queryexonappexonsgenes","INSERT INTO exonapp_exongenes (id, wef, exId, aaseq, length, codst, codend, rawst, rawend, sstype, gene_id) VALUES ('%d', '%s',                         '%s','%s',%d', '%d','%d', '%d','%d', '%s', '%d')" %                    (exonpk, str(wef), exId,aaseq, length, codst, codend, rawst,
    #                 rawend, sstype, genefk)
    # cursor.execute(sqliterations[-1])
    # print "done with care"
    if aaseq and exon.ssseq:
        for sseq in exon.ssseq:
            if sseq is not None and sseq != "NULL":
                translist = ",".join(exon.ssseq[sseq])
                for tri in trans_fkhash:
                    translist = re.sub(tri, str(trans_fkhash[tri]), translist)
                query = "INSERT INTO exonapp_exonssprediction \
                (id, list_trans_fk,exon_id, ssseq) VALUES ('%d', '%s', '%d','%s')" % \
                    (exonssspk, translist, exonpk, re.sub("C", "-", sseq))
                # print "query",query
                sqliterations += [query]
                exonssspk += 1
    if aaseq and exon.disorderseq:
        for seq in exon.disorderseq:
            if seq is not None and seq != "NULL":
                translist = ",".join(exon.disorderseq[seq])
                for tri in trans_fkhash:
                    translist = re.sub(tri, str(trans_fkhash[tri]), translist)
                query = "INSERT INTO exonapp_exonsdisorder (id, list_trans_fk,exon_id, disseq) VALUES ('%d', '%s',                     '%d','%s')" % (
                    exonsdispk, translist, exonpk, re.sub("D", "-", seq))
                sqliterations += [query]
                exonsdispk += 1

    if aaseq and exon.ID in exons_domainseq:
        for seq in exons_domainseq[exon.ID]:
            query = "INSERT INTO exonapp_exonsdomseq  (id, list_trans_fk,exon_id, domseq) VALUES ('%d', '%s',                 '%d','%s')" % (
                exonsdompk, exons_domainseq[exon.ID][seq], exonpk, seq)
            # print "domainquery", query
            sqliterations += [query]
            exonsdompk += 1
    return exonssspk, exonsdispk, exonsdompk, sqliterations


def transcriptObcreator(genefk, transpk, trans, totc, uc, exons_fk_has):
    par1_dic = {'U': 0, 'T': 0, 'M': 0, 'R': 0, 'D': 0}
    par2_dic = {'A': 0, "G": 0, "F": 0}

    pi = trans.PI
    tId = trans.ID
    length = trans.seqlen
    exonscount = len(trans.exons)
    swissprot = trans.swiss_pid
    exonsIds = ",".join([i.ID for i in trans.exons])
    for efk in exons_fk_has:
        exonsIds = re.sub(efk, str(exons_fk_has[efk]), exonsIds)
    total_domC = totc
    unique_domC = uc
    disordered_region = ''
    ssp_region = ''
    factor_length = 0

    for i in trans.exons:
        if i.seq:
            factor_length += len(i.seq)
            disp = i.out_disorderseq(trans.ID)
            ssp = i.out_secondseq(trans.ID)
            if disp and disp is not None:
                disordered_region += disp
            if ssp and ssp is not None:
                ssp_region += ssp

        if i.ID[0] != 'R':
            par1_dic[i.ID[0]] += 1
            par2_dic[i.ID.split(".")[2]] += 1
        else:
            par1_dic['R'] += 1
            par2_dic['A'] += 1
    order1 = ["U", "D", "T", "M", "R"]
    order2 = ["G", "A", "F"]
    str1 = ''
    str2 = ''
    for i in order1:
        str1 += "%s:%s, " % (i, par1_dic[i])
    for i in order2:
        str2 += "%s:%s, " % (i, par2_dic[i])
    str1 = str1[:-2]
    str2 = str2[:-2]
    if ssp_region:
        structred_regionss = round(float(ssp_region.count(
            "H")+ssp_region.count("E"))/factor_length*100, 2)
    else:
        structred_regionss = "NULL"
    if disordered_region:
        structred_regiondisp = round(
            float(disordered_region.count("S"))/factor_length*100, 2)
    else:
        structred_regiondisp = "NULL"
    exonscountUTMRD = str1
    exonscountAG = str2

    sql_iter = "INSERT INTO exonapp_transcripts \
    (id, tId, swissprot, length, pi,exonscount,exonsIds,\
    unique_domC,total_domC,exonscountUTMRD,exonscountAG,\
     geneRef_id,structured_count_ssp,structured_count_disp) \
     VALUES ('%d', '%s','%s','%d', '%d','%d', '%s','%d','%d', \
     '%s','%s','%d','%s','%s')" % (transpk, tId, swissprot, length,
                                   pi, exonscount, exonsIds, unique_domC, total_domC,
                                   exonscountUTMRD, exonscountAG, genefk, str(structred_regionss), str(structred_regiondisp))
    # print sql_iter
    return sql_iter


def domseqiterator(trans, transfk, transdoms, genecodes, domainhash):
    # transcript, trans_pk_temp[transcript.ID], dom_trans_spans,has_dom, exons_domainseq)
    # print("genecodes", genecodes)
    strdom = "-"*trans.seqlen
    # print trans.seqlen, trans.ID
    # print "".join([i.seq for i in trans.exons])
    # print strdom
    #spanset = []
    if transdoms:
        for i in transdoms:
            domcode = genecodes[transdoms[i]][0]
            # print("DomCODE", domcode)
            #spanset += [(set(range(i[0]-1, i[1])), domcode)]
            # 0 here is the code and 1 is the color
            strdom = strdom[:i[0]-1] + \
                (i[1] - i[0]+1)*domcode + strdom[i[1]:]
            # print i
            # print spanset
    # print strdom
    lengt = 0
    returnsql = []
    for exons in trans.exons:
        exspan = set(range(lengt, lengt+exons.length))
        dom_seq = strdom[lengt:lengt+exons.length]
        lengt += exons.length
        # print exons.ID
        # print dom_seq
        # print exons.seq

        exIdt = exons.ID
        if exIdt not in domainhash:
            domainhash[exIdt] = {}
            domainhash[exIdt][dom_seq] = str(transfk)
        else:
            if dom_seq not in domainhash[exIdt]:
                domainhash[exIdt][dom_seq] = str(transfk)
            else:
                domainhash[exIdt][dom_seq] += ",%s" % str(transfk)
    return domainhash


def domindObcreator(genepk, dompk, i, hasi):
    name = i[0]
    dId = i[1]
    color = hasi[1]
    code = hasi[0]
    sql_iter = "INSERT INTO exonapp_domaininfgene  (id, color, code, name, dId, gene_id)                 VALUES ('%d', '%s', '%s', '%s', '%s','%s')" % (
        dompk, color, code, name, dId, genepk)
    return sql_iter


# In[4]:


db = MySQLdb.connect("localhost", "root", "cbg_2018", "nextrapdb")
cursor = db.cursor()
cursor.execute("SET FOREIGN_KEY_CHECKS=0;")

cursor.execute("DROP TABLE IF EXISTS exonapp_gene")
sql_table_header_gene = """CREATE TABLE exonapp_gene(
         id INT NOT NULL,
         name VARCHAR(200) NOT NULL,
         entrezid INT NOT NULL,
         organism VARCHAR(50),  
         txid VARCHAR(50),
         PRIMARY KEY(id));
         """
cursor.execute(sql_table_header_gene)

cursor.execute("DROP TABLE IF EXISTS exonapp_transcripts")
sql_table_header_trans = """CREATE TABLE exonapp_transcripts(
         id INT NOT NULL,
         tId VARCHAR(20) NOT NULL,
         swissprot VARCHAR(200),
         length INT,
         pi TINYINT(1),
         exonscount INT,  
         exonsIds TEXT,
         unique_domC INT,
         total_domC INT,
         structured_count_ssp VARCHAR(7),
         structured_count_disp VARCHAR(7),
         exonscountUTMRD VARCHAR(100),
         exonscountAG VARCHAR(50),
         geneRef_id INT,
         PRIMARY KEY(id),
         FOREIGN KEY(geneRef_id) REFERENCES exonapp_gene(id));
         """
cursor.execute(sql_table_header_trans)

cursor.execute("DROP TABLE IF EXISTS exonapp_exongenes")
sql_table_header_exons = """CREATE TABLE exonapp_exongenes(
         id INT NOT NULL,
         parent VARCHAR(50),
         wef FLOAT NOT NULL,
         exId VARCHAR(50),
         length INT,
         codst INT,
         aaseq TEXT,
         codend INT,  
         rawst INT,
         rawend INT,
         sstype INT,
         gene_id INT,
         PRIMARY KEY(id),
         FOREIGN KEY(gene_id) REFERENCES exonapp_gene(id));
         """
cursor.execute(sql_table_header_exons)

cursor.execute("DROP TABLE IF EXISTS exonapp_domaininfgene")
sql_table_header_dom = """CREATE TABLE exonapp_domaininfgene(
         id INT NOT NULL,
         color VARCHAR(30) NOT NULL,
         code VARCHAR(1),
         name VARCHAR(100),
         dId VARCHAR(100),
         gene_id INT,  
         PRIMARY KEY(id),
         FOREIGN KEY(gene_id) REFERENCES exonapp_gene(id));
         """
cursor.execute(sql_table_header_dom)
cursor.execute("DROP TABLE IF EXISTS exonapp_exonssprediction")
sql_table_header_exss = """CREATE TABLE exonapp_exonssprediction(
         id INT NOT NULL,
         exon_id INT NOT NULL,
         list_trans_fk TEXT NOT NULL,
         ssseq TEXT,
         PRIMARY KEY(id),
         FOREIGN KEY(exon_id) REFERENCES exonapp_exongenes(id));
         """
cursor.execute(sql_table_header_exss)
cursor.execute("DROP TABLE IF EXISTS exonapp_exonsdomseq")
sql_table_header_exsdom = """CREATE TABLE exonapp_exonsdomseq(
         id INT NOT NULL,
         exon_id INT NOT NULL,
         list_trans_fk TEXT NOT NULL,
         domseq TEXT,
         PRIMARY KEY(id),
         FOREIGN KEY(exon_id) REFERENCES exonapp_exongenes(id));
         """
cursor.execute(sql_table_header_exsdom)
cursor.execute("DROP TABLE IF EXISTS exonapp_exonsdisorder")
sql_table_header_exdis = """CREATE TABLE exonapp_exonsdisorder(
         id INT NOT NULL,
         exon_id INT NOT NULL,
         list_trans_fk TEXT NOT NULL,
         disseq TEXT,
         PRIMARY KEY(id),
         FOREIGN KEY(exon_id) REFERENCES exonapp_exongenes(id));
         """
cursor.execute(sql_table_header_exdis)
cursor.execute("SET FOREIGN_KEY_CHECKS=1;")
db.commit()
# db.close()
genePK = 1
transPK = 1
exongenePK = 1
exontransPK = 1
domainPK = 1
exonssspk = 1
exonsdispk = 1
exonsdompk = 1
object_list = ["/home/paras/project/protein_splicing/6239/derived_data/results/objectsave_diso_6239.pick",
               "/home/paras/project/protein_splicing/9606_0/derived_data/results/objectsave_diso_9606.pick",
               "/home/paras/project/protein_splicing/7227/derived_data/results/objectsave_diso_7227.pick",
               "/home/paras/project/protein_splicing/7955/derived_data/results/objectsave_diso_7955.pick",
               "/home/paras/project/protein_splicing/10090/derived_data/results/objectsave_diso_10090.pick"
               ]


# with open("/home/paras/project/protein_splicing/9606_0/derived_data/results/objectsave_9606.pick") as fin:
#     minidata=pickle.load(fin)
#

# In[5]:


def gene_ob_populator(data, genePK, transPK, exongenePK, domainPK, exonssspk, exonsdispk, exonsdompk):
    gene_c = 0
    for gene in data:
        # if gene == 103:
        exons_pk_temp = {}
        trans_pk_temp = {}
        for i in data[gene].exons:
            exons_pk_temp[i.ID] = exongenePK
            exongenePK += 1
        for i in data[gene].transcripts:
            trans_pk_temp[i.ID] = transPK
            transPK += 1
        queriestoexecute = []
        iiter = 0
        codcol = colcode
        #print("\n\ncodcol", codcol)
        has_dom = {}
        queriestoexecute += [geneObcreator(data[gene], genePK)]
        # return sql queries

        exons_domainseq = {}
        for transcript in data[gene].transcripts:
            dom_trans = {}
            dom_trans_spans = {}
            if transcript.pfam_list:
                for domiter in transcript.pfam_list:
                    keypair = (domiter[1], domiter[3])
                    if keypair not in dom_trans:
                        dom_trans[keypair] = 1
                        # dom_trans is having pfam id,name as tuple and counting them
                    else:
                        dom_trans[keypair] += 1
                    dom_trans_spans[domiter[0]] = keypair
                    # this ill have span as key and pair of id and name as values
                    # this value will be used to retrieve the correspodnding color and dom
                    # code for the span to be populated
                    if keypair not in has_dom:
                        has_dom[keypair] = [
                            codcol[iiter][0], codcol[iiter][1]]
                        iiter += 1

            queriestoexecute += [transcriptObcreator(genePK, trans_pk_temp[transcript.ID], transcript, sum(
                dom_trans.values()), len(dom_trans.keys()), exons_pk_temp)]
            exons_domainseq = domseqiterator(
                transcript, trans_pk_temp[transcript.ID], dom_trans_spans, has_dom, exons_domainseq)
            '''
                need gene object updater here, domainexon hash will have domain sequence 
                and then the components, 
                disorder types and transcripts emclosed
                '''
        # print exons_domainseq
        for exons in data[gene].exons:
            exonssspk, exonsdispk, exonsdompk, queriessql = exonObcreator(
                genePK, exons, exons_pk_temp[exons.ID], trans_pk_temp, exons_domainseq,                          exonssspk, exonsdispk, exonsdompk)
            queriestoexecute += queriessql

        for smalldom in has_dom:
            queriestoexecute += [domindObcreator(genePK,
                                                 domainPK, smalldom, has_dom[smalldom])]
            domainPK += 1
        gene_c += 1
        # print queriestoexecute
        print(gene, gene_c, len(data), data[gene].txid)
        for que in queriestoexecute:
            # print que
            try:
                cursor.execute(que)
            except Exception as E:
                print E, que
                db.rollback()
                db.close()
        genePK += 1
        # if gene == 100:
        #   break
    return genePK, transPK, exongenePK, domainPK, exonssspk, exonsdispk, exonsdompk


# In[6]:


for objs in object_list:
    with open(objs, "rb") as fin:
        dataobj = pickle.load(fin)
    genePK, transPK, exongenePK, domainPK, exonssspk, exonsdispk, exonsdompk = gene_ob_populator(
        dataobj, genePK, transPK, exongenePK, domainPK, exonssspk, exonsdispk, exonsdompk)


# In[7]:


db.commit()
db.close()


# db = MySQLdb.connect("localhost","root","cbg_2018","Alt_spl_2018" )
# cursor = db.cursor()
# cursor.execute("DROP TABLE IF EXISTS NCBI_Feb_Apr_2018_FK_PID_ExonComb_PI_SeqLen_SP_ENS_PDB_StrucLen")
# sql_table_header = """CREATE TABLE NCBI_Feb_Apr_2018_FK_PID_ExonComb_PI_SeqLen_SP_ENS_PDB_StrucLen(
#          FK INT NOT NULL,
#          PID  VARCHAR(20),
#          ExonComb TEXT,
#          PI CHAR(1),
#          SeqLen INT,
#          SP TEXT,
#          ENS TEXT,
#          PDB VARCHAR(40),
#          StrucLen INT)"""
# cursor.execute(sql_table_header)
#
# for gene in temp:
#     for var in has_gene_2[str(gene)]:
#         fk=for_key_var_no[var]
#         coods=[]
#         with open (variants_exon_type_add+"/%s"%var) as fin:
#             temp_dat=pickle.load(fin)
#             for i in temp_dat:
#                 if i[3]=="C":
#                     coods+=[i[2]]
#                 else:
#                     coods+=[i[1]]
#             #coods=[i[2] for i in pickle.load(fin) if i[3]=="C"]
#         exoncomb=",".join([storing_var_exon_combination[gene][var][i] for i in coods])
#         #exoncomb=",".join(map(str,storing_var_exon_combination[var]))
#         if var == pi[str(gene)]:
#             pist="Y"
#         else:
#             pist="N"
#         seqlen=prot_length[var]
#         sp=swiss_h[var]
#         ens=ens_h[var]
#         pdbs1=pdb[var][0]
#         strlen=pdb[var][1]
#         if strlen=="NULL":
#             strlen=0
#         try:
#             #print fk,var,exoncomb,pist,seqlen,sp,ens,pdbs1,strlen
#            #
#             sql_iter = "INSERT INTO NCBI_Feb_Apr_2018_FK_PID_ExonComb_PI_SeqLen_SP_ENS_PDB_StrucLen \
#                     (FK, PID, ExonComb, PI, SeqLen, SP, ENS, PDB, StrucLen) \
#                     VALUES ('%d', '%s', '%s', '%s', '%d', '%s' , '%s', '%s', '%d')" % \
#                     (fk,var,exoncomb,pist,seqlen,sp,ens,pdbs1,strlen)
#             cursor.execute(sql_iter)
#         except Exception as EE:
#             print gene,var,EE
#             break
#             db.rollback()
# db.commit()
# db.close()

# In[ ]:

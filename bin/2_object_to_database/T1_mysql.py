import os
import sys
import cPickle as pickle
import MySQLdb
import re
from progress.bar import Bar
import db_modules

if len(sys.argv)!=3:
    print ("Please provide 1. dir which has the objects to be read, 2. database config")
    sys.exit()

from db_modules import colcode 

organismdict = {
    '9606': "Homo sapiens", '7227': "Drosophila melanogaster", '10090': "Mus musculus", '7955': "Danio rerio", '6239': "Caenorhabditis elegans"}

def domIndentificationSql(genepk, dompk, i, hasi):
    name = i[0]
    dId = i[1]
    color = hasi[1]
    code = hasi[0]
    sql_iter = "INSERT INTO exonapp_domaininfgene  (id, color, code, name, dId, gene_id) VALUES ('%d', '%s', '%s', '%s', '%s','%s')" % (
        dompk, color, code, name, dId, genepk)
    return sql_iter

#################################################### 1st block  #######################################################
db = MySQLdb.connect(host="localhost", user="root", passwd="cbg_2022", db="enactdb", unix_socket="/var/run/mysqld/mysqld.sock")

#db = _mysql.connect(host="localhost", user="root", passwd="cbg_2022", db="enactdb", unix_socket="/var/run/mysqld/mysqld.sock")
# if cant identify correct socket, 
# in mysql window type following "show variables like 'socket';"
cursor = db.cursor()
cursor.execute("SET FOREIGN_KEY_CHECKS=0;")
#######################################################################################################################
list_to_execute = db_modules.default_reset()
# create tables and drops them
for cmd in list_to_execute:
    cursor.execute(cmd)

cursor.execute("SET FOREIGN_KEY_CHECKS=1;")
db.commit()

genePK = 1
transPK = 1
exongenePK = 1
exontransPK = 1
domainPK = 1
exonssspk = 1
exonsdispk = 1
exonsdompk = 1
# object_list = ["/home/paras/project/protein_splicing/6239/derived_data/results/objectsave_diso_6239.pick",
#                "/home/paras/project/protein_splicing/9606_0/derived_data/results/objectsave_diso_9606.pick",
#                "/home/paras/project/protein_splicing/7227/derived_data/results/objectsave_diso_7227.pick",
#                "/home/paras/project/protein_splicing/7955/derived_data/results/objectsave_diso_7955.pick",
#                "/home/paras/project/protein_splicing/10090/derived_data/results/objectsave_diso_10090.pick"
#                ]

def gene_ob_populator(data, txid, organism,  genePK, transPK, exongenePK, domainPK, exonssspk, exonsdispk, exonsdompk):
    '''
    iterate the genes,
        1. assign prim key to exons and transcripts and update the handles asap
        "INSERT INTO exonapp_gene (id, name, entrezid, organism, txid)"
        2. generate the domainseq for exons has, and 
        ITERATE THE TRANSCRIPTS
            generate two hashes, one having spans as keys and pfamId, name as value, and other the opposite, (id, name) as key, will be used to give colors
            call to transcript ob sql's, will be used to generate 3 new vars, exons region, having information of their shapes,
            call to domseq WONT add SQL queries but update transcript sequences to domain alphabets and store than in domseqHAS with following format
                domainhash[exIdt][dom_seq] += ",%s" % str(transfk)
        3. domseqHAS will continue to update
        ITERATE the EXONS
    '''
    gene_c = 0
    bar = Bar('Processing genes:', max=len(data))
    for gene in data:
        ###################### !1 ###########################
    #   if gene == 108115:
        temptrans = [[i.PI, i] for i in  data[gene].transcripts]
        temptrans.sort(reverse=True)
        piFirstTransList= [i[1] for i in temptrans]
        # just keeping PI in the first position, 

        exons_pk_temp = {}
        trans_pk_temp = {}
        for i in data[gene].exons:
            exons_pk_temp[i.ID] = exongenePK
            exongenePK += 1
        for i in piFirstTransList:
            trans_pk_temp[i.ID] = transPK
            transPK += 1
        ###################### !1/> ###########################
        queriestoexecute = []
        queriestoexecute += [db_modules.geneObcreator(data[gene], genePK, txid, organism)]
        ###################### * 2 ############################
        handleColorUpdator = 0
        codcol = colcode
        has_color_storage = {}
        exons_domainseq = {}
        for transcript in piFirstTransList:
            # print (transcript.ID)
            dom_trans = {}
            dom_trans_spans = {}
            # print (transcript.pfam_list)
            if transcript.pfam_list:
                for domiter in transcript.pfam_list:
                    keypair = (domiter[1], domiter[3])
                    if keypair not in dom_trans:
                        dom_trans[keypair] = 0
                        # dom_trans is having pfam id,name as tuple and counting them
                    dom_trans[keypair] += 1
                    dom_trans_spans[domiter[0]] = keypair
                    # this Will have span as key and pair of id and name as values
                    # this value will be used to retrieve the correspodnding color and dom
                    # code for the span to be populated
                    if keypair not in has_color_storage:
                        has_color_storage[keypair] = (
                            codcol[handleColorUpdator][0], codcol[handleColorUpdator][1])
                        handleColorUpdator += 1
                        # colcode is list of list, spanning 0 to 25, having alphabet as 1st elem, color as 2nd
                        # has_color_storage will assign colors on the basis of the uniqueCombination of pfamName and ID
                        # means that same domain if repeated wont be given new color
                        # -> alphabets can be limitations after some time, do gave it a thought

            queriestoexecute += [db_modules.transcriptObcreator(gene, genePK, trans_pk_temp[transcript.ID], transcript, sum(
                dom_trans.values()), len(dom_trans.keys()), exons_pk_temp, dom_trans_spans, has_color_storage)]

            exons_domainseq = db_modules.domseqiterator(gene, 
                transcript, trans_pk_temp[transcript.ID], dom_trans_spans, has_color_storage, exons_domainseq)
            # transcriptOb, transPrimID, spans as Key:(pfamidName) as value, pfamidName as key: (alphbet and colname) as value, 
            # exons_domainseq will be updated

            '''
                need gene object updater here, domainexon hash will have domain sequence 
                and then the components, 
                disorder types and transcripts emclosed
            '''
        # print (dom_trans_spans)
        # print exons_domainseq
        ##############################*2/###########################################
        ##############################!3 ###########################################
        # for i in data[gene].exons:
        #     print (i.ID, i.seq)

        # for transcript in piFirstTransList:
        #     print (transcript.ID)
        #     for i in transcript.exons:
        #         print (i.ID, i.seq,'\n')

        for exons in data[gene].exons:
            exonssspk, exonsdispk, exonsdompk, queriessql = db_modules.exonObcreator(
                genePK, exons, exons_pk_temp[exons.ID], trans_pk_temp, exons_domainseq, exonssspk, exonsdispk, exonsdompk)
            queriestoexecute += queriessql
        
        for smalldom in has_color_storage:
            queriestoexecute += [domIndentificationSql(genePK,
                                                 domainPK, smalldom, has_color_storage[smalldom])]
            domainPK += 1
        gene_c += 1
        # print queriestoexecute
        # print(gene, gene_c, len(data), data[gene].txid)
        # print (queriestoexecute)
        for que in queriestoexecute:
            # print que
            try:
                cursor.execute(que)
            except Exception as E:
                print E, que, gene
                db.rollback()
                db.close()
        genePK += 1
        # if gene == 100:
        #   break
        bar.next()
    return genePK, transPK, exongenePK, domainPK, exonssspk, exonsdispk, exonsdompk

prog, inpfolder, dbparam = sys.argv
object_list = os.listdir(inpfolder)
print (object_list)
for objs in object_list:
    txid=re.match(r'objectsave_\d+',objs)
    txid=txid.group().split("_")[1]
    #print (txid)
    organism = organismdict[txid]
    #sys.exit()
    with open(os.path.join(inpfolder,objs), "rb") as fin:
        dataobj = pickle.load(fin)
    genePK, transPK, exongenePK, domainPK, exonssspk, exonsdispk, exonsdompk = gene_ob_populator(
        dataobj, txid, organism, genePK, transPK, exongenePK, domainPK, exonssspk, exonsdispk, exonsdompk)

db.commit()
db.close()

print ("DOne SuccessFully")


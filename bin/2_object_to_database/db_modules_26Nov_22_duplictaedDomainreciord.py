import re,sys

def default_reset():
    cmdtoexecute = []
    sql_table_header_gene = """CREATE TABLE exonapp_gene(
            id INT NOT NULL,
            name VARCHAR(200) NOT NULL,
            entrezid INT NOT NULL,
            organism VARCHAR(50),  
            txid VARCHAR(50),
            PRIMARY KEY(id));
            """
    # ensemblID will be varchar something
    
    sql_table_header_trans = """CREATE TABLE exonapp_transcripts(
            id INT NOT NULL,
            tId VARCHAR(20) NOT NULL,
            swissprot VARCHAR(400),
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
            domains TEXT, 
            exonsRegion TEXT,
            secondaryStructure TEXT,
            FOREIGN KEY(geneRef_id) REFERENCES exonapp_gene(id));
            """
    # -> will need to add new columns to this, 

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
    sql_table_header_exss = """CREATE TABLE exonapp_exonssprediction(
            id INT NOT NULL,
            exon_id INT NOT NULL,
            list_trans_fk TEXT NOT NULL,
            ssseq TEXT,
            PRIMARY KEY(id),
            FOREIGN KEY(exon_id) REFERENCES exonapp_exongenes(id));
            """
    sql_table_header_exsdom = """CREATE TABLE exonapp_exonsdomseq(
            id INT NOT NULL,
            exon_id INT NOT NULL,
            list_trans_fk TEXT NOT NULL,
            domseq TEXT,
            PRIMARY KEY(id),
            FOREIGN KEY(exon_id) REFERENCES exonapp_exongenes(id));
            """
    
    sql_table_header_exdis = """CREATE TABLE exonapp_exonsdisorder(
            id INT NOT NULL,
            exon_id INT NOT NULL,
            list_trans_fk TEXT NOT NULL,
            disseq TEXT,
            PRIMARY KEY(id),
            FOREIGN KEY(exon_id) REFERENCES exonapp_exongenes(id));
            """
    
    drop_geneTab = "DROP TABLE IF EXISTS exonapp_gene"
    drop_transTab = "DROP TABLE IF EXISTS exonapp_transcripts"
    drop_exonsTab = "DROP TABLE IF EXISTS exonapp_exongenes"
    drop_headerdomTab = "DROP TABLE IF EXISTS exonapp_domaininfgene"
    drop_exssTab = "DROP TABLE IF EXISTS exonapp_exonssprediction"
    drop_exonsdomTab = "DROP TABLE IF EXISTS exonapp_exonsdomseq"
    drop_exonsdisTab = "DROP TABLE IF EXISTS exonapp_exonsdisorder"
    
    cmdtoexecute += [drop_geneTab]
    cmdtoexecute += [sql_table_header_gene]
    cmdtoexecute += [drop_transTab]
    cmdtoexecute += [sql_table_header_trans]
    cmdtoexecute += [drop_exonsTab]
    cmdtoexecute += [sql_table_header_exons]
    cmdtoexecute += [drop_headerdomTab]
    cmdtoexecute += [sql_table_header_dom]
    cmdtoexecute += [drop_exssTab]
    cmdtoexecute += [sql_table_header_exss]
    cmdtoexecute += [drop_exonsdomTab]
    cmdtoexecute += [sql_table_header_exsdom]
    cmdtoexecute += [drop_exonsdisTab]
    cmdtoexecute += [sql_table_header_exdis]
    return cmdtoexecute

colorcode = {
    0: "#4d86a5",
    1: "#cf0bf1",
    2: "#12e2f1",
    3: "#3e517a",
    4: "#98da1f",
    5: "#fc9f5b",
    6: "#d60b2d",
    7: "#c3c4e9",
    8: "#9cc76d",
    9: "#2dffdf",
    10: "#d00000",
    11: "#ffba08",
    12: "#cbff8c",
    13: "#8fe388",
    14: "#1b998b",
    15: "#3185fc",
    16: "#5d2e8c",
    23: "#46237a",
    17: "#ff7b9c",
    18: "#ff9b85",
    19: "#669900",
    20: "#99cc33",
    21: "#ccee66",
    22: "#006699",
    24: "#3399cc",
    25: "#990066",
    26: "#cc3399",
    27: "#ff6600",
    28: "#ff9900",
    29: "#ffcc00",
    30: "#033270",
    31: "#1368aa",
    32: "#4091c9",
    33: "#9dcee2",
    34: "#fedfd4",
    35: "#f29479",
    36: "#f26a4f",
    37: "#ef3c2d",
    38: "#cb1b16",
    39: "#65010c"
}
alpha = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
         'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
# -> this should change to websafe colors

colcode = []
for i in range(0, 25):
    colcode += [(alpha[i], colorcode[i])]

def geneObcreator(genob, genepk, txid, organism):
    name = genob.detail
    name = re.sub(r'\'', '\'\'', name)
    entrezid = genob.ID
    sql_iter = "INSERT INTO exonapp_gene (id, name, entrezid, organism, txid)    VALUES ('%d', '%s', '%d', '%s', '%s')" % (
        genepk, name, entrezid, organism, str(txid))
    return sql_iter


def transcriptObcreator(genencbid, genefk, transpk, trans, totc, uc, exons_fk_has, dom_trans_spans, col_cood_hash):
    # transcriptObcreator(genePK, trans_pk_temp[transcript.ID], transcript, sum(dom_trans.values()), len(dom_trans.keys()), exons_pk_temp)
    # dom_trans_spans[domiter[0]] = keypair
                    # this Will have span as key and pair of id and name as values


    par1_dic = {'U': 0, 'T': 0, 'M': 0, 'R': 0, 'D': 0}
    par2_dic = {'A': 0, 'G': 0, 'F': 0}

    pi = trans.PI
    tId = trans.ID
    length = trans.seqlen
    exonscount = len(trans.exons)
    swissprot = trans.swiss_pid
    exonsIds = ",".join([i.ID for i in trans.exons])
    
    # making retention cases go first because otheroiwse normal ids, may make retention case ID's corrup and ha;f processed, 
    retentionHas={i:exons_fk_has[i] for i in exons_fk_has if i[0]=='R'}
    normalHas={i:exons_fk_has[i] for i in exons_fk_has if i[0]!='R'}
    for efk in retentionHas:
        exonsIds = re.sub(efk, str(retentionHas[efk]), exonsIds)
    for efk in normalHas:
        exonsIds = re.sub(efk, str(normalHas[efk]), exonsIds)
    
    total_domC = totc
    unique_domC = uc
    disordered_region = ''
    ssp_region = ''
    strExon = ''
    factor_length = 1
    # length of exons, it should start from 1

    for i in trans.exons:
        modExo = i.ID + '$'
        smoothDis=0
        if i.seq:
            disp = i.out_disorderseq(trans.ID)
            ssp = i.out_secondseq(trans.ID)
            #dom = i.out_domseq(trans.ID)
            if disp and disp is not None:
                disordered_region += disp
            if ssp and ssp is not None:
                ssp_region += ssp
            # -> will None be ever true, ??? in middle of sequqnce, I think No, but ssp can if length >3000 is considered

            ## ading ss and dis spans:
            
            cond1 = abs(i.coding_span[0]-i.nc_span[0])>3
            cond2 = abs(i.coding_span[1]-i.nc_span[1])>3
            if cond1 and cond2:
                smoothDis=3
            elif cond1:
                smoothDis=1
            elif cond2:
                smoothDis=2
            # else:
            #     print (cond1)
            #     print (cond2)
            #     print (i.coding_span)
            #     print (i.nc_span)
            #  # modExon will be "T.1.A.3.0.0:1,10,[0-4]" 1 to 10 its aa span and 0-4 are shapes of the boundary, 0, round both sides (default), 1 zag left (means bopundaries are longer than contributed span,), 2 zag right (3; UTR), 3 single exon case in tranrcipt and both UTRS in this. 4 absence (no aa, YET TO BE WORKED)

            modExo += '%s,%s,%s_'%(factor_length, factor_length+len(i.seq)-1,smoothDis)
            factor_length += len(i.seq)
            # if exon has length 10 in beginning , give exon coods, 1, 1+10-1(10), 
            # next from 10,20, will have, 10, 10+10-1, 10,19
        else:
            smoothDis = 4
            modExo+='%s,%s,%s,%s,%s_'%(factor_length,factor_length,smoothDis, i.nc_span[0], i.nc_span[1])
        
        if i.ID[0] != 'R':
            par1_dic[i.ID[0]] += 1
            par2_dic[i.ID.split(".")[2]] += 1
        else:
            par1_dic['R'] += 1
            par2_dic['A'] += 1

        # above two are related to exon shape rendering
        # modExon will be "T.1.A.3.0.0:1,10,[0-4]" 1 to 10 its aa span and 0-4 are shapes of the boundary, 0, round both sides (default), 1 zag left (means bopundaries are longer than contributed span,), 2 zag right (3; UTR), 3 single exon case in tranrcipt and both UTRS in this. 4 absence (no aa, YET TO BE WORKED)

        strExon += '%s'%(modExo)

    def spansGiver(prestring, needmatch, string):
        region = re.finditer(needmatch, string)
        for i in region:    
            sp=i.span()
            prestring += '%s,%s_'%(sp[0]+1,sp[1]), #from +1 deducted 1 in right coods, 
        #H:10,20_50,40_  otherwise H:
        return prestring

    strh = spansGiver('H:', 'H{1,}', ssp_region)
    stre = spansGiver('E:', 'E{1,}', ssp_region)
    strdis = spansGiver('D:', 'D{1,}', disordered_region)
    if strh=="H:" and stre=="E:":
        strh="H:0,0_"
        stre="E:"
        print ("no ss", genencbid)
    # if no ss

    secondaryStructure = strh[:-1]+'-'+stre[:-1]+'-'+strdis[:-1]
    #last element is skpped because of trailing undersocre used above
    exonsRegion = strExon[:-1]
    
    uniqueNameIdDomHash={}
    # dom_trans_spans[(1,100)]=(pfamId, pfam name)
    for spans in dom_trans_spans:
        value = dom_trans_spans[spans]
        if value not in uniqueNameIdDomHash:
            uniqueNameIdDomHash[value]=[]
        uniqueNameIdDomHash[value]+=[spans]
    
    # uniqueNameIdDomHash[(pfamId, pfamName)] = [list of spans]
    
    # sortedSpansList=[]
    listDomString=[]
    domDone_listDomString = set()
    listDomString1=[]
    domDone_listDomString1 = set()
    
    for idName in uniqueNameIdDomHash:
        pfam, domname = idName
        spansList= uniqueNameIdDomHash[idName]
        alphabet, color = col_cood_hash[idName]
        domstr = '%s,%s,%s:'%(domname,pfam, color)
        domstr1 = '%s,%s,%s:'%(domname,pfam, color)
        flag1 = 0
        for spans in spansList:
            tempRange = set(range(spans[0],spans[1]+1))
            if tempRange & domDone_listDomString:
                # means clash observed with spans in domstr, move to next data
                if not tempRange & domDone_listDomString1:
                    domstr1 += '%s,%s_'%(spans[0],spans[1])
                    flag1 = 1
                    domDone_listDomString1 |= tempRange
                else:
                    print ("clash Noted in alt track of 0omains also")
                    sys.exit()
            else:
                # default where domains are not clashing
                domstr+='%s,%s_'%(spans[0],spans[1])
                domDone_listDomString |= tempRange
            
            # -> reviist if they represent correct values, 
        listDomString += [domstr[:-1]]
        if flag1:
            listDomString1 += [domstr1[:-1]]

    if listDomString:
        domains='$'.join(listDomString)    
    else:
        print ("no domains",genencbid)
        listDomString=["NoDomainEntryPredicted,0,white:0,0"]
        domains='$'.join(listDomString)

    if listDomString1:
        domains1='$'.join(listDomString1)    
    else:
        listDomString1=["NoDomainEntryPredicted,0,white:0,0"]
        domains1='$'.join(listDomString1)

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
        structred_regionss = round((float(ssp_region.count("H")+ssp_region.count("E"))/trans.seqlen) *100, 2)
    else:
        structred_regionss = "NULL"
    if disordered_region:
        disordered_region = round((float(disordered_region.count("D"))/trans.seqlen)*100, 2)
    else:
        disordered_region = "NULL"
    exonscountUTMRD = str1
    exonscountAG = str2

    sql_iter = "INSERT INTO exonapp_transcripts \
    (id, tId, swissprot, length, pi, exonscount, exonsIds,\
    unique_domC, total_domC, exonscountUTMRD, exonscountAG,\
    geneRef_id, domains, domains1, exonsRegion, secondaryStructure, structured_count_ssp,structured_count_disp) \
     VALUES ('%d', '%s', '%s', '%d', '%d', '%d', '%s',\
        '%d', '%d', '%s','%s',\
        '%d', '%s', '%s', '%s', '%s', '%s')"\
             % (transpk, tId, swissprot, length, pi, exonscount, exonsIds, 
                unique_domC, total_domC, exonscountUTMRD, exonscountAG, 
                genefk, domains, domains1, exonsRegion, secondaryStructure, str(structred_regionss), str(disordered_region))
    # print sql_iter
    # -> need to add secondaryStructure, exonsRegion, domains aspect to it
    return sql_iter

def domseqiterator(genencbid, trans, transfk, transdoms, genecodes, domainhash):
    # transcript, trans_pk_temp[transcript.ID], dom_trans_spans,has_dom, exons_domainseq)
    # transcriptOb, transPrimID, spans as Key:(pfamidName) as value, pfamidName as key: (alphbet and colname) as value, 
    # domainHas will be updated
    # print (trans.ID, trans.seqlen)
    strdom = ["-"]*trans.seqlen
    lengthInitial=trans.seqlen
    # print ("init:",lengthInitial)
    if transdoms:
        for i in transdoms:
            # print (i)
            start=i[0]-1
            end=i[1]
            if start > lengthInitial:
                continue
            if end >lengthInitial:
                end = lengthInitial
            # above segment will help parse cases where exon intervals wont contain co mplete protein sequence, for exmaple in gene 170676, and isoform NP_570947.2
            # however partial domain assignemnt can still be done and i Ill be douing for segment covered, if start exceeds the seqlen, then I wont proceed, 
            # but if end exceeds, and start doesnt then i will modify the end
            domcode = genecodes[transdoms[i]][0]
            overlap=0
            length=end-start+1
            for j in range (start,end):
                if strdom[j]=='-':
                    strdom[j]=domcode
                else:
                    overlap+=1
                    # pass
                    # print ("Raise Flag, %s, has double dom assignment prob at index=%s"%(trans.ID, j))
            # only for case log reporting, can be hashed
            if overlap:
                print ("overlap, %s, %s, %s"%(genencbid, trans.ID, round(float(overlap)/length,2)))
    strdom = "".join(strdom)        
    lengt = 0
    for exons in trans.exons:
        dom_seq = strdom[lengt:(lengt+exons.length)]
        lengt += exons.length

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


def exonObcreator(genefk, exon, exonpk, trans_fkhash, exons_domainseq, exonssspk, exonsdispk, exonsdompk):
    # ############################# doc #######################################
    # reads PSI, exid, length, codst end, raw st end, ssaeq, parent also
    # in fitrst pass it writes generic information to the tables
    # in the sxceond third pass, it wonly writes if it has something values to be shared, like ss and disoredred seq
    # why i dont wiorte domainss equqnce like thta >
    sqliterations = []
    wif = exon.WEF
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
                      (parent, exonpk, str(wif), exId, aaseq, length, codst, codend, rawst,
                       rawend, sstype, genefk)]

    if aaseq:
        if exon.ssseq:
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
        if exon.disorderseq:
            for seq in exon.disorderseq:
                if seq is not None and seq != "NULL":
                    translist = ",".join(exon.disorderseq[seq])
                    for tri in trans_fkhash:
                        translist = re.sub(tri, str(trans_fkhash[tri]), translist)
                    query = "INSERT INTO exonapp_exonsdisorder (id, list_trans_fk,exon_id, disseq) VALUES ('%d', '%s',                     '%d','%s')" % (
                        exonsdispk, translist, exonpk, re.sub("S", "-", seq))
                    sqliterations += [query]
                    exonsdispk += 1

        if exon.ID in exons_domainseq:
            for seq in exons_domainseq[exon.ID]:
                query = "INSERT INTO exonapp_exonsdomseq  (id, list_trans_fk,exon_id, domseq) VALUES ('%d', '%s',                 '%d','%s')" % (
                    exonsdompk, exons_domainseq[exon.ID][seq], exonpk, seq)
                sqliterations += [query]
                exonsdompk += 1
    return exonssspk, exonsdispk, exonsdompk, sqliterations

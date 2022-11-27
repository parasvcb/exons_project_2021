import json,re,sys
with open ("oldSet/103_pretty.json") as fin:
    job=json.load(fin)

##########################################
###    storing data in simple dicts    ###
trans={}
exons={}
domains={}
for i in job['trans']:
    trans[i['id']]=i
    #appedning whole object
for i in job["exonsGenes"]:
    exons[i['id']]=i
for i in job['domainCod']:
    code=i['code']
    domains[code]={}
    domains[code]['color']=i['color']
    domains[code]['name']=i['dId']
    domains[code]['pfam']=i['name']
####                                  ####
##########################################

def seqRet (ex,trans,key='ss'):
    #trans here is the only foreiogn key
    if key=='ss':
        val='ssseq'
        parentCat='exonsSS'
    if key =="dom":
        val='domseq'
        parentCat='exonsDom'
    if key =="dis":
        val='disseq'
        parentCat='exonsDis'

    if len(ex[parentCat]):
        tempss={}
        for miniD in ex[parentCat]:
            if str(trans) in miniD["list_trans_fk"]:
                return miniD[val]
    return ''  
### give corresponding sequqnce of transcript property, based on its listing in exon


for t in trans:
    seqdom=''
    seqdis=''
    seqss=''
    aaseq=''
    lengthDone=1
    strExon=''
    for ex in trans[t]['exonsIds'].split(','):
        ex=int(ex)
        smoothDis=0
        modExo=exons[ex]["exId"]+'$'        
        # above two are related to exon shape rendering
        # modExon will be "T.1.A.3.0.0:1,10,[0-4]" 1 to 10 its aa span and 0-4 are shapes of the boundary, 0, round both sides (default), 1 zag left (means bopundaries are longer than contributed span,), 2 zag right (3; UTR), 3 single exon case in tranrcipt and both UTRS in this. 4 absence (no aa, YET TO BE WORKED)
        aa=exons[ex]['aaseq']
        dom=seqRet(exons[ex],t,'dom')
        dis=seqRet(exons[ex],t,'dis')
        ss=seqRet(exons[ex],t,'ss')
        if (aa):
            cond1=abs(exons[ex]['codst']-exons[ex]['rawst'])>3
            cond2=abs(exons[ex]['codend']-exons[ex]['rawend'])>3
            if cond1 and cond2:
                smoothDis=3
            elif cond1:
                smoothDis=1
            elif cond2:
                smoothDis=2
            modExo+='%s,%s,%s_'%(lengthDone,lengthDone+len(aa),smoothDis)
            lengthDone+=len(aa)
        else:
            modExo+='%s,%s,%s,%s,%s_'%(lengthDone,lengthDone,4,exons[ex]['rawst'], exons[ex]['rawend'])
            
            
        strExon+='%s'%modExo
        aaseq+=aa
        seqdom+=dom
        seqdis+=dis
        seqss+=ss
    
    strh='H:'
    stre='E:'
    strdis='S:'
    
    ssh=re.finditer('H{1,}',seqss)
    sse=re.finditer('E{1,}',seqss)
    diss=re.finditer('S{1,}',seqdis)
    
    for i in ssh:
        sp=i.span()
        strh+='%s,%s_'%(sp[0]+1,sp[1]+1)
    for i in sse:
        sp=i.span()
        stre+='%s,%s_'%(sp[0]+1,sp[1]+1)
    for i in diss:
        sp=i.span()
        strdis+='%s,%s_'%(sp[0]+1,sp[1]+1)
    #trans[t]["structuredRegion"]=strh+'-'+stre
    trans[t]["secondaryStructure"]=strh[:-1]+'-'+stre[:-1]+'-'+strdis[:-1]
    #last element is skpped because of trailing undersocre used above
    trans[t]["exonsRegion"]=strExon[:-1]

    ##########################################################
    #                         <domains>                      #
    # get list of unique alphabes represneting domains,      # 
    # iterate alphabets as domains, and serach for their     #
    # spans. one repeition entry will only give one span, 
    
    uniqDom=list(set(seqdom))
    uniqDom=[i for i in uniqDom if i  !='-']
    listDomString=[]
    for i in uniqDom:
        #print (i)
        domSearch=re.finditer(r'%s{1,}'%i,seqdom)
        strd='%s,%s,%s:'%(domains[i]['name'],domains[i]['pfam'],domains[i]['color'])
        #print (strd)
        for span in domSearch:
            sp=span.span()
            strd+='%s,%s_'%(sp[0]+1,sp[1]+1)
        #print (strd)
        #print ('*',listDomString)
        #sys.exit()
        listDomString+=[strd[:-1]]
    #print (listDomString)
    trans[t]["domains"]='$'.join(listDomString)
    #break

job["trans"]=list(trans.values())

with open("oldSet/sample.json", "w") as outfile:
    json.dump(job, outfile, indent=4)
#print (job.keys())
import cPickle as pickle
import numpy as np
import scipy,os
import common.general_modules as gm

#genes
def exonScreenerBetweenTConstitutive(geneExons,inclusive=True):
    acceptableNumericFlagList=[]
    for ex in geneExons:
        codNonCodFlag=ex.ID.split(".")[0]
        if codNonCodFlag=='T':    
            constAltFlag=ex.ID.split(".")[2]
            if constAltFlag=='G':
                numericFlag=int(ex.ID.split(".")[3])
                acceptableNumericFlagList+=[numericFlag]
    acceptableNumericFlagList.sort()
    if not inclusive:
        acceptableNumericFlagList=acceptableNumericFlagList[1:-1]
    return acceptableNumericFlagList

def exonPasser(e1,e2, consecCoding=[]):
    baseCondition=e1.length>9 and e2.length>9 and e1.ID[0]=='T' and e2.ID[0]=='T'
    if consecCoding:
        choosenExonCondition=int(e1.ID.split(".")[3]) in consecCoding and int(e2.ID.split(".")[3]) in consecCoding
        return True if (baseCondition and choosenExonCondition) else False
    else:
        if baseCondition:
            return True
        else:
            return False

def csv_writer(has, filename,fout):
    key = has.keys()
    total = sum(has.values())
    key.sort()
    fout.write("\n%s"%filename)
    with open(filename, "w") as fin:
        fin.write("SSType,Total,Freq\n")
        fout.write("\nSSType\tTotal\tFreq")
        for i in key:
            fin.write("%s,%s,%s\n" % (i, has[i], gm.div_fact(has[i],total)))
            fout.write("\n%s\t%s\t%s" % (i, has[i], gm.div_fact(has[i],total)))
            

def givemeresnoofpdb(pdb):
    with open (pdb) as fin:
        dat=fin.read().split('\n')
    seq=[]
    for i in dat:
        if len(i)>10 and i[:4]=='ATOM':
            seq+=[int(i[22:26])]
    seq=list(set(seq))
    seq.sort()
    return seq

def parseRsa(fname):
    with open (fname) as fin:
        dat = fin.readlines()
    SAR_pdb = {}  # residue number
    for lines in dat:
        if lines.startswith("RES"):
            res_num = int(lines[10:13].strip())
            res = lines[4:7].strip()
            rel_all = float(lines[22:28].strip())
            SAR_pdb[res_num] = 'E' if rel_all >5 else 'B'
    return SAR_pdb

def subsetsurfaceexposedExonWise(surfaceExposeddata,exonspanrange):
    string=''
    #print (surfaceExposeddata,'SEDATA')
    #print (exonspanrange,'EXDATA')
    for i in exonspanrange:
        string+=surfaceExposeddata[i] if i in surfaceExposeddata else 'Z'
    return string

def rsaToSeqFile(proteinPDB_seq,PI,naccessLink):
    #print (naccessLink,'naccessFile')
    #print (proteinPDB_seq,'pdbseq')
    seqrange=0
    ss_seq=''
    exonspanhas={}
    surfaceExposeddata=parseRsa(naccessLink)
    #has with resno as key and E/B as value
    for i in PI.exons:
        if i.length>0:
            exonspan=seqrange+i.length
            exonspanrange=list(range(seqrange,exonspan))
            exonspanhas[i]=exonspanrange
            surfacestring=''
            #print (exonspanrange,'exrange')
            if set(exonspanrange)<=set(proteinPDB_seq):
                #print ('*******yes\n\n')
                #means complete exon shold be resolved
                surfacestring=subsetsurfaceexposedExonWise(surfaceExposeddata,exonspanrange)
            seqrange+=exonspan
            ss_seq+=surfacestring
    return surfaceExposeddata,exonspanhas,ss_seq

#transcripts
def fastareturn(string):
    out=[]
    size=80
    out=[]
    for i in range (0,len(string),size):
        out+=[string[i:i+size]]
    res="\n".join(out)
    return res

def givemePI(gene,sortby='PI'):
    mat=[]
    index=0 if sortby=='PI' else 2 if sortby=="codingExon" else 3
    for trans in gene.transcripts:
        '''
        add a matrxi if theres a stride seq available, and then sum it all up 
        '''
        #print (trans.PI,len([exon for exon in trans.exons if exon.length>0]),trans.structure_lenAF,trans)
        mat += [[trans.PI,len([exon for exon in trans.exons if exon.length>0]),trans.structure_lenAF,trans, trans.co]]
    mat.sort()
    transInterest=0
    for i in mat[::-1]:
        if mat[-1][2]>0:
            transInterest=mat[-1][3]
            break
    return transInterest if transInterest else False

def split_cases(exonslist,out='default'):
    # splits alt and majorly into unsplitted forms and forms with variations 
    '''
    gets ids extracted from the list fo exons 
    the list recieved will be either U or T or D cases but not their intermix
    will screen from the first two letters, if their varitaions do exist (will append their n,c,b,0 case as key value)
    '''
    idlistHas={}
    idlistHasRef={}
    #idlist=[i.ID for i in exonslist]
    for ex in exonslist:
        tag=ex.ID.split('.')
        key=".".join(tag[:2])
        value=tag[4]
        if key not in idlistHas:
            idlistHas[key]=set()
            idlistHasRef[key]=[]
        idlistHas[key]|=set([value])
        idlistHasRef[key]+=[ex]
    pureAlternate=[idlistHasRef[key][0] for key in idlistHas if len(idlistHas[key])==1]
    alternateWithSS=[idlistHasRef[key][0] for key in idlistHas if len(idlistHas[key])>1]
    if out =='default':
        return len(pureAlternate),len(alternateWithSS)
    return pureAlternate,[idlistHasRef[key] for key in idlistHas if len(idlistHas[key])>1]
    #sending ouytput in 1st var as list of alternate xons with no change in splice site and n dsceond case as the total possible exons in all teh trancripts (redundant)

#exons
def ss_junctionWindow(ex1, ex2, window):
    if not (ex1 and ex2):
        return 'NA'
    # for cases in stride and other types, if not defined then do the NA
    n = ex2[:window]
    c = ex1[-(window):]
    if len(set(n)) == 1:
        nr = list(set(n[:3]))[0]
    else:
        nr = "X"

    if len(set(c)) == 1:
        cr = list(set(c))[0]
    else:
        cr = "X"
    res = cr+nr
    return res


def key_searcher(box,key):
    box.sort()
    if key in box:
        return key
    else:
        small=min(key)
        big=max(key)
        for i in box:
            if i[0]==small or i[1] == big:
                return i

def sorting_screening_junctions(has,fname_junctions,fname_junc_conservations, default_sstypes=['HH','EE','CC']):
    #has_exon_categories_combined[flagPairKey][gene,exonSequencePair,WEF]=(interestType,reprSS,reprSS_count,reprSS_Freq,temphas)
    '''
    this serves two purposes
    from parent 'has' dictionary recieved, it will iterate the exontype first (AA,AG,GG)
        from there it will start iterating the exon pairs:
            will define variables like, a) WEF of the exon junction (irresptive if SS is defined or not, but definitely coding)
            b) junctionSS, its count and freq
            -> the freq will be updated in the has called the WEFRangeHas
            the range will be later used to update another has called conservationsHas (seg2 initiated)

    '''
    f1=open(fname_junctions,'w')
    f1.write('ExonType\tWEFRange\tJunction\tCount\tFreq\n')
    conservation_has={exonType:{} for exonType in has} 
    # this will be used to get data for the cvonservation/occurence of the Junctions in the isoforms, segregated first into exonTypes, later into their occurenec range (HH is observed in 23 times out of 100 when junction type is ALTALt) and then in to HH, how oftne this was presnerved in next junctions
    for exonType in has:
        #exonType is conscons altalt altcons
        WEFRangeHas={(0,0.25):{},(0.25,0.50):{},(0.50,0.75):{},(0.75,1.01):{}}
        #seg 1 begins
        for geneExpairWefKey in has[exonType]:
            #geneExpairWefKey is tuple of (gene,Expair,Wef)
            wef=geneExpairWefKey[2]
            exonsListdetailed,junction,junctionCount,JunctionFreq,temphas=has[exonType][geneExpairWefKey]
            #print (exonsListdetailed)
            #print (junction,junctionCount,JunctionFreq,temphas)
            #has_exon_categories_combined[flagExonType2][gene,exonSequencePair,WEF]=(representative[exonSequencePair],reprSS,reprSS_count,reprSS_Freq)
            #where representative[exonSequencePair] is list of list and every list has PI,coding_exons_count,texo[i].ID,texo[i+1].ID,SS_junction+[reprSS,reprSS_count,reprSS_Freq]
            for psiExonRange in WEFRangeHas:
                #print (psiExonRange,wef)
                if psiExonRange[0]<wef<=psiExonRange[1]:
                    #print ('yes')
                    if junction not in WEFRangeHas[psiExonRange]:
                        WEFRangeHas[psiExonRange][junction]=[0,[]]
                    WEFRangeHas[psiExonRange][junction][0]+=1
                    WEFRangeHas[psiExonRange][junction][1]+=[JunctionFreq]
                    break
            # seg1, WEFRangeHas updated
            if psiExonRange not in conservation_has[exonType]:
                conservation_has[exonType][psiExonRange]={}
            # seg2, initiated, we are still in loop 2 (iterating major exjunc types) then individual pairs in them, 
            # where we are considering adding PSI value of junction to anotherhas which will calculate the occurence of major HHTYPES in them 
            for junctionAll in temphas: 
                #print (junctionAll)
                #print (temphas[junctionAll])   
                if junctionAll in default_sstypes:
                    #seg 2 now being pro[perly resumed]
                    if junctionAll not in conservation_has[exonType][psiExonRange]:
                        conservation_has[exonType][psiExonRange][junctionAll]={(0,0.10):0,(0.101,0.20):0,(0.201,0.30):0,(0.301,0.40):0,(0.401,0.50):0,(0.501,0.60):0,(0.601,0.70):0,(0.701,0.80):0,(0.801,0.90):0,(0.901,1.00):0}
                    for vals2 in conservation_has[exonType][psiExonRange][junctionAll]:
                        if vals2[0]<=temphas[junctionAll]<=vals2[1]:
                            conservation_has[exonType][psiExonRange][junctionAll][vals2]+=1
                            break
        # seg2, is now updated
        # resuming seg1 below
        list_WEF=list(WEFRangeHas.keys())
        list_WEF.sort()
        for WEFrange in list_WEF:
            #f1.write('ExonType\tWEFRange\tJunction\tCount\tFreq\n') written already
            total_junctions=sum([WEFRangeHas[WEFrange][junc][0] for junc in WEFRangeHas[WEFrange]])
            for junctions in WEFRangeHas[WEFrange]:
                count=WEFRangeHas[WEFrange][junctions][0]
                f1.write('%s\t%s\t%s\t%s\t%s\n'%(exonType,"_".join(map(str,list(WEFrange))),junctions,count,gm.div_fact(count,total_junctions)))
    #seg 1 done
    #seg 2, now lets do the second aspect of calculating the fraction
    f2=open(fname_junc_conservations,'w')
    f2.write('ExonType\tWEFRange\tJunction\tOccurenceRange\tOccureneceTotal\tOccureneceFreq\tcumulativeCount\tcumulativeFreq\n')
    a=conservation_has.keys()
    # for key in a:
    #     print (key,conservation_has[key].keys())
    
    #print (conservation_has)
    for exonType in conservation_has:
        for WEFRangeJunction in conservation_has[exonType]:
            wefrangestr="_".join(map(str,list(WEFRangeJunction)))
            for junction in conservation_has[exonType][WEFRangeJunction]:
                totalOcc=sum([conservation_has[exonType][WEFRangeJunction][junction][ti] for ti in conservation_has[exonType][WEFRangeJunction][junction]])
                occRangeList=list(conservation_has[exonType][WEFRangeJunction][junction].keys())
                occRangeList.sort()
                previousValue=0
                for occRange in occRangeList:
                    occrangestr="_".join(map(str,list(occRange)))
                    count_occrange=conservation_has[exonType][WEFRangeJunction][junction][occRange]
                    cumulative_count=previousValue+count_occrange
                    previousValue+=count_occrange
                    #print (1,exonType,WEFRangeJunction)
                    #print (2,exonType,wefrangestr,junction,occrangestr,count_occrange,gm.div_fact(count_occrange,totalOcc), cumulative_count,gm.div_fact(cumulative_count,totalOcc))
                    f2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(exonType,wefrangestr,junction,occrangestr,count_occrange,gm.div_fact(count_occrange,totalOcc), cumulative_count,gm.div_fact(cumulative_count,totalOcc)))



def junctions_from_different_isoforms_poulator(repr_has, trans,PI, window,geneExonsShortlisted =False, normalorStide='normal'):
    '''
    for a trancript:
        iterate its exons in pairs, 
        if their length is grt 9 (both):
            get their ss seqs
            if they are not retention cases and they are not empty ss's
                get their linear sequqnce tags or close to them ()
                now get their SS pair for junction
                for the linear sequqnce tag, appen informtion to the has
    return
    '''
    texo=trans.exons[:]
    coding_exons_count=len([i for i in texo if i.length>0])
    for i in range (0, len(texo)-1):
        if exonPasser(texo[i],texo[i+1], consecCoding=geneExonsShortlisted):
            if normalorStide=='normal':    
                i1seq=texo[i].out_secondseq(trans.ID)
                i2seq=texo[i+1].out_secondseq(trans.ID)
            else:
                i1seq=texo[i].out_strideseqAF(trans.ID)
                i2seq=texo[i+1].out_strideseqAF(trans.ID)
            part1=int(texo[i].ID.split(".")[3])
            part2=int(texo[i+1].ID.split(".")[3])
            #print (repr_has.keys())
            #print (part1,part2,texo[i].ID,texo[i+1].ID)
            #exon_sequence=key_searcher(repr_has.keys(),(part1,part2))
            #print (exon_sequence)
            exon_sequence=(part1,part2)
            SS_junction=ss_junctionWindow(i1seq,i2seq,window)
            if exon_sequence not in repr_has:
                repr_has[exon_sequence]=[]
            repr_has[exon_sequence]+=[[PI,coding_exons_count,texo[i].ID,texo[i+1].ID,SS_junction]]
    return repr_has

def junctions_from_different_isoforms_poulatorSURFACEEXPOSED(repr_has, trans,PI, window, pdb_dir,naccess_dir):
    '''
    for a trancript:
        iterate its exons in pairs, 
        if their length is grt 9 (both):
            get their ss seqs
            if they are not retention cases and they are not empty ss's
                get their linear sequqnce tags or close to them ()
                now get their SS pair for junction
                for the linear sequqnce tag, appen informtion to the has
    return
                exonspanhas={}
            for i in PI.exons:
                if i.length>0:
                    exonspan=seqrange+i.length
                    exonspanrange=list(range(seqrange,exonspan))
                    exonspanhas[i]=exonspanrange
                    if set(exonspanrange)<=set(proteinPDB_seq):
                        #means complete exon shold be resolved
                        surfacestring=subsetsurfaceexposedExonWise(srfaceData,exonspanrange)
                    seqrange+=exonspan
                    ss_seq+=surfacestring
    '''
    texo=trans.exons[:]
    flag=False
    if trans.structure_fileAF:
        pdbF=os.path.join(pdb_dir,trans.structure_fileAF)
        pdbresno=givemeresnoofpdb(pdbF)
        naccessFile=os.path.join(naccess_dir,trans.structure_fileAF+'.rsa')
        srfaceData,exonspanhas,ss_seq=rsaToSeqFile(pdbresno,trans,naccessFile)
        if ss_seq:
            flag=True
            #print (srfaceData)
            #print (exonspanhas)
            #print (ss_seq)

    coding_exons_count=len([i for i in texo if i.length>0])
    for i in range (0, len(texo)-1):
            if exonPasser(texo[i], texo[i+1],consecCoding=False):
                #changed !=R tp ==T
                #print ('*************no\n\n')
                i1seq=subsetsurfaceexposedExonWise(srfaceData,exonspanhas[texo[i]]) if flag else False
                i2seq=subsetsurfaceexposedExonWise(srfaceData,exonspanhas[texo[i+1]]) if flag else False
                part1=int(texo[i].ID.split(".")[3])
                part2=int(texo[i+1].ID.split(".")[3])
                exon_sequence=(part1,part2)
                SS_junction=ss_junctionWindow(i1seq,i2seq,window)
                if exon_sequence not in repr_has:
                    repr_has[exon_sequence]=[]
                repr_has[exon_sequence]+=[[PI,coding_exons_count,texo[i].ID,texo[i+1].ID,SS_junction]]
    return repr_has




def refineFcategory(has):
    #print (has.keys()[:5])
    #print (dir(has[103]))
    for gene in has:
       #if gene ==103:
        exhas={}
        #print (gene,has[gene].transcripts)
        Tcount=len(has[gene].transcripts)
        for transcripts in has[gene].transcripts:
            for exons in transcripts.exons:
                exID=exons.ID
                if exID[0]!='R':
                    ele=exID.split('.')
                    seqNumber=ele[3]
                    flag=ele[2]
                    if flag!='G':    
                        if seqNumber not in exhas:
                            exhas[seqNumber]=[]
                        exhas[seqNumber]+=[transcripts.ID]
        changeNeeded={i:0 for i in exhas if len(set(exhas[i]))==Tcount}
        # print {i:len(set(exhas[i])) for i in exhas} ,Tcount
        # print changeNeeded
        if changeNeeded:
            for transcripts in has[gene].transcripts:
                for exons in transcripts.exons:
                    exID=exons.ID
                    if exID[0]!='R':
                        ele=exID.split('.')
                        if ele[3] in changeNeeded and ele[2]!='G':
                            ele[2]='F'
                            exons.ID=".".join(ele)
            for exons in has[gene].exons:
                exID=exons.ID
                if exID[0]!='R':
                    ele=exID.split('.')
                    if ele[3] in changeNeeded and ele[2]!='G':
                        ele[2]='F'
                        exons.ID=".".join(ele)
        # for i in has[103].exons:
        #      print (i.ID)
    #one more leg of refinment needed to change such ids in the IR retention cases
    return has

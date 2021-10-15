import cPickle as pickle
import numpy as np
import scipy
def div_fact(num,denom):
        try:
                return round(float(num)/denom,3)
        except:
                return 0

def readPickle(fname):
    with open (fname) as fin:
        has=pickle.load(fin)
    return has

def writePickle (fname,has):
    with open(fname, "w") as fin:
        pickle.dump(has, fin)

def stats(lis):
    lis=[i for i in lis if i>0]
    if len(lis):
        string="\t".join(map(str,[len(lis),sum(lis),np.max(lis),np.min(lis),round(np.mean(lis),3),round(np.median(lis),3),scipy.stats.mode(lis)[0][0],scipy.stats.mode(lis)[1][0],round(np.std(lis),3)]))
    else:
        string="\t".join(map(str,[len(lis),sum(lis),0,0,0,0,0,0,0]))
    return string

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


def ss_junctionWindow(ex1, ex2, window):
    if not (ex1 and ex2):
        return 'NA'
    # for cases in stride and other types, if not defined then do the NA
    n = ex2[:window-1]
    c = ex1[-(window-1):]
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

def csv_writer(has, filename,fout):
    key = has.keys()
    total = sum(has.values())
    key.sort()
    fout.write("\n%s"%filename)
    with open(filename, "w") as fin:
        fin.write("SSType,Total,Freq\n")
        fout.write("\nSSType\tTotal\tFreq")
        for i in key:
            fin.write("%s,%s,%s\n" % (i, has[i], div_fact(has[i],total)))
            fout.write("\n%s\t%s\t%s" % (i, has[i], div_fact(has[i],total)))


def sorting_screening_junctions(has,fname_junctions,fname_junc_conservations):
    #has_exon_categories_combined[flagPairKey][gene,exonSequencePair,WEF]=(interestType,reprSS,reprSS_count,reprSS_Freq)
    f1=open(fname_junctions,'w')
    f1.write('ExonType\tWEFRange\tJunction\tCount\tFreq\n')
    conservation_has={exonType:{} for exonType in has} 
    # this will be used to get data for the cvonservation/occurence of the Junctions in the isoforms, segregated first into exonTypes, later into their occurenec range (HH is observed in 23 times out of 100 when junction type is ALTALt) and then in to HH, how oftne this was presnerved in next junctions
    for exonType in has:
        #exonType is conscons altalt altcons
        WEFRangeHas={(0,0.25):{},(0.25,0.50):{},(0.50,0.75):{},(0.75,1.01):{}}
        for geneExpairWefKey in has[exonType]:
            #geneExpairWefKey is tuple of (gene,Expair,Wef)
            wef=geneExpairWefKey[2]
            exonsListdetailed,junction,junctionCount,JunctionFreq=has[exonType][geneExpairWefKey]
            #has_exon_categories_combined[flagExonType2][gene,exonSequencePair,WEF]=(representative[exonSequencePair],reprSS,reprSS_count,reprSS_Freq)
            #where representative[exonSequencePair] is list of list and every list has PI,coding_exons_count,texo[i].ID,texo[i+1].ID,SS_junction+[reprSS,reprSS_count,reprSS_Freq]
            for tupleRange in WEFRangeHas:
                #print (tupleRange,wef)
                if tupleRange[0]<wef<=tupleRange[1]:
                    #print ('yes')
                    if junction not in WEFRangeHas[tupleRange]:
                        WEFRangeHas[tupleRange][junction]=[0,[]]
                    WEFRangeHas[tupleRange][junction][0]+=1
                    WEFRangeHas[tupleRange][junction][1]+=[JunctionFreq]
                    break
            if tupleRange not in conservation_has[exonType]:
                conservation_has[exonType][tupleRange]={}
            if junction in ['HH','EE','CC']:
                if junction not in conservation_has[exonType][tupleRange]:
                    conservation_has[exonType][tupleRange][junction]={(0,0.10):0,(0.101,0.20):0,(0.201,0.30):0,(0.301,0.40):0,(0.401,0.50):0,(0.501,0.60):0,(0.601,0.70):0,(0.701,0.80):0,(0.801,0.90):0,(0.901,1.00):0}
                for vals2 in conservation_has[exonType][tupleRange][junction]:
                    if vals2[0]<=JunctionFreq<=vals2[1]:
                        conservation_has[exonType][tupleRange][junction][vals2]+=1
                        break
        #Here The WEFRangeHas has been filled
        list_WEF=list(WEFRangeHas.keys())
        list_WEF.sort()
        for WEFrange in list_WEF:
            total_junctions=sum([WEFRangeHas[WEFrange][junc][0] for junc in WEFRangeHas[WEFrange]])
            for junctions in WEFRangeHas[WEFrange]:
                count=WEFRangeHas[WEFrange][junctions][0]
                f1.write('%s\t%s\t%s\t%s\t%s\n'%(exonType,"_".join(map(str,list(WEFrange))),junctions,count,div_fact(count,total_junctions)))
    #1st aspect is written
    #now lets do the second aspect of calculating the fraction
    f2=open(fname_junc_conservations,'w')
    f2.write('ExonType\tWEFRange\tJunction\tOccurenceRange\tOccureneceTotal\tOccureneceFreq\n')
    a=conservation_has.keys()
    for key in a:
        print (key,conservation_has[key].keys())
    
    for exonType in conservation_has:
        for WEFRangeJunction in conservation_has[exonType]:
            for junction in conservation_has[exonType][WEFRangeJunction]:
                totalOcc=sum([conservation_has[exonType][WEFRangeJunction][junction][ti] for ti in conservation_has[exonType][WEFRangeJunction][junction]])
                for occRange in conservation_has[exonType][WEFRangeJunction][junction]:
                    print (1,exonType,WEFRangeJunction)
                    print (2,exonType,"_".join(map(str,list(WEFRangeJunction))),junction,"_".join(map(str,list(occRange))),conservation_has[exonType][WEFRangeJunction][junction][occRange],div_fact(conservation_has[exonType][WEFRangeJunction][junction][occRange],totalOcc),'\n')
                    f2.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(exonType,"_".join(map(str,list(WEFRangeJunction))),junction,"_".join(map(str,list(occRange))),conservation_has[exonType][WEFRangeJunction][junction][occRange],div_fact(conservation_has[exonType][WEFRangeJunction][junction][occRange],totalOcc)))


def junctions_from_different_isoforms_poulator(repr_has, trans,PI, window,normalorStide='normal'):
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
        if texo[i].length>9 and texo[i+1].length>9:
            if normalorStide=='normal':    
                i1seq=texo[i].out_secondseq(trans.ID)
                i2seq=texo[i+1].out_secondseq(trans.ID)
            else:
                i1seq=texo[i].out_strideseqAF(trans.ID)
                i2seq=texo[i+1].out_strideseqAF(trans.ID)
            if texo[i].ID[0]!='R' and texo[i+1].ID[0]!='R' and i1seq != False and i2seq != False:
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

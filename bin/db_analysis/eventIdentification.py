'''
Before creation of object Builder type 2, the object needs pre filtering for the crrect depiction of the panels having keywords as AGF in exon nomenclature 
'''

import sys
import re, os
from constructing_data import Classes_exons
import common.general_modules as cm
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
from progress.bar import Bar

matrix = matlist.blosum62

if len(sys.argv)!=3:
        print ("Please give 1 object 2 outalnDir")
        sys.exit()
prog,source_gene_object, outdir = sys.argv

def stats (aln):
    # identity, cov1, cov2
    seq1 = aln[0]
    seq2 = aln[1]
    match = 0
    nongaps = 0 
    for ind, val in enumerate (seq1):
        if seq2[ind] == val:
            match +=1
        if seq2[ind]!='-' and val!='-':
            nongaps +=1
    seq1Length = len(seq1) - seq1.count('-')
    seq2Length = len(seq2) - seq2.count('-')
    identity1 = round(float(match)/seq1Length,3)
    identity2 = round(float(match)/seq2Length,3)
    cov1 = round(float(nongaps)/seq1Length,3)
    cov2 = round(float(nongaps)/seq2Length,3)
    return seq1Length, seq2Length, identity1, identity2, cov1, cov2
    
def refinestandard(s):
    # A CDEFGHI KLMN PQRST VW Y
    #  B       J    O     U  X Z
    lischange = ["B","J","O","U","X","Z"]
    for i in lischange:
        if i in s:
                s = re.sub(r'%s'%i,'X',s)
    return s
humanGeneObject=cm.readPickle(source_gene_object)


filestats = "name\tgene\tisoform\tisoformCount\texonCount\tACase\tDcase\tncbCase\taaChangeCases\tfirstOrlast\taaChangeInExon\tfractionLength\tid\tcov\tstatLine\texonList\n"

bar = Bar('Processing genes:', max=len(humanGeneObject))
for gene in humanGeneObject:
    first = ''
    last = ''
    seqlist = []
    hasaaseq= {} # pos and default (1) var seq
    hasaaChangeLog = {} #ID and seq
    aalisExons = ""
    for ex in humanGeneObject[gene].exons:
        if ex.ID[0]!='R':
            utmrd, code, cons, pos, ncb, occ = ex.ID.split('.')
            aalisExons += "%s: %s\n"%(ex.ID, ex.seq)
            seqlist += [int(pos)]
            if int(code)==1:
                hasaaseq[(int(pos), ncb, int(occ))]=ex.seq
                # 4,n,0, 1 | 4,0,0, 1| 4,c,0, 1| 4,c,1, 1 || 1st default aa sequence is given, last diogit is only imaginary
                # hasaaseq[int(pos)]=ex.seq
            if int(code) >=2 and len(ex.seq)>9: #-> shall i add length condition here
                tuplevar = (int(pos), ncb, int(occ), int(code))
                # 4,n,0,1, 4,n,0,2, 4,0,0,0  4,0,0,2 | 4,c,0,2
                if tuplevar not in hasaaChangeLog:
                    hasaaChangeLog[tuplevar] = ex.seq
    seqlist.sort()
    first = seqlist[0]
    last = seqlist[-1]
    '''
    Above phase has so far only recorded the aaseq with keys 4,n,0 | 4,0,0, | 4,c,0 considering on;y the default 1st var occurence, no n1n2n3 forms have been considered, though they may be recoerded like 4,n,1 with aaseq but for the 1st ocuurenecce 
    '''
    nameAddition = ""
    if hasaaChangeLog:
        flagWrite = False
        removeEntry = []
        stringAln = ''
        for exon in hasaaChangeLog:
                pos, ncb, occ, code = exon 
                seq2 = hasaaChangeLog[exon]
                try:
                        seq1 = hasaaseq [(pos, ncb, occ)]
                except Exception as E:
                        removeEntry += [exon]
                        print (E)
                        print (gene)
                        print (aalisExons)
                        # sys.exit()
                        continue
                lengthFraction = round(float(len(seq2))/len(seq1)*100, 3) if len(seq1) and len(seq2) else 0
                seq1 = refinestandard (seq1)
                seq2 = refinestandard (seq2)
                
                align = pairwise2.align.globaldx(seq1, seq2, matrix)
                count = 0
                stringAln += 'exon_%s_upperloc_exonvar_%s_lowerloc\n'%(pos,code)
                # 
                for a in align:
                        flagWrite = True
                        if count == 0:
                                seq1Length, seq2Length, identity1, identity2, cov1, cov2 = stats(a)
                                if 0.8<=identity2<=0.9:
                                        nameAddition ='1' 
                                elif 0.7<=identity2<=0.8:
                                        nameAddition ='2'
                                elif 0.6<=identity2<=0.7:
                                        nameAddition ='3'
                                elif 0.5<=identity2<=0.4:
                                        nameAddition ='4'
                                elif 0.4<=identity2<=0.3:
                                        nameAddition ='5'
                                

                        seq1Le, seq2Le, ide1, ide2, cv1, cv2 = stats(a)
                        stringAln += "type %s/%s\n"%(count+1,len(align))
                        stringAln += format_alignment(*a, full_sequences=True)
                        stringAln += "\nseq1Length:%s, seq2Length:%s, identity1:%s, identity2:%s, cov1:%s, cov2:%s\n\n"%(seq1Le, seq2Le, ide1, ide2, cv1, cv2)         
                        count += 1       
                hasaaChangeLog[exon] = [seq2, lengthFraction, seq1Length, seq2Length, identity1, identity2, cov1, cov2]
        if removeEntry:
                hasaaChangeLog = {i: hasaaChangeLog[i] for i in hasaaChangeLog if i not in removeEntry}
        if flagWrite:
                with open (os.path.join(outdir,"%s_exons.aln"%(gene)),'w') as fout:
                        fout.write("%s"%(stringAln))

    for trans in humanGeneObject[gene].transcripts:
        string=''
        events= [0,0,0]
        aaChange = []
        aaChangeCases = 0
        if trans.ID[0]=='N' and trans.seqlen < 1000:
                string = " | ".join([ex.ID for ex in trans.exons])
                for ex in trans.exons:
                    if ex.ID[0]!='R':
                        utmrd, code, cons, pos, ncb, occ = ex.ID.split('.')
                        tuplevar = (int(pos), ncb, int(occ), int(code))
                        if tuplevar in hasaaChangeLog:
                                aaChangeCases += 1
        dCount = string.count('D')
        Asps = string.count('n') + string.count('c') + string.count('b')
        acount = string.count('F')
        flag = False
        
        if dCount:
            events[0] += dCount
            flag = True
        if Asps:
            events[1] += Asps
            flag = True
        if acount:
            events[2] += acount
            flag = True
        if flag or aaChangeCases:
                for ex in trans.exons:
                        if ex.ID[0]!='R':
                                utmrd, code, cons, pos, ncb, occ = ex.ID.split('.')
                                aaChangeInExon = "F"
                                identity = "F"
                                coverage = "F"
                                statLine = "F"
                                lfraction = "F"
                                firstlast =  "1st" if int(pos)==first else "last" if int(pos)==last else int (pos)
                                tuplevar = (int(pos), ncb, int(occ), int(code))
                                if tuplevar in hasaaChangeLog:
                                        seq2, lfraction, seq1Length, seq2Length, identity1, identity2, cov1, cov2 = hasaaChangeLog[tuplevar]
                                        identity = identity1
                                        coverage = cov1
                                        statLine = "s1:%s, s2:%s, id1:%s, id2:%s, cv1:%s, cv2:%s" %(seq1Length, seq2Length, identity1, identity2, cov1, cov2)
                                # filestats = "name\tgene\tisoform\tisoformCount\texonCount\tACase\tDcase\tncbCase\taaChangeCases\tfirstOrlast\taaChangeInExon\tfractionLength\tid\tcov\tstatLine\texonList\n"
                                filestats += "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(humanGeneObject[gene].detail,gene,trans.ID, len(humanGeneObject[gene].transcripts), len(trans.exons), acount, dCount, Asps, aaChangeCases,  firstlast, pos+','+ncb+','+occ+','+code, lfraction, identity, coverage, statLine)
    bar.next()
    if nameAddition:
        print (gene, nameAddition)
with open ("filestatsEvenetsDb",'w') as fout:
        fout.write(filestats)


#python db_analysis/eventIdentification.py ../outdir/oldSet_human19_orgs_18/objectsToDB/objectsave_9606__noStructure.pick ../outdir/oldSet_human19_orgs_18/9606_1/derived_data/db_alignments/
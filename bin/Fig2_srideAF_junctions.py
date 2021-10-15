import cPickle as pickle
import modules_common as cm
import modules_analysis as am
import sys
if len(sys.argv)!=4:
        print ("Please give 1 object 2 conditionFile 3 results dir")
        sys.exit()
prog,source_gene_object,gene_conFile,results_dir_csv = sys.argv

fileswriter=open(results_dir_csv+"ss.log","w")
window=3

humanGeneObject=cm.readPickle(source_gene_object)

CONDITION_GENES=cm.readPickle(gene_conFile)
'''
1. for the PI, get the exon types and the junctions and just calculate the basic types like before
2. Screen all the isoforms for pairiwse junctions poulation, and keep adding the junctions, 
    Add their WEF, and junction types te junctions, in sequqnce, popultate the window
'''

def sstype_ends_indi(gene_source,condition, filename, window,fout):
    cons_cons = {}
    alt_alt = {}
    cons_alt = {}
    all_junct = {}
    backgroud = {}
    irrelavant=0
    tupleWritten=0
    desiredlistofgenes=[i for i in gene_source if i in condition]
    for gene in desiredlistofgenes:
        PI=am.givemePI(gene_source[gene])
        if PI:    
            pi_seq=''.join([i.seq for i in PI.exons])
            ss_seq=''
            for i in PI.exons:
                if i.length>0:
                    ssres=i.out_strideseqAF(PI.ID)
                    if ssres:
                        ss_seq+=ssres
                    else:
                        pass

            tuples=[]
            if len(pi_seq)==len(ss_seq):
                for i in range (0 ,len(ss_seq)-((window*2)-1)):
                    tuples+=[(ss_seq[i:i+(window*2)])]

            PI_representative=[]
            for i in range(0,len(PI.exons)-1):
                if PI.exons[i].length>9 and PI.exons[i+1].length>9:
                    PI_representative+=[(PI.exons[i],PI.exons[i+1])]
                else:
                    irrelavant+=1

            if PI_representative:
                for i in PI_representative:
                    i1seq=i[0].out_strideseqAF(PI.ID)
                    i2seq=i[1].out_strideseqAF(PI.ID)
                    if i[0].ID[0]!='R' and i[1].ID[0]!='R': # and i1seq != False and i2seq != False:
                        c_af1=i[0].ID.split(".")[2]
                        c_af2=i[1].ID.split(".")[2]
                        valueJunction=am.ss_junctionWindow(i1seq, i2seq, window)
                        if valueJunction not in all_junct:
                            all_junct[valueJunction]=0
                        all_junct[valueJunction]+=1
                        if c_af1 + c_af2  in ['GG','FF','FG','GF']:
                            if valueJunction not in cons_cons:
                                cons_cons[valueJunction]=0
                            cons_cons[valueJunction]+=1
                            
                        elif c_af1 + c_af2  in ['AA']:
                            if valueJunction not in alt_alt:
                                alt_alt[valueJunction]=0
                            alt_alt[valueJunction]+=1
                        else:
                            if valueJunction not in cons_alt:
                                cons_alt[valueJunction]=0
                            cons_alt[valueJunction]+=1

            if tuples:
                tupleWritten+=1
                #crearing background if the length matches with PI, full isoform
                for i in tuples:
                    sliceRange=int(len(i)/window)
                    temp1=i[-window:]
                    temp2=i[:window]
                    valueJunction=am.ss_junctionWindow(temp1,temp2, window)
                    if valueJunction not in backgroud:
                        backgroud[valueJunction]=0
                    backgroud[valueJunction]+=1
    #print irrelavant
    print ('tuplles background wriottein from %s genes'%tupleWritten)
    cm.csv_writer(backgroud, filename+"background_window%s.csv" % window,fout)
    cm.csv_writer(all_junct, filename+"alljunct_window%s.csv" % window,fout)
    cm.csv_writer(cons_cons, filename+"ConsCons_window%s.csv" % window,fout)
    cm.csv_writer(alt_alt, filename+"AltAlt_window%s.csv" % window,fout)
    cm.csv_writer(cons_alt, filename+"ConsAlt_window%s.csv" % window,fout)
    
def sstype_ends_coservation(gene_source,condition, filename, window, fout):
    total_junctions={}
    fout.write("Helix junction freq\tGene\tExon_pair\tGene\n")
    '''
    get exons in gene, make pairs of them (linear pairs not from exon ids)
    screen transcripts, and for every trancript, get its aa sequence, 
    and fro unique protein coding sequqnces only (going for every transcript, byoassing previuos condition), populate the repr has and its duplicate like below:
        iterate its exons in pairs, 
        if their length is grt 9 (both):
            get their ss seqs
            if they are not retention cases and they are not empty ss's
                get their linear sequqnce tags or close to them ()
                now get their junction sequqnces, both amino acid and SS (4tup), [however not essential anymore now. (will be doing th WEF count of them)]
                if that junction 4tup pair is not there in dup(linear key), add it [bypassed]
                and make similar addition to repr_has version but add there in value (exon id of pairs, amino acid seq and PI tag if any)
    now get to every pair in the linear seq fashion:
        reember a pir is choosen to be diffrenet in either AAseq or SSseq sequqnce, hence being varied
        templis
        if atleast two variations:
            for every variation:
                store their N and C ter Flags as SS + PI flag in list _templis
        sort the templis if it is:
            get PI form as refernec, and value as SS flag
            add SS flag as refernce count
            if its HH EE CC:
                chcek how often this was maintained
    '''
    has_exon_categories_combined={'cons_cons':{},'cons_alt':{},'alt_alt':{}}
    desiredlistofgenes=[i for i in gene_source if i in condition]
    for gene in desiredlistofgenes:
            representative={}            
            for trans in gene_source[gene].transcripts:
                representative= am.junctions_from_different_isoforms_poulator(representative, trans, trans.PI,window, normalorStide='stride')
                #repr_has[exon_sequence]+=[[PI,coding_exons_count,texo[i].ID,texo[i+1].ID,SS_junction]]
            '''
            this will for pairwise exonjunctions and for every such pair will include the ioccurence of the seuqqnces
            # -> need a a way to tell that to either get the 
            '''
            
            for exonSequencePair in representative:
                #getting for every exon
                temp_lis=[]
                #print (exonSequencePair)
                numberofjunctionrpeats=len(representative[exonSequencePair])
                numberofssNotfalse=len([i[-1] for i in representative[exonSequencePair] if i[-1]!='NA'])
                if numberofjunctionrpeats>2 and numberofssNotfalse:
                    # -> means SS should be there, in prediction keep later and part condition equqla to the number of junction repeats 
                    # means junction now has representatives more than 2, not caring about the structure or other definitions are presnet or NA's                    #print ('begin',representative[exonSequencePair])
                    # values per exon in this fashion, PI,coding_exons_count,texo[i].ID,texo[i+1].ID,SS_junction
                    representative[exonSequencePair].sort()
                    interestType=''
                    for temp_pair in representative[exonSequencePair][::-1]:
                        if temp_pair[-1]!='NA':
                            interestType=temp_pair
                            break
                    if interestType:   # its hould always be true 
                        #print (interestType)
                        flagExonType1 = interestType[2].split(".")[2]
                        flagExonType2 = interestType[3].split(".")[2]
                        flagPairKey = 'cons_cons' if flagExonType1+flagExonType2 in ['GG','FF','FG','GF'] else 'alt_alt' if flagExonType1+flagExonType2 =='AA' else 'cons_alt'
                        WEF=round(float(numberofjunctionrpeats)/float(len(gene_source[gene].transcripts)),2)
                        #print (gene,flagPairKey,numberofjunctionrpeats,len(aaseq),WEF)
                        reprSS=interestType[-1]
                        reprSS_count=sum([1 for ex in representative[exonSequencePair] if ex [-1]==reprSS])
                        reprSS_Freq=round(float(reprSS_count)/float(numberofjunctionrpeats),2)
                        temphas={}
                        for tempIteration in representative[exonSequencePair]:
                            tempPIFlag,temp_coding_exons_count,tempex1ID,tempex2ID,tempSS_junction=tempIteration
                            if tempSS_junction not in temphas:
                                temphas[tempSS_junction]=0
                            temphas[tempSS_junction]+=1
                        #print (temphas, 'temphas')
                        for tempjunction in temphas:
                            temphas[tempjunction]=float(temphas[tempjunction]/len(representative[exonSequencePair]))
                        #print (gene,exonSequencePair,interestType,flagExonType1,flagExonType2,flagPairKey,numberofjunctionrpeats,len(gene_source[gene].transcripts),WEF,reprSS,reprSS_count,reprSS_Freq)
                        has_exon_categories_combined[flagPairKey][gene,exonSequencePair,WEF]=(interestType,reprSS,reprSS_count,reprSS_Freq, temphas)
    am.sorting_screening_junctions(has_exon_categories_combined,filename+'WEF_and_exonsJunction.csv',filename+'WEF_and_majorSS_exonsJunction_conservation.csv')
        
#sstype_ends_coservation(humanGeneObject, CONDITION_GENES, results_dir_csv +'ExonsJunctionStrideAF_', window,fileswriter)
sstype_ends_indi(humanGeneObject, CONDITION_GENES, results_dir_csv +'IndividualStrideAF_', window,fileswriter)


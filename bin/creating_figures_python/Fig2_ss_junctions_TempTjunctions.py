import cPickle as pickle
import common.general_modules as gm
import figPanels.modules_analysis as am
import sys
if len(sys.argv)!=4:
        print ("Please give 1 object 2 conditionFile 3 results dir")
        sys.exit()
prog,source_gene_object,gene_conFile,results_dir_csv = sys.argv

fileswriter=open(results_dir_csv+"ss.log","w")
window=3

humanGeneObject=gm.readPickle(source_gene_object)

CONDITION_GENES=gm.readPickle(gene_conFile)
'''
1. for the PI, get the exon types and the junctions and just calculate the basic types like before
2. Screen all the isoforms for pairiwse junctions poulation, and keep adding the junctions, 
    Add their WEF, and junction types te junctions, in sequqnce, popultate the window
'''

def updateHas_increment(has,key):
    if key not in has:
        has[key]=0
    has[key]+=1
    return has


def sstype_ends_indi(gene_source,condition, filename, window,fout):
    cons_cons = {}
    alt_alt = {}
    cons_alt = {}
    all_junct = {}
    backgroud = {}
    irrelavant=0
    tupleWritten=0
    
    desiredlistofgenes=[i for i in gene_source if i in condition]
    desiredlistofgenes=[103]
    for gene in desiredlistofgenes:
        PI=[i for i in gene_source[gene].transcripts if i.PI][0]
        if PI:
            tExonList=am.exonScreenerBetweenTConstitutive(gene_source[gene].exons)
            #print (tExonList)
            pi_seq=''.join([i.seq for i in PI.exons])
            ss_seq=''
            for i in PI.exons:
                if i.length>0:
                    ssres=i.out_secondseq(PI.ID)
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
                print (PI.exons[i].ID, PI.exons[i].length, PI.exons[i+1].ID, PI.exons[i+1].length)
                if am.exonPasser(PI.exons[i], PI.exons[i+1], consecCoding=tExonList):
                    PI_representative+=[(PI.exons[i],PI.exons[i+1])]
                    print ('[passed]')
                else:
                    print ('failed')
                    irrelavant+=1

            print (PI_representative)
            print (irrelavant)
            if PI_representative:
                for i in PI_representative:
                    i1seq=i[0].out_secondseq(PI.ID)
                    i2seq=i[1].out_secondseq(PI.ID)
                    if 1:#i1seq != False and i2seq != False:
                        c_af1=i[0].ID.split(".")[2]
                        c_af2=i[1].ID.split(".")[2]
                        valueJunction=am.ss_junctionWindow(i1seq, i2seq, window)
                        all_junct=updateHas_increment(all_junct,valueJunction)
                        if c_af1 + c_af2  in ['GG','FF','FG','GF']:
                            cons_cons=updateHas_increment(cons_cons,valueJunction)
                        elif c_af1 + c_af2  in ['AA']:
                            alt_alt=updateHas_increment(alt_alt,valueJunction)
                        else:
                            cons_alt=updateHas_increment(cons_alt,valueJunction)
            print (cons_cons)
            sys.exit()
            if tuples:
                tupleWritten+=1
                #crearing background if the length matches with PI, full isoform
                for i in tuples:
                    sliceRange=int(len(i)/window)
                    temp1=i[-window:]
                    temp2=i[:window]
                    valueJunction=am.ss_junctionWindow(temp1,temp2, window)
                    backgroud=updateHas_increment(backgroud,valueJunction)
    #print irrelavant
    am.csv_writer(backgroud, filename+"background_window%s.csv" % window,fout)
    am.csv_writer(all_junct, filename+"alljunct_window%s.csv" % window,fout)
    am.csv_writer(cons_cons, filename+"ConsCons_window%s.csv" % window,fout)
    am.csv_writer(alt_alt, filename+"AltAlt_window%s.csv" % window,fout)
    am.csv_writer(cons_alt, filename+"ConsAlt_window%s.csv" % window,fout)
    
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
            tExonsBtweenFirstConstitutiveToLast=am.exonScreenerBetweenTConstitutive(gene_source[gene].exons)
            representative={}
            for trans in gene_source[gene].transcripts:
                representative= am.junctions_from_different_isoforms_poulator(representative, trans, trans.PI, window, tExonsBtweenFirstConstitutiveToLast)

            #get their WEF's also listed
            for exonSequencePair in representative:
                #getting for every exon
                temp_lis=[]
                #print (exonSequencePair)
                numberofjunctionrpeats=len(representative[exonSequencePair])
                numberofssNotfalse=len([i[-1] for i in representative[exonSequencePair] if i[-1]!='NA'])
                if len(representative[exonSequencePair])>2 and numberofssNotfalse >2:
                    #chnaged suh that we will only go for the scenarios such that 
                    #print ('begin',representative[exonSequencePair])
                    #atleast two combinations:
                    # values per exon in this fashion, PI,coding_exons_count,texo[i].ID,texo[i+1].ID,SS_junction
                    representative[exonSequencePair].sort()
                    interestType=''
                    for temp_pair in representative[exonSequencePair][::-1]:
                        if temp_pair[-1]!='NA':
                            interestType=temp_pair
                            break
                    if interestType:    
                        flagExonType1=interestType[2].split(".")[2]
                        flagExonType2=interestType[3].split(".")[2]
                        flagPairKey='cons_cons' if flagExonType1+flagExonType2 in ['GG','FF','FG','GF'] else 'alt_alt' if flagExonType1+flagExonType2 =='AA' else 'cons_alt'
                        WEF=round(float(len(representative[exonSequencePair]))/float(len(gene_source[gene].transcripts)),2)
                        #print (gene,flagPairKey,len(representative[exonSequencePair]),len(aaseq),WEF)
                        reprSS=interestType[-1]
                        reprSS_count=sum([1 for ex in representative[exonSequencePair] if ex [-1]==reprSS])
                        reprSS_Freq=round(float(reprSS_count)/float(len(representative[exonSequencePair])),3)
                        temphas={}
                        for tempIteration in representative[exonSequencePair]:
                            tempPIFlag,temp_coding_exons_count,tempex1ID,tempex2ID,tempSS_junction=tempIteration
                            if tempSS_junction not in temphas:
                                temphas[tempSS_junction]=0
                            temphas[tempSS_junction]+=1
                        for tempjunction in temphas:
                            temphas[tempjunction]=round(float(temphas[tempjunction])/len(representative[exonSequencePair]),3)
                        
                        #print (gene,exonSequencePair,interestType,flagExonType1,flagExonType2,flagPairKey,len(representative[exonSequencePair]),len(gene_source[gene].transcripts),WEF,reprSS,reprSS_count,reprSS_Freq)
                        has_exon_categories_combined[flagPairKey][gene,exonSequencePair,WEF]=(representative[exonSequencePair],reprSS,reprSS_count,reprSS_Freq, temphas)
    am.sorting_screening_junctions(has_exon_categories_combined,filename+'WEF_and_exonsJunction.csv',filename+'WEF_and_majorSS_exonsJunction_conservation.csv')
        
#sstype_ends_coservation(humanGeneObject, CONDITION_GENES, results_dir_csv +'F2_ExonsJunction', window,fileswriter)
sstype_ends_indi(humanGeneObject, CONDITION_GENES, results_dir_csv +'F2_Individual_', window,fileswriter)





#projectDir/output/derived_data/structure_data/surfaceExposedAF
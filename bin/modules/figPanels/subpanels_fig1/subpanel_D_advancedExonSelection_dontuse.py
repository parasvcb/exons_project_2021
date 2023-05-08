# uncompyle6 version 3.7.5.dev0
# Python bytecode 2.7 (62211)
# Decompiled from: Python 3.7.4 (default, Aug 13 2019, 20:35:49) 
# [GCC 7.3.0]
# Embedded file name: /media/paras/254266be-73cf-476c-b872-c471c9a9dba0/paras/project/protein_splicing/projectDir/bin/subpanels_fig1/subpanel_D.py
# Compiled at: 2021-08-30 15:58:04
import re, os, sys
# import figPanels.modules_common as cm
import common.general_modules as cm
from progress.bar import Bar
import figPanels.modules_analysis as ca

def constitutive_alternate_with_freq(geneob, gene):
    exonMatrix = ca.positionalExonMatrix_forExCharacterization(geneob)
            #{1: ['U', 'F', [2, 0, 2], 0, 1, 0, 0], 2: ['U', 'A', [0, 0, 0], 0, 1, 0, 0], 3: ['T', 'A', [0, 0, 0], 0, 0, 0, 0], 4: ['T', 'F', [0, 1, 0], 0, 0, 0, 0]}

    # print (exonMatrix)
    # temp1, temp2, temp3, temp4, temp5 = ca.giveExonSelectionBasic(exonMatrix,"D")
    temp11, temp22, temp33, temp44, temp55 = ca.giveExonSelectionBasic(exonMatrix,"T")
    # pos_CodingStrict = temp1 + temp11
    # pos_CodingMajorlyConstitutive = temp2 + temp22
    # pos_CodingAlternate = temp3 + temp33
    # pos_CodingAlternateWithSS = temp4 + temp44
    # pos_CodingConstitutive = temp5 + temp55
    
    pos_CodingStrict = temp11
    pos_CodingMajorlyConstitutive = temp22
    pos_CodingAlternate = temp33
    pos_CodingAlternateWithSS = temp44
    pos_CodingConstitutive = temp55
    


    exSpliceSiteChange = []
    exCommon = []
    pos_CodingAlternate_varnoncod = []
    pos_CodingAlternateWithSS_varnoncod = []
    
    def calculatevaluePerPos(pos, geneob):
        def recalcWEF(pos, geneob):
            # print (' in ()', pos, geneob)
            codingVarCount=[]
            for trans in geneob.transcripts:
                for exon in trans.exons:
                    if exon.ID[0]!='R':
                        newpos=int(exon.ID.split('.')[3])
                        # print (newpos,pos==newpos, exon.length>0, pos==newpos and exon.length>0)
                        if pos==newpos and exon.length>0:
                            codingVarCount+=[exon.ID]
                            break
                            # break here is essential as if we 6c1 6c2 6c3 are all presnet in single isoform, then three repeats of theirs will be calucltaed
            # print (len(codingVarCount), len(geneob.transcripts), round(float(len(codingVarCount))/len(geneob.transcripts),3))
            # print (codingVarCount)
            return round(float(len(codingVarCount))/len(geneob.transcripts),3)

        tlis = []
        flag = []
        # print ('POS',pos)
        for i in geneob.exons:
            # if i.ID[0]!='R' and i.ID[0] in ['T','D'] and int(i.ID.split('.')[3])  == pos and i.length >0:
            if i.ID[0]!='R' and int(i.ID.split('.')[3]) == pos:
                if i.ID[0] == 'T' and i.length >0:
                    tlis+=[(i.WEF,pos)]
                else: 
                    flag += [i.ID]
                    # means this position has some var for the exon position whihc are neither coding nor alternate 

        if tlis:
            # print ("YES tlis")
            if len(tlis)!=1:
                value = recalcWEF(pos,geneob)
                # print (value, flag)
                return [value, flag]
            else:
                # print (tlis[0][0], flag)
                return [tlis[0][0], flag]
            
        return [False, False]
    
    # print ('gene ', gene)
    # print (pos_CodingAlternate)
    # print (pos_CodingAlternateWithSS)

    for position in pos_CodingAlternate:
        valueCommon, flagA = calculatevaluePerPos(position, geneob)
        if valueCommon:
            exCommon += [valueCommon]
            if flagA:
                pos_CodingAlternate_varnoncod += [(flagA)]
        else:
            if pos_CodingAlternate:
                print (gene, "Alt", pos_CodingAlternate, "lacks evena  singel entity")
                # sys.exit()
    
    for position in pos_CodingAlternateWithSS:
        valueSS, flagAss= calculatevaluePerPos(position, geneob)
        if valueSS:
            exSpliceSiteChange += [valueSS]
            if flagAss:
                pos_CodingAlternateWithSS_varnoncod += [(flagAss)]
        else:
            if pos_CodingAlternateWithSS:
                print (gene, "AltSS", pos_CodingAlternateWithSS, "lacks evena  singel entity")
                # sys.exit()
    # flag A and flag Ass are the exonic loci positions which does have the positions variatsions which are not (T or coding) 
    return (exCommon, exSpliceSiteChange, pos_CodingAlternate_varnoncod, pos_CodingAlternateWithSS_varnoncod)


def alternate_WEF(has, condition, res_dir, fout):
    """
    Here this will get the information for the every exon in two different ways, instead of calculating all the instances.
    we would like to calcute the fraction a) only for the latrenate exons with nos plice site vanges
    b) combuned instances of the alt exons toghether in all the transcripts.
    # this is idffernt because we are now focusiing on only teh fraction of exons which remaiins coding all the times
    """
    has_bins_AltOnly = {(0, 0.2): 0, (0.201, 0.4): 0, (0.401, 0.6): 0, (0.601, 0.8): 0, (0.801, 1.001): 0}
    has_bins_AltWithSS = {(0, 0.2): 0, (0.201, 0.4): 0, (0.401, 0.6): 0, (0.601, 0.8): 0, (0.801, 1.001): 0}
    valsAltOnly = []
    ValsAltWithSSChange = []
    bar = Bar('Processing genes D:', max=len(condition))
    genesAlt = 0
    genesAltwSS = 0
    flagA_noncodvarcases =[]
    flagAss_noncodvarcases =[]
    

    for gene in has:
        if gene in condition:
            commonvals, valscumulativefromSSchanges, flagA, flagAss = constitutive_alternate_with_freq(has[gene],gene)
            valsAltOnly += commonvals
            ValsAltWithSSChange += valscumulativefromSSchanges
            if commonvals:
                genesAlt +=1
            if valscumulativefromSSchanges:
                genesAltwSS +=1
            if flagA:
                flagA_noncodvarcases += [(gene,flagA)]
            if flagAss:
                flagAss_noncodvarcases += [(gene,flagAss)]                
            bar.next()
    bar.finish()
    flagA_noncodvarcases_exonic_entities = sum([1 for i in flagA_noncodvarcases for j in i[1]])
    flagAss_noncodvarcases_exonic_entities = sum([1 for i in flagAss_noncodvarcases for j in i[1]])

    print ("AlyOnlyFraction: %s exonic entities, from %s genes, sample vals %s" %(len(valsAltOnly), genesAlt , valsAltOnly[:5]))
    print ("*AlyOnlyFraction where exonic loci had variations: %s exonic entities, from %s genes" %(flagA_noncodvarcases_exonic_entities, len(flagA_noncodvarcases)))
    print ("AlySSFraction: %s exonic entities, from %s genes, sample vals %s" %(len(ValsAltWithSSChange), genesAltwSS , ValsAltWithSSChange[:5]))
    print ("*AltSSFraction where exonic loci had variations: %s exonic entities, from %s genes" %(flagAss_noncodvarcases_exonic_entities, len(flagAss_noncodvarcases)))

    with open(os.path.join(res_dir, 'Alternate_exons_and_SS_where_pos_has_non_cod_var.tsv'), 'w') as fin:
        fin.write("Category\tgene\tpos\n")
        for i in flagA_noncodvarcases:
            for j in i[1]:
                fin.write("A\t%s\t%s"%(i,j))
        for i in flagAss_noncodvarcases:
            for j in i[1]:
                fin.write("Ass\t%s\t%s"%(i,j))
        
    for wef in valsAltOnly:
        for key in has_bins_AltOnly:
            if key[0] <= wef <= key[1]:
                has_bins_AltOnly[key] += 1
                break

    for wef in ValsAltWithSSChange:
        for key in has_bins_AltWithSS:
            if key[0] <= wef <= key[1]:
                has_bins_AltWithSS[key] += 1
                break

    fout.write('\n\n=============>\n')
    fout.write('Inclusion Frequency of Alternate exons\n')
    with open(os.path.join(res_dir, 'Alternate_exons_WEF.csv'), 'w') as fin:
        fout.write('\tInclusion_Range\tType\tTotal\tFrequency\n')
        fin.write('Inclusion_Range,Type,Total,Frequency\n')
        suAlTonly = sum(has_bins_AltOnly.values())
        suAlTwithSS = sum(has_bins_AltWithSS.values())
        keys = list(has_bins_AltOnly.keys())
        keys.sort()
        for i in keys:
            fin.write('%s,%s,%s,%s\n' % (('-').join(map(str, i)), 'AltOnly', has_bins_AltOnly[i], cm.div_fact(has_bins_AltOnly[i], suAlTonly)))
            fout.write('\n\t%s\t%s\t%s\t%s' % (('-').join(map(str, i)), 'AltOnly', has_bins_AltOnly[i], cm.div_fact(has_bins_AltOnly[i], suAlTonly)))

        for i in keys:
            fin.write('%s,%s,%s,%s\n' % (('-').join(map(str, i)), 'AltWithSS', has_bins_AltWithSS[i], cm.div_fact(has_bins_AltWithSS[i], suAlTwithSS)))
            fout.write('\n\t%s\t%s\t%s\t%s' % (('-').join(map(str, i)), 'AltWithSS', has_bins_AltWithSS[i], cm.div_fact(has_bins_AltWithSS[i], suAlTwithSS)))

    with open(os.path.join(res_dir, 'Alternate_exons_WEF_raw.csv'), 'w') as fin:
        fin.write('Type,Value\n')
        for wef in valsAltOnly:
            fin.write ("%s,%s\n"%("AltOnly",wef))
        for wef in ValsAltWithSSChange:
            fin.write ("%s,%s\n"%("AltWithSS",wef))

   
    fout.write('\nCSV_File:%s' % res_dir + 'RAW_Alternate_exons_inclusion_freq.csv')
# okay decompiling /media/paras/254266be-73cf-476c-b872-c471c9a9dba0/paras/project/protein_splicing/projectDir/bin/modules/figPanels/subpanels_fig1/subpanel_D.pyc

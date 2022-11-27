# uncompyle6 version 3.7.5.dev0
# Python bytecode 2.7 (62211)
# Decompiled from: Python 3.7.4 (default, Aug 13 2019, 20:35:49) 
# [GCC 7.3.0]
# Embedded file name: /media/paras/254266be-73cf-476c-b872-c471c9a9dba0/paras/project/protein_splicing/projectDir/bin/subpanels_fig1/subpanel_D.py
# Compiled at: 2021-08-30 15:58:04
import re, os
# import figPanels.modules_common as cm
import common.general_modules as cm

import figPanels.modules_analysis as ca

def constitutive_alternate_with_freq(gene):
    tempExAltCoding = [ i for i in gene.exons if re.match('^[T]\\.[1-n]+\\.A', i.ID) ]
    res1, res2 = ca.split_cases(tempExAltCoding, out='exons')
    UTRAlternate = res1
    UTRAlternateWithSS = res2
    mat = []
    for i in gene.transcripts:
        mat += [i.exons]

    exCommon = []
    exSpliceSiteChange = []
    for i in UTRAlternate:
        exCommon += [i.WEF]

    for i in UTRAlternateWithSS:
        count = len([ 1 for transcriptExonList in mat for j in i if j in transcriptExonList ])
        exSpliceSiteChange += [round(float(count) / len(mat), 3)]

    return (
     exCommon, exSpliceSiteChange)


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
    for gene in has:
        if gene in condition:
            commonvals, valscumulativefromSSchanges = constitutive_alternate_with_freq(has[gene])
            valsAltOnly += commonvals
            ValsAltWithSSChange += valscumulativefromSSchanges

    print valsAltOnly[:5]
    print ValsAltWithSSChange[:5]
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
    with open(os.path.join(res_dir, 'Alternate_exons_WEF.csv'), 'w') as (fin):
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

    fout.write('\nCSV_File:%s' % res_dir + 'RAW_Alternate_exons_inclusion_freq.csv')
# okay decompiling /media/paras/254266be-73cf-476c-b872-c471c9a9dba0/paras/project/protein_splicing/projectDir/bin/modules/figPanels/subpanels_fig1/subpanel_D.pyc

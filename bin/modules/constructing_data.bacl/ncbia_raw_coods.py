import cPickle as pickle


def aachangecases(key, exon_pack, cod_status_pre, coding_cood, leftCoodExoni, rightCoodExoni, tupRangeExoni, exons_aa_change, parent_or_alt_exons, parent_exons_eve):

    temp_l1 = {i[0]: (i[1], i[2]) for i in key}
    element = (exon_pack[0], temp_l1[exon_pack[0]]
               [0], temp_l1[exon_pack[0]][1])
    # print element
    # print parent_or_alt_exons[element]

    event_occurence = parent_exons_eve[element][3]
    cod_status = event_occurence if cod_status_pre > 0 else cod_status_pre
    newidt = parent_or_alt_exons[element][0][:]
    newidt[1] = cod_status
    exons_aa_change[exon_pack] = [newidt, coding_cood,
                                  leftCoodExoni, rightCoodExoni, tupRangeExoni]
    parent_exons_eve[element][3] += 1
    parent_exons_eve[exon_pack] = [
        {}, {}, {}, 2 if cod_status_pre == 1 else 1, {}]
    return exons_aa_change, parent_exons_eve


def hidden_events1(exon_pack,  coding_cood, leftCoodExoni, rightCoodExoni, tupRangeExoni, parent_exons, alt_exons, parent_exons_eve, exons_aa_change):
    cod_status_pre = -2 if exon_pack[2] == "N" else 0 if exon_pack[2] == "M" else - \
        1 if exon_pack[2] == "C" and exon_pack[1] == '' else 1
    key = parent_exons.keys()
    key.sort()
    if exon_pack[0] in [i[0] for i in key]:
        exons_aa_change, parent_exons_eve = aachangecases(
            key, exon_pack, cod_status_pre, coding_cood, leftCoodExoni, rightCoodExoni, tupRangeExoni, exons_aa_change, parent_exons, parent_exons_eve)
        return parent_exons, alt_exons, parent_exons_eve, exons_aa_change

    key_alt = alt_exons.keys()
    key_alt.sort()
    if exon_pack[0] in [i[0] for i in key_alt]:
        # print "true"
        exons_aa_change, parent_exons_eve = aachangecases(
            key_alt,  exon_pack, cod_status_pre, coding_cood, leftCoodExoni, rightCoodExoni, tupRangeExoni, exons_aa_change, alt_exons, parent_exons_eve)
        return parent_exons, alt_exons, parent_exons_eve, exons_aa_change

    overlap = [i for i in key if set(tupRangeExoni) & set(parent_exons[i][4])]

    if len(overlap) == 1:
        parent_exon1 = overlap[0]
        parent_exon_id = parent_exons[parent_exon1][0][3]
        parent_exon_cod_type = parent_exons[parent_exon1][0][0]
        nside = leftCoodExoni - parent_exon1[0][0]
        cside = rightCoodExoni - parent_exon1[0][1]

        if cside == 0 and abs(nside):
            event_occurence = len(parent_exons_eve[parent_exon1][0])
            newid = [parent_exon_cod_type, cod_status_pre,
                     'A', parent_exon_id, 'n', event_occurence]
            parent_exons_eve[parent_exon1][0][nside] = newid

        elif nside == 0 and abs(cside):
            event_occurence = len(parent_exons_eve[parent_exon1][1])
            newid = [parent_exon_cod_type, cod_status_pre,
                     'A', parent_exon_id, 'c', event_occurence]
            parent_exons_eve[parent_exon1][1][cside] = newid

        elif abs(cside) and abs(nside):
            event_occurence = len(parent_exons_eve[parent_exon1][2])
            newid = [parent_exon_cod_type, cod_status_pre,
                     'A', parent_exon_id, 'b', event_occurence]
            parent_exons_eve[parent_exon1][2][(nside, cside)] = newid

        alt_exons[exon_pack] = [newid, coding_cood,
                                leftCoodExoni, rightCoodExoni, tupRangeExoni]
        parent_exons_eve[exon_pack] = [
            {}, {}, {}, 2 if cod_status_pre == 1 else 1, {}]
    else:
        # print "tru"
        overlap.sort()
        '''
        for i in [[i11[:3], alt_exons[i11][2], alt_exons[i11]
                   [3], alt_exons[i11][0]] for i11 in alt_exons]:
            print i, "alt"
            print "out"
            print parent_exons_eve
        '''
        hasleftb = {}
        hasrightb = {}
        for ti1 in parent_exons:
            # print ti1, "ti1"
            # print parent_exons_eve[ti1]
            for tkeys in parent_exons_eve[ti1][2].keys():

                hasleftb[tkeys[0]] = parent_exons_eve[ti1][2][tkeys]
                hasrightb[tkeys[1]] = parent_exons_eve[ti1][2][tkeys]

        # print "error"
        # print parent
        first_parent_exon = overlap[0]
        last_parent_exon = overlap[-1]
        '''
        print hasleftb, "left"
        print hasrightb, "right"
        print exon_pack, "exon_pack"
        print parent_exons[first_parent_exon][0], "first"
        print parent_exons[last_parent_exon][0], "last"
        print nside, cside
        '''
        nside = leftCoodExoni - first_parent_exon[0][0]
        cside = rightCoodExoni - last_parent_exon[0][1]
        newid = []

        if nside in parent_exons_eve[first_parent_exon][0]:
            newid += [parent_exons_eve[first_parent_exon][0][nside]]
            # print "t2"
        elif nside in hasleftb:
            newid += [hasleftb[nside]]
            # print "t4"
        else:
            # print "tru2"
            event_occurence = len(parent_exons_eve[first_parent_exon][0])
            parent_exon_id = parent_exons[first_parent_exon][0][3]
            parent_exon_cod_type = parent_exons[first_parent_exon][0][0]
            parent_exon_nature = parent_exons[first_parent_exon][0][2]
            parent_exons_status = parent_exons[first_parent_exon][0][1]
            tnewidn = [parent_exon_cod_type, parent_exons_status,
                       parent_exon_nature, parent_exon_id, 'n', event_occurence]
            # print tnewidn
            newid += [tnewidn]
            parent_exons_eve[first_parent_exon][0][nside] = tnewidn
        ''' c side'''
        newidinner = str(len(parent_exons_eve[last_parent_exon][4]))
        if cside in parent_exons_eve[last_parent_exon][1]:
            # print "t21"
            newid += [parent_exons_eve[last_parent_exon][1][cside]]
        elif cside in hasrightb:
            # print "t22"
            newid += [hasrightb[cside]]
        else:
            # print "t23"
            event_occurence = len(parent_exons_eve[last_parent_exon][1])
            parent_exon_id = parent_exons[last_parent_exon][0][3]
            parent_exon_cod_type = parent_exons[last_parent_exon][0][0]
            parent_exon_nature = parent_exons[last_parent_exon][0][2]
            parent_exons_status = parent_exons[last_parent_exon][0][1]
            tnewidc = [parent_exon_cod_type, parent_exons_status,
                       parent_exon_nature, parent_exon_id, 'c', event_occurence]
            parent_exons_eve[last_parent_exon][1][cside] = tnewidc
            newid += [tnewidc]
        newid = ['R', cod_status_pre, ".".join(
            map(str, newid[0])), newidinner, ".".join(map(str, newid[1]))]
        # print newid, "newid()"
        alt_exons[exon_pack] = [newid, coding_cood,
                                leftCoodExoni, rightCoodExoni, tupRangeExoni]
        parent_exons_eve[exon_pack] = [
            {}, {}, {}, 2 if cod_status_pre == 1 else 1, {}]

    return parent_exons, alt_exons, parent_exons_eve, exons_aa_change


def template_feeder(parent_list, child_list, mega_set=None):
    flag = 1 if mega_set is not None else 0
    for ind, val in enumerate(child_list):
        rawCood = val[1]
        codingCood = val[2]
        aaSeq = val[0]
        exStatus = val[3]
        tupleExonKey = (rawCood, aaSeq, exStatus)
        # print rawCood
        parent_list += [(tupleExonKey, codingCood, rawCood[0],
                         rawCood[1], tuple(range(rawCood[0], rawCood[1]+1)))]
        if flag:
            mega_set |= set(range(rawCood[0], rawCood[1]+1))
    if flag:
        return parent_list, mega_set
    else:
        return parent_list


def insertion_filter(insertion_cases, extension_cases, parent_exons_template):
    '''
    for i in [i[:3] for i in insertion_cases]:
        print i, "old_insertion"
    for i in parent_exons_template:
        print i[:3], "old"
    '''
    temp_lis_temp = [[len(i[4]), i] for i in insertion_cases]
    temp_lis_temp.sort()
    insertion_cases = [i[1] for i in temp_lis_temp]

    while insertion_cases:
        temp_lis = [(len(insertion_cases[0][4]), insertion_cases[0])]
        ele_remove = [insertion_cases[0]]
        for j in range(1, len(insertion_cases)):
            if set(insertion_cases[j][4]) & set(insertion_cases[0][4]):
                temp_lis += [(len(insertion_cases[j][4]), insertion_cases[j])]
                ele_remove += [insertion_cases[j]]
        insertion_cases = [i for i in insertion_cases if i not in ele_remove]

        element_of_interest = False
        if len(temp_lis) > 1:
            temp_lis.sort()
            if temp_lis[-1][0] > 29:
                for i in temp_lis:
                    if i[0] > 29:
                        element_of_interest = i
                        break
            else:
                element_of_interest = temp_lis[-1]
            for i in temp_lis:
                if i != element_of_interest:
                    extension_cases += [i[1]]
                else:
                    parent_exons_template += [i[1]]
        else:
            parent_exons_template += [temp_lis[0][1]]
        '''
        for i in [i[:3] for i in insertion_cases]:
            print i, "new_insertion"
        for i in temp_lis:
            print i[0], i[1][:3], "temp_lis"
        for i in parent_exons_template:
            print i[:3], "pet"
        '''
    return insertion_cases, extension_cases, parent_exons_template


def exon_manipulator(raw_var, parent_exons, parent_exons_eve, isoforms):
    #print "\nrawvar:\n%s\n" % raw_var
    # print "\nparent exons: \n %s\n" % parent_exons
    #print "\nparent_Exons_eve:\n%s\n" % parent_exons_eve
    #print "in exon manipulator"
    #print isoforms
    #for i in set(raw_var):
        #print raw_var.count(i), "count", i
    #print "\n\n"
    # print "parent exons,,", parent_exons
    #for i in [i[0] for i in parent_exons]:
        #print i, "i"
    #print "parent_exons_eve:", parent_exons_eve
    for exon in parent_exons:
        consAlt = 'G' if raw_var.count(exon[0]) == isoforms else 'A'
        parent_exons[exon][0][2] = consAlt
        #print consAlt
        #print exon
        #print parent_exons[exon]
        #print raw_var
        # break
    template_temp = []
    template_temp_repeat_count = []
    #print "templating"
    for exon in parent_exons:
        valtemp = parent_exons[exon]
        if valtemp[0][2] == 'A' and exon[0] not in template_temp_repeat_count:
            #print exon[0], valtemp[0], raw_var.count(exon[0])
            for ji in range(0, raw_var.count(exon[0])):
                template_temp += [valtemp[0][3]]
            template_temp_repeat_count += [exon[0]]
    #template_temp = [parent_exons[exon][0][3] for exon in parent_exons if parent_exons[exon][0][2] == 'A']
    #print "parentexonsi", [parent_exons[i][0] for i in parent_exons]
    #print "template_temp ", template_temp
    overlap_t = {}
    for val in template_temp:
        if val not in overlap_t:
            overlap_t[val] = 1
        else:
            overlap_t[val] += 1
    #print overlap_t, "ovt"
    '''
    adding category F
    '''
    for exon in parent_exons:
        first_digit = parent_exons[exon][0][3]
        if parent_exons[exon][0][2] == 'A':
            #print parent_exons[exon][0], overlap_t[first_digit], isoforms
            if overlap_t[first_digit] == isoforms:
                # print "inthi"
                parent_exons[exon][0][2] = 'F'
                # print parent_exons[exon][0]
    has = {tuple(parent_exons[i][0][3:]): [] for i in parent_exons}
    for i in parent_exons:
        if tuple(parent_exons[i][0][3:]) in has:
            has[tuple(parent_exons[i][0][3:])] += [(parent_exons[i][0][3],
                                                    parent_exons[i][0][4], parent_exons[i][0][5], parent_exons[i][0][1])]
    has1 = {}
    # print has
    for i in has:
        t = []
        for j in has[i]:
            t += [j[-1]]
        if len(set(t)) == 1:
            if t[0] == -2:
                has1[i] = 'U'
            elif t[0] == 0:
                has1[i] = 'M'
            else:
                has1[i] = 'T'
        else:
            if -2 in t:
                has1[i] = 'D'
            else:
                has1[i] = 'T'
    for i in parent_exons:
        val_upd = has1[tuple(parent_exons[i][0][3:])]
        parent_exons[i][0][0] = val_upd
    for i in parent_exons_eve:
        for ind in [0, 1, 2, 4]:
            for j in parent_exons_eve[i][ind]:
                val_upd = has1[tuple(parent_exons_eve[i][ind][j][3:])]
                parent_exons_eve[i][ind][j][0] = val_upd

    return parent_exons, parent_exons_eve


def ncbia(var_list, pi, exons_add):

    with open(exons_add+"%s" % (pi)) as fin:
        dat = pickle.load(fin)

    parent_exons_template = []
    mega_Set = set()
    # print dat, "PI template_raw"
    parent_exons_template, mega_Set = template_feeder(
        parent_list=parent_exons_template, child_list=dat, mega_set=mega_Set)
    raw_var = [i[0][0] for i in parent_exons_template]
    # print len(raw_var), "lwen1"
    '''
    for i in [i[:3] for i in parent_exons_template]:
        print i, "parent_Exons_template"
    '''
    sorted_vars = var_list[:]
    sorted_exons_o = []
    isoforms = len(var_list)
    for var in sorted_vars:
        # iterating the variable now for the gene group
        if var != pi:
            with open(exons_add+"%s" % var) as fin:
                exons = pickle.load(fin)
            sorted_exons_o = template_feeder(
                parent_list=sorted_exons_o, child_list=exons)

    sorted_exons = list(set(sorted_exons_o))
    # sorted_exons.sort()
    # print sorted_exons,"sortedexons"
    '''
    for i in [i[:3] for i in sorted_exons]:
        print i, "sorted_exons_o"
    '''
    insertion_cases = []
    extension_cases = []
    for i in sorted_exons:
        if set(i[4]) & mega_Set:
            extension_cases += [i]
        else:
            insertion_cases += [i]
    '''
    for i in [i[:3] for i in insertion_cases]:
        print i, "insertion_cases"
    for i in [i[:3] for i in extension_cases]:
        print i, "extensioncases"
    # print [i[:-1] for i in extension_cases],"exten"
    '''
    insertion_cases.sort()
    insertion_cases, extension_cases, parent_exons_template\
        = insertion_filter(insertion_cases, extension_cases, parent_exons_template)
    # print [i[:-1] for i in parent_exons_template], "parent2"
    # print [i[:-1] for i in extension_cases], "extension2"
    '''
    for i in [i[:3] for i in parent_exons_template]:
        print i, "parent_exons_template"
    for i in [i[:3] for i in extension_cases]:
        print i, "extensioncases_ref"
    '''
    parent_exons_template.sort()
    parent_exons = {}
    alt_exons = {}
    parent_exons_eve = {}
    exons_aa_change = {}
    for ind, val in enumerate(parent_exons_template):
        idexon = ind+1
        unique_key = val[0]
        cod_status = -2 if unique_key[2] == "N" else 0 if unique_key[2] == "M" else - \
            1 if unique_key[2] == "C" and unique_key[1] == '' else 1
        codCood = val[1]
        leftpos = val[2]
        rightpos = val[3]
        range_tup = val[4]
        uniqueIdLis = ['E', cod_status, '', idexon, 0, 0]
        parent_exons[unique_key] = [uniqueIdLis,
                                    codCood, leftpos, rightpos, range_tup]
        nidentity = leftpos-leftpos
        cidentity = rightpos-rightpos
        parent_exons_eve[unique_key] = [{nidentity: uniqueIdLis}, {cidentity: uniqueIdLis},
                                        {(nidentity, cidentity): uniqueIdLis}, 2 if cod_status == 1 else 1, {}]
    intron_retention_cases = []
    temp_list_refine = []
    for i in extension_cases:
        ic = 0
        for j in parent_exons_template:
            if set(i[4]) & set(j[4]):
                ic += 1
        temp_list_refine += [(ic, i)]
    extension_cases = []
    for i in temp_list_refine:
        # print i[0], i[1][:3], "temp_ref"
        if i[0] == 1:
            extension_cases += [i[1]]
        else:
            intron_retention_cases += [i[1]]
    '''
    for i in [i[:3] for i in intron_retention_cases]:
        print i, "intron_ret"
    for i in [i[:3] for i in extension_cases]:
        print i, "extensioncases_ref3"
    for i in parent_exons:
        print i, parent_exons[i][0], "indparent"
    '''
    parent_exons_upd = {}
    extension_cases.sort()
    for i in extension_cases:
        # print i[:3], "extensioncases"
        # i[0] is the exonic unique key
        if not (i[0] in parent_exons or i[0] in alt_exons or i[0] in exons_aa_change):
            parent_exons, alt_exons, parent_exons_eve, exons_aa_change = hidden_events1(
                i[0], i[1], i[2], i[3], i[4], parent_exons, alt_exons, parent_exons_eve, exons_aa_change)

    parent_exons_upd.update(parent_exons)
    parent_exons_upd.update(alt_exons)
    parent_exons_upd.update(exons_aa_change)
    raw_var += [i[0][0] for i in sorted_exons_o]
    # print len(raw_var), "len2"
    parent_exons_upd, parent_exons_eve = exon_manipulator(
        raw_var, parent_exons_upd, parent_exons_eve, isoforms)
    for i in parent_exons:
        parent_exons[i][0] = parent_exons_upd[i][0]
    for i in alt_exons:
        alt_exons[i][0] = parent_exons_upd[i][0]
    for i in exons_aa_change:
        exons_aa_change[i][0] = parent_exons_upd[i][0]
    '''
    keynew = parent_exons_upd.keys()
    keynew.sort()
    for i in keynew:
        print i, parent_exons_upd[i][:3], "this"
    '''
    intron_retention_cases.sort()
    for i in intron_retention_cases:
        # if i[0][0] == (4817, 6739):
        # print i[:4], "int"
        if not (i[0] in parent_exons or i[0] in alt_exons or i[0] in exons_aa_change):
            parent_exons, alt_exons, parent_exons_eve, exons_aa_change = hidden_events1(
                i[0], i[1], i[2], i[3], i[4], parent_exons, alt_exons, parent_exons_eve, exons_aa_change)
    parent_exons.update(alt_exons)
    parent_exons.update(exons_aa_change)

    '''
    keynew = parent_exons.keys()
    keynew.sort()
    for i in keynew:
        print i, parent_exons[i][:3], "this"
    '''
    parent_exons_id_format = {}
    for i in parent_exons:
        if parent_exons[i][0][0] == "R":
            parent_exons_id_format[i] = [
                ":".join(map(str, parent_exons[i][0])), parent_exons[i][1]]
        else:
            parent_exons_id_format[i] = [
                ".".join(map(str, parent_exons[i][0])), parent_exons[i][1]]
    '''
    keynew = parent_exons_id_format.keys()
    keynew.sort()
    for i in keynew:
        print i, parent_exons_id_format[i][:3], "this"
    '''
    return parent_exons_id_format

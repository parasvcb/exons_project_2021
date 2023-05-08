import cPickle as pickle
import sys, os

tupleSpanIndex = 2


def sub_aachangecases(key, exon_pack, cod_status_pre, coding_cood, tupRangeExoni, exons_aa_change, parent_or_alt_exons, parent_exons_eve):
    # 2-5 are args of exon in question
    # key are list of parent exons (3 ele tuples)
    temp_exonListParentOrAlt = {i[0]: (i[1], i[2]) for i in key}
    parentImitation = (exon_pack[0], temp_exonListParentOrAlt[exon_pack[0]]
                       [0], temp_exonListParentOrAlt[exon_pack[0]][1])
    # rawCoods tuple (st,ed)(curr exon), aaseq(par exon), codingStatus (par exon)
    # pseudo parent form reconstruction

    # format: parent_exons_eve[unique_key] = [{nidentity: uniqueIdLis}, {cidentity: uniqueIdLis}, {bidentity: uniqueIdLis}, 2 if cod_status == 1 else 1, {}]

    event_occurence = parent_exons_eve[parentImitation][3]
    # it can be 2 if parent was coding , else 1
    cod_status = event_occurence if cod_status_pre > 0 else cod_status_pre
    # for new exon, if its curr cod status is coding with aaseq, then transfer value from parent (transferParentStatus (>1 any)), else keep its form (-2,0,-1)
    # why to transfer, because i have already increremented it +1 during inititalization (events keep record)
    # chnaging forms is only nenecessary it its comong with aaseq
    newidt = parent_or_alt_exons[parentImitation][0][:]
    # 6 letter list copied from the parent
    newidt[1] = cod_status
    exons_aa_change[exon_pack] = [newidt, coding_cood, tupRangeExoni]
    parent_exons_eve[parentImitation][3] += 1
    # parent_exons_eve[exon_pack] = [{}, {}, {}, 2 if cod_status_pre == 1 else 1, {}]
    return exons_aa_change, parent_exons_eve


def cod_status_primitive(value, aaseq):
    return -2 if value == "N" \
        else 0 if value == "M" \
        else -1 if value == "C" and aaseq == '' \
        else 1


def hidden_events1(exon_pack, coding_cood, tupRangeExoni, parent_exons, alt_exons, parent_exons_eve, exons_aa_change):
    # first 3 args are for incoming exon
    try:
        leftCoodExoni = exon_pack[0][0]
        rightCoodExoni = exon_pack[0][1]
        cod_status_pre = cod_status_primitive(exon_pack[2], exon_pack[1])
        '''
        This routine checks if the incomingExon entry has identical raw span coods, if yes, then must be a change in the aa, part 1, seond part looks at if it shares
        '''
        parExonsSign = parent_exons.keys()
        parExonsSign.sort()
        altExonSign = alt_exons.keys()
        altExonSign.sort()

        # print ("(),1")
        # TODO: certain notions arent accepted inthe python 3 format, like above
        # ################################### PART1 ################################################################################################################
        if exon_pack[0] in [i[0] for i in parExonsSign]:
            exons_aa_change, parent_exons_eve = sub_aachangecases(
                parExonsSign, exon_pack, cod_status_pre, coding_cood, tupRangeExoni, exons_aa_change, parent_exons, parent_exons_eve)
            return parent_exons, alt_exons, parent_exons_eve, exons_aa_change

        # print ("(),2")
        if exon_pack[0] in [i[0] for i in altExonSign]:
            exons_aa_change, parent_exons_eve = sub_aachangecases(
                altExonSign, exon_pack, cod_status_pre, coding_cood, tupRangeExoni, exons_aa_change, alt_exons, parent_exons_eve)
            return parent_exons, alt_exons, parent_exons_eve, exons_aa_change
            # removed entry of aa change in event trackers (parent_exons_eve), i think that wasnt needed
        # print ("(),3")
        # ###################################################     PART 1 ends here     #############################################################################

        # print (parExonsSign)
        # print (tupRangeExoni)
        # print ("0", parent_exons[(
        #     (6629, 6938), 'MSFKDPPTLQQLARRSLLKDEALTISALPNLPVQLFPPLFKDAFTSRQRKILSLMVATWPFPVLPVGALCGIDHLETLKAVLDGLDLLMSQKDRPS', 'C')])
        # print ("1", parent_exons[(
        #     (6629, 6938), 'MSFKDPPTLQQLARRSLLKDEALTISALPNLPVQLFPPLFKDAFTSRQRKILSLMVATWPFPVLPVGALCGIDHLETLKAVLDGLDLLMSQKDRPS', 'C')][tupleSpanIndex])
        # print ("2", set(parent_exons[(
        #     (6629, 6938), 'MSFKDPPTLQQLARRSLLKDEALTISALPNLPVQLFPPLFKDAFTSRQRKILSLMVATWPFPVLPVGALCGIDHLETLKAVLDGLDLLMSQKDRPS', 'C')][tupleSpanIndex]))

        overlap = [i for i in parExonsSign if set(
            tupRangeExoni) & set(parent_exons[i][tupleSpanIndex])]
        # print ("(),4")
        if len(overlap) == 1:
            parent_exon_curr = overlap[0]
            parent_exon_rank = parent_exons[parent_exon_curr][0][3]
            # FORMAT: [uniqueIdLis, codCood, leftpos, rightpos, range_tup]
            #  uniqueIdLis = ['E', cod_status, '', idexon, 0, 0]
            parent_exon_cod_type = parent_exons[parent_exon_curr][0][0]
            # its defaults to 'E'
            nside = leftCoodExoni - parent_exon_curr[0][0]
            cside = rightCoodExoni - parent_exon_curr[0][1]
            # print ("(),5")
            if cside == 0 and abs(nside):
                event_occurence = len(parent_exons_eve[parent_exon_curr][0])
                # n,c,b identities are hashes, new events will be on basis of length update.
                newid = [parent_exon_cod_type, cod_status_pre,
                         'A', parent_exon_rank, 'n', event_occurence]
                parent_exons_eve[parent_exon_curr][0][nside] = newid
                # why occurence of same length is not being accounted for is, beacuse on 2nd repeat of same form, it woould have neen alreay dealt (we are dealing set, or it will go to aa chnage case)
                # print ("(),6")
            elif nside == 0 and abs(cside):
                event_occurence = len(parent_exons_eve[parent_exon_curr][1])
                newid = [parent_exon_cod_type, cod_status_pre,
                         'A', parent_exon_rank, 'c', event_occurence]
                parent_exons_eve[parent_exon_curr][1][cside] = newid
                # print ("(),7")
            elif abs(cside) and abs(nside):
                event_occurence = len(parent_exons_eve[parent_exon_curr][2])
                newid = [parent_exon_cod_type, cod_status_pre,
                         'A', parent_exon_rank, 'b', event_occurence]
                parent_exons_eve[parent_exon_curr][2][(nside, cside)] = newid
            # print ("(),8")
            alt_exons[exon_pack] = [newid, coding_cood, tupRangeExoni]
            parent_exons_eve[exon_pack] = [
                {}, {}, {}, 2 if cod_status_pre == 1 else 1, {}]
        else:
            # print ("overlap", overlap)
            overlap.sort()
            first_parent_exon = overlap[0]
            last_parent_exon = overlap[-1]

            hasLeftBformer = {}
            hasRightBlatter = {}
            for bspan in parent_exons_eve[first_parent_exon][2].keys():
                # this is b case
                leftChange = bspan[0]
                if leftChange not in hasLeftBformer:
                    # nside change as key, and complete ID as value, ID is []
                    hasLeftBformer[leftChange] = []
                hasLeftBformer[leftChange] += [parent_exons_eve[first_parent_exon][2][bspan]]
            for bspan in parent_exons_eve[last_parent_exon][2].keys():
                # this is b case
                rightChange = bspan[1]
                if rightChange not in hasRightBlatter:
                    # cside change as key, and complete ID as value, ID is []
                    hasRightBlatter[rightChange] = []
                hasRightBlatter[rightChange] += [
                    parent_exons_eve[last_parent_exon][2][bspan]]

            # b changes can be many, (1-10, 1-11,), LOL bug
            # technically retenetion may begin from any of those changed values, but precise selection of any of them is tricky and hence not specifed, (first come, first selected basis)

            nside = leftCoodExoni - first_parent_exon[0][0]
            cside = rightCoodExoni - last_parent_exon[0][1]
            newid = []

            '''n side'''
            newidinner = str(len(parent_exons_eve[first_parent_exon][4]))
            if nside in parent_exons_eve[first_parent_exon][0]:
                newid += [parent_exons_eve[first_parent_exon][0][nside]]
            elif nside in hasLeftBformer:
                newid += [hasLeftBformer[nside][0]]
            else:
                event_occurence = len(parent_exons_eve[first_parent_exon][0])
                parent_exon_rank = parent_exons[first_parent_exon][0][3]
                parent_exon_cod_type = parent_exons[first_parent_exon][0][0]
                parent_exon_nature = parent_exons[first_parent_exon][0][2]
                parent_exons_status = parent_exons[first_parent_exon][0][1]
                tnewidn = [parent_exon_cod_type, parent_exons_status,
                           parent_exon_nature, parent_exon_rank, 'n', event_occurence]
                newid += [tnewidn]
                parent_exons_eve[first_parent_exon][0][nside] = tnewidn

            '''c side'''
            if cside in parent_exons_eve[last_parent_exon][1]:
                newid += [parent_exons_eve[last_parent_exon][1][cside]]
            elif cside in hasRightBlatter:
                newid += [hasRightBlatter[cside][0]]
            else:
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
            # TODO: newidinner, when will this be updated ???
            alt_exons[exon_pack] = [newid, coding_cood, tupRangeExoni]
            parent_exons_eve[exon_pack] = [
                {}, {}, {}, 2 if cod_status_pre == 1 else 1, {}]
        return parent_exons, alt_exons, parent_exons_eve, exons_aa_change
    except Exception as E:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        print ("ERROR: ncbia_raw..()/hidden_events", E, exc_tb.tb_lineno)


def template_feeder(parent_list, child_list, mega_set=None):
    '''
    the first and foremsot call will have parent list empty, chjild list as PI, and mega set empty
    '''
    flag = 1 if mega_set is not None else 0
    for ind, val in enumerate(child_list):
        rawCood = val[1]
        codingCood = val[2]
        aaSeq = val[0]
        exStatus = val[3]
        tupleExonKey = (rawCood, aaSeq, exStatus)

        parent_list += [(tupleExonKey, codingCood,
                         tuple(range(rawCood[0], rawCood[1] + 1)))]
        if flag:
            mega_set |= set(range(rawCood[0], rawCood[1] + 1))
    if flag:
        return parent_list, mega_set
    else:
        return parent_list
    '''
    It iterate the child list and records their rawcood, aaseq and exStatus as tuple(key) and add them to parent list, it also keep records of their rawcood spans in megaset (if its non empty)
    '''


def insertion_filter(insertion_cases, extension_cases, parent_exons_template):
    temp_lis_temp = [[len(i[2]), i] for i in insertion_cases]
    temp_lis_temp.sort()
    insertion_cases = [i[1] for i in temp_lis_temp]
    # sorting insertion cases on the basis of their tuiple spans which is the 4the element

    '''
    construct below:
    1. till insertion cases list exists,
        2. record the length of insertion case (from its tuple span) and ob, 0th case [temp lis]
        3. record it in the ele_remove lis
        4. iterate the 1st to last element in insertion list,
            5. if it intersects with 0th element
                6. record this element also in temp lis and ele_remove lis
        [by now we have all elements in insertion list whose spans intersects wuth 0th element]
        7. update insertion cases by remove those lemenst which intersects with the 0th elemenst [remove elemenst in step 3 from 1. ]

        check if more than 1 elem in step 2 [these are the group of intersecting exons]
            osrt thenm on length basis,
            if smallest one is larger than the 29 nucleotides then proceed:
                iterate sorted list, and once exon of span greater than 29 nucletides is captured, break the loop
            else:
                record 1st as element if interest
            all elemnets whcih aret of ineterst are the extension cases hence and elemnt of ineterste in the novel insertion case and is being added to the parent exons template
    return empty insertion cases, added exnetsion cases and parent exons emplate
    '''
    while insertion_cases:
        temp_lis = [(len(insertion_cases[0][2]), insertion_cases[0])]
        ele_remove = [insertion_cases[0]]
        for j in range(1, len(insertion_cases)):
            if set(insertion_cases[j][2]) & set(insertion_cases[0][2]):
                temp_lis += [(len(insertion_cases[j][2]), insertion_cases[j])]
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
    return insertion_cases, extension_cases, parent_exons_template


def exon_manipulator(storerawExons, raw_var, parent_exons, parent_exons_eve):
    '''
    difficult to comprehend the code
    Give A and G tags on basis of the count of raw coods and the isoforms
    '''
    try:
        for exon in parent_exons:
            # exon is key (rawcood, aaseq, exStatus), value = parent_exons[(rawcood, aaseq, exStatus)]=[['E', cod_status, '', idexon, 0, 0], codCood, tupRange]
            consAlt = 'G' if raw_var.count(
                exon[0]) == len(storerawExons) else 'A'
            # counting instances of the rawCood
            parent_exons[exon][0][2] = consAlt

        interMediateIds = {}  # store the rank/number/sequqncePositionSpecifier/1 or 2/or 3
        for exon in parent_exons:
            rawcoods = exon[0]
            ID = parent_exons[exon][0][3]
            interMediateIds[rawcoods] = ID

        # print (interMediateIds)
        exon_trans = {}
        for var in storerawExons:
            for exon in storerawExons[var]:
                # print(exon, exon[1], "temp")
                # exon[1] is rawCoods, pos is the rank/number/sequqncePositionSpecifier/1 or 2/or 3
                if exon[1] in interMediateIds:
                    # ommiting cases of retention, which will be captured later
                    pos = interMediateIds[exon[1]]
                    if pos not in exon_trans:
                        exon_trans[pos] = []
                    exon_trans[pos] += [var]
                # if exon3c1 and 3n1 occured in same trans T, value will be added twice and in set, will be conted only once
        '''
        Template temp repeat count(A) and template temp(B) are related lists
        if that is Alt exon and if its id isnt present in the A, add its rank occurences in B and record its tuple in B
        Below: ADDING F cATEGORY
        '''
        for exon in parent_exons:
            first_digit = parent_exons[exon][0][3]
            # but what about the retenstion cases, it will be former exon id
            if parent_exons[exon][0][2] == 'A':
                if len(set(exon_trans[first_digit])) == len(storerawExons):
                    # if the number of entries (transcripts) is equaal to total number of transcripts,
                    # (26 sep 22) because 4c1, 4b1 can be presne in 1 transcript and in 2 transripts hence, 4 can be presnet 4 times, so canging == to >=
                    # concern (27 sep 2022): changed, codes have been changed, and its almost optimal
                    # (27 sep 22), techinaclly, ity needs to be sorted and compared with transcript IDS, involving them overhere would be quite cumbersome, and hence being avoided.
                    # if best of correction is being seeked, look at the refineF categories routine in themain ()
                    parent_exons[exon][0][2] = 'F'
        has = {tuple(parent_exons[i][0][3:]): [] for i in parent_exons}
        for i in parent_exons:
            if tuple(parent_exons[i][0][3:]) in has:
                has[tuple(parent_exons[i][0][3:])] += [(parent_exons[i][0][3],
                                                        parent_exons[i][0][4], parent_exons[i][0][5], parent_exons[i][0][1])]
                # has[(4.0.0)]=[4.0.0, 2, 4.0.0, -1, 4.0.0, 0]
        has1 = {}
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
    except Exception as E:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        print ("ERROR: ncbia_raw..()/manipuator()", E, exc_tb.tb_lineno)


def ncbia(var_list, pi, exons_add):
    try:
        storerawExons = {}
        with open(os.path.join(exons_add, pi)) as fin:
            dat = pickle.load(fin)
        storerawExons[pi] = dat
        parent_exons_template = []
        mega_Set = set()
        parent_exons_template, mega_Set = template_feeder(
            parent_list=parent_exons_template, child_list=dat, mega_set=mega_Set)
        '''
        parent list is empty, child list is principal isoform and megaset is empty,
        output: parent_exons will have record of tuple and coding coods, megaset will have record of the raw cood span

        tupleExonKey = (rawCood, aaSeq, exStatus)
        parent_list += [(tupleExonKey, codingCood, tuple(range(rawCood[0], rawCood[1]+1)))]
        '''
        # print ("2nd part")
        raw_var = [i[0][0] for i in parent_exons_template]
        # list of raw_coods (genespan)
        sorted_vars = var_list[:]
        sorted_exons_o = []
        isoforms = len(var_list)
        for var in sorted_vars:
            # iterating the variable now for the gene group
            if var != pi:
                with open(os.path.join(exons_add,  var)) as fin:
                    exons = pickle.load(fin)
                storerawExons[var] = exons
                sorted_exons_o = template_feeder(
                    parent_list=sorted_exons_o, child_list=exons)
        sorted_exons = list(set(sorted_exons_o))
        # create 'parent_exons_template' equivalent for isoforms other than PI, which will have same child elements, no reference to megaset
        insertion_cases = []
        extension_cases = []
        for i in sorted_exons:
            if set(i[2]) & mega_Set:
                extension_cases += [i]
            else:
                insertion_cases += [i]
        # if raw coding spans intersect, they must be extension cases to parent_exons_template otherwise they will be insertion cases.

        insertion_cases.sort()
        insertion_cases, extension_cases, parent_exons_template\
            = insertion_filter(insertion_cases, extension_cases, parent_exons_template)

        # insertion cases would now be null, (appended into the parent exons template,)
        parent_exons_template.sort()
        parent_exons = {}
        alt_exons = {}

        parent_exons_eve = {}
        exons_aa_change = {}
        # above two categories will be populated later, will depict what they mean actually.

        for ind, val in enumerate(parent_exons_template):
            # val is [(tupleExonKey, codingCood, tuple(range(rawCood[0], rawCood[1]+1)) ) ]
            # tupleExonKey is (rawCood, aaSeq, exStatus)
            idexon = ind + 1
            unique_key = val[0]  # (rawCood, aaSeq, exStatus)
            cod_status = cod_status_primitive(unique_key[2], unique_key[1])
            codCood = val[1]
            leftpos = unique_key[0][0]
            rightpos = unique_key[0][1]
            range_tup = val[tupleSpanIndex]
            uniqueIdLis = ['E', cod_status, '', idexon, 0, 0]
            # 6 letter maturation begins
            # print ("2.1")
            parent_exons[unique_key] = [uniqueIdLis,
                                        codCood, range_tup]
            nidentity = leftpos - leftpos
            cidentity = rightpos - rightpos
            bidentity = (nidentity, cidentity)
            # print ("2.3")
            parent_exons_eve[unique_key] = [{nidentity: uniqueIdLis}, {cidentity: uniqueIdLis},
                                            {bidentity: uniqueIdLis}, 2 if cod_status == 1 else 1, {}]
            # print ("2.4")
        intron_retention_cases = []
        temp_list_refine = []

        # print ("3rd part")
        for i in extension_cases:
            ic = 0
            for j in parent_exons_template:
                if set(i[tupleSpanIndex]) & set(j[tupleSpanIndex]):
                    ic += 1
            temp_list_refine += [(ic, i)]

        # print ("4th")
        extension_cases = []
        for i in temp_list_refine:
            if i[0] == 1:
                extension_cases += [i[1]]
            else:
                intron_retention_cases += [i[1]]
        # print ("5th")
        parent_exons_upd = {}
        extension_cases.sort()
        for i in extension_cases:
            # print (i)
            if not (i[0] in parent_exons or i[0] in alt_exons or i[0] in exons_aa_change):
                # seems like a redundant condition, but not changing it now (PV: sep 19, 2022)
                # if key tuple is not already present
                # dont know if it eevr will be false
                # print ("5.0.1")
                parent_exons, alt_exons, parent_exons_eve, exons_aa_change = hidden_events1(
                    i[0], i[1], i[2], parent_exons, alt_exons, parent_exons_eve, exons_aa_change)
                # i[0] is tuple key (3 elem), i[1] are coding coods,
            else:
                pass
                # TODO: write to file o check if this ever turns true

        # ####################### DOne Till above ######################
        parent_exons_upd.update(parent_exons)
        parent_exons_upd.update(alt_exons)
        parent_exons_upd.update(exons_aa_change)
        raw_var += [i[0][0] for i in sorted_exons_o]
        # it had raw coods for the exons of principal template
        # used to count numbe rof exona dn tracripts
        parent_exons_upd, parent_exons_eve = exon_manipulator(
            storerawExons, raw_var, parent_exons_upd, parent_exons_eve)
        for i in parent_exons:
            parent_exons[i][0] = parent_exons_upd[i][0]
        for i in alt_exons:
            alt_exons[i][0] = parent_exons_upd[i][0]
        for i in exons_aa_change:
            exons_aa_change[i][0] = parent_exons_upd[i][0]
        intron_retention_cases.sort()
        for i in intron_retention_cases:
            # if i[0][0] == (4817, 6739):
            # print i[:4], "int"
            # print (i)
            if not (i[0] in parent_exons or i[0] in alt_exons or i[0] in exons_aa_change):
                parent_exons, alt_exons, parent_exons_eve, exons_aa_change = hidden_events1(
                    i[0], i[1], i[2],  parent_exons, alt_exons, parent_exons_eve, exons_aa_change)
        parent_exons.update(alt_exons)
        parent_exons.update(exons_aa_change)

        parent_exons_id_format = {}
        for i in parent_exons:
            if parent_exons[i][0][0] == "R":
                parent_exons_id_format[i] = [
                    ":".join(map(str, parent_exons[i][0])), parent_exons[i][1]]
            else:
                parent_exons_id_format[i] = [
                    ".".join(map(str, parent_exons[i][0])), parent_exons[i][1]]
        return parent_exons_id_format

    except Exception as E:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        print ("ERROR: ncbia_raw..()", E, exc_tb.tb_lineno)

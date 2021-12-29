#!/usr/bin/python3
##inp: ./name np_notA_oneID_self_xyzw_LR_mdl_pdelta.txt out
##find BP

import sys

def bp_check(LR, sline):
    if LR == "L":
        bp = sline[6]
    if LR == "R":
        bp = sline[5]
    return (bp)

def d_creation(d, chrn, gr, rn):
    rn = [rn]
    if chrn in d:
        d[chrn] = d_creation1(d[chrn], gr, rn)
    else:
        d[chrn] = {}
        d[chrn] = d_creation1(d[chrn], gr, rn)
    return(d)

def d_creation1(d, ID, args):
    if ID not in d:
        d[ID] = []
        d[ID].extend(i for i in args)
    else:
        d[ID].extend(i for i in args)
    return (d)

def first_check(chr, intrsctn, d, gr):
    LR = set()
    for rn in intrsctn:
        temp_id = chr + '+' + rn + '+' + gr
        LR.add(d[temp_id][0])
    if len(LR) == 2:
        return True
    else:
        return False

def full_check(chr1, chr2, rn, d, gr1, gr2):
    temp_id1 = chr1 + '+' + rn + '+' + gr1
    temp_id2 = chr2 + '+' + rn + '+' + gr2
    if d[temp_id1][0] == d[temp_id2][0]: #check RL are same
        if d[temp_id1][1] != d[temp_id2][1]: #check flag is different
            line_type = [chr1, chr2, rn, d[temp_id1][0], d[temp_id1][2], gr1, d[temp_id2][2], gr2]
            return(line_type)

with open(sys.argv[1], 'rt') as inp:
    with open(sys.argv[2], 'wt') as out:
#with open('synthr24ins3k_self_zw_sort3_p0.txt', 'rt') as inp:
#    with open('test.txt', 'wt') as out:
        
        header = inp.readline()  ##skip header line in inp
        out.write("chr1\tchr2\trn\tL|R\tbp1\tgr1\tbp2\tgr2\n")
        main_d = {}  # chr1 : {gr1 : (rn1, rn2, rn3, ...)
                     #         gr2 : (rn1, rn2, rn3, ...)
                     #         ... }
                     #chr2 : {...}
        full_d = {}  # chr_rn_gr : [flag, L|R, bp]
        gr = 1
        for line in inp:
            if line == ".\n":
                gr += 1
            else:
                sline = line.strip().split()
                rn = sline[0]
                chrn = sline[2]
                LR = sline[8]
                flag = sline[4]
                bp = bp_check(LR, sline)
                full_id = chrn + '+' + rn + '+' + str(gr)
                full = [LR, flag, bp]
                main_d = d_creation(main_d, chrn, str(gr), rn)
                full_d = d_creation1(full_d, full_id, full)
        
        
        for chr1 in main_d:
            for chr2 in main_d:
                if chr1 != chr2:
                    for gr1 in main_d[chr1]:
                        for gr2 in main_d[chr2]:
                            intrsctn = list(set(main_d[chr1][gr1]).intersection(set(main_d[chr2][gr2])))
                            if len(intrsctn) > 1:
                                lines = []
                                if first_check(chr1, intrsctn, full_d, gr1) and first_check(chr2, intrsctn, full_d, gr2):
                                    for rn in intrsctn:
                                        tempLine = full_check(chr1, chr2, rn, full_d, gr1, gr2)
                                        if type(tempLine) == list:
                                            lines.append(tempLine)
                                    if len(lines) == len(intrsctn):
                                        for line_type in lines:
                                            out.write('%s\n' % '\t'.join(line_type))

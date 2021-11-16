#!/usr/bin/python3
##inp: ./name np_vs_np_notA_0_16.sam np_vs_np_notA_0_16_xyzw.txt 300 np_vs_np_notA_0_16_xyzw_LR.txt 2>123noG.log &

import sys
import math

from decimal import Decimal
import datetime
import re
from multiprocessing import Pool


def CIGAR_translator(name, CIGAR, Z, chrL, U):
    HS = ['H', 'S']
    ##['D','N','H','P'] absent in read
    ##['I','S','H','P'] absent in chr

    full_s = re.findall('([0-9]+)([A-Z]+)', CIGAR)
    # list of pairs full_s[0][0] - first number, full_s[0][1] - first letter.

    ##check unclipped ends
    S1 = 0
    S2 = 0
    if full_s[-1][1] in HS:
        S2 = int(full_s[-1][0])
        full_s = full_s[:-1]
        ##del right S or H if exist
    if full_s[0][1] in HS:
        S1 = int(full_s[0][0])

    DI = 0
    W = Z  ##end pos on chr
    # X = 0  ##start pos in read
    # Y = 0  ##end pos in read
    M = 0
    for MatchState in full_s:
        Amount = int(MatchState[0])  #int function
        # if MatchState[1] in HS:
        #     X += Amount + 1
        #     Y += Amount
        if MatchState[1] == 'M':
            W += Amount
            #Y += Amount
            M += Amount
        elif MatchState[1] == 'D':
            W += Amount
            DI += Amount
        elif MatchState[1] == 'I':
            #Y += Amount
            DI += Amount
    qual = (Decimal(DI) + 1) / Decimal(M)
    log_qual = math.log10(qual)
    
    #print('Z,S1,S2,W,U', Z,S1,S2,W,U)
    if (Z > U and S1 > U) and (S2 < 100 or (chrL-W) < 100):
        LR = "R"
    elif ((chrL-W) > U and S2 > U) and (Z < 100 or S1 < 100): 
        LR = "L"
    else:
        LR = "N/A"
        
    return (Z, W, log_qual, LR)


def ReadOtherLine(Dictionary):
    rL = Dictionary['rL']
    min_rl = Dictionary['min_rl']
    if rL <= min_rl:
        return
    sline = Dictionary['line']
    chrL = Dictionary['chrL']
    name = sline[0]
    chr_name = sline[2]
    flag = int(sline[1])
    CIGAR = sline[5]
    Z = int(sline[3])  # start pos on chr (left)
    Z, W, log_qual, LR = CIGAR_translator(name, CIGAR, Z, chrL, U)
    if LR == "N/A":
        return
    line_type = [str(name), str(rL), str(chr_name), str(chrL), str(flag),
                 str(Z), str(W), str(log_qual), str(LR)]
    return (log_qual, line_type)


def ReadSQLine(line):  
    name = line[1][3:]
    if name not in np_read_len_d:
        read_len = int(line[2][3:])
        np_read_len_d[name] = read_len


if __name__ == '__main__':  #essential line
    ThreadCount = 8
    time = datetime.datetime.now()
    print("Start")
    with open(sys.argv[1], 'rt') as inp_sam:
        with open(sys.argv[2], 'wt') as out:
            U = int(sys.argv[3])

            ##pc_read less10more90
            np_read_len_d = {}
            logqual_list = []
            pc_read_th = 10
            ##pc_read_th < pc_read < (100 - pc_read_th)
            min_rl = 5000
            out.write('##read_name\tread_len\tchr_name\tchr_len\tflag\t' +
                      'Z\tW\tlog_qual\tL|R\n')

            SQLines = []
            OtherLines = []
            for line in inp_sam:
                sline = line.strip().split()
                if sline[0] == '@SQ':
                    SQLines.append(sline)

                elif sline[0] != '@PG':
                    if sline[0] != sline[2]:
                        if sline[2] != '*':
                            OtherLines.append(sline)

            print("Lines Sorted, SQLinesCount=" + str(len(SQLines)) + " OtherLines=" + str(len(OtherLines)))

            for OneSQLine in SQLines:
                ReadSQLine(OneSQLine)
            print("SQ Lines Calculated, it took=" + str(datetime.datetime.now() - time))

            ArgumentList = map(lambda x: {'line': x, 'rL': np_read_len_d[x[0]], 'chrL': np_read_len_d[x[2]], 'min_rl': min_rl},
                    OtherLines)
            
            with Pool(processes=ThreadCount) as executor:
                Result = executor.map(ReadOtherLine, ArgumentList)
            #Pool(processes=ThreadCount)
            #Result = map(ReadOtherLine, ArgumentList)
            
            print("All Lines Calculated, it took=" + str(datetime.datetime.now() - time))
            
            for OneResult in Result:
                if not OneResult:
                    continue
                logqual_list.append(float(OneResult[0]))
                out.write('%s\n' % '\t'.join(OneResult[1]))
                    ##instead of %s writes elements of line_type list divided by \t
                

            print("Lines Written, it took=" + str(datetime.datetime.now() - time))
            import numpy

            with open(sys.argv[2], 'rt') as inp:
                with open(sys.argv[4], 'wt') as out3:
                    SDq = numpy.std(logqual_list)
                    print(SDq)
                    mean_logqual = numpy.mean(logqual_list)
                    R = 1  # logqual filtration power
                    header = inp.readline()  ##skip header line in inp
                    out3.write(header)
                    for line in inp:
                        strip_line = line.strip()
                        sline = line.strip().split()
                        log_qual = float(sline[-2])

                        if log_qual < mean_logqual - R * SDq:
                            out3.write(line)

    time = datetime.datetime.now() - time
    print("Total Execution Time="+str(time))

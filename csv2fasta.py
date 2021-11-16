#!/usr/bin/python3
##inp: ./IDT2fasta.py IDT_PrimerQuest_Export.csv out.fasta
##IDT primer quest tool parser

import sys

with open(sys.argv[1], 'rt') as inp:
    try:
        out = open(sys.argv[2], 'a')
    except:
        out = open(sys.argv[2], 'wt')
    with out:
        header = inp.readline()
        for line in inp:
            if "Product" not in line:
                sline = line.strip().split('\t')
                out.write(">%s %s\n%s\n" % (sline[0], sline[1], sline[2]))

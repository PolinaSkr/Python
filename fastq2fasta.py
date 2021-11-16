#!/usr/bin/python3
##inp: ./fastq2fasta.py np_notA_0_16_oneID.fastq np_notA_0_16_oneID.fasta 2>log1.log &

import sys

with open(sys.argv[1], 'rt') as input:
    with open(sys.argv[2], 'wt') as out:
        line = input.readline()
        while line:
            out.write('>' + line[1:])
            out.write(input.readline())
            line = input.readline()
            line = input.readline()
            line = input.readline()

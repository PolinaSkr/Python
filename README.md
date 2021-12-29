# Python code examples

- **csv2fasta.py** and **fastq2fasta.py** are small scripts that I often use in my work, they are format converters.
- **bp_search_part1.py** is a part of my current project to find chromosome rearrangements 
using NGS datasets sequenced with the Oxford Nanopore platform. This program is parallelized because 
the size of the input data is quite large (41 GB). This part of code convert 
alignment from SAM format into ZW format, select only alignments which satisfies the boundary 
conditions and min read length, classify alignments into L and R, filter by quality. The output 
data is processed by another code part, I do not provide the full project code, 
since it has not been published yet.
- **bp_search_part5.py** also a part of chromosome rearrangements search project

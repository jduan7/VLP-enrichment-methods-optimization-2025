#!/usr/bin/python3

import pysam
import numpy as np
from Bio import SeqIO
import sys

# Check if correct number of arguments are passed
if len(sys.argv) != 3:
    print("Usage: python script.py <BAM file> <FASTA file>")
    sys.exit(1)

# Open the BAM file
samfile = pysam.AlignmentFile(sys.argv[1], "rb")

# Initialize containers
contigs = {}
mapped_set = set()
total_set = set()

# Define functions to filter reads
def read_good(bam):
    min_mapq = 10
    return bam.mapping_quality >= min_mapq

def read_bad(bam):
    max_mapq = 10
    return bam.mapping_quality < max_mapq

# Process reads to populate total_set and mapped_set
for read in samfile.fetch(until_eof=True):
    total_set.add(read.query_name)
    if read_good(read):
        mapped_set.add(read.query_name)

# Process FASTA file to extract contig information
for read in SeqIO.parse(open(sys.argv[2]), 'fasta'):
    contigs[read.id] = str(read.seq)

# Print header
print('\t'.join(["contig.id", "contig.length", "mapped", "mapped.low.qc", "total.mapped", "total"]))

# Analyze each contig
mapped_global_set = set()
for contig_id in contigs.keys():
    contig_length = len(contigs[contig_id])
    count = 0
    countbad = 0
    for read in samfile.fetch(contig_id):
        if read.is_secondary or read.is_supplementary:
            continue
        if read.query_name in mapped_global_set:
            continue
        mapped_global_set.add(read.query_name)
        
        if read_good(read):
            count +=1
        elif read_bad(read):
            countbad +=1
    total_mapped = count + countbad
    #count = samfile.count(contig_id, read_callback=read_good)
    #countbad = samfile.count(contig_id, read_callback=read_bad)
    #total_mapped = samfile.count(contig_id)
    print('\t'.join([contig_id, str(contig_length), str(count), str(countbad), str(total_mapped), str(len(total_set))]))

# Close the BAM file
samfile.close()


#%%
import pysam
import re
import gzip
import os
import sys
from collections import defaultdict, Counter
import argparse


def parse_args():

    parser = argparse.ArgumentParser(description='Add tags to bam file')
    parser.add_argument('-i', '--input', dest='input_bam', type=str, required=True,
                        help='BAM path for which to deduplicate')
    parser.add_argument('-o', '--out_bam', dest='output_bam', type=str, required=True,
                        help='Out BAM path')
    parser.add_argument('-n', '--num_tags', dest='num_tags', type=int, required=True, 
                        help="Number of tags")

    args = parser.parse_args()
    return args


#%%
def main():
    opts = parse_args()
    get_clusters(opts.input_bam, opts.output_bam, opts.num_tags)

#%%
def get_clusters(bam, output_bam, num_tags):
    """Parses a BAM file, groups reads into clusters according to their
    barcodes, filter for unique reads in each cluster, write reads to new BAM file.

    Each BAM record must have the barcode stored in the query name like so:

    ORIGINAL_READ_NAME::[Tag1][Tag2][Tag3][DPM/RPM]

    The tags should be enclosed in brackets and separated from
    the original read name with a double-colon.
    """

    clusters = defaultdict(set)
    pattern = re.compile('::' + num_tags * '\[([a-zA-Z0-9_\-]+)\]')
    total_reads = 0
    reads_written = 0
    
    with pysam.AlignmentFile(bam, "rb") as in_bam, \
    pysam.AlignmentFile(output_bam, "wb", template=in_bam) as out_bam:
        for read in in_bam.fetch(until_eof = True):
            total_reads +=1
            name = read.query_name
            match = pattern.search(name)
            barcode = ".".join(list(match.groups()))
            sequence = read.query_sequence
            currClus = clusters[barcode]
            if sequence in currClus:
               continue
            clusters[barcode].add(sequence)
            reads_written +=1
            out_bam.write(read)
    print('Total reads:', total_reads)
    print('Reads written out:', reads_written)

#%%
if __name__ == "__main__":
    main()

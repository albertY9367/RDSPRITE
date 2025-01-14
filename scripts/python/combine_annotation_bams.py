#%%
import pysam
from collections import defaultdict
import argparse
import re


def parse_args():

    parser = argparse.ArgumentParser(description='Add tags to bam file')
    parser.add_argument('-i', '--input', dest='input_bams', nargs='+', type=str, required=True,
                        help='BAM(s) path(s) for which to consolidate the XS/XT tag')
    parser.add_argument('-i2', '--input_2', dest='master_bam', type=str, required=True,
                        help='Master BAM to which to write consolidated tags')
    parser.add_argument('-o', '--out_bam', dest='output_bam', type=str, required=True,
                        help='Out BAM path')

    args = parser.parse_args()
    return args


#%%
def main():
    
    opts = parse_args()

    anno_dict = get_annotation(opts.input_bams)
    anno_clean = clean_annotation(anno_dict)
    add_annotation(opts.master_bam, opts.output_bam, anno_clean)




#%%
def get_annotation(bams):
    '''Read through bam files and extract featureCounts annotation for each read
    and combine intron, exon and repeat annotation

    Args:
        bams(list): A list of bam file paths

    '''
    read_annotation = defaultdict(set)

    for bam in bams:
        try:
            with pysam.AlignmentFile(bam, 'rb') as f:
                for read in f.fetch(until_eof = True):
                    name = read.query_name
                    #get featureCounts annotation
                    if read.has_tag('XT'):
                        anno = read.get_tag('XT')
                        for single in anno.split(','):
                            read_annotation[name].add(single)
                    #only XS flag present if no feature identified
                    elif read.has_tag('XS'):
                        anno = read.get_tag('XS')
                        #XS field exists for some other purpose with an integer value i type vs Z type
                        if anno == 'Unassigned_Ambiguity':
                            read_annotation[name].add('Ambiguous.none')
                        else: 
                            read_annotation[name].add('NoFeatures.none')
                    else:
                        raise Exception('XS tag missing, was featureCounts run?')
                    
        except ValueError:
            print('BAM file provided is not a BAM or is empty!')

    return read_annotation

#%%
def add_annotation(bam, output_bam, anno_dict):
    '''Add read annotation to a singular bam file

    Args:
        bam(str): Input bam to which annotation will be added
        output_bam(str): Output of input bam with annotation
        anno_dict(dict): read (key) and annotation (value) dictionary
    '''
    count = 0
    skipped = 0
    with pysam.AlignmentFile(bam, "rb") as in_bam, \
    pysam.AlignmentFile(output_bam, "wb", template=in_bam) as out_bam:

        for read in in_bam.fetch(until_eof = True):
            count += 1
            name = read.query_name
            anno = anno_dict.get(name, 'NoFeature.none')
            try:
                read.tags += [('XT', ''.join(anno))] #should be a single string (could be multiple annotations)
                out_bam.write(read)
            except KeyError:
                skipped += 1
    
    print('Total reads:', count)
    print('Reads with an error not written out:', skipped)

#%%
def clean_annotation(anno_dict):
    '''For a read with annotation, remove NoFeatures.none and concatenate the rest
    '''
    anno_out = defaultdict(str)
    for k, v in anno_dict.items():
        if len(v) > 1:
            new_anno = []
            for i in v:
                if i != 'NoFeatures.none':
                    new_anno.append(i)
            anno_out[k] = ';'.join(new_anno)
        else:
            anno_out[k] = ''.join(v)

    return anno_out

#%%
if __name__ == "__main__":
    main()

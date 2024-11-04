import pandas as pd
import numpy as np
import re
import tqdm


def parse_barcodes(clusters):
    barcodes, rpm, dpm, size = list(), list(), list(), list()
    with open(clusters, "r") as c:
        for line in tqdm.tqdm(c):
            barcode, *reads = line.rstrip('\n').split('\t')
            barcodes.append(barcode)
            total_size = len(reads)
            dpm_reads = [r for r in reads if r.startswith("DPM")]
            dpm_size = len(dpm_reads)
            size.append(total_size)
            dpm.append(dpm_size)
            rpm.append(total_size - dpm_size)
    df = pd.DataFrame({'Name': pd.Series(barcodes, dtype=str),
                       'RPM': pd.Series(rpm, dtype=int),
                       'DPM': pd.Series(dpm, dtype=int),
                       'Size': pd.Series(size, dtype=int)})
    return df

def parse_dpm(clusters):
    barcodes, chroms, starts, ends, strands, size, dpm_sizes = list(), list(), list(), list(), list(), list(), list()
    pattern = re.compile('([a-zA-Z0-9]+)\[(.+)\]_(.+):([0-9]+)\-([0-9]+)')
    with open(clusters, "r") as c:
        for line in tqdm.tqdm(c):
            barcode, *reads = line.rstrip('\n').split('\t')
            total_size = len(reads)
            dpm_reads = [r for r in reads if r.startswith('DPM')]
            dpm_size = len(dpm_reads)
            for read in dpm_reads:
                match = pattern.search(read)
                read_type, feature, chrom, start, end = match.groups()
                size.append(total_size)
                chroms.append(chrom)
                strands.append(feature)
                starts.append(int(start))
                ends.append(int(end))
                dpm_sizes.append(dpm_size)
                barcodes.append(barcode)
    df = pd.DataFrame({'Chromosome': pd.Series(chroms, dtype=str),
                       'Start': pd.Series(starts, dtype=int),
                       'End': pd.Series(ends, dtype=int),
                       'Name': pd.Series(barcodes, dtype=str),
                       'Strand': pd.Series(strands, dtype=str),
                       'Size': pd.Series(size, dtype=int),
                       'DPMSize': pd.Series(dpm_sizes, dtype=int)})
    return df

def parse_rpm(clusters):
    barcodes, chroms, starts, ends, strands, features, size, dpm_sizes = list(), list(), list(), list(), list(), list(), list(), list()
    pattern = re.compile('([a-zA-Z0-9]+)\[(.+)\]_(.+):([0-9]+)\-([0-9]+)')
    with open(clusters, "r") as c:
        for line in tqdm.tqdm(c):
            barcode, *reads = line.rstrip('\n').split('\t')
            total_size = len(reads)
            rpm_reads = [r for r in reads if r.startswith("RPM")]
            dpm_size = total_size- len(rpm_reads)
            for read in rpm_reads:
                match = pattern.search(read)
                read_type, feature, chrom, start, end = match.groups()
                if chrom.startswith('chr'):
                    size.append(total_size)
                    chroms.append(chrom)
                    anno, strand = feature.rsplit(';', 1)
                    strands.append(strand)
                    starts.append(int(start))
                    ends.append(int(end))
                    dpm_sizes.append(dpm_size)
                    barcodes.append(barcode)
                    features.append(anno)
    df = pd.DataFrame({'Chromosome': pd.Series(chroms, dtype=str),
                       'Start': pd.Series(starts, dtype=int),
                       'End': pd.Series(ends, dtype=int),
                       'Name': pd.Series(barcodes, dtype=str),
                       'Strand': pd.Series(strands, dtype=str),
                       'Feature': pd.Series(features, dtype=str),
                       'Size': pd.Series(size, dtype=int),
                       'DPMSize': pd.Series(dpm_sizes, dtype=int)})
    return df


def parse_repeats(clusters):
    barcodes, repeats, size, dpm_sizes  = list(), list(), list(), list()
    pattern = re.compile('([a-zA-Z0-9]+)\[(.+)\]_(.+):([0-9]+)\-([0-9]+)') 
    with open(clusters, "r") as c:
        for line in tqdm.tqdm(c):
            barcode, *reads = line.rstrip('\n').split('\t')
            total_size = len(reads)
            rpm = [r for r in reads if r.startswith("RPM")]
            dpm_size =total_size- len(rpm)
            for read in rpm:
                match = pattern.search(read)
                read_type, feature, chrom, start, end = match.groups()
                if not chrom.startswith('chr'):
                    size.append(total_size)
                    dpm_sizes.append(dpm_size)
                    barcodes.append(barcode)
                    repeats.append(chrom)
    df = pd.DataFrame({'Name': pd.Series(barcodes, dtype=str),
                       'Repeat': pd.Series(repeats, dtype=str),
                       'Size': pd.Series(size, dtype=int),
                       'DPMSize': pd.Series(dpm_sizes, dtype=int)})
    return df


def parse_feature(feature_string):
    feature_dict = {}
    for item in feature_string.split(';'):
        anno, tag = item.rsplit('.',1)
        if tag in feature_dict:
            feature_dict[tag] = 'AMB'
        else:
            feature_dict[tag] = anno
    return (feature_dict.get('exon', str(np.nan)),
            feature_dict.get('intron', str(np.nan)),
            feature_dict.get('repeat', str(np.nan)))

def split_feature(df):
   df['Split'] = df['Feature'].apply(parse_feature)
   df[['Exon', 'Intron', 'Repeat']] =  pd.DataFrame(df['Split'].tolist(), index=df.index)
   del df['Split']
   del df['Feature']
   return df

def simplify_barcodes(barcode_dict, df, repeat=False):
    df.loc[:,'Name'] = df['Name'].map(barcode_dict)
    df = df.astype({'Name':'int64'}, copy=False)
    if repeat == True:
       df.sort_values(by=['Name', 'Repeat'], inplace=True)
    else:
       df.sort_values(by=['Name', 'Chromosome'], inplace=True)
    return df

import pandas as pd
import conversion as cv
import numpy as np
import re
import argparse

def main():
    args = parse_arguments()

    if args.conversion_type == 'DPM':
        convert_clusters_dpm(args.clusters, args.output, args.min_size, args.max_size)
    elif args.conversion_type == 'RPM':
        convert_clusters_rpm(args.clusters, args.output, args.min_size, args.max_size) 
    elif args.conversion_type == 'DPM_RPM':
        convert_clusters_rpm_dpm(args.clusters, args.output, args.min_size, args.max_size)
 
def convert_clusters_dpm(clusters, output, min_size, max_size):

    store = pd.HDFStore(output, 'w', complevel=0, complib='blosc') #lzo

    df = cv.parse_dpm(clusters)
    barcodes = cv.parse_barcodes(clusters)
    barcodes = barcodes.query('DPM!=0')

    # Filter for size
    barcodes  = barcodes.query('Size>=@min_size and Size<=@max_size')
    df = df.query('Size>=@min_size and Size<=@max_size')

    # Reassign the barcodes
    barcodes = barcodes.reset_index(drop=True)
    barcode_dict = barcodes.reset_index().set_index('Name')['index'].to_dict()
    barcodes = barcodes.reset_index()
    barcodes.rename(columns = {'Name':'Barcode', 'index':'Name'}, inplace = True)
    print(barcodes.head())

    df = cv.simplify_barcodes(barcode_dict,df)
    print(df.head())

    # Split the dpm reads
    dpm_only = df.query('Size==DPMSize')
    dpm_rpm = df.query('Size!=DPMSize')
 
    # Add to the store
    store.append('DPM_only', dpm_only, format='table', append=True, data_columns=True, index=False)
    store.append('DPM_RPM', dpm_rpm, format='table', append=True, data_columns=True, index=False)    
    store.append('Names', barcodes, format='table', append=True, data_columns=True, index=False)

    store.create_table_index('DPM_only', columns=list(dpm_only.columns), optlevel=9, kind='full')
    store.create_table_index('DPM_RPM', columns=list(dpm_rpm.columns), optlevel=9, kind='full')
    store.create_table_index('Names', columns=list(barcodes.columns), optlevel=9, kind='full')

    store.close()

def convert_clusters_rpm(clusters, output, min_size, max_size):

    store = pd.HDFStore(output, 'w', complevel=0, complib='blosc') #lzo

    df = cv.parse_rpm(clusters)
    repeats = cv.parse_repeats(clusters)
    barcodes = cv.parse_barcodes(clusters)

    # Filter for rna barcodes
    barcodes = barcodes.query('RPM!=0')

    # Filter for size
    barcodes  = barcodes.query('Size>=@min_size and Size<=@max_size')
    df = df.query('Size>=@min_size and Size<=@max_size')
    repeats = repeats.query('Size>=@min_size and Size<=@max_size')

    # Reassign the barcodes
    barcodes = barcodes.reset_index(drop=True)
    barcode_dict = barcodes.reset_index().set_index('Name')['index'].to_dict()
    barcodes = barcodes.reset_index()
    barcodes.rename(columns = {'Name':'Barcode', 'index':'Name'}, inplace = True)
    print(barcodes.head())

    df = cv.simplify_barcodes(barcode_dict,df)
    repeats = cv.simplify_barcodes(barcode_dict, repeats, True)
   
    # Parse features of RPM
    df = cv.split_feature(df)
    print(df.head())

    # Split the rpm
    df_only = df.query('DPMSize==0')
    df_dpm = df.query('DPMSize!=0')
    repeats_only = repeats.query('DPMSize==0')
    repeats_dpm = repeats.query('DPMSize!=0')
    
    print('RPM only: ', len(df_only)) 
    print('RPM DPM: ', len(df_dpm))
    print('Repeats only: ', len(repeats_only))
    print('Repeats DPM: ', len(repeats_dpm))

    # Add to the store
    store.append('RPM', df_dpm, format='table', append=True, data_columns=True, index=False)
    store.append('RepeatRPM', repeats_dpm, format='table', append=True, data_columns=True, index=False)    
    store.append('RPM_only', df_only, format='table', append=True, data_columns=True, index=False)
    store.append('RepeatRPM_only', repeats_only, format='table', append=True, data_columns=True, index=False)    
    store.append('Names', barcodes, format='table', append=True, data_columns=True, index=False)

    store.create_table_index('RPM', columns=list(df_dpm.columns), optlevel=9, kind='full')
    store.create_table_index('RepeatRPM', columns=list(repeats_dpm.columns), optlevel=9, kind='full') 
    store.create_table_index('RPM_only', columns=list(df_only.columns), optlevel=9, kind='full')
    store.create_table_index('RepeatRPM_only', columns=list(repeats_only.columns), optlevel=9, kind='full')
    store.create_table_index('Names', columns=list(barcodes.columns), optlevel=9, kind='full')

    store.close()


def convert_clusters_rpm_dpm(clusters, output, min_size, max_size):

    store = pd.HDFStore(output, 'w', complevel=0, complib='blosc') #lzo

    rpm = cv.parse_rpm(clusters)
    dpm = cv.parse_dpm(clusters)
    repeats = cv.parse_repeats(clusters)
    barcodes = cv.parse_barcodes(clusters)

    # Filter for rna-dna
    barcodes = barcodes.query('RPM!=0 and DPM!=0')
    rpm = rpm.query('DPMSize!=0')
    repeats = repeats.query('DPMSize!=0')
    dpm = dpm.query('DPMSize!=Size')

    # Filter for size
    barcodes  = barcodes.query('Size>=@min_size and Size<=@max_size')
    dpm = dpm.query('Size>=@min_size and Size<=@max_size')
    rpm = rpm.query('Size>=@min_size and Size<=@max_size') 
    repeats = repeats.query('Size>=@min_size and Size<=@max_size')

    # Reassign the barcodes
    barcodes = barcodes.reset_index(drop=True)
    barcode_dict = barcodes.reset_index().set_index('Name')['index'].to_dict()
    barcodes = barcodes.reset_index()
    barcodes.rename(columns = {'Name':'Barcode', 'index':'Name'}, inplace = True)
    print(barcodes.head())

    dpm = cv.simplify_barcodes(barcode_dict,dpm)
    rpm = cv.simplify_barcodes(barcode_dict,rpm)
    repeats = cv.simplify_barcodes(barcode_dict, repeats, True)
    print(dpm.head())
    print(repeats.head())
   
    # Annotate RNAs
    rpm = cv.split_feature(rpm) 
    print(rpm.head())

    # Add to the store
    store.append('RPM', rpm, format='table', append=True, data_columns=True, index=False)
    store.append('DPM', dpm, format='table', append=True, data_columns=True, index=False)
    store.append('RepeatRPM', repeats, format='table', append=True, data_columns=True, index=False)    
    store.append('Names', barcodes, format='table', append=True, data_columns=True, index=False)

    store.create_table_index('RPM', columns=list(rpm.columns), optlevel=9, kind='full')
    store.create_table_index('DPM', columns=list(dpm.columns), optlevel=9, kind='full')
    store.create_table_index('RepeatRPM', columns=list(repeats.columns), optlevel=9, kind='full')
    store.create_table_index('Names', columns=list(barcodes.columns), optlevel=9, kind='full')

    store.close()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Generates a cooler file from a cluster file.')
    parser.add_argument('-c', '--clusters',
                        metavar = "FILE",
                        action = "store",
                        required=True,
                        help = "The input cluster file.")
    parser.add_argument('-o', '--output',
                        metavar = "FILE",
                        action = "store",
                        required=True,
                        help = "The output cooler file.")
    parser.add_argument('--min_size',
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        required=True,
                        help = "The minimum cluster size to use.")
    parser.add_argument('--max_size',
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        required=True,
                        help = "The maximum cluster size to use.")
    parser.add_argument('--conversion_type',
                        metavar = 'TYPE',
                        action = 'store',
                        required=True,
                        help = "The type of output file to generate")
    return parser.parse_args()

if __name__ == "__main__":
    main()





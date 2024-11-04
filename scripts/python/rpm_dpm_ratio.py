import pandas as pd
import numpy as np
import argparse
import glob
import matplotlib.pyplot as plt
import tqdm
import matplotlib
import matplotlib.font_manager
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.size']=12

def main():
    args = parse_arguments()
    pattern = args.directory + "/*" + args.pattern
    files = glob.glob(pattern)
    results_binned = []
    results_sizes = []
    results_types = []
    names = []
    for f in files:
       print('Working on file: ', f)
       counted = count_rpm_dpm(f)
       binned = bin_counts(counted)
       results_binned.append(binned)
       results_sizes.append(normalize_sizes(binned))
       results_types.append(normalize_types(binned))
       names.append(f.split('/')[-1])
    save_dataframes(results_binned, names, args.directory)
    plot_size_profile(results_sizes, names, args.directory)
    plot_type_profile(results_types, names, args.directory)

def save_dataframes(dfs, names, output):
    for i in range(len(dfs)):
        df = dfs[i]
        name = names[i]
        out = output + '/' + name + '.csv'
        df.to_csv(out, sep = '\t')
  

def count_rpm_dpm(filename):
    counts = []
    rpm = []
    dpm = []
    repeats = []
    with open(filename, 'r') as f:
        for line in tqdm.tqdm(f):
            reads = line.split('\t')[1:]
            total = len(reads)
            rpm_reads = [r for r in reads if r.startswith('RPM')]
            dpm_counts = total - len(rpm_reads)
            if len(rpm_reads) != 0:
                rpm_aligned = len([r for r in rpm_reads if '_chr' in r])
                rpm.append(rpm_aligned)
                repeats.append(len(rpm_reads) - rpm_aligned)
            else:
                rpm.append(0)
                repeats.append(0)
            counts.append(total)
            dpm.append(dpm_counts)
    df = pd.DataFrame.from_dict({'Total':counts,'RNA': rpm,'DNA': dpm, 'Repeats':repeats})
    return df    

def bin_counts(df):   
    bins = pd.IntervalIndex.from_tuples([(0, 1), (1, 10), (10, 100), (100,1000), (1000,1000000)])
    df['Groups'] = pd.cut(df['Total'],bins)
    df_summed = df.groupby('Groups').sum()
    return df_summed

def normalize_sizes(df_summed):
    df_trans = df_summed.T
    df_trans.columns = df_trans.columns.astype(str)
    rename_mapping = {"(0, 1]": "Single", "(1, 10]":"2-10", "(10, 100]":"11-100", "(100, 1000]":"101-1000", "(1000, 1000000]": ">1000"}
    df_trans = df_trans.rename(columns=rename_mapping)
    df_trans['Total'] = df_trans.sum(axis=1)
    df_normed = df_trans.T
    df_normed = df_normed.divide(df_normed['Total'], axis=0)[['DNA', 'RNA', 'Repeats']]
    return df_normed

def normalize_types(df_summed):
    df_normed = df_summed.div(df_summed.sum(axis=0), axis=1)[['Total', 'DNA', 'RNA', 'Repeats']]
    df_trans = df_normed.T
    df_trans.columns = df_trans.columns.astype(str)
    rename_mapping = {"(0, 1]": "Single", "(1, 10]":"2-10", "(10, 100]":"11-100", "(100, 1000]":"101-1000", "(1000, 1000000]": ">1000"}
    df_trans = df_trans.rename(columns=rename_mapping)
    df_trans = df_trans[['Single','2-10','11-100','101-1000','>1000']]
    return df_trans


def plot_size_profile(dfs, names, output):
    fig, axes = plt.subplots(nrows=1, ncols=len(dfs), figsize=(16,8))
    if len(dfs) == 1:
        axes = [axes]
    ax_position = 0
    for i in range(len(dfs)):
        subset = dfs[i]
        print(subset)
        ax = subset.plot(kind="bar", stacked=True, ax=axes[ax_position])
        #ax.set_title(names[i])
        ax.set_title('Replicate ' + str(i+1))
        ax.set_ylabel("Proportion"),
        ax.set_xlabel("Cluster Size")
        ax.set_ylim(0, 1)
        ax.set_xticklabels(labels=['Single', '2-10', '11-100', '101-1000', '>1000', 'Total'], rotation=90, minor=False)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend().set_visible(False)
        ax_position += 1
    axes[len(dfs)-1].legend(['DNA', 'RNA', 'Repeats'],loc='center left', bbox_to_anchor = [1.0,0.5])
    plt.tight_layout()
    fig.savefig(output + '/RNA_DNA_Ratios.pdf', bbox_inches = 'tight')

def plot_type_profile(dfs, names, output):
    fig, axes = plt.subplots(nrows=1, ncols=len(dfs), figsize=(16,8))
    if len(dfs) == 1:
        axes = [axes]
    ax_position = 0
    for i in range(len(dfs)):
        subset = dfs[i]
        print(subset)
        ax = subset.plot(kind="bar", stacked=True, ax=axes[ax_position])
        #ax.set_title(names[i])
        ax.set_title('Replicate '+ str(i+1))
        ax.set_ylabel("Proportion"),
        ax.set_xlabel("Read Type")
        ax.set_ylim(0, 1)
        ax.set_xticklabels(labels=['Total', 'DNA', 'RNA', 'Repeat'], rotation=90, minor=False)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend().set_visible(False)
        ax_position += 1
    axes[len(dfs)-1].legend(['Single','2-10','11-100', '101-1000', '>1000'],loc='center left', bbox_to_anchor = [1.0,0.5])
    plt.tight_layout()
    fig.savefig(output + '/cluster_size_profiles.pdf', bbox_inches = 'tight')


def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Generate the cluster size distribution plot.')

    parser.add_argument('--directory',
                        metavar = "FILE",
                        action = "store",
                        help = "The directory with SPRITE clusters")
    parser.add_argument('--pattern',
                        metavar = "FILE",
                        action = "store",
                        help = "The pattern for file names")
    return parser.parse_args()

if __name__ == "__main__":
    main()
            













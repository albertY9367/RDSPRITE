import glob
import argparse
import pysam
import numpy as np
from collections import Counter, defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import pandas as pd

"""
Gnerate all cluster statistics from master bam file
Merging the function of generate_cluster_statistics.py, max_representation_ecdfs_perlib.py, get_bead_size_distribution.py

Specifically,
    Count 1) number of clusters, 2) number of DPM reads (aligned) and 3) number of  BPM reads within each clusterfile for a directory of clusterfiles.
    Generate maximum representation ecdfs for bead type representation within clusters.
    Profiles the proportion of clusters and proportion of reads within various cluster size categories. Considers DPM and BPM reads seperately (based on the input parameter readtype).
"""


def main():
    args = parse_arguments()
    search = args.directory + "/*" + args.pattern
    files = glob.glob(search)
    
    cluster_counts_dpm = []
    read_counts_dpm = []
    cluster_counts_rpm = []
    read_counts_rpm = []

    for f in files:
        
        df1_dpm, df2_dpm, df1_rpm, df2_rpm = generate_statistics(f, args.directory)
        cluster_counts_dpm.append(df1_dpm)
        read_counts_dpm.append(df2_dpm)
        cluster_counts_rpm.append(df1_rpm)
        read_counts_rpm.append(df2_rpm)

    # DNA statistics
    cluster_df_dpm = pd.concat(cluster_counts_dpm, axis=1).transpose()
    read_df_dpm = pd.concat(read_counts_dpm, axis=1).transpose()
    cluster_fig = plot_profile(cluster_df_dpm)
    cluster_fig.savefig(
        args.directory + "/" + "DPM" + "_cluster_distribution.pdf",
        bbox_inches="tight",
    )
    read_fig = plot_profile(read_df_dpm)
    read_fig.savefig(
        args.directory + "/" + "DPM" + "_read_distribution.pdf",
        bbox_inches="tight",
    )
    
    # RNA statistics
    cluster_df_rpm = pd.concat(cluster_counts_rpm, axis=1).transpose()
    read_df_rpm = pd.concat(read_counts_rpm, axis=1).transpose()
    cluster_fig = plot_profile(cluster_df_rpm)
    cluster_fig.savefig(
        args.directory + "/" + "RPM" + "_cluster_distribution.pdf",
        bbox_inches="tight",
    )
    read_fig = plot_profile(read_df_rpm)
    read_fig.savefig(
        args.directory + "/" + "RPM" + "_read_distribution.pdf",
        bbox_inches="tight",
    )

def generate_statistics(bamfile, dir):
    """
    Loop through all clusters within a bamfile, counting DPM and BPM reads

    Args:
        bamfile(str): Path to bamfile
    """
    # count statistics
    cluster = 0
    dpm = 0
    rpm = 0
    barcodes = set()

    # get bead size distribution
    # count rpm and dpm at the same time
    bins = np.array([0, 1, 5, 10, 20, 50, 100, 200])  # bins with data are 1-8

    # DPM
    cluster_counts_dpm = defaultdict(int)
    read_counts_dpm = defaultdict(int)
    # RPM
    cluster_counts_rpm = defaultdict(int)
    read_counts_rpm = defaultdict(int)

    # the cluster dictionaries
    dpm_dict = defaultdict(int)
    rpm_dict = defaultdict(int)

    count = 0
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            count += 1
            if count % 100000 == 0:
                print(count)
            read_type = read.get_tag("RT")
            barcode = read.get_tag("RC")
            try:
                if "DPM" in read_type:
                    dpm_dict[barcode] += 1 # get bead size distribution
                elif "RPM" in read_type:
                    rpm_dict[barcode] += 1
            except KeyError:
                pass

    for bin in np.arange(
        1, len(bins) + 1
    ):  # initialize all bins in case any end up being empty categories
        # DPM
        cluster_counts_dpm[bin] = 0
        read_counts_dpm[bin] = 0
        # RPM
        cluster_counts_rpm[bin] = 0
        read_counts_rpm[bin] = 0
    
    # DPM cluster size distribution
    for bc in dpm_dict.keys():
        count_dpm = dpm_dict[bc]
        if bc not in barcodes:
            cluster += 1
        dpm += count_dpm
        bin_dpm = np.digitize(count_dpm, bins, right=True)
        cluster_counts_dpm[bin_dpm] += 1
        read_counts_dpm[bin_dpm] += count_dpm
    # RPM cluster size distribution
    for bc in rpm_dict.keys():
        count_rpm = rpm_dict[bc]
        if bc not in barcodes:
            cluster += 1
        rpm += count_rpm
        bin_rpm = np.digitize(count_rpm, bins, right=True)
        cluster_counts_rpm[bin_rpm] += 1
        read_counts_rpm[bin_rpm] += count_rpm


    # generate cluster statistics
    print("For bamfile ", bamfile)
    print("Total number of clusters: ", cluster)
    print("Total number of DPM: ", dpm)
    print("Total number of RPM: ", rpm)

    filepath = os.path.join(dir, "cluster_statistics.txt")
    f = open(filepath, "a")
    f.write("For bamfile " + str(bamfile) + "\n")
    f.write("Total number of clusters: " + str(cluster) + "\n")
    f.write("Total number of DPM: " + str(dpm) + "\n")
    f.write("Total number of RPM: " + str(rpm) + "\n")
    f.close()


    # get bead size distribution
    df_cluster_counts_dpm = pd.DataFrame.from_dict(cluster_counts_dpm, orient="index")
    df_cluster_counts_dpm.columns = [bamfile.rsplit("/", 1)[-1].rsplit(".bam")[0]]
    df_read_counts_dpm = pd.DataFrame.from_dict(read_counts_dpm, orient="index")
    df_read_counts_dpm.columns = [bamfile.rsplit("/", 1)[-1].rsplit(".bam")[0]]

    df_cluster_counts_rpm = pd.DataFrame.from_dict(cluster_counts_rpm, orient="index")
    df_cluster_counts_rpm.columns = [bamfile.rsplit("/", 1)[-1].rsplit(".bam")[0]]
    df_read_counts_rpm = pd.DataFrame.from_dict(read_counts_rpm, orient="index")
    df_read_counts_rpm.columns = [bamfile.rsplit("/", 1)[-1].rsplit(".bam")[0]]

    return df_cluster_counts_dpm, df_read_counts_dpm, df_cluster_counts_rpm, df_read_counts_rpm

def plot_profile(df):
    """
    Plot the proportion of reads within each size category as a stacked bar graph
    Args:
        df(dataframe): binned counts
    """
    columns = ["1", "2-5", "6-10", "11-20", "21-50", "51-100", "101-200", "201+"]
    df.columns = columns
    df = df.div(df.sum(axis=1), axis=0)
    plot = df.plot(
        kind="bar", stacked=True, ylabel="Proportion", cmap=cm.get_cmap("Dark2")
    )
    plot.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
    plot.set_ylabel = "Proportion"
    return plot.get_figure()

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate the statistics for all clusterfiles in a directory"
    )
    parser.add_argument(
        "--directory",
        metavar="FILE",
        action="store",
        help="The directory of clusters file",
    )
    parser.add_argument(
        "--pattern", action="store", help="The pattern of cluster file names"
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
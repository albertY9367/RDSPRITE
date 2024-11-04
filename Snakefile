'''
Aim: A Snakemake workflow to process RNA-DNA and DNA-DNA SPRITE-seq data
'''

import os 
import sys
import json
import numpy as np
import datetime
from pathlib import Path


##############################################################################
# Initialize settings
##############################################################################


# Copy config file into logs
v = datetime.datetime.now()
run_date = v.strftime('%Y.%m.%d.')

try:
    config_path = config["config_path"]
except:
    config_path = 'config.yaml'

configfile: config_path

try:
    email = config['email']
except:
    email = None
    print("Will not send email on error")


##############################################################################
# Location of scripts
##############################################################################

barcode_id_jar = "scripts/java/BarcodeIdentification_v1.2.0.jar"
lig_eff = "scripts/python/get_ligation_efficiency.py"
split_fq = "scripts/python/split_dpm_rpm_fq_barcodes.py"
atttb = "scripts/python/add_tnx_tag_to_bam.py"
add_chr = "scripts/python/ensembl2ucsc.py"
get_clusters = "scripts/python/get_clusters.py"
comb_anno = "scripts/python/combine_annotation_bams.py"
remove_dup = "scripts/python/remove_duplicate_rRNA.py"
tag_bam = "scripts/python/tag_bam.py"
split_fastq = "scripts/split_fastq.sh"
generate_statistics = "scripts/python/generate_statistics.py"

##############################################################################
# Load other general settings
##############################################################################

try:
    bid_config = config['bID']
    print('Using BarcodeID config', bid_config)
except:
    bid_config = 'config.txt'
    print('Config "bID" not specified, looking for config at:', bid_config)

try:
    formatfile = config['format']
    print('Using Barcode Formatfile', formatfile)
except:
    formatfile = 'format.txt'
    print('Formatfile not specified, looking for file at:', formatfile)

try:
    num_tags = config['num_tags']
    print('Using', num_tags, 'tags')
except:
    num_tags = "5"
    print('Config "num_tags" not specified, using:', num_tags)

try:
    adaptors = "-g file:" + config['cutadapt']
    print('Using cutadapt sequence file', adaptors)
except:
    adaptors = "-g GGTGGTCTTT -g GCCTCTTGTT \
        -g CCAGGTATTT -g TAAGAGAGTT -g TTCTCCTCTT -g ACCCTCGATT"
    print("No file provided for cutadapt. Using standard cutadapt sequences")

try:
    assembly = config['assembly']
    assert assembly in ['mm10', 'hg38'], 'Only "mm10" or "hg38" currently supported'
    print('Using', assembly)
except:
    print('Config "assembly" not specified, defaulting to "mm10"')
    assembly = 'mm10'

try:
    samples = config['samples']
    print('Using samples file:', samples)
except:
    samples = './samples.json'
    print('Defaulting to working directory for samples json file')

try:
    out_dir = config['output_dir']
    print('All data will be written to:', out_dir)
except:
    out_dir = ''
    print('Defaulting to working directory as output directory')


try:
    num_chunks = config['num_chunks']
    print('Samples will be chunked into:', num_chunks)
except:
    num_chunks = 2
    print('Num chunks not specified. Defaulting to 2')

try:
    conda_env = config['conda_env']
    print('Using conda environment:', conda_env)
except:
    conda_env = "envs/sprite.yaml"
    print('Conda environment not specified. Defaulting to envs/sprite.yaml')
################################################################################
# Annotations
################################################################################

try:
    anno_repeats_saf = config['anno_repeats_saf'][config['assembly']]
    mask = config['mask'][config['assembly']]
    exon_intron_saf = config['exon_intron_saf'][config['assembly']]
    antisense_saf = config['antisense_saf'][config['assembly']]
except:
    print('Annotation or mask path not specified in config.yaml')
    sys.exit() #no default, exit

try:
    anno_rename = config['anno_rename'][config['assembly']]
except:
    anno_rename = ''
    print('Annotations will not be renamed')

################################################################################
# Aligner Indexes
################################################################################

# DNA aligner
try:
    bowtie2_index = config['bowtie2_index'][config['assembly']]
    bowtie2_repeat_index = config['bowtie2_repeat_index'][config['assembly']]
except:
    print('Bowtie2 index not specified in config.yaml')
    sys.exit() #no default, exit


# RNA aligner
try:
    hisat2_index = config['hisat2_index'][config['assembly']]
    hisat2_ss = config['hisat2_splice_sites'][config['assembly']]
except:
    print('RNA aligner not specified in config.yaml')
    sys.exit()


################################################################################
# Make output directories (aren't created automatically on cluster)
################################################################################

Path(out_dir + "workup/logs/cluster").mkdir(parents=True, exist_ok=True)
out_created = os.path.exists(out_dir + "workup/logs/cluster")
print('Output logs path created:', out_created)

################################################################################
# Get sample files
###############################################################################

# Prep samples from fastq directory using fastq2json_updated.py, now load json file
FILES = json.load(open(samples))
ALL_SAMPLES = sorted(FILES.keys())

ALL_FASTQ = []
for SAMPLE, file in FILES.items():
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R1')])
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R2')])

CONFIG = [out_dir + "workup/logs/config_" + run_date + "yaml"]

NUM_CHUNKS = [f"{i:03}" for i in np.arange(0, num_chunks)]

################################################################################
# Trimming
################################################################################
SPLIT_FQ = expand(out_dir + "workup/splitfq/{sample}_{read}.part_{splitid}.fastq.gz", sample=ALL_SAMPLES, read = ["R1", "R2"], splitid=NUM_CHUNKS)

TRIM = expand([out_dir + "workup/trimmed/{sample}_R1.part_{splitid}_val_1.fq.gz", out_dir + "workup/trimmed/{sample}_R2.part_{splitid}_val_2.fq.gz"], sample = ALL_SAMPLES, splitid=NUM_CHUNKS)
TRIM_LOG = expand(out_dir + "workup/trimmed/{sample}_{read}.part_{splitid}.fastq.gz_trimming_report.txt", sample = ALL_SAMPLES, read = ["R1", "R2"], splitid = NUM_CHUNKS)
TRIM_RD = expand([out_dir + "workup/trimmed/{sample}_R1.part_{splitid}_val_1_RDtrim.fq.gz", out_dir + "workup/trimmed/{sample}_R2.part_{splitid}_val_2_RDtrim.fq.gz"], sample = ALL_SAMPLES, splitid=NUM_CHUNKS)

################################################################################
# Logging
################################################################################

LE_LOG_ALL = [out_dir + "workup/ligation_efficiency.txt"]

MULTI_QC = [out_dir + "workup/qc/multiqc_report.html"]

################################################################################
# Barcoding
################################################################################
BARCODEID = expand(out_dir + "workup/fastqs/{sample}_{read}.part_{splitid}.barcoded.fastq.gz", sample = ALL_SAMPLES, read = ["R1", "R2"], splitid=NUM_CHUNKS)

SPLIT_ALL = expand([out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_rpm.fastq.gz", out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_dpm.fastq.gz"], sample=ALL_SAMPLES, splitid=NUM_CHUNKS)


################################################################################
# Bowtie2 alignment + workup
################################################################################

Bt2_DNA_ALIGN = expand(out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.DNA.bowtie2.mapq20.bam", sample=ALL_SAMPLES, splitid=NUM_CHUNKS)

Bt2_RNAr = expand([out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNAr.bowtie2.bam", out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNAr.bowtie2.unmapped.bam"], sample=ALL_SAMPLES, splitid=NUM_CHUNKS)

Bt2_TAG_ALL = expand(out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNAr.bowtie2.tag.bam",sample=ALL_SAMPLES, splitid=NUM_CHUNKS)

################################################################################
# HISAT2 Alignment
################################################################################

Ht2_RNA_ALIGN = expand([out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNA.hisat2.mapq20.bam", out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNA.hisat2.unmapped.lowmq.fq.gz"], sample=ALL_SAMPLES, splitid=NUM_CHUNKS)

Ht2_ANNO_RNA = expand([out_dir + "workup/annotations/{sample}.RNAexin.hisat2.mapq20.bam.featureCounts.bam", out_dir + "workup/annotations/{sample}.RNAr.hisat2.mapq20.bam.featureCounts.bam", out_dir + "workup/annotations/{sample}.RNAanti.hisat2.mapq20.bam.featureCounts.bam"], sample=ALL_SAMPLES, splitid=NUM_CHUNKS)

RNA_COMBINE = expand(out_dir + "workup/annotations/{sample}.RNA.merged.anno.bam", sample=ALL_SAMPLES, splitid=NUM_CHUNKS)

################################################################################
# Filtering
################################################################################
CHR_ALL = expand([out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.DNA.chr.bam", out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNA.chr.bam", out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNAr.chr.bam"], sample=ALL_SAMPLES, splitid=NUM_CHUNKS)


MASKED = expand(out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.DNA.chr.masked.bam", sample=ALL_SAMPLES, splitid=NUM_CHUNKS)

################################################################################
# Merge Alignments
################################################################################
MERGE_DNA = expand(out_dir + "workup/alignments/{sample}.DNA.merged.tagged.bam", sample = ALL_SAMPLES)

MERGE_RNA = expand(out_dir + "workup/annotations/{sample}.RNA.merged.anno.tagged.bam", sample = ALL_SAMPLES)

MERGE_RNAr = expand(out_dir + "workup/alignments/{sample}.RNAr.merged.tagged.bam", sample = ALL_SAMPLES)

# Bt2_DEDUP = expand(out_dir + "workup/alignments/{sample}.RNAr.merged.dedup.tagged.bam", sample=ALL_SAMPLES)


################################################################################
# Clustering
################################################################################

CLUSTERS = expand(out_dir + "workup/clusters/{sample}.bam", sample=ALL_SAMPLES)

CLUSTER_RENAME = expand(out_dir + "workup/clusters/{sample}.renamed.bam", sample=ALL_SAMPLES)

################################################################################
# Statistics
################################################################################

CLUSTER_STATISTICS = [out_dir + "workup/clusters/cluster_statistics.txt"]

CLUSTER_SIZES = [out_dir + "workup/clusters/DPM_read_distribution.pdf",
            out_dir + "workup/clusters/DPM_cluster_distribution.pdf",
            out_dir + "workup/clusters/RPM_read_distribution.pdf",
            out_dir + "workup/clusters/RPM_cluster_distribution.pdf"]

################################################################################
################################################################################
# RULE ALL
################################################################################
################################################################################

rule all:
    input: ALL_FASTQ + TRIM + TRIM_LOG + TRIM_RD + BARCODEID + LE_LOG_ALL + SPLIT_ALL + Bt2_DNA_ALIGN + Ht2_RNA_ALIGN + Bt2_TAG_ALL + CHR_ALL + MASKED + MERGE_DNA + MERGE_RNA + MERGE_RNAr + Ht2_ANNO_RNA + RNA_COMBINE + CLUSTERS + MULTI_QC + CLUSTER_RENAME + CLUSTER_STATISTICS + CLUSTER_SIZES


# Send and email if an error occurs during execution
onerror:
    shell('mail -s "an error occurred" ' + email + ' < {log}')

wildcard_constraints:
    sample = "[^\.]+"


####################################################################################################
# Trimming and barcode identification
####################################################################################################

rule splitfq_R1:
    input:
        lambda wildcards: FILES[wildcards.sample]['R1']
    output:
        expand(out_dir + "workup/splitfq/{{sample}}_R1.part_{splitid}.fastq.gz",
             splitid=NUM_CHUNKS)
    params:
        dir = out_dir + "workup/splitfq/",
        prefix = "{sample}_R1.part_0",
    log:
        out_dir + "workup/logs/{sample}.splitfq_r1.log"
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        {{
            mkdir -p "{params.dir}"
            bash "{split_fastq}" "{input}" {num_chunks} "{params.dir}" "{params.prefix}" {threads}
        }} &> "{log}"
        '''

# Split fastq files into chunks to processes in parallel
rule splitfq_R2:
    input:
        lambda wildcards: FILES[wildcards.sample]['R2']
    output:
        expand(out_dir + "workup/splitfq/{{sample}}_R2.part_{splitid}.fastq.gz",
             splitid=NUM_CHUNKS)
    params:
        dir = out_dir + "workup/splitfq/",
        prefix = "{sample}_R2.part_0",
    log:
        out_dir + "workup/logs/{sample}.splitfq_r2.log"
    conda:
        conda_env
    threads:
        4
    shell:
        '''
        {{
            mkdir -p "{params.dir}"
            bash "{split_fastq}" "{input}" {num_chunks} "{params.dir}" "{params.prefix}" {threads}
        }} &> "{log}"
        '''

# Trim adaptors
# multiple cores requires pigz to be installed on the system
rule adaptor_trimming_pe:
    input:
        [out_dir + "workup/splitfq/{sample}_R1.part_{splitid}.fastq.gz",
         out_dir + "workup/splitfq/{sample}_R2.part_{splitid}.fastq.gz"]
    output:
         out_dir + "workup/trimmed/{sample}_R1.part_{splitid}_val_1.fq.gz",
         out_dir + "workup/trimmed/{sample}_R1.part_{splitid}.fastq.gz_trimming_report.txt",
         out_dir + "workup/trimmed/{sample}_R2.part_{splitid}_val_2.fq.gz",
         out_dir + "workup/trimmed/{sample}_R2.part_{splitid}.fastq.gz_trimming_report.txt"
    threads:
        10
    log:
        out_dir + "workup/logs/{sample}.{splitid}.trim_galore.log"
    conda:
        conda_env
    shell:
        '''
        if [[ {threads} -gt 8 ]]
        then
            cores=2
        else
            cores=1
        fi

        trim_galore \
        --paired \
        --gzip \
        --cores $cores \
        --quality 20 \
        --fastqc \
        --fastqc_args "-t 4" \
        -o {out_dir}workup/trimmed/ \
        {input} &> {log}
        '''

# Cutadapt
rule cutadapt:
    '''
    Trim DPM RPM if read through reads
    RPM from right ATCAGCACTTA
    DPM from right GATCGGAAGAG
    DPM from left GGTGGTCTT ^ anchored (only appears at the start of read)
    DPM5bot2-B1    /5Phos/TGACTTGTCATGTCTTCCGATCTGGTGGTCTTT
    DPM5bot3-C1    /5Phos/TGACTTGTCATGTCTTCCGATCTGCCTCTTGTT
    DPM5bot26-B4   /5Phos/TGACTTGTCATGTCTTCCGATCTCCAGGTATTT
    DPM5bot44-D6   /5Phos/TGACTTGTCATGTCTTCCGATCTTAAGAGAGTT
    DPM5bot85-E11  /5Phos/TGACTTGTCATGTCTTCCGATCTTTCTCCTCTT
    DPM5bot95-G12  /5Phos/TGACTTGTCATGTCTTCCGATCTACCCTCGATT
    '''
    input:
        [out_dir + "workup/trimmed/{sample}_R1.part_{splitid}_val_1.fq.gz", 
        out_dir + "workup/trimmed/{sample}_R2.part_{splitid}_val_2.fq.gz"]
    output:
        fastq1=out_dir + "workup/trimmed/{sample}_R1.part_{splitid}_val_1_RDtrim.fq.gz",
        fastq2=out_dir + "workup/trimmed/{sample}_R2.part_{splitid}_val_2_RDtrim.fq.gz",
        qc=out_dir + "workup/trimmed/{sample}.part_{splitid}.RDtrim.qc.txt"
    threads: 10
    params:
        adapters_r1 = "-a GATCGGAAGAG -a ATCAGCACTTA " + adaptors,
        adapters_r2 = "",
        others = "--minimum-length 20"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.cutadapt.log"
    conda:
        conda_env
    shell:
        '''
        (cutadapt \
         {params.adapters_r1} \
         {params.adapters_r2} \
         {params.others} \
         -o {output.fastq1} \
         -p {output.fastq2} \
         -j {threads} \
         {input} > {output.qc}) &> {log}

        fastqc {output.fastq1}
        '''

# Identify barcodes using BarcodeIdentification_v1.2.0.jar
rule barcode_id:
    input:
        r1 = out_dir + "workup/trimmed/{sample}_R1.part_{splitid}_val_1_RDtrim.fq.gz",
        r2 = out_dir + "workup/trimmed/{sample}_R2.part_{splitid}_val_2_RDtrim.fq.gz"
    output:
        r1_barcoded = out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded.fastq.gz",
        r2_barcoded = out_dir + "workup/fastqs/{sample}_R2.part_{splitid}.barcoded.fastq.gz"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.bID.log"
    shell:
        "java -jar {barcode_id_jar} \
        --input1 {input.r1} --input2 {input.r2} \
        --output1 {output.r1_barcoded} --output2 {output.r2_barcoded} \
        --config {bid_config} &> {log}"


# Get ligation efficiency
rule get_ligation_efficiency:
    input:
        r1 = out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded.fastq.gz" 
    output:
        temp(out_dir + "workup/{sample}.part_{splitid}.ligation_efficiency.txt")
    conda:
        conda_env
    shell:
        "python {lig_eff} {input.r1} > {output}"


rule cat_ligation_efficiency:
    input:
        expand(out_dir + "workup/{sample}.part_{splitid}.ligation_efficiency.txt", sample=ALL_SAMPLES, splitid=NUM_CHUNKS)
    output:
        out_dir + "workup/ligation_efficiency.txt"
    shell:
        "tail -n +1 {input} > {output}"



rule split_rpm_dpm:
    '''
    split rpm and dpm will also remove incomplete barcodes
    '''
    input:
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded.fastq.gz"
    output:
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_dpm.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_rpm.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_other.fastq.gz",
        out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_short.fastq.gz"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.RPM_DPM.log"
    shell:
        "python {split_fq} --r1 {input} --format {formatfile} &> {log}"



################################################################################
# DNA workup
################################################################################

rule bowtie2_align:
    '''
    MapQ filter 20, -F 4 only mapped reads, -F 256 remove not primary alignment reads
    -F: Do not output alignments with any bits set in INT present in the FLAG field
    '''
    input:
        fq=out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_dpm.fastq.gz"
    output:
        sorted = out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.DNA.bowtie2.mapq20.bam",
        bam = temp(out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.unsorted.bam")
    threads: 10
    log:
        out_dir + "workup/logs/{sample}.{splitid}.bowtie2.log"
    conda:
        conda_env
    shell:
        '''
        (bowtie2 \
        -p 10 \
        -t \
        --phred33 \
        -x {bowtie2_index} \
        -U {input.fq} | \
        samtools view -bq 20 -F 4 -F 256 - > {output.bam}) &> {log}
        samtools sort -@ {threads} -o {output.sorted} {output.bam}
        '''

################################################################################
# RNA alignment
################################################################################

# RNA alignment to genome, no soft clipping and high penalty
rule hisat2_align:
    input:
        fq=out_dir + "workup/fastqs/{sample}_R1.part_{splitid}.barcoded_rpm.fastq.gz"
    output:
        all_reads=temp(out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNA.hisat2.bam"),
        mapped=temp(out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNA.hisat2.mapq20.unsorted.bam"),
        sorted = out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNA.hisat2.mapq20.bam",
        merged=out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNA.hisat2.unmapped.lowmq.bam",
        fq_gz=out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNA.hisat2.unmapped.lowmq.fq.gz"
    params:
        fq=out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNA.hisat2.unmapped.lowmq.fq",
        dir=out_dir + "workup/alignments/"
    threads: 10
    conda:
        conda_env
    log:
        out_dir + "workup/logs/{sample}.{splitid}.hisat2.log"
    shell:
        '''
        (hisat2 --end-to-end \
        -p 10 \
        -t \
        --phred33 \
        --known-splicesite-infile {hisat2_ss} \
        -x {hisat2_index} \
        -U {input.fq} | \
        samtools view -b -F 256 - > {output.all_reads}) &> {log}
        #split out unmapped and low mapq reads for realignment to repeats
        samtools view -bq 20 -U {output.merged} -F 4 {output.all_reads} > {output.mapped}
        samtools sort -@ {threads} -o {output.sorted} {output.mapped}
        mkdir -p {params.dir}
        samtools bam2fq -@ {threads} {output.merged} > {params.fq}
        pigz {params.fq}
        '''

################################################################################
# Repeats alignment
################################################################################

rule bowtie2_align_rna_repeats:
    '''
    Align RNA with Bowtie2 to repeats
    MapQ filter 20, -F 4 only mapped reads, -F 256 remove not primary alignment reads
    -F: Do not output alignments with any bits set in INT present in the FLAG field
    NO mapq20 because of repeat reads -> mapq20 only to get unique reads
    '''
    input:
        fq = out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNA.hisat2.unmapped.lowmq.fq.gz"
    output:
        all_reads = temp(out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNA.bowtie2.all.bam"),
        mapped = temp(out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNAr.bowtie2.unsorted.bam"),
        sorted = out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNAr.bowtie2.bam",
        unmapped = out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNAr.bowtie2.unmapped.bam"
    log:
        out_dir + "workup/logs/{sample}.{splitid}.RNAr.bowtie2.log"
    threads: 10
    conda:
        conda_env
    shell:
        '''
        (bowtie2 \
        -p 10 \
        -t \
        --phred33 \
        --very-sensitive \
        -x {bowtie2_repeat_index} \
        -U {input.fq} | \
        samtools view -bS -F 256 - > {output.all_reads}) &> {log}
        #Split out unmapped to keep just in case
        samtools view -b -U {output.unmapped} -F 4 {output.all_reads} > {output.mapped}
        samtools sort -@ {threads} -o {output.sorted} {output.mapped}
        '''

# Add chromosome of repeats as tag the featureCounts uses
rule add_tags_bowtie2:
    input:
        out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNAr.bowtie2.bam"
    output:
        out=out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNAr.bowtie2.tag.bam",
        idx=temp(out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNAr.bowtie2.bam.bai"),
    log:
        out_dir + "workup/logs/{sample}.{splitid}.RNAr.bowtie2.tag.log"
    threads:
        10
    conda:
        conda_env
    shell:
        '''
        samtools index {input}

        python {atttb} -i {input} -o {output.out}
        '''

###############################################################################
# Filter and chr
###############################################################################

rule add_chr_RNA_DNA:
    input:
        dpm=out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.DNA.bowtie2.mapq20.bam",
        rpm=out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNA.hisat2.mapq20.bam",
        rpm_repeat= out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNAr.bowtie2.tag.bam"
    output:
        dpm=out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.DNA.chr.bam",
        rpm=out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNA.chr.bam",
        rpm_repeat=out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.RNAr.chr.bam"
    log:
        dpm=out_dir + "workup/logs/{sample}.{splitid}.DNA_bcs.log",
        rpm=out_dir + "workup/logs/{sample}.{splitid}.RNA_bcs.log",
        rpm_repeat=out_dir + "workup/logs/{sample}.{splitid}.RNAr_bcs.log"
    conda:
        conda_env
    shell:
        '''
        python {add_chr} -i {input.dpm} -o {output.dpm} --assembly {assembly} &> {log.dpm}
        python {add_chr} -i {input.rpm} -o {output.rpm} --assembly {assembly} &> {log.rpm}
        python {add_chr} -i {input.rpm_repeat} -o {output.rpm_repeat} --assembly none &> {log.rpm_repeat}
        '''


rule repeat_mask:
    input:
        out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.DNA.chr.bam"
    output:
        out_dir + "workup/alignments_pieces/{sample}.part_{splitid}.DNA.chr.masked.bam"
    conda:
        conda_env
    shell:
        '''
        bedtools intersect -v -a {input} -b {mask} > {output}
        '''

################################################################################
# Merge and Deduplicate
################################################################################

rule merge_dna:
    input:
        expand(out_dir + "workup/alignments_pieces/{{sample}}.part_{splitid}.DNA.chr.masked.bam", splitid=NUM_CHUNKS)
    output:
        out = temp(out_dir + "workup/alignments/{sample}.DNA.merged.bam"),
        tagged = out_dir + "workup/alignments/{sample}.DNA.merged.tagged.bam"
    conda:
        conda_env
    threads:
        8
    log:
        out_dir + "workup/logs/{sample}.merge_DNA.log"
    shell:
        '''
        (samtools merge -@ {threads} {output.out} {input}) &> {log}
        python "{tag_bam}" --input_bam "{output.out}" --output_bam "{output.tagged}" --num_tags "{num_tags}"
        '''

rule merge_rna:
    input:
        expand(out_dir + "workup/alignments_pieces/{{sample}}.part_{splitid}.RNA.chr.bam", splitid=NUM_CHUNKS)
    output:
        out_dir + "workup/alignments/{sample}.RNA.merged.bam"
    conda:
        conda_env
    threads:
        8
    log:
        out_dir + "workup/logs/{sample}.merge_RNA.log"
    shell:
        '''
        (samtools merge -@ {threads} {output} {input}) &> {log}
        '''

# deduplication is included in tag_bam step
rule merge_rrna:
    input:
        expand(out_dir + "workup/alignments_pieces/{{sample}}.part_{splitid}.RNAr.chr.bam", splitid=NUM_CHUNKS)
    output:
        out = temp(out_dir + "workup/alignments/{sample}.RNAr.merged.bam"),
        tagged = out_dir + "workup/alignments/{sample}.RNAr.merged.tagged.bam"
    conda:
        conda_env
    threads:
        8
    log:
        out_dir + "workup/logs/{sample}.merge_RNAr.log"
    shell:
        '''
        (samtools merge -@ {threads} {output.out} {input}) &> {log}
        python "{tag_bam}" --input_bam "{output.out}" --output_bam "{output.tagged}" --num_tags "{num_tags}"
        '''

################################################################################
# RNA annotation
################################################################################

# Annotate rna in three rounds, allow multimapping, stranded and antisense annotations
rule annotate_rna:
    input:
        out_dir + "workup/alignments/{sample}.RNA.merged.bam"
    threads: 10
    output:
        bam_exon_intron=out_dir + "workup/annotations/{sample}.RNAexin.hisat2.mapq20.bam.featureCounts.bam",
        counts_exon_intron=out_dir + "workup/annotations/{sample}.RNAexin.hisat2.mapq20.bam.featureCounts.txt",
        bam_rrna=out_dir + "workup/annotations/{sample}.RNAr.hisat2.mapq20.bam.featureCounts.bam",
        counts_rrna=out_dir + "workup/annotations/{sample}.RNAr.hisat2.mapq20.bam.featureCounts.txt",
        bam_antisense=out_dir + "workup/annotations/{sample}.RNAanti.hisat2.mapq20.bam.featureCounts.bam",
        counts_antisense=out_dir + "workup/annotations/{sample}.RNAanti.hisate2.mapq20.bam.featureCounts.txt"
    log:
        out_dir + "workup/logs/{sample}.anno.log"
    params:
        ex_rename = out_dir + "workup/annotations/{sample}.RNAexin.hisat2.mapq20.bam",
        r_rename = out_dir + "workup/annotations/{sample}.RNAr.hisat2.mapq20.bam",
        anti_rename  = out_dir + "workup/annotations/{sample}.RNAanti.hisat2.mapq20.bam"
    conda:
        conda_env
    shell:
        '''
        mv {input} {params.ex_rename}
        featureCounts -T {threads}  \
        -R BAM -M -O -s 1 \
        -F SAF -a {exon_intron_saf} \
        -o {output.counts_exon_intron} \
        {params.ex_rename}

        mv {params.ex_rename} {params.r_rename}
        featureCounts -T {threads} \
        -R BAM -M -O -s 1 --fracOverlap 0.85 \
        -F SAF -a {anno_repeats_saf} \
        -o {output.counts_rrna} \
        {params.r_rename}

        mv {params.r_rename} {params.anti_rename}
        featureCounts -T {threads} \
        -R BAM -M -O -s 2 \
        -F SAF -a {antisense_saf} \
        -o {output.counts_antisense} \
        {params.anti_rename}

        mv {params.anti_rename} {input}

        '''

rule combine_annotations_rna:
    input:
        bam=out_dir + "workup/alignments/{sample}.RNA.merged.bam",
        bam_exon=out_dir + "workup/annotations/{sample}.RNAexin.hisat2.mapq20.bam.featureCounts.bam",
        bam_rrna=out_dir + "workup/annotations/{sample}.RNAr.hisat2.mapq20.bam.featureCounts.bam",
        bam_anti=out_dir + "workup/annotations/{sample}.RNAanti.hisat2.mapq20.bam.featureCounts.bam"
    output:
        out_dir + "workup/annotations/{sample}.RNA.merged.anno.bam"
    log:
        out_dir + "workup/logs/{sample}.RNA_anno_combine.log"
    conda:
        conda_env
    shell:
        "python {comb_anno} -i {input.bam_exon} {input.bam_rrna} {input.bam_anti} \
            -i2 {input.bam} \
            -o {output} &> {log}"

rule tag_annotated_rna:
    input:
        out_dir + "workup/annotations/{sample}.RNA.merged.anno.bam",
    output:
        out_dir + "workup/annotations/{sample}.RNA.merged.anno.tagged.bam",
    log:
        out_dir + "workup/logs/{sample}.RNA_anno_tag.log"
    conda:
        conda_env
    shell:
        """
        python "{tag_bam}" --input_bam "{input}" --output_bam "{output}" --num_tags "{num_tags}"
        """

################################################################################
# Make clusters
################################################################################

rule merge_samp:
    input:
        [out_dir + "workup/alignments/{sample}.DNA.merged.tagged.bam",
        out_dir + "workup/annotations/{sample}.RNA.merged.anno.tagged.bam",
        out_dir + "workup/alignments/{sample}.RNAr.merged.tagged.bam"]
    output:
        out_dir + "workup/clusters/{sample}.bam"
    log:
        out_dir + "workup/clusters/{sample}.make_clusters.log"
    conda:
        conda_env
    shell:
        """
        samtools merge -@ {threads} "{output}" {input} &> "{log}"
        """

rule rename_clusters:
    input:
       out_dir + "workup/clusters/{sample}.bam"
    output:
       out_dir + "workup/clusters/{sample}.renamed.bam"
    shell:
       "sed -f {anno_rename} {input} > {output}"

################################################################################
# MultiQC
################################################################################

rule multiqc:
    input:
        expand([out_dir + "workup/clusters/{sample}.bam"], sample=ALL_SAMPLES) 
    output:
        out_dir + "workup/qc/multiqc_report.html"
    log:
        out_dir + "workup/logs/multiqc.log"
    conda:
        conda_env
    shell:
        "multiqc {out_dir}workup -o {out_dir}workup/qc"

################################################################################
# Generate Statistics
################################################################################

rule generate_statistics:
    input:
        expand([out_dir + "workup/clusters/{sample}.bam"], sample=ALL_SAMPLES)
    output:
        CLUSTER_STATISTICS + CLUSTER_SIZES
    log:
        out_dir + "workup/logs/generate_statistics.log"
    conda:
        conda_env
    shell: 
        """ 
        python {generate_statistics} --directory '{out_dir}workup/clusters' --pattern '.bam'
        """
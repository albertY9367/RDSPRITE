import glob
import os
from tag_bam import *

"""
loc = "/central/scratchio/mblanco/albert/MegaSprite/*.bam"
bams = glob.glob(loc)
print("Total bams: " + str(len(bams)))
counter = 0

new_dir = "/central/scratchio/mblanco/albert/MegaSpriteClustersNew/"

for bam in bams:

    counter += 1
    out_bam = bam.split(".")[:-1]
    out_bam.append("labeled.bam")
    out_bam = ".".join(out_bam)

    out_bam = os.path.join(new_dir, os.path.basename(out_bam))

    if os.path.isfile(out_bam):
        print("Already done " + str(counter))
        pass
    else:
        print("Doing " + str(counter))
        label_bam_file(bam, out_bam, 5)
"""


dna_loc = "/central/scratchio/mblanco/albert/MegaSpriteClustersNew/*.DNA.merged.labeled.bam"
dnas = glob.glob(dna_loc)

counter = 0

for dna in dnas:
    
    basename = dna.split(".")[0]
    rna = basename + ".RNA.merged.anno.labeled.bam"
    rnar = basename + ".RNAr.merged.dedup.labeled.bam"

    output = basename + ".bam"

    if os.path.isfile(output):
        pass
    else:
        os.system(f"samtools merge -o '{output}' {dna} {rna} {rnar}")

    counter += 1

    print(f"finished {counter}")
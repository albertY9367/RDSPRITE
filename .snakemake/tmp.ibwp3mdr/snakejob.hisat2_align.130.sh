#!/bin/sh
# properties = {"type": "single", "rule": "hisat2_align", "local": false, "input": ["/central/scratchio/mblanco/albert/sprite_test_output/workup/fastqs/DPM_first_R1.part_006.barcoded_rpm.fastq.gz"], "output": ["/central/scratchio/mblanco/albert/sprite_test_output/workup/alignments_pieces/DPM_first.part_006.RNA.hisat2.bam", "/central/scratchio/mblanco/albert/sprite_test_output/workup/alignments_pieces/DPM_first.part_006.RNA.hisat2.mapq20.unsorted.bam", "/central/scratchio/mblanco/albert/sprite_test_output/workup/alignments_pieces/DPM_first.part_006.RNA.hisat2.mapq20.bam", "/central/scratchio/mblanco/albert/sprite_test_output/workup/alignments_pieces/DPM_first.part_006.RNA.hisat2.unmapped.lowmq.bam", "/central/scratchio/mblanco/albert/sprite_test_output/workup/alignments_pieces/DPM_first.part_006.RNA.hisat2.unmapped.lowmq.fq.gz"], "wildcards": {"sample": "DPM_first", "splitid": "006"}, "params": {"fq": "/central/scratchio/mblanco/albert/sprite_test_output/workup/alignments_pieces/DPM_first.part_006.RNA.hisat2.unmapped.lowmq.fq", "dir": "/central/scratchio/mblanco/albert/sprite_test_output/workup/alignments/"}, "log": ["/central/scratchio/mblanco/albert/sprite_test_output/workup/logs/DPM_first.006.hisat2.log"], "threads": 10, "resources": {"mem_mb": 1000, "mem_mib": 954, "disk_mb": 1000, "disk_mib": 954, "tmpdir": "<TBD>"}, "jobid": 130, "cluster": {"time": "12:00:00", "mem": "50g", "cpus": 10, "nodes": 1, "output": "workup/logs/cluster/hisat2_align.DPM_first.006.out", "error": "workup/logs/cluster/hisat2_align.DPM_first.006.err"}}
cd '/central/home/zyang4/RDSPRITE_v4' && /home/zyang4/anaconda3/envs/snakemake/bin/python3.11 -m snakemake --snakefile '/central/home/zyang4/RDSPRITE_v4/Snakefile' --target-jobs 'hisat2_align:sample=DPM_first,splitid=006' --allowed-rules 'hisat2_align' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=1000' 'mem_mib=954' 'disk_mb=1000' 'disk_mib=954' --wait-for-files '/central/home/zyang4/RDSPRITE_v4/.snakemake/tmp.ibwp3mdr' '/central/scratchio/mblanco/albert/sprite_test_output/workup/fastqs/DPM_first_R1.part_006.barcoded_rpm.fastq.gz' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'input' 'code' 'params' 'software-env' 'mtime' --skip-script-cleanup  --use-conda  --conda-frontend 'mamba' --conda-base-path '/home/zyang4/anaconda3' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --configfiles '/central/home/zyang4/RDSPRITE_v4/config.yaml' --latency-wait 5 --scheduler 'ilp' --scheduler-solver-path '/home/zyang4/anaconda3/envs/snakemake/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && touch '/central/home/zyang4/RDSPRITE_v4/.snakemake/tmp.ibwp3mdr/130.jobfinished' || (touch '/central/home/zyang4/RDSPRITE_v4/.snakemake/tmp.ibwp3mdr/130.jobfailed'; exit 1)


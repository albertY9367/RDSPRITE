#!/bin/sh
# properties = {"type": "single", "rule": "add_tags_bowtie2", "local": false, "input": ["/central/scratchio/mblanco/albert/sprite_test_output/workup/alignments_pieces/DPM_first.part_001.RNAr.bowtie2.bam"], "output": ["/central/scratchio/mblanco/albert/sprite_test_output/workup/alignments_pieces/DPM_first.part_001.RNAr.bowtie2.tag.bam", "/central/scratchio/mblanco/albert/sprite_test_output/workup/alignments_pieces/DPM_first.part_001.RNAr.bowtie2.bam.bai"], "wildcards": {"sample": "DPM_first", "splitid": "001"}, "params": {}, "log": ["/central/scratchio/mblanco/albert/sprite_test_output/workup/logs/DPM_first.001.RNAr.bowtie2.tag.log"], "threads": 10, "resources": {"mem_mb": 1000, "mem_mib": 954, "disk_mb": 1000, "disk_mib": 954, "tmpdir": "<TBD>"}, "jobid": 146, "cluster": {"time": "12:00:00", "mem": "50g", "cpus": 10, "nodes": 1, "output": "workup/logs/cluster/add_tags_bowtie2.DPM_first.001.out", "error": "workup/logs/cluster/add_tags_bowtie2.DPM_first.001.err"}}
cd '/central/home/zyang4/RDSPRITE_v4' && /home/zyang4/anaconda3/envs/snakemake/bin/python3.11 -m snakemake --snakefile '/central/home/zyang4/RDSPRITE_v4/Snakefile' --target-jobs 'add_tags_bowtie2:sample=DPM_first,splitid=001' --allowed-rules 'add_tags_bowtie2' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=1000' 'mem_mib=954' 'disk_mb=1000' 'disk_mib=954' --wait-for-files '/central/home/zyang4/RDSPRITE_v4/.snakemake/tmp.ibwp3mdr' '/central/scratchio/mblanco/albert/sprite_test_output/workup/alignments_pieces/DPM_first.part_001.RNAr.bowtie2.bam' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'input' 'code' 'params' 'software-env' 'mtime' --skip-script-cleanup  --use-conda  --conda-frontend 'mamba' --conda-base-path '/home/zyang4/anaconda3' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --configfiles '/central/home/zyang4/RDSPRITE_v4/config.yaml' --latency-wait 5 --scheduler 'ilp' --scheduler-solver-path '/home/zyang4/anaconda3/envs/snakemake/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && touch '/central/home/zyang4/RDSPRITE_v4/.snakemake/tmp.ibwp3mdr/146.jobfinished' || (touch '/central/home/zyang4/RDSPRITE_v4/.snakemake/tmp.ibwp3mdr/146.jobfailed'; exit 1)


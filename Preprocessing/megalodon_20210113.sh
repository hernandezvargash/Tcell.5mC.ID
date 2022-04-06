 #! /bin/sh.

######################################
# Mar 2022

# NVIDIA-SMI 510.39.01
# CUDA Version: 11.6
# re-installed tensorflow with pip
# check version with: pip list | grep tensorflow
# tensorflow                         2.7.0
# tensorflow-estimator               2.7.0
# tensorflow-io-gcs-filesystem       0.23.1

# python                    3.7.4
# numpy                     1.21.4

# guppy_basecaller -v : Version 6.1.1
# minimap2 version 2.22

# megalodon                 2.5.0
# ont-pyguppy-client-lib    6.1.1
# ont-remora		    1.0.0

# fresh Rerio install


# demultiplex after re-basecalling with guppy
######################################

# re-basecall

cd '/mnt/sdb/projects/nano/Tcell_ID/20210113/'

guppy_basecaller --input_path fast5 --save_path demultiplexed_fast5s --flowcell FLO-MIN106 --kit SQK-LSK109 --barcode_kits EXP-NBD104 --fast5_out --nested_output_folder --device cuda:0
#--fast5_out parameter is not necessary


# demux_fast5 (using ont_fast5_api)

cd '/mnt/sdb/projects/nano/Tcell_ID/20210113/'

demux_fast5 --input /mnt/sdb/projects/nano/Tcell_ID/20210113/fast5 --save_path /mnt/sdb/projects/nano/Tcell_ID/20210113/demultiplexed_fast5s/demultiplexed_reads --summary_file /mnt/sdb/projects/nano/Tcell_ID/20210113/demultiplexed_fast5s/sequencing_summary.txt



# 5hmC & 5mC in CpG context (all context still not available for Remora)
######################################


conda activate megalodon

OUTDIR=/mnt/sdb/results/nano/Tcell_ID/20210113
FAST5PATH=/mnt/sdb/projects/nano/Tcell_ID/20210113/demultiplexed_fast5s/demultiplexed_reads/
GUPPY=/opt/ont/guppy/bin/guppy_basecall_server
REF=/mnt/sdb/refs/GRCm38.p6.genome.fa
CONFIG=dna_r9.4.1_450bps_fast.cfg

#!/bin/sh
for barcode in $(find $FAST5PATH -mindepth 1 -maxdepth 1 -type d -printf '%f\n')
do 
	megalodon $FAST5PATH/$barcode/ \
        --guppy-server-path $GUPPY \
        --guppy-config $CONFIG \
        --remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5hmc_5mc CG 0 \
        --outputs basecalls mods per_read_mods mod_basecalls mappings mod_mappings mods \
        --reference $REF \
        --devices 0 \
        --processes 8 \
        --output-directory ${OUTDIR}/${barcode}
done

conda deactivate


######################################






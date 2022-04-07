 #! /bin/sh.

######################################
# Apr 2022


# re-basecall
######################################

cd '/mnt/sdb/projects/nano/Tcell_ID/20201210/no_sample/20201210_1408_MN32885_FAO95644_2de3ca5f/'

guppy_basecaller --input_path fast5 --save_path demultiplexed_fast5s --flowcell FLO-MIN106 --kit SQK-LSK109 --barcode_kits EXP-NBD104 --nested_output_folder --device cuda:0
#--fast5_out parameter is not necessary


# demux_fast5 (using ont_fast5_api)
######################################

demux_fast5 --input /mnt/sdb/projects/nano/Tcell_ID/20201210/no_sample/20201210_1408_MN32885_FAO95644_2de3ca5f/fast5 --save_path /mnt/sdb/projects/nano/Tcell_ID/20201210/no_sample/20201210_1408_MN32885_FAO95644_2de3ca5f/demux --summary_file /mnt/sdb/projects/nano/Tcell_ID/20201210/no_sample/20201210_1408_MN32885_FAO95644_2de3ca5f/demultiplexed_fast5s/sequencing_summary.txt



# 5hmC & 5mC in CpG context (all context still not available for Remora)
######################################

conda activate megalodon

OUTDIR=/mnt/sdb/results/nano/Tcell_ID/20201210
FAST5PATH=/mnt/sdb/projects/nano/Tcell_ID/20201210/no_sample/20201210_1408_MN32885_FAO95644_2de3ca5f/demux/
GUPPY=/opt/ont/guppy/bin/guppy_basecall_server
REF=/mnt/sdb/refs/GRCm38.p6.genome.fa
CONFIG=dna_r9.4.1_450bps_fast.cfg

for barcode in $(basename -a $(ls -d $FAST5PATH/*/))
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










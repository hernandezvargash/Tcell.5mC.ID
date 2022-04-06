
######################################
# Apr 2022


# fast5 already demultiplexed


# 5hmC & 5mC in CpG context (all context still not available for Remora)
######################################


conda activate megalodon

OUTDIR=/mnt/sdb/results/nano/Tcell_ID/20210312
FAST5PATH=/mnt/sdb/projects/nano/Tcell_ID/20210312/Tcells/20210312_1659_MN32885_FAP75660_eeee68e1/fast5_pass/
GUPPY=/opt/ont/guppy/bin/guppy_basecall_server
REF=/mnt/sdb/refs/GRCm38.p6.genome.fa
CONFIG=dna_r9.4.1_450bps_fast.cfg

#!/bin/sh
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










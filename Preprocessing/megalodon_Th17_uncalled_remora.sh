
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



# 5hmC & 5mC in CpG context (all context still not available for Remora)

######################################
# 20210419
######################################

conda activate megalodon

OUTDIR=/mnt/sdb/results/nano/Tcell_ID/Th17/WT_remora
FAST5PATH=/mnt/sdb/projects/nano/Tcell_ID/Th17/uncalled/20210419/fast5
GUPPY=/opt/ont/guppy/bin/guppy_basecall_server
REF=/mnt/sdb/refs/GRCm38.p6.genome.fa
CONFIG=dna_r9.4.1_450bps_fast.cfg

megalodon \
    $FAST5PATH \
    --guppy-server-path $GUPPY \
    --guppy-config $CONFIG \
    --remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5hmc_5mc CG 0 \
    --outputs basecalls mappings mod_mappings mods \
    --reference $REF \
    --output-directory $OUTDIR \
    --overwrite \
    --devices 0 \
    --processes 8

conda deactivate


######################################
# 20210402
######################################

conda activate megalodon

OUTDIR=/mnt/sdb/results/nano/Tcell_ID/Th17/KO_remora
FAST5PATH=/mnt/sdb/projects/nano/Tcell_ID/Th17/uncalled/20210402/fast5
GUPPY=/opt/ont/guppy/bin/guppy_basecall_server
REF=/mnt/sdb/refs/GRCm38.p6.genome.fa
CONFIG=dna_r9.4.1_450bps_fast.cfg

megalodon \
    $FAST5PATH \
    --guppy-server-path $GUPPY \
    --guppy-config $CONFIG \
    --remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5hmc_5mc CG 0 \
    --outputs basecalls mappings mod_mappings mods \
    --reference $REF \
    --output-directory $OUTDIR \
    --overwrite \
    --devices 0 \
    --processes 8

conda deactivate


######################################





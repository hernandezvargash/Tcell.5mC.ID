 #! /bin/sh.

######################################

# System Info:
# Memory: 125.5 GiB
# Processor: Intel Xeon(R)W-2175 CPU@2.50GHzx28
# Graphics: Quadro RTX 6000/PCIe/SSE2
# Ubuntu 18.04.6 LTS

# NVIDIA-SMI 510.39.01
# CUDA Version: 11.6
# re-installed tensorflow with pip
# check version with: pip list | grep tensorflow
# tensorflow                         2.7.0
# tensorflow-estimator               2.7.0
# tensorflow-io-gcs-filesystem       0.23.1

# python                    3.7.4
# numpy                     1.21.4

# guppy_basecaller -v : Version 6.1.7
# minimap2 version 2.22

# megalodon                 2.5.0 (installed in conda environment)
# ont-pyguppy-client-lib    6.1.1
# ont-remora		    1.0.0

# check available remora models:
# conda activate megalodon
# remora model list_pretrained


# 5hmC & 5mC in CpG context (all context still not available for Remora)
######################################

cd '/mnt/sdb/projects/nano/mouse_brain/'

conda activate megalodon

OUTDIR=/mnt/sdb/results/nano/mouse_brain/20221215
FAST5PATH=/mnt/sdb/projects/nano/mouse_brain/20210115_1044_MN32885_FAO34395_4f44cddc/fast5
GUPPY=/opt/ont/guppy/bin/guppy_basecall_server
REF=/mnt/sdb/refs/GRCm38.p6.genome.fa
CONFIG=dna_r9.4.1_450bps_fast.cfg

megalodon $FAST5PATH \
        --guppy-server-path $GUPPY \
        --guppy-config $CONFIG \
        --remora-modified-bases dna_r9.4.1_e8 hac 0.0.0 5hmc_5mc CG 0 \
        --outputs basecalls mods per_read_mods mod_basecalls mappings mod_mappings mods \
        --reference $REF \
        --devices 0 \
        --processes 12 \
        --output-directory ${OUTDIR}

conda deactivate


######################################


























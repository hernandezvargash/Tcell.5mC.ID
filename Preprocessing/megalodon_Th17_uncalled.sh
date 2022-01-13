
######################################
# 180521
######################################

conda activate megalodon

DIR=/mnt/sdb/results/Th17/WT
FAST5PATH=/mnt/sdb/projects/Th17/uncalled/20210419/fast5
REFERENCE=/mnt/sdb/refs/GRCm38.p6.genome.fa
RERIO_PATH=/mnt/sdb/scripts/rerio/basecall_models/
CONFIG="res_dna_r941_min_modbases_5mC_v001.cfg"

megalodon $FAST5PATH \
    --guppy-params "-d $RERIO_PATH" \
    --guppy-config $CONFIG \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server  \
    --outputs basecalls mods per_read_mods mod_basecalls mappings mod_mappings  \
    --write-mod-log-probs \
    --write-mods-text \
    --mappings-format bam \
    --mod-motif m CN 0 \
    --sort-mappings \
    --mod-map-emulate-bisulfite \
    --mappings-format bam \
    --mod-map-base-conv C T \
    --mod-map-base-conv m C \
    --mod-binary-threshold 0.8 \
    --reference ${REFERENCE} \
    --mod-output-formats bedmethyl modvcf wiggle \
    --output-directory ${DIR} \
    --devices "cuda:all:100%" \
    --processes 8

conda deactivate


######################################


conda activate megalodon

DIR=/mnt/sdb/results/Th17/KO
FAST5PATH=/mnt/sdb/projects/Th17/uncalled/20210402/fast5
REFERENCE=/mnt/sdb/refs/GRCm38.p6.genome.fa
RERIO_PATH=/mnt/sdb/scripts/rerio/basecall_models/
CONFIG="res_dna_r941_min_modbases_5mC_v001.cfg"

megalodon $FAST5PATH \
    --guppy-params "-d $RERIO_PATH" \
    --guppy-config $CONFIG \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server  \
    --outputs basecalls mods per_read_mods mod_basecalls mappings mod_mappings  \
    --write-mod-log-probs \
    --write-mods-text \
    --mappings-format bam \
    --mod-motif m CN 0 \
    --sort-mappings \
    --mod-map-emulate-bisulfite \
    --mappings-format bam \
    --mod-map-base-conv C T \
    --mod-map-base-conv m C \
    --mod-binary-threshold 0.8 \
    --reference ${REFERENCE} \
    --mod-output-formats bedmethyl modvcf wiggle \
    --output-directory ${DIR} \
    --devices "cuda:all:100%" \
    --processes 8

conda deactivate


######################################






######################################
# Jan 2022
# megalodon v.2.4.1
# guppy_basecaller v.6.0.1
######################################

conda activate megalodon

DIR=/mnt/sdb/results/nano/Tcell_ID/20210113
FAST5PATH=/mnt/sdb/projects/nano/Tcell_ID/20210113/fast5
REFERENCE=/mnt/sdb/refs/GRCm38.p6.genome.fa
RERIO_PATH=/mnt/sdb/scripts/rerio/basecall_models/
CONFIG="res_dna_r941_min_modbases_5mC_5hmC_v001.cfg"

megalodon $FAST5PATH \
    --guppy-params "-d $RERIO_PATH" \
    --guppy-config $CONFIG \
    --guppy-server-path /opt/ont/guppy/bin/guppy_basecall_server  \
    --outputs basecalls mods per_read_mods mod_basecalls mappings mod_mappings  \
    --write-mod-log-probs \
    --write-mods-text \
    --mappings-format bam \
    --sort-mappings \
    --mod-map-emulate-bisulfite \
    --mappings-format bam \
    --mod-map-base-conv C T \
    --mod-map-base-conv m C \
    --reference ${REFERENCE} \
    --mod-output-formats bedmethyl modvcf wiggle \
    --output-directory ${DIR} \
    --devices "cuda:all:100%" \
    --processes 16

conda deactivate


######################################





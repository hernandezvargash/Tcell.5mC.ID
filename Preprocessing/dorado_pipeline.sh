#! /bin/sh.

# fresh Ubuntu 22.04 install
# dorado installation 061123
# dorado exectutable (dorado-0.4.2-linux-x64.tar.gz) was downloaded directly from https://github.com/nanoporetech/dorado#platforms

# NVIDIA Drivers for Quadro RTX6000
# 24GB
# CUDA v12
# monitor GPU with 'watch nvidia-smi'

# useful websites:
# https://github.com/nanoporetech/dorado
# https://github.com/nanoporetech/modkit
# https://github.com/nanoporetech/modkit/blob/master/book/src/advanced_usage.md
# https://labs.epi2me.io/how-to-mod/

##############################

# export paths
export PATH=/home/hh/Software/dorado-0.4.2-linux-x64/bin:$PATH
export PATH=/home/hh/Software/modkit_v0.2.2_centos7_x86_64/dist:$PATH

# definitions
FAST5PATH=fast5/
BAMPATH=demux/
BEDPATH=bed/
REF=/home/hh/Documents/references/GRCm38.p6.genome.fa
MODEL=/home/hh/Software/dorado-0.4.2-linux-x64/models/dna_r9.4.1_e8_hac@v3.3
KIT=EXP-NBD104

# convert to pod5 format
pod5 convert fast5 $FAST5PATH -r --output converted.pod5

# align during basecalling:
dorado basecaller $MODEL converted.pod5 --modified-bases 5mCG_5hmCG > calls.bam --kit-name $KIT --reference $REF

# demultiplex bam file
dorado demux --output-dir demux --no-classify calls.bam

dorado demux --output-dir demux_bam_aligned --no-classify calls_demux_aligned.bam

# sort and index with samtools (for IGV)
for bamfile in demux/*.bam  
 do  
   samtools sort ${bamfile} > ${bamfile}.sorted.bam   
   samtools index ${bamfile}.sorted.bam
   rm ${bamfile}
 done

# pileup with ModKit
# add the --only-tabs parameter to improve importing into R
mkdir bed
for bamfile in demux/*.bam  
 do  
   modkit pileup ${bamfile} $BEDPATH$(basename ${bamfile}).bed --log-filepath $BEDPATH$(basename ${bamfile}).log --cpg --ref $REF --combine-strands --only-tabs
 done


##############################




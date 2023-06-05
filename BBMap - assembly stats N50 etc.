#!/bin/bash
# BBMap Analysis with Hifiasm/Miniasm Output (primary contig assembly)
# This stats.sh allows you to calculate:
## G/C/A/T percentage
## GC %
## number of scaffold/contig total
## total scaffold/contig sequence
## N50/N90 of scaffold/contig
## max scaffold/contig length
## number of scaffolds
## %genome in scaffolds >50kB

# load BBMap
module spider bbmap
# select the latest version and load.
module load bbmap/38.86

# convert .gfa (output from Hifiasm) to .fa (fasta file for BBMap input)
awk '/^S/{print ">"$2"\n"$3}' filename.asm.bp.p_ctg.gfa | fold > filename.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' miniasm_Lingonberry_RedCandy_assembly.gfa | fold > miniasm_assembly.fasta

# Run stats.sh on your Hifiasm fasta file
stats.sh in=file.asm.bp.p_ctg.fa k=13 #addition of k somehow made it work.
stats.sh miniasm_assembly.fasta -Xmx4g > assembly.n50.txt

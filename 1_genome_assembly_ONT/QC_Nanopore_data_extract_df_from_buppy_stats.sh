#!/bin/bash

#Quality check of basecalled MinION data

#From guppy output: sequencing_summary.txt file
for file in `ls`;
do
  less -S $file | cut -f 14-15 > ${file}_basecalled_reads_stats.txt
done

#From guppy-duplex output: sequencing_summary.txt file
for file in `ls *duplex_by_guppy*`;
do
  less -S $file | cut -f 15-16 > ${file}_basecalled_reads_stats.txt
done

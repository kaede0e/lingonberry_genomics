## Genome assembly pipeline 

### sequence ###
Sequence raw reads (Nanopore MinION) \
Basecall with super accurate mode (Guppy) \
Duplex-basecall (Guppy_duplex) \
Combine the duplexed reads with sup singleplex reads --> combined_all_reads.fa

### de novo assenble ###
Assemble draft genome (smartdenovo) --> asm.fa \
Polish the draft assembly with Nanopore long-reads x3 (Nextpolish) \
Polish the draft assembly with Illumina short-reads x3 ([Pilon](https://github.com/broadinstitute/pilon/wiki)) --> polished_asm.fa

### assemble using reference ###
Scaffold draft genome with cranberry reference genome ([RagTag](https://github.com/malonge/RagTag/wiki)) \
Clean haplotigs with Illumina short-reads ([purge_haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/)), or \
Clean haplotigs with Nanopore long-reads (map raw reads onto draft genome and calculate freq kmer [Winnowmap](https://github.com/marbl/Winnowmap), then purge haplotigs based on read depth [purge_dups](https://github.com/dfguan/purge_dups))

### annotate ###
Gene annotation on the assembly using my RNA seq ...?
1. Run transcript assembler (Stringtie)
2. Align transcript to genome (Hisat or STAR)

Gene annotation on the assembly using published data (SRP110973; berry development transcriptome)
1. Download available data (paired Illumina reads) --> SRA__.fastq
2. Quality check public libraries (fastqc) 
3. Assemble transcript from multiple individuals (Stringtie)
4. Align transcript to genome (Hisat or STAR, make sure above 80% alignment rates)



## Genome assembly pipeline 

### sequence ###
Sequence raw reads (Nanopore MinION) \
Basecall with super accurate mode (Guppy) \
Duplex-basecall (Guppy_duplex) \
Combine the duplexed reads with sup singleplex reads --> combined_all_reads.fa

### de novo assenble ###
Assemble draft genome (smartdenovo) --> asm.fa \
Polish the draft assembly with Nanopore long-reads x3 (Nextpolish) \
Polish the draft assembly with Illumina short-reads x3 ([Pilon](https://github.com/broadinstitute/pilon/wiki) or [ntEdit](https://github.com/bcgsc/ntEdit)) --> polished_asm.fa

### assemble using reference ###
Scaffold draft genome with cranberry reference genome ([RagTag](https://github.com/malonge/RagTag/wiki)) \
Clean haplotigs with Illumina short-reads ([purge_haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/)), or \
Clean haplotigs with Nanopore long-reads (map raw reads onto draft genome and calculate freq kmer [Winnowmap](https://github.com/marbl/Winnowmap), then purge haplotigs based on read depth [purge_dups](https://github.com/dfguan/purge_dups))

### annotate ###
Gene annotation on the assembly using my RNA seq ...?
1. Adapter trimming of raw Illumina outputs. 
2. Run fastqc to make sure the data quality is good for downstream analysis - what do I do if the quality of reads are bad? 
3. Align RNA seq = cDNA library seq data, to the genome ([Hisat](http://daehwankimlab.github.io/hisat2/manual/) or STAR)
4. Run transcript assembler ([Stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)) --> .gtf which is basically .gff
5. 

Gene annotation on the assembly using published data (SRP110973; berry development transcriptome)
1. Download available data (paired Illumina reads) --> SRA__.fastq 
2. Quality check public libraries (fastqc) --> all good to use. 
3. Align transcript to genome (Hisat or STAR, make sure above 80% alignment rates) --> 95.5% alignment rate on Hisat2! 
4. Assemble transcript from multiple individuals (Stringtie) --> .gtf - this annotation contains "transcript" and "exon" 
5. Gene prediction: Cui et al. 2022 used [TransDecoder](https://github.com/TransDecoder/TransDecoder) and [GeneMark-ET](http://exon.gatech.edu/GeneMark/gmes_instructions.html)

Useful protein search resources (used in corporation with the BlastP search: [uniProt](https://www.uniprot.org/uniprotkb?facets=model_organism%3A3702&query=arabidopsis), it has Arabidopsis 136,000+ genes and 1,200+ Vaccinium genes). \
I downloaded them and created the blastdb. 


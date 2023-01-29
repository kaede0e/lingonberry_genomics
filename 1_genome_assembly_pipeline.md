## Genome assembly pipeline 

### sequence ###
Sequence raw reads (Nanopore MinION) \
Basecall with super accurate mode (Guppy) \
Duplex-basecall (Guppy_duplex) \
Combine the duplexed reads with sup singleplex reads --> combined_all_reads.fa

Sequence raw reads (Illumina Novaseq for both DNA & RNA) \
Data download from GenomeBC \
Clean reads \
Quality check the reads (fastqc)

### de novo assenble ###
Assemble draft genome (smartdenovo) --> asm.fa \
Correct reads and assemble draft genome (canu_smartdenovo) --> asm.fa \
Polish the draft assembly with Nanopore long-reads x3 (Nextpolish) \
Polish the draft assembly with Illumina short-reads x3 ([Pilon](https://github.com/broadinstitute/pilon/wiki) or [ntEdit](https://github.com/bcgsc/ntEdit)) --> polished_asm.fa

### assemble using reference ###
Clean haplotigs with Illumina short-reads or Nanopore long-reads ([purge_haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/); apparently it works best with the reads that were used to assemble the draft assembly), or \
Clean haplotigs with Nanopore long-reads (map raw reads onto draft genome and calculate freq kmer [Winnowmap](https://github.com/marbl/Winnowmap), then purge haplotigs based on read depth [purge_dups](https://github.com/dfguan/purge_dups)), or do both (purge_haplotigs and purge_dups curated better) \
Remove any contaminants or organeller sequences from assembly at this stage. \
Scaffold draft genome with cranberry reference genome ([RagTag](https://github.com/malonge/RagTag/wiki))

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

Alternative option is [BRAKER pipeline](https://github.com/Gaius-Augustus/BRAKER)
1. Map RNAseq to genome --> .bam
2. Sort the bam file by read names (samtools -sort) --> sorted.bam
3. Run BRAKER pipeline
  1. GeneMark-ES will use unsupervised gene prediction based on RNAseq mapped onto genome (based on eukaryotic coding genes features)
  2. Train AUGUSTUS based on GeneMark-ES gene prediction
  3. AUGUSTUS performs ab initio gene prediction
4. At the end of BRAKER, it produces genes annotated to the genome (?) --> .gff3?
5. 

Because I am having trouble installing BRAKER, alternative pipeline is [EVidenceModeler](https://github.com/EVidenceModeler/EVidenceModeler/wiki) which was actually what was used in Cui et al. 2022 V. darrowii paper and also Greg's colleague's sunflower genome annotation. 


Useful protein search resources (used in corporation with the BlastP search: [uniProt](https://www.uniprot.org/uniprotkb?facets=model_organism%3A3702&query=arabidopsis), it has Arabidopsis 136,000+ genes and 1,200+ Vaccinium genes). \
I downloaded them and created the blastdb. 


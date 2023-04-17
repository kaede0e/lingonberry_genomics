## Population genetics and demographic history of lingonberry subspecies 
Questions to address:
* How do the two subspecies genomes differ? Are they different enough to be considered subspecies? (ssp. minus vs vitis-idaea) 
* How much are their genomes differentiated? 
* When did the two subspecies split? What are their common ancestors? 
* What is their demograhpic history? How did their populations evolve/change over time?

### paired genome comparison ### 
Whole genome alignment can be achieved (minimap2) \
Variants may be detected (SyRI); however, large SVs will likely be missed because the assemblies are not de novo assembled \
What is the best two-genome comparison here? What questions am I addressing? 

### population models ### 
Useful link that lists all possible software packages here: https://www.duckdna.org/softwares/

Heterogeneity estimation using Illumina reads mapped to genome ([BWA_mem](https://github.com/lh3/bwa)) \
SNP calling within the single genome can be performed using [GATK](https://github.com/broadinstitute/gatk/blob/master/README.md), nucleotide diversity can be calculated using VCFtools \
The genome wide mask file should be created based on repetitive reads using [GenMap](https://github.com/cpockrandt/genmap) \
Use multiple phased genome sequence simultaneously, by separating into maternal and paternal haplotypes ([MSMC](http://www.github.com/stschiff/msmc-tools)) \
This computes the most probable ancestral recombination histories based on genealogical trees. \
It approximates coalescence rates using the Hidden Markov Model (HMM). 

Use pairwise Markovian coalescent ([PSMC](https://github.com/lh3/psmc)) \
This computes population models similar to MSMC, but computes based on pairwise coalescence. Instead of multiple individuals it can only perform one to one at a time.

For unphased genomes with lower coverage, [SMC++](https://github.com/popgenmethods/smcpp) is an option. Note that this was developed for population-level sampling and generally recommended with >20 samples. It is described to perform worse than PSMC when sample size is small. \
The model is based on sample frequency spectrum (SFS) instead of the SMC, so if the sample size is small there is just the minimal number of haplotypes to do this (if n=1 for a diploid individual, it's a sampling freq based on just two haplotypes). That's how accurate it can get. 

For computing runs of homozygosity to estimate the inbreeding coefficient, [ROHan](https://github.com/grenaud/ROHan) is a useful programme Jordan recommends.

Justin says it's useful to check the stats on multihetsep files (https://github.com/stschiff/msmc-tools/blob/master/getStats.d) but I need to compile D to run it. - if I can't figure out he said he can get stats from my multihetsep files! 


### signals of natural selection? ### 
This is kind of ambitious to do with only two genomes sequencing data... but here is idea: 
* population branch statistic (PBS) = summary statistic of pairwise FST values among three populations


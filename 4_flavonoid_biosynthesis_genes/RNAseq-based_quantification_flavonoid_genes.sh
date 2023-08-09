#!/bin/bash
#Genome annotation with genes
export PATH=$PATH:/~/reference_data/RNAseq/
module load hisat2 stringtie samtools bedtools

## Use annotated genes to count read depth and normalize = gene expression 

#6. Find flavonoid biosynthesis pathway genes from UniProt hits
#- I searched up the Colle, et al. 2019 paper's mentioend enzymes in UniProt
#- downloaded the FASTA for most from V.corymbosum, 2 from Arabidopsis and 1 from Vitis vinifera
#- use BLASTP to find which genes annotated in the V.corymbosum protein set is the one I'm looking for
#- find the orthologous pair from the Orthofinder output N0.tsv. 

#7. Determine expression levels of flavonoid biosynthesis pathway genes in lingonberry fruits and flowers
#StringTie can find expression estimate, DESeq2 seems classic for performing any sort of statistical analysis on the expression levels. 
hisat2-build $input_refgenome Lingonberry_RedCandy_hisat2_index #took about 20min

#1.1 Red berry
RNAseq_m1_list_of_files=/~/raw_fastq_files/redberry_SRR5799279_1.fastq.gz
RNAseq_m2_list_of_files=/~/raw_fastq_files/redberry_SRR5799279_2.fastq.gz
hisat2 \
-q -x Lingonberry_RedCandy_hisat2_index -1 `echo $RNAseq_m1_list_of_files` -2 `echo $RNAseq_m2_list_of_files` -S RedCandy_redberry_RNAseq_on_genome_hisat2_aln.sam
samtools sort -o RedCandy_redberry_RNAseq_on_genome_hisat2_sorted.aln.bam RedCandy_redberry_RNAseq_on_genome_hisat2_aln.sam

#1.2 White berry
RNAseq_m1_list_of_files=/~/raw_fastq_files/whiteberry_SRR5799277_1.fastq.gz
RNAseq_m2_list_of_files=/~/raw_fastq_files/whiteberry_SRR5799277_2.fastq.gz
hisat2 \
-q -x Lingonberry_RedCandy_hisat2_index -1 `echo $RNAseq_m1_list_of_files` -2 `echo $RNAseq_m2_list_of_files` -S RedCandy_whiteberry_RNAseq_on_genome_hisat2_aln.sam
samtools sort -o RedCandy_whiteberry_RNAseq_on_genome_hisat2_sorted.aln.bam RedCandy_whiteberry_RNAseq_on_genome_hisat2_aln.sam

#1.3 Green berry
RNAseq_m1_list_of_files=/~/raw_fastq_files/greenberry_SRR5799278_1.fastq.gz
RNAseq_m2_list_of_files=/~/raw_fastq_files/greenberry_SRR5799278_2.fastq.gz
hisat2 \
-q -x Lingonberry_RedCandy_hisat2_index -1 `echo $RNAseq_m1_list_of_files` -2 `echo $RNAseq_m2_list_of_files` -S RedCandy_greenberry_RNAseq_on_genome_hisat2_aln.sam
samtools sort -o RedCandy_greenberry_RNAseq_on_genome_hisat2_sorted.aln.bam RedCandy_greenberry_RNAseq_on_genome_hisat2_aln.sam

#1.4 my berry
RNAseq_m1_list_of_files=/~/raw_fastq_files/Lingonberry_RedCandy_berry_F123133_1_paired_trimmomatic.fastq.gz
RNAseq_m2_list_of_files=/~/raw_fastq_files/Lingonberry_RedCandy_berry_F123133_2_paired_trimmomatic.fastq.gz
hisat2 \
-q -x Lingonberry_RedCandy_hisat2_index -1 `echo $RNAseq_m1_list_of_files` -2 `echo $RNAseq_m2_list_of_files` -S RedCandy_myberry_RNAseq_on_genome_hisat2_aln.sam
samtools sort -o RedCandy_myberry_RNAseq_on_genome_hisat2_sorted.aln.bam RedCandy_myberry_RNAseq_on_genome_hisat2_aln.sam

#1.5 Flower
RNAseq_m1_list_of_files=/~/raw_fastq_files/Lingonberry_RedCandy_flower_F123132_1_paired_trimmomatic.fastq.gz
RNAseq_m2_list_of_files=/~/raw_fastq_files/Lingonberry_RedCandy_flower_F123132_2_paired_trimmomatic.fastq.gz
hisat2 \
-q -x Lingonberry_RedCandy_hisat2_index -1 `echo $RNAseq_m1_list_of_files` -2 `echo $RNAseq_m2_list_of_files` -S RedCandy_flower_RNAseq_on_genome_hisat2_aln.sam
samtools sort -o RedCandy_flower_RNAseq_on_genome_hisat2_sorted.aln.bam RedCandy_flower_RNAseq_on_genome_hisat2_aln.sam

#1.6 lastly, perform transcription level estimate
for files in `ls RedCandy_*_RNAseq_on_genome_hisat2_sorted.aln.bam`;
do
       stringtie -G /~/output/gene_annotation_pipeline/~/Lingonberry_RedCandy_genes_anno.gff3 \
       -A ${files}_gene_abund.tab -o ${files}_output.gtf \
       -e $files;
done

#Use Gene abundance -A as a proxy for expression level using StringTie. 
## I need to get the corresponding genes in blueberry. 
for enzyme in `cat enzyme_names.txt`;
do
       cat flavonoid_enzymes_Colle2019.txt | grep "$enzyme" > stringnames_${enzyme}.txt
done
mkdir orthologous_strings
for enzyme in `cat enzyme_names.txt`;
do
       for gene in `cat stringnames_${enzyme}.txt`;
       do
              cat N0.tsv | grep "$gene" > orthologous_strings/${enzyme}_${gene}.gff3
       done
done
for enzyme in `cat ../enzyme_names.txt`;
do
       cat ${enzyme}_* > ${enzyme}_orthologs.txt
done
for enzyme in `cat ../enzyme_names.txt`;
do
       cat ${enzyme}_orthologs.txt | cut -f 7 | sed 's/\./    /2' | cut -f 1 | sort | uniq > lingonberry_${enzyme}_orthologs.txt;
done

#extract the orthologous genes in StringTie output
#1 PAL
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.19798\|STRG.2345\|STRG.36377" | awk '{print $0 "   " "PAL"}' > RedCandy_${tissue}_RNAseq_on_genome_PAL.tab
done
#2 HQT
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.14993" | awk '{print $0 "   " "HQT"}' > RedCandy_${tissue}_RNAseq_on_genome_HQT.tab
done
#3 HCT
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.14980\|STRG.14993\|STRG.15008\|STRG.29314\|STRG.3188" | awk '{print $0 "   " "HCT"}' > RedCandy_${tissue}_RNAseq_on_genome_HCT.tab
done
#4 CoAl_4CL
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.15185\|STRG.29359\|STRG.35013\|STRG.37315" | awk '{print $0 "   " "4CL"}' > RedCandy_${tissue}_RNAseq_on_genome_CoAl_4CL.tab
done
#5 CHS
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.11010\|STRG.25785\|STRG.4773" | awk '{print $0 "   " "CHS"}' > RedCandy_${tissue}_RNAseq_on_genome_CHS.tab
done
#6 CHI
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.36288" | awk '{print $0 "   " "CHI"}' > RedCandy_${tissue}_RNAseq_on_genome_CHI.tab
done
#7 C4H
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.19324\|STRG.37808" | awk '{print $0 "   " "C4H"}' > RedCandy_${tissue}_RNAseq_on_genome_C4H.tab
done
#8 C3H
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.957" | awk '{print $0 "   " "C3H"}' > RedCandy_${tissue}_RNAseq_on_genome_C3H.tab
done
#9 FLS
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.13136\|STRG.13258\|STRG.5589" | awk '{print $0 "   " "FLS"}' > RedCandy_${tissue}_RNAseq_on_genome_FLS.tab
done
#10 FHT
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.15352\|STRG.25711" | awk '{print $0 "   " "FHT"}' > RedCandy_${tissue}_RNAseq_on_genome_FHT.tab
done
#11 UFGT
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.15158\|STRG.15162\|STRG.34161\|STRG.34247" | awk '{print $0 "   " "UFGT"}' > RedCandy_${tissue}_RNAseq_on_genome_UFGT.tab
done
#12 TT19
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.25529\|STRG.5915" | awk '{print $0 "   " "TT19"}' > RedCandy_${tissue}_RNAseq_on_genome_TT19.tab
done
#13 TT12
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.12763\|STRG.24688\|STRG.34849\|STRG.5930\|STRG.7941\|STRG.9615\|STRG.9628\|STRG.9629" | awk '{print $0 "   " "TT12"}' > RedCandy_${tissue}_RNAseq_on_genome_TT12.tab
done
#14 OMT
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.30037\|STRG.42261" | awk '{print $0 "   " "OMT"}' > RedCandy_${tissue}_RNAseq_on_genome_OMT.tab
done
#15 LAR
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.33735\|STRG.6628\|STRG.9223" | awk '{print $0 "   " "LAR"}' > RedCandy_${tissue}_RNAseq_on_genome_LAR.tab
done
#16 F3pH
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.27277" | awk '{print $0 "   " "F3pH"}' > RedCandy_${tissue}_RNAseq_on_genome_F3pH.tab
done
#17 F3p5pH
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.39736\|STRG.39737" | awk '{print $0 "   " "F3pH5pH"}' > RedCandy_${tissue}_RNAseq_on_genome_F3p5pH.tab
done
#18 DFR
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.24785\|STRG.26863" | awk '{print $0 "   " "DFR"}' > RedCandy_${tissue}_RNAseq_on_genome_DFR.tab
done
#19 ANS
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.39690" | awk '{print $0 "   " "ANS"}' > RedCandy_${tissue}_RNAseq_on_genome_ANS.tab
done
#20 ANR
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_hisat2_sorted.aln.bam_gene_abund.tab | grep -i -w "STRG.20771" | awk '{print $0 "   " "ANR"}' > RedCandy_${tissue}_RNAseq_on_genome_ANR.tab
done

#combine finally!
for tissue in `cat tissue_types.txt`;
do
       cat RedCandy_${tissue}_RNAseq_on_genome_*.tab > RedCandy_${tissue}_RNAseq_on_genome_gene_abund.txt
done







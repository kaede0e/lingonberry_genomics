## 4. Anthocyanin production and phenylpropanoid pathway in lingonberry
As the last part of my thesis, it would be nice to highlight the presence of important enzymatic pathways in anthocyanin biosynthesis (phenylpropanoid pathway, flavonoid pathway) in berries. \
My goal is then to show the presence of those genes in my assembly, hopefully quantify them in different tissue types (various berry and flower development stages).

How do you do it? 
I have: 
- Annotation file in .gff3 
- Annotated genes in .fasta, .faa where the feature name corresponds to each gene
- EggNog mapper results for orthologous genes in model organisms, labelled with KEGG, GO terms, etc. 
- RNAseq data from different tissue types
- Orthofinder is great at finding orthologs
My method overview: 
1. Run Orthofinder with tetraploid V. corymbosum Draper included. Orthofinder is run best with more genome inputs, so I chose to include V. macrocarpon Stevens, V. vitis-idaea, and Rhododendron williamsianum as the closely related outgroup. (found 174,171 genes in orthogroups, 1,403 single-copy orthogroups)
2. On the N0.tsv, identify the orthologs in lingonberry that correspond to the blueberry enzyme. 
3. Run HiSAT + StringTie (-A) to get gene abundance estimate - use this as a proxy for gene expression level. 
4. Plot the FPKM by tissue type & enzyme, and map on flavonoid biosynthesis pathway. 


Example papers: 
- Blueberry paper (Colle, et al. 2019) is a good example showing the panel of flavonoid biosynthetic pathways in berries at different developmental stages (Figure 3). 
- Subtropical blueberry paper (Cui et al. 2022) performed GO enrichment analysis and tissue-specific expression profiles (Figure 6). Cuticle formation pathway is highlighted for similar type of analysis (Figure 7). 
- Another blueberry paper (Yu et al. 2019) is a good example showing flavonoid biosynthesis pathways in berries, directly highlighting the differentially expressed genes (DEGs) identified by different berry developmental stages on KEGG map (Figure 5).

Methods (use Colle, et al. 2019 resource): 
1. Look to see in Orthofinder if I have flavonoid biosynthesis genes categorized as orthologous in lingonberry and start from there. 
2. If not, then get the NCBI sequences from this V. corymbosum genome and download .fasta sequences. 
3. BLAST enzymatic genes to lingonberry genome. 
4. Annotate those genes with obvious "genes" annotation. 
5. Align RNAseq to the genome, try to quantify expression levels. 
6. Cui et al. used [StringTie](https://ccb.jhu.edu/software/stringtie/) to calculate the gene expression levels, and [DESeq2](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) to perform differentially expressed genes analysis - more of a statistical implications of whether to call something significalty more/less expressed. \
DESeq2 seems to appear in other papers too so maybe I should give it a try. 



## Vaccinium phylogenetics 
This aims to compare lingonberry (V. vitis-idaea) with other Vaccinium species and ultimately re-define their phylogenetic relationships using whole-genome assembly. 
Questions to address: 
* How are they related? What is the correct phylogeny? 
* When did lingonberry diverge from the closest relative? 
* Are there any genes undergoing selection or adaptation particularly for lingonberry? some known traits like cold adapted, nutritional requirement, berry colors, etc.

Genome assembly available for: 

| Species | Assembly level  | ploidy | Annotated |
| :------------------- |:-----:| :------------------:| :-------------: |
| V. macrocarpon Stevens    | chr | 2x | yes |
| V. macrocarpon Ben Lear | chr | 2x | yes |
| V. oxycoccos | scaf | 2x | yes |
| V. microcarpum | scaf | 2x | yes |
| V. corymbosum Draper | chr | 4x | yes |
| V. darrowii | chr | 2x | yes |
| V. caesariense | scaf | 2x | yes |
| V. myrtillus | chr | 2x | yes |
| V. bacteatum | chr | 2x | yes |
| The pangenome data availability is questionable; it is still restricted access due to pre-publication. | scaf | var. | yes |
V. bacteatum genome is published, but the annotation and protein seq data is not available for download. 

### Phylogenetic tree construction 
How is a tree constructed? \
There is genetic tree and species tree, I'm interested in doing both? \
Phylogeny of species in the same genus (relatively closely related species) seems to be built from: 
* Protein sequences (Diaz-Garcia, et al. 2021 used 25 phylogenetically close species to Vmac and identified orthologs with [Orthofinder](https://github.com/davidemms/OrthoFinder); Cui et al. 2022 used 10 related species proteins and did pretty much the same thing except that they used only single-copy genes in generating tree using [MAFFT](https://mafft.cbrc.jp/alignment/software/manual/manual.html) for MSA and [IQTREE](http://www.iqtree.org/doc/Concordance-Factor) for tree construction), which is also what Xiao et al. 2023 Greg sent me did with salsugineum genome. 
* Plastid genome sequence, like those encoding chloroplast and mitochondria (Polashock et al. 2014) 
* A single selected gene of interest (Wu et al. 2021 was interested in MYB genes evolution so they used BLAST to identify orthologs, aligned with MUSCLE, and constructed tree with PhyML)
* To date the evolutionary tree, you need fossil calibration from one or two Vaccinium ancestors. (Schwery, et al. 2015 used Vaccinium creedensis fossil data, dating back to 26.5 myr min constraint) 
* The closest relative to Vaccinium would be (Kron, et al. 2002): Chamaedaphne and Gaultheria, both in Vaccinioideae (Schwery, et al. 2015 used Chamaedaphne calyculata fossil data, dating 2.58 myr min constraint) 

Multiple sequence comparisons ([MUSCLE](https://drive5.com/muscle5/)) seems good, superior to Clustral and MAFFT for multiple sequence alignment \
Alternativey, phylogenetic estimation usingn maximum likelihood ([PHYML](https://github.com/stephaneguindon/phyml)) seems to produce a robust tree, but might need a fossil calibration which I don't have. \
Once alignments data is compiled, [IQTREE](http://www.iqtree.org/doc/Concordance-Factor) should be good to make a phylogenetic tree. \
Based on specific genes of interest (single copy genes for example), [ASTRAL](https://github.com/smirarab/ASTRAL) makes an accurate species tree. 

### What can I do with the Orthofinder results? 
Species_Tree/SpeciesTree_rooted.txt provides the overall species tree inferred from the whole sets of gene tree. \
The STAG support value (at the node, looks like boostrap value but is not) represents the proportion of gene trees that supports this node topology. It is cocnsidered more robust than the conventional bootstrap support values. \
Gene_Trees/OG000XXXX_tree.txt provides individual gene tree for each orthogroup. \
If I can pick a gene I'm interested in, I should pull out the specific gene tree and see if the pattern matches the species tree, or supported by the earlier works. \
If I search up in "Arabidopsis thaliana" "one2one" in the Eggnog_mapper output, I can get a rough sense of orthologous matches in model organisms, but not all genes are there in the Orthofinder output. I need to know which genes to search up.

### gene evolution in Vaccinium species
It would be interesting to investigate further into specific gene evolution within this genus; considering the diverse morphology and biochemical applications. \
Transcripts and chemical properties are generally available for many Vaccinium berries. 

How can I find a support and good explanations for the seen species tree pattern? \
Why is there a discrepancy from previous knowledge (i.e., not clustering with cranberry, but its closest relative is bilberry)? 
* Maybe the preliminary dataset did not include bilberry in it so assumed it was a branch off from cranberry? 
* Maybe the preliminary studies did not thoroughly investigate all the genes? 
* Maybe gene duplication events for those targeted genes before had topologies different from the overall species tree? 


### Flavonoid biosynthesis pathway 
Goal here is to look for evolutionary questions regarding flavonoid biosynthesis in Vaccinium. How does lingonberry develop the red pigment colouration that is quite divergent from blueberry and bilberry? 



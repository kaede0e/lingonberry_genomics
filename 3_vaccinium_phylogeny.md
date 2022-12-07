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
| V. corymbosum Draper | chr | 4x | yes |
| V. darrowii | chr | 2x | yes |
| V. caesariense | scaf | 2x | yes |
| V. myrtillus | chr | 2x | yes |
| The pangenome data availability is questionable; it is still restricted access due to pre-publication. | scaf | var. | yes

### Phylogenetic tree construction 
How is a tree constructed? \
There is genetic tree and species tree, I'm interested in doing both? \
Multiple sequence comparisons ([MUSCLE](https://drive5.com/muscle5/)) seems good, superior to Clustral and MAFFT for multiple sequence alignment \
Alternativey, phylogenetic estimation usingn maximum likelihood ([PHYML](https://github.com/stephaneguindon/phyml)) seems to produce a robust tree, but might need a fossil calibration which I don't have. \
Once alignments data is compiled, [IQTREE](http://www.iqtree.org/doc/Concordance-Factor) should be good to make a phylogenetic tree. 

### gene evolution in Vaccinium species
It would be interesting to investigate further into specific gene evolution within this genus; considering the diverse morphology and biochemical applications. \
Transcripts and chemical properties are generally available for many Vaccinium berries. \


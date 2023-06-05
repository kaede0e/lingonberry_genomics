#!/bin/bash
## maffilter - statistical tools for analyzing MAF (multiple alignment file) ##

# ---------------------------------------------------------------------
export PATH=$PATH:/project/ctb-grego/khirabayashi/bin/minimap2
export PATH=$PATH:/project/ctb-grego/khirabayashi/bin/minimap2/misc
export PATH=$PATH:/project/ctb-grego/khirabayashi/bin/MafFilter #stand-alone version v1.3.1
# ---------------------------------------------------------------------

#For computing pairwise divergence, align two species genomes with identifier.chr# and convert to maf. 
minimap2 -ax asm5 --cs=long -t 40 Lingonberry_minus_asm_7.ragtag.scaffold.chr.fasta Lingonberry_RedCandy_asm_7.2.ragtag.scaffold.chr.fasta > minimap2_Lingonberry_minus_RedCandy.chr.aln.sam
paftools.js sam2paf minimap2_Lingonberry_minus_RedCandy.chr.aln.sam >  minimap2_Lingonberry_minus_RedCandy.chr.aln.paf
paftools.js view -f maf minimap2_Lingonberry_minus_RedCandy.chr.aln.paf > minimap2_Lingonberry_minus_RedCandy.chr.aln.maf
#you can replace chrnames after the fact to species.chr# format (not sure if this is necessary)
cat minimap2_Lingonberry_minus_RedCandy.chr.aln.maf | sed s/'Chr01_Vaccinium_vitis-idaea_ssp_minus'/'Vvitisminus.1'/g | sed s/'Chr02_Vaccinium_vitis-idaea_ssp_minus'/'Vvitisminus.2'/g | sed s/'Chr03_Vaccinium_vitis-idaea_ssp_minus'/'Vvitisminus.3'/g | sed s/'Chr04_Vaccinium_vitis-idaea_ssp_minus'/'Vvitisminus.4'/g | sed s/'Chr05_Vaccinium_vitis-idaea_ssp_minus'/'Vvitisminus.5'/g | sed s/'Chr06_Vaccinium_vitis-idaea_ssp_minus'/'Vvitisminus.6'/g | sed s/'Chr07_Vaccinium_vitis-idaea_ssp_minus'/'Vvitisminus.7'/g | sed s/'Chr08_Vaccinium_vitis-idaea_ssp_minus'/'Vvitisminus.8'/g | sed s/'Chr09_Vaccinium_vitis-idaea_ssp_minus'/'Vvitisminus.9'/g | sed s/'Chr10_Vaccinium_vitis-idaea_ssp_minus'/'Vvitisminus.10'/g | sed s/'Chr11_Vaccinium_vitis-idaea_ssp_minus'/'Vvitisminus.11'/g | sed s/'Chr12_Vaccinium_vitis-idaea_ssp_minus'/'Vvitisminus.12'/g \
| sed s/'Chr01_Vaccinium_vitis-idaea_var_RedCandy'/'Vvitisvitis.1'/g | sed s/'Chr02_Vaccinium_vitis-idaea_var_RedCandy'/'Vvitisvitis.2'/g | sed s/'Chr03_Vaccinium_vitis-idaea_var_RedCandy'/'Vvitisvitis.3'/g | sed s/'Chr04_Vaccinium_vitis-idaea_var_RedCandy'/'Vvitisvitis.4'/g | sed s/'Chr05_Vaccinium_vitis-idaea_var_RedCandy'/'Vvitisvitis.5'/g | sed s/'Chr06_Vaccinium_vitis-idaea_var_RedCandy'/'Vvitisvitis.6'/g | sed s/'Chr07_Vaccinium_vitis-idaea_var_RedCandy'/'Vvitisvitis.7'/g | sed s/'Chr08_Vaccinium_vitis-idaea_var_RedCandy'/'Vvitisvitis.8'/g | sed s/'Chr09_Vaccinium_vitis-idaea_var_RedCandy'/'Vvitisvitis.9'/g | sed s/'Chr10_Vaccinium_vitis-idaea_var_RedCandy'/'Vvitisvitis.10'/g | sed s/'Chr11_Vaccinium_vitis-idaea_var_RedCandy'/'Vvitisvitis.11'/g | sed s/'Chr12_Vaccinium_vitis-idaea_var_RedCandy'/'Vvitisvitis.12'/g \
> minimap2_Lingonberry_minus_RedCandy.maffilter.chr.aln.clean.maf 
#and make sure to clean up your MAF or else it will cause problems. 
#examples are... 
cat minimap2_Lingonberry_minus_RedCandy.maffilter.chr.aln.clean.maf | sed s/'a 0'/'a'/g > minimap2_Lingonberry_minus_RedCandy.maffilter.chr.aln.clean2.maf

#maffilter needs an option file(option_file.bpp) that defines what you want to do. 
maffilter param=option_file.bpp
#below is how I should make the option file for doing this analysis in 10kb windows, 1kbp minimum alignment block size. 

#1. Remove alignments where there are duplicated regions aligned to the same part
DATA=minimap2_Lingonberry_minus_RedCandy.maffilter.chr.aln.clean2
input.file=$(DATA).maf
input.file.compression=none
input.format=Maf
output.log=$(DATA).maffilter.dup.log
maf.filter=\
    Subset(species=(Vvitisminus,Vvitisvitis),\
    strict=yes,\
    keep=no,\
    remove_duplicates=yes),\
    Output(file=$(DATA).filt.uniq.maf)

#2. Break the MAF into sliding windows & calculate pairwise divergence in one step
DATA=minimap2_Lingonberry_minus_RedCandy.maffilter.chr.aln.clean2.filt.uniq
input.file=$(DATA).maf
input.file.compression=none
input.format=Maf
output.log=$(DATA).maffilter.dup.log
maf.filter=\
    WindowSplit(preferred_size=10000, align=ragged_left),\
    MinBlockLength(min_length=1000),\
    Output(file=$(DATA).windowed.maf),\
    SequenceStatistics(\
        statistics=(\
            PairwiseDivergence(\
                species1=Vvitisminus,\
                species2=Vvitisvitis)),\
        ref_species=Vvitisminus,\
        file=$(DATA).windowed.divergence.stats.txt)

###The above doesn't filter min_length 1000bp as expected, only looks for those in ref_species. 
cat minimap2_Lingonberry_minus_RedCandy.maffilter.chr.aln.clean2.filt.uniq.windowed.maf | tr -s " " | perl maf_grab_lengths.pl | sed s/' '/'   '/g | awk '{if ($4>=1000) {print}}' > minimap2_Lingonberry_minus_RedCandy.maffilter.chr.aln.clean2.filt.uniq.windowed.above1konvitis.maf
cat minimap2_Lingonberry_minus_RedCandy.maffilter.chr.aln.clean2.filt.uniq.windowed.above1konvitis.maf | awk '{ if ($8>=1000) {print}}' > minimap2_Lingonberry_minus_RedCandy_ragtag.scaff.maffilter.chr.aln.clean2.filt.uniq.above1k.maf

#this should remove all maf alignments with <1000bp aligned. 




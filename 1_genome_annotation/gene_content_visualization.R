library(tidyverse)
library(dplyr)
library(ggplot2)
#install.packages("devtools")
#library(devtools)
#devtools::install_github('TickingClock1992/RIdeogram')
#install.packages("RIdeogram") #not sure if I need to reload again?

getwd()
setwd("/Users/hirabayashikaede/Documents/UVic/Masters/Lingonberry data/Annotation")

#RIdeogram: drawing SVG graphics to visualize and map genome-wide data on idiograms

require(RIdeogram)
data(human_karyotype, package="RIdeogram")
data(gene_density, package="RIdeogram")
data(Random_RNAs_500, package="RIdeogram")


lingonberry_karyotype <- read.table("lingonberry_karyotype.txt", header = TRUE)
lingonberry_gene_density <- GFFex(input = "Lingonberry_RedCandy_genes_only_anno.gff3", karyotype = "lingonberry_karyotype.txt", feature = "gene", window = 1000000)
lingonberry_LTR_retrotransposon_density <- GFFex(input = "Lingonberry_RedCandy_asm_7.2.ragtag.scaffold.chr.fasta.chrnames.replaced.mod.EDTA.TEanno.gff3", karyotype = "lingonberry_karyotype.txt", feature = "LTR_retrotransposon", window = 1000000)%>% 
  mutate(Color="fc8d62")
lingonberry_hAT_TIR_transposon_density <- GFFex(input = "Lingonberry_RedCandy_asm_7.2.ragtag.scaffold.chr.fasta.chrnames.replaced.mod.EDTA.TEanno.gff3", karyotype = "lingonberry_karyotype.txt", feature = "hAT_TIR_transposon", window = 1000000)%>% 
  mutate(Color="fc8d10")
lingonberry_TE_density <- GFFex(input = "Lingonberry_RedCandy_EDTA_all_TE_anno.gff3", karyotype = "lingonberry_karyotype.txt", feature = "TE", window = 1000000)%>% 
  mutate(Color="#99d8c9")

ideogram(karyotype = lingonberry_karyotype, overlaid = lingonberry_gene_density, label = lingonberry_TE_density, label_type = "heatmap", colorset1 = c("#ffe8e8", "#ff99a8", "#d24140"), colorset2 = c("#9CB8D3", "#7486B5", "#432371"))
convertSVG("chromosome.svg", device = "pdf")
ideogram(karyotype = lingonberry_karyotype, overlaid = lingonberry_gene_density, label = lingonberry_LTR_retrotransposon_density, label_type = "polygon", colorset1 = c("#ffe8e8", "#ff99a8", "#d24140"))
convertSVG("chromosome.svg", device = "pdf")
ideogram(karyotype = lingonberry_karyotype, overlaid = lingonberry_gene_density, label = lingonberry_hAT_TIR_transposon_density, label_type = "polygon", colorset1 = c("#ffe8e8", "#ff99a8", "#d24140"))
convertSVG("chromosome.svg", device = "pdf")

#---- example -----#
data(liriodendron_karyotype, package="RIdeogram") #load the karyotype data
data(Fst_between_CE_and_CW, package="RIdeogram") #load the Fst data for overlaid heatmap
data(Pi_for_CE, package="RIdeogram") #load the Pi data for one-polygon label
ideogram(karyotype = liriodendron_karyotype, overlaid = Fst_between_CE_and_CW, label = Pi_for_CE, label_type = "polygon", colorset1 = c("#e5f5f9", "#99d8c9", "#2ca25f"))
convertSVG("chromosome.svg", device = "png")
ideogram(karyotype = human_karyotype, overlaid = gene_density, label = LTR_density, label_type = "heatmap", colorset1 = c("#f7f7f7", "#e34a33"), colorset2 = c("#f7f7f7", "#2c7fb8")) #use the arguments 'colorset1' and 'colorset2' to set the colors for gene and LTR heatmaps, separately.





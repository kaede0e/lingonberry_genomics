library(tidyverse)
library(dplyr)
library(ggplot2)
#install.packages("devtools")
#library(devtools)
#devtools::install_github('TickingClock1992/RIdeogram')
#install.packages("RIdeogram") #not sure if I need to reload again?

#RIdeogram: drawing SVG graphics to visualize and map genome-wide data on idiograms
require(RIdeogram)

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


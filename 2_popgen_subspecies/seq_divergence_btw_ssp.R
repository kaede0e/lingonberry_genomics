library(ggplot2)
library(compiler)
library(tidyverse)
library(dplyr)

#Computing sequence divergence between lingonberry subspecies with pairwise nucleotide alignment 

#------------------ plotting ------------------------#

#maybe something like coordinate plot
pdf("lingonberry_ssp_aln_visualization.pdf", height=24,width=24)
ssp_paf %>% 
  filter(aln_block_length <3000) %>%
  ggplot(., aes(x=ref_start, y=qry_start)) +
  facet_grid(ref ~ qry) +
  geom_point(aes(colour = strand))
dev.off()

#------------------ plotting by windows ------------------------#
#maffilter did calculation by 10kbpp windows for this. 
valid_maffilter_output %>% #10x8
  ggplot()+
  theme_minimal()+
#  coord_cartesian(ylim=c(0, 3))+ #this will keep all datapoints
  geom_line(aes(x=Start, y=Div.Vvitisminus.Vvitisvitis, group = Chr), color=c("#777e33")) +
  facet_grid(rows = vars(Chr))+
  xlab("ssp. minus chromosome position (bp)")+ ylab("Basepair mismatch in 10kb window (%)")

valid_maffilter_output %>% #10x7
  ggplot()+
  geom_histogram(aes(x=Div.Vvitisminus.Vvitisvitis, fill = above1k))+
  theme_minimal()+
  scale_fill_manual(values = c("#777e33"))+
  theme(legend.position = "none")+
  xlab("Pairwise divergence between ssp. minus and ssp. vitis-idaea (%)")+ ylab("Number of 10 kbp windows")


#--------combined plots genome-wide for presenting -------#
library("stringr")
rep_str = c('1'='Chr01','2'='Chr02','3'='Chr03','4'='Chr04','5'='Chr05','6'='Chr06','7'='Chr07','8'='Chr08','9'='Chr09','10'='Chr10','11'='Chr11','12'='Chr12') #do this sequentially
valid_maffilter_output$Chr <- str_replace_all(valid_maffilter_output$Chr, rep_str)
rep_str = c('Chr01Chr01'='Chr11','Chr01Chr02'='Chr12','Chr0Chr10'='Chr10')
valid_maffilter_output$Chr <- str_replace_all(valid_maffilter_output$Chr, rep_str)
valid_maffilter_output_Chr <- valid_maffilter_output %>% 
  mutate(window = floor(Start/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>% 
  select(Chr, window, window_middle, Div.Vvitisminus.Vvitisvitis, above1k)%>% 
  rename(pairwise_div = Div.Vvitisminus.Vvitisvitis, 
         CHROM = Chr)

lingonberry_genome_windowed <- left_join(valid_maffilter_output_Chr, LC1_allele_freq_with_zeros_mask_for_plot)
LW1_allele_freq_with_zeros_mask_for_plot_2 <- LW1_allele_freq_with_zeros_mask_for_plot%>% 
  rename(minus_percent_het = percent_het)%>% 
  select(CHROM, window, window_middle, minus_percent_het)
lingonberry_genome_windowed_both_ssp <- left_join(lingonberry_genome_windowed, LW1_allele_freq_with_zeros_mask_for_plot_2)

colors <- c("Pairwise divergence (%)" = "#777e33", "ssp.vitis-idaea %het" = "#d24140", "ssp.minus %het" = "#5E8DB8")
lingonberry_genome_windowed_both_ssp %>% #15x12 on portrait
  ggplot()+
  theme_minimal()+
  geom_line(aes(x=window,y=percent_het, group = CHROM, color="ssp.vitis-idaea %het")) +
  geom_line(aes(x=window,y=minus_percent_het, group = CHROM, color="ssp.minus %het")) +
  geom_line(aes(x=window,y=pairwise_div, group = CHROM, color="Pairwise divergence (%)"),size = 0.2) +
  coord_cartesian(ylim=c(0, 3))+
  facet_grid(rows = vars(CHROM))+
  xlab("Chromosome position (bp)")+ ylab("")+
  labs(color = "")+
  scale_color_manual(values = colors)





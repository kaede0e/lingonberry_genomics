library(tidyverse)
library(dplyr)
library(ggplot2)

getwd()
setwd("/Users/hirabayashikaede/Documents/UVic/Masters/Lingonberry data/Lingonberry_genomics/popgen_results")

#Troubleshooting my popgen simulation - is it inbreeding?

LC1_allele_freq <- read.table("Lingonberry_RedCandy.vcf.allele.frq.filtered.bed", header = FALSE) %>%
  select(V1, V2, V4)%>% 
  rename(CHROM = V1, 
         POS = V2,
         allele_freq = V4)
LW1_allele_freq <- read.table("Lingonberry_minus.vcf.allele.frq.filtered.bed", header = FALSE) %>%
  select(V1, V2, V4)%>% 
  rename(CHROM = V1, 
         POS = V2,
         allele_freq = V4)

#Is there chunks of heterozygosity/homozygosity? 

###### Makes sure all windows are in the dataframe, with the 0 hits ##### 
window_size = 100000 #100kb window

empty_windows <- read.table("Lingonberry_minus_asm_7.ragtag.scaffold.chr.bed") %>%
  as.tibble(.)%>% 
  select(V1, V3) %>% 
  rename(CHROM = V1, length = V3) %>% 
  mutate(max_window_size = floor(length/window_size))

## Do this per chromosome ##
# refer to the dataframe_prep_for_het_plots.R
data_frame = data.frame(
    window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:378) {   #the max. number = empty_windows$max_window_size+1                   
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Chr12_Vaccinium_myrtillus_NK2018_v1_RagTag")    
    # assigning this vector to ith row
  data_frame[i, ] <- vec
}
chr12_10kb_windows <- data_frame %>% 
  mutate(window = window_1 - window_size,
         window_middle = window + (0.5*window_size)) %>%
  select(CHROM, window, window_middle)

LW1_allele_counts_per_10kb_window <- LW1_allele_freq %>% 
  filter(CHROM == "Chr12_Vaccinium_myrtillus_NK2018_v1_RagTag")%>% 
  mutate(window = floor(POS/window_size)*window_size,
                           window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(n_alleles = n())
LC1_allele_counts_per_10kb_window <- LC1_allele_freq %>% 
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(n_alleles = n())

chr12_allelef_freq <- left_join(chr12_10kb_windows, LW1_allele_counts_per_10kb_window) %>% 
  replace(is.na(.), 0)

### Combine chromosomes ###
LW1_allele_freq_with_zeros_mask <- rbind(chr1_allelef_freq, chr2_allelef_freq, chr3_allelef_freq, chr4_allelef_freq, chr5_allelef_freq, chr6_allelef_freq,
      chr7_allelef_freq, chr8_allelef_freq, chr9_allelef_freq, chr10_allelef_freq, chr11_allelef_freq, chr12_allelef_freq)
LC1_allele_freq_with_zeros_mask <- LW1_allele_freq_with_zeros_mask %>%
  select(CHROM, window, window_middle) %>% 
  left_join(., LC1_allele_counts_per_10kb_window)%>% 
  replace(is.na(.), 0)

##### Visualize heterozygosity #####
LW1_het_counts_per_10kb_window_mask <- LW1_allele_freq %>% as_tibble(.)%>%
  filter(allele_freq == 0.5) %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(n_het_alleles = n())
LC1_het_counts_per_10kb_window_mask <- LC1_allele_freq %>% as_tibble(.)%>%
  filter(allele_freq == 0.5) %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(n_het_alleles = n())

#masked regions considered 
LC1_callable_sites <- read.table("Lingonberry_RedCandy.callable_sites.mapped.filtered.bed") %>% as_tibble(.)%>%
  rename(CHROM = V1, 
         chr_start = V2, 
         chr_end = V3)
LW1_callable_sites <- read.table("Lingonberry_minus.callable_sites.mapped.filtered.bed") %>% as_tibble(.) %>%
  rename(CHROM = V1, 
         chr_start = V2, 
         chr_end = V3)
#---------This ensures all sites are called, including those that cross windows ---------#
crossing_window <- LC1_callable_sites %>% 
  mutate(n_callable_sites = chr_end-chr_start, 
         window = floor(chr_start/window_size)*window_size,
         window_middle = window + (0.5*window_size), 
         site_window_end = floor(chr_end/window_size)*window_size)%>%
  mutate(crs = case_when(window != site_window_end ~ "cross_window", 
                         window == site_window_end ~ "within_window"))%>% 
  filter(crs == "cross_window")
start <- crossing_window %>% 
  select(CHROM, chr_start, site_window_end)
start_fixed <- as.integer(start$site_window_end)%>% as.tibble(.)
start_fixed_df <- cbind(start, start_fixed)%>% as.tibble(.)%>% 
  select(CHROM, chr_start, value) %>% 
  rename(chr_end = value)
end <- crossing_window %>% 
  select(CHROM, site_window_end, chr_end)
end_fixed <- as.integer(end$site_window_end)%>% as.tibble(.)
end_fixed_df <- cbind(end, end_fixed)%>% as.tibble(.)%>% 
  select(CHROM, value, chr_end) %>% 
  rename(chr_start = value)
LC1_crossing_window_df <- rbind(start_fixed_df, end_fixed_df)%>% 
  mutate(n_callable_sites = chr_end-chr_start, 
         window = floor(chr_start/window_size)*window_size,
         window_middle = window + (0.5*window_size), 
         site_window_end = floor(chr_end/window_size)*window_size)

crossing_window <- LW1_callable_sites %>% 
  mutate(n_callable_sites = chr_end-chr_start, 
         window = floor(chr_start/window_size)*window_size,
         window_middle = window + (0.5*window_size), 
         site_window_end = floor(chr_end/window_size)*window_size)%>%
  mutate(crs = case_when(window != site_window_end ~ "cross_window", 
                         window == site_window_end ~ "within_window"))%>% 
  filter(crs == "cross_window")
start <- crossing_window %>% 
  select(CHROM, chr_start, site_window_end)
start_fixed <- as.integer(start$site_window_end)%>% as.tibble(.)
start_fixed_df <- cbind(start, start_fixed)%>% as.tibble(.)%>% 
  select(CHROM, chr_start, value) %>% 
  rename(chr_end = value)
end <- crossing_window %>% 
  select(CHROM, site_window_end, chr_end)
end_fixed <- as.integer(end$site_window_end)%>% as.tibble(.)
end_fixed_df <- cbind(end, end_fixed)%>% as.tibble(.)%>% 
  select(CHROM, value, chr_end) %>% 
  rename(chr_start = value)
LW1_crossing_window_df <- rbind(start_fixed_df, end_fixed_df)%>% 
  mutate(n_callable_sites = chr_end-chr_start, 
         window = floor(chr_start/window_size)*window_size,
         window_middle = window + (0.5*window_size), 
         site_window_end = floor(chr_end/window_size)*window_size)
#------------------#

LC1_callable_sites_within_windows <- LC1_callable_sites %>% 
  mutate(n_callable_sites = chr_end-chr_start, 
         window = floor(chr_start/window_size)*window_size,
         window_middle = window + (0.5*window_size), 
         site_window_end = floor(chr_end/window_size)*window_size)%>% 
  mutate(crs = case_when(window != site_window_end ~ "cross_window", 
                         window == site_window_end ~ "within_window"))%>% 
  filter(crs == "within_window")%>% 
  select(-crs)%>% 
  rbind(., LC1_crossing_window_df)

LW1_callable_sites_within_windows <- LW1_callable_sites %>% 
  mutate(n_callable_sites = chr_end-chr_start, 
         window = floor(chr_start/window_size)*window_size,
         window_middle = window + (0.5*window_size), 
         site_window_end = floor(chr_end/window_size)*window_size)%>% 
  mutate(crs = case_when(window != site_window_end ~ "cross_window", 
                         window == site_window_end ~ "within_window"))%>% 
  filter(crs == "within_window")%>% 
  select(-crs)%>% 
  rbind(., LW1_crossing_window_df)

LC1_callable_sites_per_window <- LC1_callable_sites_within_windows %>%
  group_by(CHROM, window, window_middle) %>% 
  summarise(n_callable_sites_total = sum(n_callable_sites))

LW1_callable_sites_per_window <- LW1_callable_sites_within_windows%>%
  group_by(CHROM, window, window_middle) %>% 
  summarise(n_callable_sites_total = sum(n_callable_sites))


#--------------------------------------------------#
#plotting genome wide heterozygosity
LC1_callable_sites_per_window <- read_csv("LC1_callable_sites_per_window.csv")
LW1_callable_sites_per_window <- read_csv("LW1_callable_sites_per_window.csv")
LC1_het_counts_per_10kb_window_mask <- read.csv("LC1_het_counts_per_10kb_window_mask.csv")
LW1_het_counts_per_10kb_window_mask <- read.csv("LW1_het_counts_per_10kb_window_mask.csv")
LC1_allele_freq_with_zeros_mask <- read_csv("LC1_allele_freq_with_zeros_mask.csv")
LW1_allele_freq_with_zeros_mask <- read_csv("LW1_allele_freq_with_zeros_mask.csv")

LW1_allele_freq_with_zeros_mask_for_plot <- LW1_allele_freq_with_zeros_mask %>%#(export pdf as 10x10, portrait)
  left_join(., LW1_het_counts_per_10kb_window_mask)%>% 
  left_join(., LW1_callable_sites_per_window)%>%
  mutate(percent_het = 100 * (n_het_alleles/n_callable_sites_total),
         percent_calls = 100 * (n_alleles/n_callable_sites_total)) %>%
  filter(n_callable_sites_total >= 1000)%>%
  separate(CHROM, into = c("CHROM", "Vaccinium", "vitis-idaea", "ssp", "minus"), sep = "_") %>%
  select(CHROM, window, window_middle, percent_het, percent_calls)
  
LW1_allele_freq_with_zeros_mask_for_plot %>% 
  ggplot() +
  theme_minimal()+
  geom_line(aes(x=window,y=percent_het, group = CHROM), color=c("#5E8DB8")) +
  facet_grid(rows = vars(CHROM))+
  xlab("ssp. minus Chromosome position (bp)")+ ylab("Percent heterozygosity within 100 kb window (%)")
  
LC1_allele_freq_with_zeros_mask_for_plot <- LC1_allele_freq_with_zeros_mask %>%
  left_join(., LC1_het_counts_per_10kb_window_mask)%>% 
  left_join(., LC1_callable_sites_per_window)%>%
  mutate(percent_het = 100 * (n_het_alleles/n_callable_sites_total), 
         percent_calls = 100 * (n_alleles/n_callable_sites_total)) %>%
  filter(n_callable_sites_total >= 1000)%>%
  separate(CHROM, into = c("CHROM", "Vaccinium", "vitis-idaea", "ssp", "minus"), sep = "_") %>%
  select(CHROM, window, window_middle, percent_het, percent_calls)

LC1_allele_freq_with_zeros_mask_for_plot %>% 
  ggplot() +
  theme_minimal()+
  geom_line(aes(x=window,y=percent_het, group = CHROM), color=c("#d24140")) +
  facet_grid(rows = vars(CHROM))+
  xlab("ssp. vitis-idaea Chromosome position (bp)")+ ylab("Percent heterozygosity within 100 kb window (%)")



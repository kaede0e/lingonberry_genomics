#Gene expression analysis with heatmap
library(tidyverse)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)


lingonberry_flavonoid_exp <- read.table("flavonoid_biosynthesis_transcript_count_matrix.txt", header = FALSE)%>%
  as.tibble()%>% 
  rename(enzyme = V1, 
         transcript_id = V2, 
         RedCandy_flower = V3, 
         RedCandy_greenberry = V4, 
         RedCandy_leaf = V5, 
         RedCandy_myberry = V6, 
         RedCandy_redberry = V7,
         RedCandy_root = V8, 
         RedCandy_whiteberry = V9)

flower <- read.table("RedCandy_flower_RNAseq_on_genome_gene_abund.txt", header = FALSE) %>%
  as.tibble()%>% 
  rename(Gene_ID = V1,
         Gene = V2, 
         Reference = V3, 
         Strand = V4, 
         Start = V5, 
         End =V6, 
         Coverage =V7, 
         FPKM = V8, 
         TPM = V9, 
         enzyme = V10)%>% 
  mutate(sample = "RedCandy_flower")%>% 
  select(Gene_ID, enzyme, sample, FPKM)
greenberry <- read.table("RedCandy_greenberry_RNAseq_on_genome_gene_abund.txt", header = FALSE) %>%
  as.tibble()%>% 
  rename(Gene_ID = V1,
         Gene = V2, 
         Reference = V3, 
         Strand = V4, 
         Start = V5, 
         End =V6, 
         Coverage =V7, 
         FPKM = V8, 
         TPM = V9, 
         enzyme = V10)%>% 
  mutate(sample = "Sunna_greenberry")%>% 
  select(Gene_ID, enzyme, sample, FPKM)
whiteberry <- read.table("RedCandy_whiteberry_RNAseq_on_genome_gene_abund.txt", header = FALSE) %>%
  as.tibble()%>% 
  rename(Gene_ID = V1,
         Gene = V2, 
         Reference = V3, 
         Strand = V4, 
         Start = V5, 
         End =V6, 
         Coverage =V7, 
         FPKM = V8, 
         TPM = V9, 
         enzyme = V10)%>% 
  mutate(sample = "Sunna_whiteberry")%>% 
  select(Gene_ID, enzyme, sample, FPKM)
redberry <- read.table("RedCandy_redberry_RNAseq_on_genome_gene_abund.txt", header = FALSE) %>%
  as.tibble()%>% 
  rename(Gene_ID = V1,
         Gene = V2, 
         Reference = V3, 
         Strand = V4, 
         Start = V5, 
         End =V6, 
         Coverage =V7, 
         FPKM = V8, 
         TPM = V9, 
         enzyme = V10)%>% 
  mutate(sample = "Sunna_redberry")%>% 
  select(Gene_ID, enzyme, sample, FPKM)
myberry <- read.table("RedCandy_myberry_RNAseq_on_genome_gene_abund.txt", header = FALSE) %>%
  as.tibble()%>% 
  rename(Gene_ID = V1,
         Gene = V2, 
         Reference = V3, 
         Strand = V4, 
         Start = V5, 
         End =V6, 
         Coverage =V7, 
         FPKM = V8, 
         TPM = V9, 
         enzyme = V10)%>% 
  mutate(sample = "RedCandy_berry")%>% 
  select(Gene_ID, enzyme, sample, FPKM)
leaf <- read.table("RedCandy_leaf_RNAseq_on_genome_gene_abund.txt", header = FALSE) %>%
  as.tibble()%>% 
  rename(Gene_ID = V1,
         Gene = V2, 
         Reference = V3, 
         Strand = V4, 
         Start = V5, 
         End =V6, 
         Coverage =V7, 
         FPKM = V8, 
         TPM = V9, 
         enzyme = V10)%>% 
  mutate(sample = "RedCandy_leaf")%>% 
  select(Gene_ID, enzyme, sample, FPKM)
root <- read.table("RedCandy_root_RNAseq_on_genome_gene_abund.txt", header = FALSE) %>%
  as.tibble()%>% 
  rename(Gene_ID = V1,
         Gene = V2, 
         Reference = V3, 
         Strand = V4, 
         Start = V5, 
         End =V6, 
         Coverage =V7, 
         FPKM = V8, 
         TPM = V9, 
         enzyme = V10)%>% 
  mutate(sample = "RedCandy_root")%>% 
  select(Gene_ID, enzyme, sample, FPKM)

lingonberry_flavonoid_exp_lev <- rbind(leaf, root, flower, myberry, greenberry, whiteberry, redberry)
#----------------------------------

#Check genes in lingonberry are algined well to blueberry orthologs from BLAST alignment# 
blastp_align_df <- read.table("RedCandy_Vcorymbosum_flavonoid_genes_prot_blast.txt") %>% 
  as_tibble()
blueberry_gene_length <- read.table("blueberry_orthologs_length.txt", header = TRUE)%>% 
  as_tibble()

blastp_align <- blastp_align_df %>% 
  rename(percent_ID = V3, 
         align_length = V4) %>%
  select(V1, V2, percent_ID, align_length) %>%
  left_join(blueberry_gene_length) %>%
  filter(percent_ID >= 95) %>% #% identity of sequences 
  mutate(percent_align = align_length/gene_length) %>%
  filter(percent_align >= 0.8) #% alignment length 

### Check that alignments are >95% in positions identity between bluebery and lingonberry 
blastp_align_df_PAL <- blastp_align %>% filter(V1 == "STRG.19798.1.p1" | V1 == "STRG.2345.1.p1" | V1 == "STRG.36377.1.p1") %>% 
  mutate(enzyme = "PAL")
blastp_align_df_HQT <- blastp_align %>% filter(V1 == "STRG.14993.1.p1" | V1 == "STRG.14983.1.p1" | V1 == "STRG.14996.1.p1") %>% 
  mutate(enzyme = "HQT")
blastp_align_df_HCT <- blastp_align %>% filter(V1 == "STRG.14980.1.p1" | V1 == "STRG.14993.1.p1" | V1 == "STRG.15008.1.p1" | V1 == "STRG.29314.1.p1" | V1 == "STRG.14983.1.p1" | V1 == "STRG.3188.1.p1" | V1 == "STRG.14996.1.p1" | V1 == "STRG.15009.1.p1" | V1 == "STRG.15013.1.p1" | V1 == "STRG.15005.1.p1") %>% 
  mutate(enzyme = "HCT")
blastp_align_df_4CL <- blastp_align %>% filter(V1 == "STRG.15185.1.p1" | V1 == "STRG.29359.1.p1" | V1 == "STRG.35013.1.p1" | V1 == "STRG.37315.1.p1") %>% 
  mutate(enzyme = "4CL")
blastp_align_df_CHS <- blastp_align %>% filter(V1 == "STRG.11010.1.p1" | V1 == "STRG.25785.1.p1" | V1 == "STRG.4773.1.p1" | V1 == "STRG.6031.1.p1" | V1 == "STRG.6030.1.p1") %>% 
  mutate(enzyme = "CHS")
blastp_align_df_CHI <- blastp_align %>% filter(V1 == "STRG.36288.1.p1") %>% 
  mutate(enzyme = "CHI")
blastp_align_df_C4H <- blastp_align %>% filter(V1 == "STRG.19324.1.p1" | V1 == "STRG.37808.1.p1" | V1 == "STRG.40911.1.p1") %>% 
  mutate(enzyme = "C4H")
blastp_align_df_C3H <- blastp_align %>% filter(V1 == "STRG.957.6.p1" | V1 == "STRG.14981.2.p1") %>% 
  mutate(enzyme = "C3H")
blastp_align_df_FLS <- blastp_align %>% filter(V1 == "STRG.13136.1.p1" | V1 == "STRG.13258.1.p1" | V1 == "STRG.5589.1.p1" | V1 == "STRG.20701.1.p1") %>% 
  mutate(enzyme = "FLS")
blastp_align_df_FHT <- blastp_align %>% filter(V1 == "STRG.15352.1.p1" | V1 == "STRG.25711.1.p1") %>% 
  mutate(enzyme = "FHT")
blastp_align_df_UFGT <- blastp_align %>% filter(V1 == "STRG.31622.1.p1" | V1 == "STRG.15158.1.p1" | V1 == "STRG.15162.1.p1" | V1 == "STRG.34161.1.p1" | V1 == "STRG.34247.1.p1" | V1 == "STRG.41819.1.p1" | V1 == "STRG.34246.1.p1" | V1 == "STRG.34249.1.p1"| V1 == "STRG.34165.1.p1" | V1 == "STRG.34164.1.p1") %>% 
  mutate(enzyme = "UFGT")
blastp_align_df_TT19 <- blastp_align %>% filter(V1 == "STRG.15311.1.p1" | V1 == "STRG.25529.1.p1" | V1 == "STRG.5915.1.p1") %>% 
  mutate(enzyme = "TT19")
blastp_align_df_TT12 <- blastp_align %>% filter(V1 == "STRG.21368.1.p1" | V1 == "STRG.40478.1.p1" | V1 == "STRG.23147.1.p1" | V1 == "STRG.12763.1.p1" | V1 == "STRG.24688.1.p1" | V1 == "STRG.34849.1.p1" | V1 == "STRG.5930.1.p1" | V1 == "STRG.7941.1.p1"| V1 == "STRG.9615.1.p1" | V1 == "STRG.9628.1.p1" | V1 == "STRG.9629.1.p1"| V1 == "STRG.9625.1.p1"| V1 == "STRG.9355.1.p1"| V1 == "STRG.9630.1.p1"| V1 == "STRG.10493.1.p1"| V1 == "STRG.12768.1.p1"| V1 == "STRG.32226.1.p1") %>% 
  mutate(enzyme = "TT12")
blastp_align_df_OMT <- blastp_align %>% filter(V1 == "STRG.30037.1.p1" | V1 == "STRG.42261.1.p1" | V1 == "STRG.28075.1.p1") %>% 
  mutate(enzyme = "OMT")
blastp_align_df_LAR <- blastp_align %>% filter(V1 == "STRG.33735.1.p1" | V1 == "STRG.6628.1.p1" | V1 == "STRG.9223.1.p1") %>% 
  mutate(enzyme = "LAR")
blastp_align_df_F3pH <- blastp_align %>% filter(V1 == "STRG.27277.1.p1") %>% 
  mutate(enzyme = "F3pH")
blastp_align_df_F3p5pH <- blastp_align %>% filter(V1 == "STRG.27544.1.p1" | V1 == "STRG.39736.1.p1" | V1 == "STRG.39737.1.p1" | V1 == "STRG.39731.1.p1" | V1 == "STRG.39733.1.p1" | V1 == "STRG.39734.1.p1" | V1 == "STRG.39738.1.p1" ) %>% 
  mutate(enzyme = "F3p5pH")
blastp_align_df_DFR <- blastp_align %>% filter(V1 == "STRG.24785.1.p1" | V1 == "STRG.26863.5.p1" | V1 == "STRG.26865.1.p1" ) %>% 
  mutate(enzyme = "DFR")
blastp_align_df_ANS <- blastp_align %>% filter(V1 == "STRG.39690.1.p1") %>% 
  mutate(enzyme = "ANS")
blastp_align_df_ANR <- blastp_align %>% filter(V1 == "STRG.20771.1.p1") %>% 
  mutate(enzyme = "ANR")

blastp_align_comb <- read_csv("blastp_align_combined_filtered.csv")
blastp_align_comb <- rbind(blastp_align_df_PAL, blastp_align_df_HQT, blastp_align_df_HCT, blastp_align_df_4CL, blastp_align_df_CHS, 
                           blastp_align_df_CHI, blastp_align_df_C4H, blastp_align_df_C3H, blastp_align_df_FLS, blastp_align_df_FHT, 
                           blastp_align_df_UFGT, blastp_align_df_TT19, blastp_align_df_TT12, blastp_align_df_OMT, blastp_align_df_LAR, 
                           blastp_align_df_F3pH, blastp_align_df_F3p5pH, blastp_align_df_DFR, blastp_align_df_ANS, blastp_align_df_ANR)

filtered_gene_hits <- blastp_align_comb %>% 
  select(V1, enzyme) %>% 
  distinct() %>% 
  separate(V1, into = c("Gene_ID"), sep = ".1.p1", 
           convert = TRUE, remove = TRUE) #and I manually trimmed the non-1.p1 genes
## Those are the only genes we should consider including lingonberry gene expression analysis!!


#----------------------------------
#visualize in heatmap# 

lingonberry_flavonoid_exp_lev <- read_csv("lingonberry_flavonoid_exp_lev_rev.csv")
lingonberry_flavonoid_exp_lev_filtered <- filtered_gene_hits %>% 
  left_join(., lingonberry_flavonoid_exp_lev)

lingonberry_flavonoid_exp_lev_filtered$sample_f = factor(lingonberry_flavonoid_exp_lev_filtered$sample, levels=c('RedCandy_rhizome','RedCandy_leaf','RedCandy_flower','RedCandy_berry', 'Sunna_greenberry', 'Sunna_whiteberry', 'Sunna_redberry'))
lingonberry_flavonoid_exp_lev_filtered %>%
#  filter(enzyme == "PAL")%>%
#  filter(enzyme == "C4H")%>%
#  filter(enzyme == "4CL")%>%
#  filter(enzyme == "HCT")%>%
#  filter(enzyme == "CHS")%>%
#  filter(enzyme == "CHI")%>%
#  filter(enzyme == "FHT")%>%
#  filter(enzyme == "F3pH")%>%
#  filter(enzyme == "F3p5pH")%>%
#  filter(enzyme == "FLS")%>%
#  filter(enzyme == "ANS")%>%
#  filter(enzyme == "ANR")%>%
#  filter(enzyme == "LAR")%>%
#  filter(enzyme == "UFGT")%>%
#  filter(enzyme == "OMT")%>%
#  filter(enzyme == "C3H")%>%
#  filter(enzyme == "DFR")%>% 
#  filter(enzyme == "HQT")%>% 
#  filter(enzyme == "TT12")%>% #not on my diagram 
  filter(enzyme == "TT19")%>%
  ggplot()+
  geom_tile(aes(x=sample_f, y=Gene_ID, fill=FPKM))+
  scale_fill_gradient(low="#ffe8e8", high="#d24140")+
  theme_classic()+
  xlab("")+ ylab("")+
  facet_grid(cols = vars(enzyme))

lingonberry_flavonoid_exp_lev_filtered$enzyme_f = factor(lingonberry_flavonoid_exp_lev_filtered$enzyme, levels=c("PAL","C4H","4CL","HCT","CHS","CHI","FHT","F3pH","F3pH5pH","FLS","ANS","ANR","LAR","DFR","UFGT","OMT","C3H","HQT","TT12","TT19"))
lingonberry_flavonoid_exp_lev_filtered %>% #5x10 pixel on PDF landscape
  ggplot()+
  geom_tile(aes(x=sample_f, y=Gene_ID, fill=FPKM))+
  scale_fill_gradient(low="#FFFFFF", high="#d24140")+
  theme_classic()+
  xlab("Enzyme")+ ylab("")
lingonberry_flavonoid_exp_lev_filtered %>% filter(sample == "RedCandy_berry")%>% 
  group_by(enzyme, sample) %>% count() %>% 
  ungroup() %>% summarise(total_genes = sum(n)) #how many flavonoid genes were found in total?

lingonberry_flavonoid_exp_lev_filtered_zeros_replaced <- lingonberry_flavonoid_exp_lev_filtered %>%
  filter(FPKM == 0) %>% 
  mutate(FPKM_non_zero = 0.0001)%>% 
  select(Gene_ID, enzyme, sample, FPKM_non_zero, sample_f, enzyme_f)%>% 
  rename(FPKM = FPKM_non_zero)

lingonberry_flavonoid_exp_lev_filtered_no_zeros <- lingonberry_flavonoid_exp_lev_filtered %>%
  filter(FPKM >0) %>% 
  rbind(., lingonberry_flavonoid_exp_lev_filtered_zeros_replaced)
lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>%
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>%
  ggplot()+
  geom_tile(aes(x=sample_f, y=enzyme_gene, fill=z_FPKM))+
  scale_fill_gradientn(colours =c("#FFFFFF","#ffe8e8","#d24140"))+
  theme_classic()+
  xlab("")+ ylab("")
    
#make individually z-scaled dataframe and combine 
norm_4CL <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "4CL") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_ANS <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "ANS") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_ANR <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "ANR") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_C3H <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "C3H") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_C4H <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "C4H") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_CHI <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "CHI") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_CHS <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "CHS") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_DFR <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "DFR") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_F3pH <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "F3pH") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_F3pH5pH <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "F3p5pH") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_FHT <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "FHT") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_FLS <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "FLS") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_HCT <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "HCT") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_HQT <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "HQT") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_LAR <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "LAR") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_OMT <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "OMT") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_PAL <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "PAL") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_TT12 <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "TT12") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_TT19 <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "TT19") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_UFGT <- lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "UFGT") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_lingonberry_flavonoid_exp_lev_filtered_no_zeros<- rbind(norm_4CL, norm_ANR, norm_ANS, norm_C3H, norm_C4H, norm_CHI, norm_CHS, norm_DFR, norm_F3pH, norm_F3pH5pH, norm_FHT, norm_FLS, norm_HCT,  norm_HQT, norm_LAR, norm_OMT, norm_PAL, norm_TT12, norm_TT19, norm_UFGT)

norm_lingonberry_flavonoid_exp_lev_filtered_no_zeros$sample_f = factor(norm_lingonberry_flavonoid_exp_lev_filtered_no_zeros$sample, levels=c('RedCandy_rhizome','RedCandy_leaf','RedCandy_flower','RedCandy_berry', 'Sunna_greenberry', 'Sunna_whiteberry', 'Sunna_redberry'))
norm_lingonberry_flavonoid_exp_lev_filtered_no_zeros %>% #10x13 PDF portrait
  ggplot()+
  geom_tile(aes(x=sample_f, y=enzyme_gene, fill=z_FPKM))+
  scale_fill_gradientn(colours =c("#FFFFFF","#d24140"))+
  theme_classic()+
  xlab("")+ ylab("")



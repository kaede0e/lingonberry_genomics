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
#visualize in heatmap# 
lingonberry_flavonoid_exp_lev <- read_csv("lingonberry_flavonoid_exp_lev.csv")

lingonberry_flavonoid_exp_lev$sample_f = factor(lingonberry_flavonoid_exp_lev$sample, levels=c('RedCandy_root','RedCandy_leaf','RedCandy_flower','RedCandy_berry', 'Sunna_greenberry', 'Sunna_whiteberry', 'Sunna_redberry'))
lingonberry_flavonoid_exp_lev %>%
#  filter(enzyme == "PAL")%>%
#  filter(enzyme == "C4H")%>%
#  filter(enzyme == "4CL")%>%
#  filter(enzyme == "HCT")%>%
#  filter(enzyme == "CHS")%>%
#  filter(enzyme == "CHI")%>%
#  filter(enzyme == "FHT")%>%
#  filter(enzyme == "F3pH")%>%
#  filter(enzyme == "F3pH5pH")%>%
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

lingonberry_flavonoid_exp_lev$enzyme_f = factor(lingonberry_flavonoid_exp_lev$enzyme, levels=c("PAL","C4H","4CL","HCT","CHS","CHI","FHT","F3pH","F3pH5pH","FLS","ANS","ANR","LAR","DFR","UFGT","OMT","C3H","HQT","TT12","TT19"))
lingonberry_flavonoid_exp_lev %>% #5x10 pixel on PDF landscape
  ggplot()+
  geom_tile(aes(x=sample_f, y=Gene_ID, fill=FPKM))+
  scale_fill_gradient(low="#FFFFFF", high="#d24140")+
  theme_classic()+
  xlab("Enzyme")+ ylab("")
  facet_grid(rows = vars(enzyme))

#make individually z-scaled dataframe and combine (Supplementary Figure 
norm_4CL <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "4CL") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_ANS <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "ANS") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_ANR <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "ANR") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_C3H <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "C3H") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_C4H <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "C4H") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_CHI <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "CHI") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_CHS <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "CHS") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_DFR <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "DFR") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_F3pH <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "F3pH") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_F3pH5pH <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "F3pH5pH") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_FHT <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "FHT") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_FLS <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "FLS") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_HCT <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "HCT") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_HQT <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "HQT") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_LAR <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "LAR") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_OMT <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "OMT") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_PAL <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "PAL") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_TT12 <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "TT12") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_TT19 <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "TT19") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_UFGT <- lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  mutate(enzyme_gene = paste0(enzyme,"_",Gene_ID)) %>% filter(enzyme == "UFGT") %>% 
  mutate(log_FPKM = log10(FPKM))%>% 
  mutate(z_FPKM = scale(FPKM)) %>% select(sample, enzyme_gene, enzyme, FPKM, z_FPKM)
norm_lingonberry_flavonoid_exp_lev_no_zeros<- rbind(norm_4CL, norm_ANR, norm_ANS, norm_C3H, norm_C4H, norm_CHI, norm_CHS, norm_DFR, norm_F3pH, norm_F3pH5pH, norm_FHT, norm_FLS, norm_HCT,  norm_HQT, norm_LAR, norm_OMT, norm_PAL, norm_TT12, norm_TT19, norm_UFGT)

norm_lingonberry_flavonoid_exp_lev_no_zeros$sample_f = factor(norm_lingonberry_flavonoid_exp_lev_no_zeros$sample, levels=c('RedCandy_rhizome','RedCandy_leaf','RedCandy_flower','RedCandy_berry', 'Sunna_greenberry', 'Sunna_whiteberry', 'Sunna_redberry'))
norm_lingonberry_flavonoid_exp_lev_no_zeros %>% #10x13 PDF portrait
  ggplot()+
  geom_tile(aes(x=sample_f, y=enzyme_gene, fill=z_FPKM))+
  scale_fill_gradientn(colours =c("#FFFFFF","#d24140"))+
  theme_classic()+
  xlab("")+ ylab("")

    

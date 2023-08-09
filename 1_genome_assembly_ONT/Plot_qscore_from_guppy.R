library(tidyverse)
library(ggplot2)

qscore <- read_table("Lingonberry_minus_Takara_02_nsSRE_01_4_Dec16_2022_duplex_basecalled_stats.txt")
qscore %>%
  ggplot()+
  geom_density(aes(x=mean_qscore_template))

qscore %>%
  filter(sequence_length_template >= 10000) %>%
  ggplot()+
  geom_histogram(aes(x=sequence_length_template))

#Get N50 and Qscore mean for raw reads that guppy basecalled
qscore %>%
  summarise(mean(mean_qscore_template))
qscore <- qscore[order(qscore$sequence_length_template, decreasing= TRUE),]
qscore$cum_sum<- cumsum(qscore$sequence_length_template)
qscore[qscore$cum_sum > max(qscore$cum_sum)/2,][1,] #length corresponds to N50

#How much data produced 
qscore %>% filter(mean_qscore_template>=10)%>%
  summarise(data = sum(sequence_length_template))%>% 
  mutate(`data(Gb)`= data / 1000000000)
  
#Guppy basecaller pass/fail
qscore <- read_table("minus_Takara_02_nsSRE_01_2_Dec12_2022_sup_qc_data.txt")
qscore %>% 
  group_by(filename)%>% 
  count
  ggplot()+
  geom_histogram(aes(x=mean_qscore_template, 
                     fill = passes_filtering, colour = passes_filtering))
  
  

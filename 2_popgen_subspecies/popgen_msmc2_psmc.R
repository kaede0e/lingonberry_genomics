library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)
getwd()
setwd("/Users/hirabayashikaede/Documents/UVic/Masters/Lingonberry data/Lingonberry_genomics")

#--Github MSMC2 examples with human dataset------------------------------
mu <- 1.25e-8
gen <- 30
afrDat<-read.table("popgen_results/NA12878.chr1.msmc2.final.txt", header=TRUE)#single individual 
eurDat<-read.table("popgen_results/NA19240.chr1.msmc2.final.txt", header=TRUE)#single individual
afrDat<-read.table("popgen_results/AFR.msmc2.final.txt", header=TRUE)
eurDat<-read.table("popgen_results/EUR.msmc2.final.txt", header=TRUE)
plot(afrDat$left_time_boundary/mu*gen, (1/afrDat$lambda)/(mu), log="x",ylim=c(0,100000),
     type="n", xlab="Years ago", ylab="effective population size")
lines(afrDat$left_time_boundary/mu*gen, (1/afrDat$lambda)/(mu), type="s", col="red")#the tutorial is apparently wrong; with 2*mu 
lines(eurDat$left_time_boundary/mu*gen, (1/eurDat$lambda)/(mu), type="s", col="blue")
legend("topright",legend=c("African", "European"), col=c("red", "blue"), lty=c(1,1))
crossPopDat<-read.table("EUR_AFR.combined.msmc2.final.txt", header=TRUE)
plot(crossPopDat$left_time_boundary/mu*gen, 2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11),
     xlim=c(1000,500000),ylim=c(0,1), type="n", xlab="Years ago", ylab="relative cross-coalescence rate")
lines(crossPopDat$left_time_boundary/mu*gen, 2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11), type="s")

#--Lingonberry dataset------------------------------
mu <- 3e-9 #mutation rate in Arabidopsis thaliana
gen <- 10 #generation time of 10 years
minusDat<-read.table("popgen_results/Lingonberry_minus.msmc2.final.txt", header=TRUE) %>%
  as_tibble()%>%
  mutate(species = "minus")%>%
  mutate(Ne = (1/lambda)/(2*mu)) %>% 
  mutate(Years_ago = left_time_boundary/mu*gen)
RedCandyDat<-read.table("popgen_results/Lingonberry_RedCandy.msmc2.final.txt", header=TRUE)%>% 
  as_tibble()%>%
  mutate(species = "vitis-idaea")%>%
  mutate(Ne = (1/lambda)/(2*mu)) %>% 
  mutate(Years_ago = left_time_boundary/mu*gen)
LingonberryDat <- rbind(minusDat,RedCandyDat)

LingonberryDat%>% #save as 7x4.26 pdf
  mutate(Ne_10000 = Ne/10000)%>%
  ggplot()+
  geom_step(aes(Years_ago, Ne_10000, colour = species))+
  scale_x_log10()+
  theme_classic()+
  theme(legend.position = "na")+
  scale_colour_manual(values = c(LGpal[1], LGpal[11]))+
  xlab("Years ago (g=10, u=0.3e10^-8)")+ ylab("Effective population size (x10^4)")
  
  
plot(minusDat$left_time_boundary/mu*gen, (1/minusDat$lambda)/(2*mu)/10000, log="x",ylim=c(0,85),
     type="n", xlab="Years ago", ylab="Effective population size (x10^4)")
lines(minusDat$left_time_boundary/mu*gen, (1/minusDat$lambda)/(2*mu)/10000, type="s", col=LGpal[1])
lines(RedCandyDat$left_time_boundary/mu*gen, (1/RedCandyDat$lambda)/(2*mu)/10000, type="s", col=LGpal[11])
legend("topright",legend=c("V. vitis-idaea ssp. minus", "V. vitis-idaea ssp. vitis-idaea"), col=c(LGpal[1], LGpal[11]), lty=c(1,1))

crossPopDat<-read.table("popgen_results/Lingonberry_minus_RedCandy.msmc2.combined.msmc2.final.txt", header=TRUE)
plot(crossPopDat$left_time_boundary/mu*gen, 2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11),
     xlim=c(1000,30000000),ylim=c(0,1), type="n", xlab="Years ago", ylab="relative cross-coalescence rate")
lines(crossPopDat$left_time_boundary/mu*gen, 2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11), type="s")
#the single-individual (only using offspring haplotypes) resulted in x10 older estimate

#PSMC results 
minusDat_psmc<-read.table("popgen_results/Lingonberry_minus_5yr_psmc.0.txt", header=FALSE) %>%
  as_tibble()%>%
  rename(Time = V1, 
         Ne_tenthousands = V2)%>% 
  mutate(Ne = Ne_tenthousands * 10000, 
         ssp = "minus")%>% 
  select(Time, Ne, ssp)
RedCandyDat_psmc<-read.table("popgen_results/Lingonberry_RedCandy_5yr_psmc_2.0.txt", header=FALSE) %>%
  as_tibble()%>%
  rename(Time = V1, 
         Ne_tenthousands = V2)%>% 
  mutate(Ne = Ne_tenthousands * 10000, 
         ssp = "RedCandy")%>% 
  select(Time, Ne, ssp)
LingonberryDat_psmc <- rbind(minusDat_psmc, RedCandyDat_psmc)

LingonberryDat_psmc %>% 
  ggplot()+
  geom_step(aes(Time, Ne, colour = ssp))+
  scale_x_log10(limits = c(50000,11000000))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.position = "na")+
  scale_colour_manual(values = c(LGpal[1], LGpal[11]))+
  xlab("Years ago (g=5, u=0.3e10^-8)")+ ylab("Effective population size")



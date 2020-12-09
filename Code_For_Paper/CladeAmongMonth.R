library(stringr)
library(lubridate)
library(ggplot2)
library(vroom)
library(dplyr)
library(ggthemes)
library(tidyr)
​
#snps=read.delim("10SNPsSubset.txt")
snps=vroom("/Users/wenhui/Desktop/❤/Fall\ 2020/PUBH6885/Project/raw\ data/fullvcfOfInterest.tsv")
metadata=vroom('/Users/wenhui/Desktop/❤/Fall\ 2020/PUBH6885/Project/raw\ data/metadata_2020-11-04_09-37.tsv')
ID=str_extract(colnames(snps),"EPI_ISL_[0-9]+")
ID=ID[-c(1:9)]
genomes=snps[,10:length(colnames(snps))]
colnames(genomes)=ID
metadataOfInterest=metadata[which(metadata$gisaid_epi_isl %in% ID),]
genomes=genomes[,order(colnames(genomes))]
metadataOfInterest=metadataOfInterest[order(metadataOfInterest$gisaid_epi_isl),]
start=length(colnames(genomes))-3
end=length(colnames(genomes))
genomes=genomes[,-c(start:end)]
#genomes=rbind(genomes,metadataOfInterest$country)
genomes=genomes[,which(colnames(genomes) %in% metadataOfInterest$gisaid_epi_isl)]
monthList=month(as.POSIXlt(as.character(metadataOfInterest$date), format="%Y-%m-%d"))
genomes=rbind(genomes,monthList)
​
#with clade
genomesclade=rbind(genomes,metadataOfInterest$pangolin_lineage)
clade <- substr(genomesclade[12,], 1, 1)
genomesclade=rbind(genomesclade,clade)
genomesclade=genomesclade[-12,]
#12months
December <- genomesclade %>% select(which(genomesclade[11,]==12))
January <- genomesclade %>% select(which(genomesclade[11,]==1))
Feburary <- genomesclade %>% select(which(genomesclade[11,]==2))
March <- genomesclade %>% select(which(genomesclade[11,]==3))
April <- genomesclade %>% select(which(genomesclade[11,]==4))
May <- genomesclade %>% select(which(genomesclade[11,]==5))
June <- genomesclade %>% select(which(genomesclade[11,]==6))
July <- genomesclade %>% select(which(genomesclade[11,]==7))
August <- genomesclade %>% select(which(genomesclade[11,]==8))
September <- genomesclade %>% select(which(genomesclade[11,]==9))
October <- genomesclade %>% select(which(genomesclade[11,]==10))
November <- genomesclade %>% select(which(genomesclade[11,]==11))
#proprotion of A clade
prop_12A <- length(which(December[12,]=='A'))/length(December[12,])
prop_1A <- length(which(January[12,]=='A'))/length(January[12,])
prop_2A <- length(which(Feburary[12,]=='A'))/length(Feburary[12,])
prop_3A <- length(which(March[12,]=='A'))/length(March[12,])
prop_4A <- length(which(April[12,]=='A'))/length(April[12,])
prop_5A <- length(which(May[12,]=='A'))/length(May[12,])
prop_6A <- length(which(June[12,]=='A'))/length(June[12,])
prop_7A <- length(which(July[12,]=='A'))/length(July[12,])
prop_8A <- length(which(August[12,]=='A'))/length(August[12,])
prop_9A <- length(which(September[12,]=='A'))/length(September[12,])
prop_10A <- length(which(October[12,]=='A'))/length(October[12,])
prop_11A <- length(which(November[12,]=='A'))/length(November[12,])
#proprotion of B clade
prop_12B <- length(which(December[12,]=='B'))/length(December[12,])
prop_1B <- length(which(January[12,]=='B'))/length(January[12,])
prop_2B <- length(which(Feburary[12,]=='B'))/length(Feburary[12,])
prop_3B <- length(which(March[12,]=='B'))/length(March[12,])
prop_4B <- length(which(April[12,]=='B'))/length(April[12,])
prop_5B <- length(which(May[12,]=='B'))/length(May[12,])
prop_6B <- length(which(June[12,]=='B'))/length(June[12,])
prop_7B <- length(which(July[12,]=='B'))/length(July[12,])
prop_8B <- length(which(August[12,]=='B'))/length(August[12,])
prop_9B <- length(which(September[12,]=='B'))/length(September[12,])
prop_10B <- length(which(October[12,]=='B'))/length(October[12,])
prop_11B <- length(which(November[12,]=='B'))/length(November[12,])
#putting together
Clade_A <- c(prop_12A,prop_1A,prop_2A,prop_3A,prop_4A,prop_5A,prop_6A,prop_7A,prop_8A,prop_9A,prop_10A,prop_11A)
Clade_B <- c(prop_12B,prop_1B,prop_2B,prop_3B,prop_4B,prop_5B,prop_6B,prop_7B,prop_8B,prop_9B,prop_10B,prop_11B)
month <- c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov")
CladeData <- as.data.frame(cbind(month,Clade_A,Clade_B))
#tidy
NewcladeData <- CladeData %>%
  pivot_longer(cols = c(`Clade_A`, `Clade_B`),
               names_to = "Clade", values_to = "proportions")
NewcladeData$month=factor(NewcladeData$month,levels=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))
NewcladeData=NewcladeData[-23,]
NewcladeData=NewcladeData[-23,]
#plot
ggplot(NewcladeData,aes(month,proportions,fill=Clade))+
  geom_bar(stat="identity",position="stack")+
  ggtitle("Clade among Month")+
  theme_wsj()+
  scale_fill_wsj("rgby", "")+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  theme(axis.text.y = element_blank())
                                 
#with region
genomes=genomes[-11,]
genomesOfRegion=rbind(genomes,metadataOfInterest$region)
genomesOfRegion=rbind(genomesOfRegion,metadataOfInterest$pangolin_lineage)
clade2 <- substr(genomesOfRegion[12,], 1, 1)
genomesOfRegion=rbind(genomesOfRegion,clade2)
genomesOfRegion=genomesOfRegion[-12,]
#continents
Asia <- genomesOfRegion %>% select(which(genomesOfRegion[11,]=='Asia'))
NorthAmerica <- genomesOfRegion %>% select(which(genomesOfRegion[11,]=='North America'))
Europe <- genomesOfRegion %>% select(which(genomesOfRegion[11,]=='Europe'))
Oceania <- genomesOfRegion %>% select(which(genomesOfRegion[11,]=='Oceania'))
SouthAmerica <- genomesOfRegion %>% select(which(genomesOfRegion[11,]=='South America'))
Africa <- genomesOfRegion %>% select(which(genomesOfRegion[11,]=='Africa'))
#proprotion of A clade
prop_AsiaA <- length(which(Asia[12,]=='A'))/length(Asia[12,])
prop_NorthAmericaA <- length(which(NorthAmerica[12,]=='A'))/length(NorthAmerica[12,])
prop_EuropeA <- length(which(Europe[12,]=='A'))/length(Europe[12,])
prop_OceaniaA <- length(which(Oceania[12,]=='A'))/length(Oceania[12,])
prop_SouthAmericaA <- length(which(SouthAmerica[12,]=='A'))/length(SouthAmerica[12,])
prop_AfricaA <- length(which(Africa[12,]=='A'))/length(Africa[12,])
#proprotion of B clade
prop_AsiaB <- length(which(Asia[12,]=='B'))/length(Asia[12,])
prop_NorthAmericaB <- length(which(NorthAmerica[12,]=='B'))/length(NorthAmerica[12,])
prop_EuropeB <- length(which(Europe[12,]=='B'))/length(Europe[12,])
prop_OceaniaB <- length(which(Oceania[12,]=='B'))/length(Oceania[12,])
prop_SouthAmericaB <- length(which(SouthAmerica[12,]=='B'))/length(SouthAmerica[12,])
prop_AfricaB <- length(which(Africa[12,]=='B'))/length(Africa[12,])
#putting together
continentname <- c("Asia","NorthAmerica","Europe","Oceania","SouthAmerica","Africa")
A_Clade <- c(prop_AsiaA,prop_NorthAmericaA,prop_EuropeA,prop_OceaniaA,prop_SouthAmericaA,prop_AfricaA)
B_Clade <- c(prop_AsiaB,prop_NorthAmericaB,prop_EuropeB,prop_OceaniaB,prop_SouthAmericaB,prop_AfricaB)
continents <- as.data.frame((cbind(continentname,A_Clade,B_Clade)))
#tidy
continents <- continents %>%
  pivot_longer(cols = c(`A_Clade`, `B_Clade`),
               names_to = "Clade_conti", values_to = "proportion_conti")
#plots
ggplot(continents,aes(continentname,proportion_conti,fill=Clade_conti))+
  geom_bar(stat="identity",position="dodge")+
  ggtitle("Clades by continents")+
  theme_wsj()+
  scale_fill_wsj("rgby", "")+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  coord_flip()+
  theme(axis.text.x = element_blank())

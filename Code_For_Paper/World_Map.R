library("ggplot2")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(rgeos)
library(vroom)
library(dplyr)
library(plyr)
setwd("C:/Users/tyson/OneDrive/Desktop/Coronavirus Proteins/November4MachineLearning")
metadata=vroom("metadata_2020-11-04_09-37.tsv")
metadata$country[which(metadata$country=="USA")]="United States of America"
countryCases=as.data.frame(table(metadata$country))
colnames(countryCases)=c("admin","caseCount")
joined=full_join(world,countryCases)
theme_set(theme_bw())
countryAndClade=metadata[,c(7,19)]
countryAndClade$pangolin_lineage=substr(countryAndClade$pangolin_lineage,1,1)
BClade=countryAndClade[which(countryAndClade$pangolin_lineage=="B"),]
AClade=countryAndClade[which(countryAndClade$pangolin_lineage=="A"),]
ACladeTable=as.data.frame(table(AClade$country))
colnames(ACladeTable)=c("Country","A_Count")
BCladeTable=as.data.frame(table(BClade$country))
colnames(BCladeTable)=c("Country","B_Count")
cladejoin=full_join(ACladeTable,BCladeTable)
cladejoin[is.na(cladejoin)] <- 0
cladejoin$A_Proportion=cladejoin$A_Count/(cladejoin$A_Count+cladejoin$B_Count)
z=data.frame(ddply(countryAndClade,.(pangolin_lineage),summarise))
colnames(cladejoin)=c("admin","A_Count","B_Count", "A_Proportion")
joined=full_join(joined,cladejoin)
world <- ne_countries(scale = "medium", returnclass = "sf")

png(filename="World_A_Proportion.png",height=2000,width=4000,res=500)
ggplot(data = joined) +
  geom_sf(aes(fill = A_Proportion)) +
  #scale_color_gradient(low="#FF0000",high="#008000")
  scale_fill_gradientn(colours = rev(hcl.colors(7,palette="YlOrRd")),
                       #breaks = c(2, 4, 10, 100, 1000, 10000),
                       trans = "sqrt")
  #scale_fill_discrete()
  #scale_fill_brewer(palette = "RdYlGn")
  #scale_fill_viridis_c(option = "magma", trans = "sqrt")
dev.off()

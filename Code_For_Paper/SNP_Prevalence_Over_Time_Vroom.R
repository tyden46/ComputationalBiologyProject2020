install.packages("vroom", repos='http://cran.us.r-project.org')
library(stringr)
library(lubridate)
library(ggplot2)
library(vroom)
#snps=read.delim("10SNPsSubset.txt")
snps=vroom("fullvcfOfInterest.tsv")
metadata=vroom('../metadata_2020-11-04_09-37.tsv')
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
genomes=genomes[,which(colnames(genomes) %in% metadataOfInterest$gisaid_epi_isl)]
monthList=month(as.POSIXlt(as.character(metadataOfInterest$date), format="%Y-%m-%d"))
genomes=rbind(genomes,monthList)

for(x in 1:10){
  mutation=as.data.frame(genomes[,which(genomes[x,]==1)])
  decProp=length(which(mutation[11,]==12))/length(mutation[11,])
  janProp=length(which(mutation[11,]==1))/length(mutation[11,])
  febProp=length(which(mutation[11,]==2))/length(mutation[11,])
  marProp=length(which(mutation[11,]==3))/length(mutation[11,])
  aprProp=length(which(mutation[11,]==4))/length(mutation[11,])
  mayProp=length(which(mutation[11,]==5))/length(mutation[11,])
  junProp=length(which(mutation[11,]==6))/length(mutation[11,])
  julProp=length(which(mutation[11,]==7))/length(mutation[11,])
  augProp=length(which(mutation[11,]==8))/length(mutation[11,])
  sepProp=length(which(mutation[11,]==9))/length(mutation[11,])
  octProp=length(which(mutation[11,]==10))/length(mutation[11,])
  novProp=length(which(mutation[11,]==11))/length(mutation[11,])
  myDF=as.data.frame(cbind(c("Dec",
                             "Jan",
                             "Feb",
                             "Mar",
                             "Apr",
                             "May",
                             "Jun",
                             "Jul",
                             "Aug",
                             "Sep",
                             "Oct",
                             "Nov"),c(decProp,
                                      janProp,
                                      febProp,
                                      marProp,
                                      aprProp,
                                      mayProp,
                                      junProp,
                                      julProp,
                                      augProp,
                                      sepProp,
                                      octProp,
                                      novProp)),stringsAsFactors = FALSE)
  colnames(myDF)=c("Month","Proportion")
  myDF$Month=factor(myDF$Month,levels=c("Dec",
                                        "Jan",
                                        "Feb",
                                        "Mar",
                                        "Apr",
                                        "May",
                                        "Jun",
                                        "Jul",
                                        "Aug",
                                        "Sep",
                                        "Oct",
                                        "Nov"))
  assign(paste("ggplot",x,sep=""),ggplot(data=myDF, aes(x=Month, y=as.numeric(Proportion), group=1)) +
           geom_line()+
           geom_point()+
           ggtitle(paste("Frequency Over Time of","Mutations at Basepair",snps$POS[x],sep='\n'))+
           theme(plot.title = element_text(hjust = 0.5)))
}
library(cowplot)
png("frequencyPlots.png",width=2000,height=1000)
plot_grid(ggplot1,
          ggplot2,
          ggplot3,
          ggplot4,
          ggplot5,
          ggplot6,
          ggplot7,
          ggplot8,
          ggplot9,
          ggplot10,nrow=2,ncol=5)
dev.off()
save.image(file="newPrevalence2.RData")

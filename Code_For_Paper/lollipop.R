library(g3viz)
library(vroom)
library(Rfast)
library(ggplot2)
setwd("C:\\Users\\tyson\\OneDrive\\Desktop\\Coronavirus Proteins\\November4MachineLearning")
table=vroom("11_4_vcf_0.00005maf.tsv")
myData=as.data.frame(table$POS)
myData$count=rowSums(table[,11:163342])
colnames(myData)=c("Position","Count")
myData=myData[-c(1,2,3),]
myData=myData[-217,]
myData=myData[-214,]
png(filename="Lollipop.png",height=2000,width=10000,res=800)
ggplot(myData, aes(x=Position, y=Count)) +
  geom_segment( aes(x=Position, xend=Position, y=0, yend=Count)) +
  geom_point( color="orange", size=2) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("Base Pair Position With Mutation") +
  ylab("Number of Genomes With Mutation")
dev.off()
mutTable=as.data.frame(rep("symbol",571))
mutTable$Hugo_Symbol="SPIKE_SARS2"
mutTable$Chromosome="1"
mutTable$Start_Position=table$POS
mutTable$End_Position=table$POS
mutTable$Strand="+"
mutTable$Variant_Classification="SNP"
mutTable$Variant_Type="SNP"
mutTable$Reference_Allele=table$REF 
mutTable$Tumor_Seq_Allele1=table$REF      
mutTable$Tumor_Seq_Allele2=table$REF      

mutTable$HGVSp="x"
mutTable$HGVSp_Short="x"           
mutTable$COSMIC="x"
mutTable$Mutation_Class="x"        
mutTable$AA_Position="x"  

# generate chart
g3Lollipop(mutTable,
           gene.symbol = "SPIKE_SARS2",
           plot.options = chart.options,
           btn.style = "blue",
           output.filename = "default_theme")

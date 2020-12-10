library(caret)
library(kernlab)
library(stringr)
library(pls)
library(vroom)
#setwd("C:\\Users\\tyson\\OneDrive\\Desktop\\Coronavirus Proteins\\November4MachineLearning\\Tutorial")
table=vroom("ncbivcf0.05MAF.tsv")
#Read in the table and clean the data
snps=as.data.frame(t(table))
colnames(snps)=snps[2,]
snps=snps[-c(1:9),]
#Add a lineages column (this info is from GISAID)
lineages=c(rep("A",23),rep("B",19))
#Additional data cleaning
colnames(snps)=table$POS
DF=snps
DF[sapply(DF, is.character)] <- lapply(DF[sapply(DF, is.character)], 
                                       as.factor)
snps=DF
snps$lineages=factor(lineages)

#Determine which SNPs come from which clade
numSNPS=snps
numSNPS[] <- lapply(snps[,1:length(colnames(snps))-1], as.numeric)
ASums=colSums(numSNPS[which(snps$lineages=="A"),])
BSums=colSums(numSNPS[which(snps$lineages=="B"),])
ASums=length(which(snps$lineages=="A"))/ASums
BSums=length(which(snps$lineages=="B"))/BSums

library(ComplexHeatmap)
library(RColorBrewer)

#Set up for our heatmap
split = rep(1:3, each = 6)
ha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE),
  foo = anno_block(gp = gpar(fill = 2:4), labels = c("Alpha","Beta","Coil"))
)
#Create proportions table with A clade, B clade,
#and ratio between the two for each SNP
proportions=as.data.frame(rbind(ASums[-length(colnames(snps))],BSums[-length(colnames(snps))]))
#Calculate ratio
newVec=c()
for(x in 1:length(proportions[1,])){
  prop=max(proportions[1,x],proportions[1,x])/min(proportions[1,x],proportions[2,x])
  newVec=append(newVec,prop)
}
proportions[3,]=newVec
#Find top 10 SNPs that differ the most between the clades
#ie the 10 that have the highest ratio
top10=which(proportions[3,] %in% sort(proportions[3,],decreasing = TRUE)[1:10])
proportions=proportions[,order(proportions[3,], decreasing=TRUE)]
top10=which(proportions[3,] %in% sort(proportions[3,],decreasing = TRUE)[1:10])
#Generate and save the heatmap
png(filename="HeatmapNewProp.png", height=1000, width=6000, res=750)
Heatmap(as.matrix(proportions[1:2,top10]), name = "mat2",
                   #column_split = split,
                   #top_annotation = ha,
                   column_title = NULL,col=colorRampPalette(brewer.pal(9,"Blues"))(100),cluster_rows = FALSE,cluster_columns = TRUE,
                   column_names_gp = gpar(fontsize = 8),)
dev.off()
#Create a dataframe with the top10 snps
top10SNPs=snps[,c(top10,length(colnames(snps)))]
#Split into train and test set
test_index <- createDataPartition(top10SNPs$lineages, times = 1, p = 0.2, list = FALSE)    # create a 20% test set
testSet <- top10SNPs[test_index,]
trainingSet <- top10SNPs[-test_index,]
#Train our model
myfit=train(lineages ~ ., data=trainingSet,
              method = 'multinom')
#Predict on the test set
prediction=predict(object=myfit, newdata=testSet)
testSet$lineages=factor(testSet$lineages)
#Generate confusion matrix
confusionMatrix(prediction, reference = testSet$lineages)$overall["Accuracy"]
xtab <- table(prediction, testSet$lineages)
confusionMatrix(xtab)

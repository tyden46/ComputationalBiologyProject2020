#install.packages("ComplexHeatmap", repos='http://cran.us.r-project.org')
#install.packages("RColorBrewer", repos='http://cran.us.r-project.org')
library(caret)
library(fastAdaboost)
library(kernlab)
library(nnet)
library(stringr)
library(pls)
library(doParallel)
library(vroom)
setwd("C:\\Users\\tyson\\OneDrive\\Desktop\\Coronavirus Proteins\\November4MachineLearning")
table=vroom("2000vcf.tsv")
snps=as.data.frame(t(table))
colnames(snps)=snps[2,]
snps=snps[-c(1:9),]

snps=snps[-which(is.na(str_extract(row.names(snps),"EPI_ISL_[0-9]*"))),]
row.names(snps)=str_extract(row.names(snps),"EPI_ISL_[0-9]*")
metadata=vroom("metadataIDAndClade.tsv")
metadata=metadata[which(metadata$gisaid_epi_isl %in% row.names(snps)),]
snps=snps[order(row.names(snps)),]
metadata=metadata[order(metadata$gisaid_epi_isl),]
lineages=metadata$pangolin_lineage
# sort variance, grab index of the top 50

DF=snps
DF[sapply(DF, is.character)] <- lapply(DF[sapply(DF, is.character)], 
                                       as.factor)
savesnps=snps
snps=DF
variances <- apply(X=snps, MARGIN=2, FUN=var)
variances=as.data.frame(variances)
snps=snps[,which(variances$variances!=0)]
snps=snps[which(rownames(snps) %in% metadata$gisaid_epi_isl),]
snps$lineages=substring(lineages,1,1)
snps=snps[sort(c(which(str_detect(lineages,"B")),which(str_detect(lineages,"A")))),]
snps$lineages=factor(snps$lineages)

#test_index <- createDataPartition(snps$lineages, times = 1, p = 0.2, list = FALSE)    # create a 20% test set
#testSet <- snps[test_index,]
#trainingSet <- snps[-test_index,]
#set.seed(2)
#cl <- makePSOCKcluster(24)
#registerDoParallel(cl)
#lvq_fit=train(lineages ~ ., data=trainingSet,
#                method = 'lvq')
#stopCluster(cl)
#prediction=predict(object=lvq_fit, newdata=testSet)
#testSet$lineages=factor(testSet$lineages)
#confusionMatrix(prediction, reference = testSet$lineages)$overall["Accuracy"]
#importance=varImp(lvq_fit)
#png(file="importance.png",height=1500,width=500)
#plot(importance)
#dev.off()
#xtab <- table(prediction, testSet$lineages)
#confusionMatrix(xtab)
#save.image(file="myData.RData") 

numSNPS=snps
numSNPS[] <- lapply(snps[,1:length(colnames(snps))-1], as.numeric)
ASums=colSums(numSNPS[which(snps$lineages=="A"),])
BSums=colSums(numSNPS[which(snps$lineages=="B"),])
ASums=length(which(snps$lineages=="A"))/ASums
BSums=length(which(snps$lineages=="B"))/BSums
save.image(file="64.RData")
# M <- cor(as.data.frame(rbind(ASums[-length(colnames(snps))],BSums[-length(colnames(snps))])))
# png("corplotTes2.png",height=2000,width=2000)
# corrplot(M, method = "circle")
# dev.off()
library(ComplexHeatmap)
library(RColorBrewer)


split = rep(1:3, each = 6)
ha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE),
  foo = anno_block(gp = gpar(fill = 2:4), labels = c("Alpha","Beta","Coil"))
)
proportions=as.data.frame(rbind(ASums[-length(colnames(snps))],BSums[-length(colnames(snps))]))
newVec=c()
for(x in 1:length(proportions[1,])){
  prop=max(proportions[1,x],proportions[1,x])/min(proportions[1,x],proportions[2,x])
  newVec=append(newVec,prop)
}
proportions[3,]=newVec
write.csv(proportions, "proportions.csv")
top10=which(proportions[3,] %in% sort(proportions[3,],decreasing = TRUE)[1:10])
save.image(file="81.RData")
pos=read.delim("pos.txt")
colnames(proportions)[2:572]=pos$POS
proportions=proportions[,order(proportions[3,], decreasing=TRUE)]
top10=which(proportions[3,] %in% sort(proportions[3,],decreasing = TRUE)[1:10])
png(filename="HeatmapNewProp.png", height=1000, width=6000, res=750)
Heatmap(as.matrix(proportions[1:2,top10]), name = "mat2",
                   #column_split = split,
                   #top_annotation = ha,
                   column_title = NULL,col=colorRampPalette(brewer.pal(9,"Blues"))(100),cluster_rows = FALSE,cluster_columns = TRUE,
                   column_names_gp = gpar(fontsize = 8),)
dev.off()
save.image(file="89.RData")
write.csv(snps,"snps.csv")
top10SNPs=snps[,c(top10,length(colnames(snps)))]
write.csv(top10SNPs, "top10SNPs.csv")
test_index <- createDataPartition(top10SNPs$lineages, times = 1, p = 0.2, list = FALSE)    # create a 20% test set
testSet <- top10SNPs[test_index,]
trainingSet <- top10SNPs[-test_index,]

knn_fit=train(lineages ~ ., data=trainingSet,
              method = 'multinom')
prediction=predict(object=knn_fit, newdata=testSet)
testSet$lineages=factor(testSet$lineages)
confusionMatrix(prediction, reference = testSet$lineages)$overall["Accuracy"]
xtab <- table(prediction, testSet$lineages)
confusionMatrix(xtab)
save.image(file="102.RData")
colnames(proportions)=table$POS
top10=which(proportions[3,] %in% sort(proportions[3,],decreasing = TRUE)[1:60])
top10SNPs=snps[,c(top10,length(colnames(snps)))]
colnames(top10SNPs)=c(table$POS,"Lineage")
proportions=proportions[,order(proportions[3,],decreasing = TRUE)]
write.csv(top10SNPs,file="top100.csv")

#install.packages("randomForest")

library(randomForest)
set.seed(123456)
setwd("../10/")
getwd()
inputFile="../08/normalize.txt"

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))

# Use make.names to create valid column names
colnames(data) <- make.names(colnames(data))

rf=randomForest(as.factor(group)~., data=data, ntree=500)
pdf(file="forest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()

optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)

importance=importance(x=rf2)
importance

pdf(file="geneImportance.pdf", width=6.2, height=7)
varImpPlot(rf2, main="TOP30")
dev.off()

rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
#rfGenes=names(rfGenes[rfGenes>2])
rfGenes=names(rfGenes[1:50])
write.table(rfGenes, file="rfGenes.txt", sep="\t", quote=F, col.names=F, row.names=F)

sigExp=t(data[,rfGenes])
sigExpOut=rbind(ID=colnames(sigExp),sigExp)
write.table(sigExpOut, file="rfGeneExp.txt", sep="\t", quote=F, col.names=F)
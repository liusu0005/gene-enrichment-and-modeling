#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")

library(limma)
library(ggpubr)
library(pROC)
setwd("../12/")
getwd()

expFile="../sample_data/geneMatrix.txt" 
conFile="../sample_data/sample1.txt"
treatFile="../sample_data/sample2.txt"
geneFile="../11/intersect.txt"

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
rt=avereps(data)

qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
	rt[rt<0]=0
	rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

con=read.table(conFile, header=F, sep="\t", check.names=F)
treat=read.table(treatFile, header=F, sep="\t", check.names=F)
conData=data[,as.vector(con[,1])]
treatData=data[,as.vector(treat[,1])]
data=cbind(conData, treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

Type=c(rep("con",conNum), rep("treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="test.normalize.txt", sep="\t", quote=F, col.names=F)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),,drop=F]

Type=c(rep("Con",conNum), rep("Treat",treatNum))
my_comparisons=list()
my_comparisons[[1]]=levels(factor(Type))

newGeneLists=c()
outTab=data.frame()
for(i in row.names(data)){
	#data[i,][data[i,]>quantile(data[i,], 0.99)]=quantile(data[i,], 0.99)
	rt1=data.frame(expression=data[i,], Type=Type)

	boxplot=ggboxplot(rt1, x="Type", y="expression", color="Type",
				      xlab="",
				      ylab=paste(i, "expression"),
				      legend.title="",
				      palette = c("blue", "red"),
				      add = "jitter")+ 
		stat_compare_means(comparisons = my_comparisons)

	pdf(file=paste0("boxplot.",i,".pdf"), width=5, height=4.5)
	print(boxplot)
	dev.off()
}
expFile="test.normalize.txt"  
geneFile="intersect.txt" 

rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
y=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y=ifelse(y=="con", 0, 1)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

for(x in as.vector(geneRT[,1])){
	roc1=roc(y, as.numeric(rt[x,]))
	ci1=ci.auc(roc1, method="bootstrap")
	ciVec=as.numeric(ci1)
	pdf(file=paste0("ROC.",x,".pdf"), width=5, height=5)
	plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=x)
	text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
	dev.off()
}
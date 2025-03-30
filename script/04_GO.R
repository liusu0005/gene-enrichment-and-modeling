#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

library("org.Hs.eg.db")
setwd("../04/")
getwd()
inputFile="../02/GEO_diff.txt" 
rt=read.table(inputFile,sep="\t",check.names=F,header=F)  
genes=as.vector(rt[,1])      
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)

out=cbind(rt,entrezID=entrezIDs)
colnames(out)[1]="Gene"
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)

pvalueFilter=0.05
qvalueFilter=1

rt=read.table("id.txt",sep="\t",header=T,check.names=F)  
rt=rt[is.na(rt[,"entrezID"])==F,]                        
gene=rt$entrezID

colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)

showNum=10
if(nrow(GO)<30){
	showNum=nrow(GO)
}

pdf(file="barplot.pdf",width = 9,height =7)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
		
pdf(file="bubble.pdf",width = 9,height =7)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()
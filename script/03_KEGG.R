
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
setwd("../03")
getwd()


library("org.Hs.eg.db")                                          
rt=read.table("../02/GEO_diff.txt",sep="\t",check.names=F,header=T)    
rt <- rt[!is.na(rownames(rt)), ]
genes=as.vector(rt[,1])
genes <- as.character(genes)  # Ensure character type
genes <- genes[!is.na(genes) & genes != ""]  # Remove NA and empty strings
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)      
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)       
pvalueFilter=0.05
qvalueFilter=1

rt=read.table("id.txt",sep="\t",header=T,check.names=F)          
#rt=rt[is.na(rt[,"entrezID"])==F,]                                
gene=rt$entrezID
geneFC=rt$logFC
names(geneFC)=gene

colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)   #????????
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$symbol[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)                          #???渻??????

showNum=30
if(nrow(KEGG)<30){
	showNum=nrow(KEGG)
}

pdf(file="barplot.pdf",width = 10,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
dev.off()

pdf(file="bubble.pdf",width = 10,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)
dev.off()

#Ȧͼ
pdf(file="circos.pdf",width = 8,height = 5.5)
kkx=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(kkx, foldChange=geneFC, circular = TRUE, colorEdge = TRUE)
dev.off()

#Enrichment Map
#pdf(file="emapplot.pdf",width = 7,height = 6)
#emapplot(kk, showCategory = nrow(KEGG), color = colorSel)
#dev.off()

###################above are original codes, below from GPT##########

# Generate the term similarity matrix
kk <- pairwise_termsim(kk)

# Now create the enrichment map plot with optimized parameters
pdf(file = "emapplot.pdf", width = 10, height = 8) # Increased the plot size
emapplot(kk, 
         showCategory = nrow(KEGG), 
         color = colorSel, 
         layout = "nicely",          # Use a layout that minimizes overlaps
         min_edge = 0.2)            # Set a minimum edge value to reduce clutter
         #node_label_size = 3)        # Adjust the label size for better readability
dev.off()
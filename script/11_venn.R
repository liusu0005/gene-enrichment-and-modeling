#install.packages("VennDiagram")

library(VennDiagram)
setwd("../11/")
files=dir()
files=grep("txt",files,value=T)
geneList=list()

for(i in 1:length(files)){
    inputFile=files[i]
	if(inputFile=="intersect.txt"){next}
    rt=read.table(inputFile,header=F)
    header=unlist(strsplit(inputFile,"\\.|\\-"))
    geneList[[header[1]]]=as.vector(rt[,1])
    uniqLength=length(unique(as.vector(rt[,1])))
    print(paste(header[1],uniqLength,sep=" "))
}

venn.plot=venn.diagram(geneList,filename=NULL,fill=rainbow(length(geneList)))
pdf(file="venn.pdf",width=6,height=6)
grid.draw(venn.plot)
dev.off()

intersectGenes=Reduce(intersect,geneList)
write.table(file="intersect.txt",intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)
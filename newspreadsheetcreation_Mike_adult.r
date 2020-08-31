source("https://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#install.packages("BiocManager")
#BiocManager::install("reshape2")
#BiocManager::install("rlist")
library(reshape2)
library(rlist)
library(qlcMatrix)
library(Seurat)
library(ggplot2)

filedir<-"~/Desktop/Astrocyte_Adult_scRNAseq_SmartSeq2"
filenames <- list.files(filedir, pattern="*.gz", full.names=TRUE)
filenamessmall<- list.files(filedir,pattern="*.gz",full.names=F)
am <- as.sparse(read.table(filenames[[35]],header = F, row.names=1,col.names=c("Gene",strsplit(filenames[[35]],"/")[[1]][6])), stringsAsFactors = F)

sm <- as.sparse(read.table(filenames[[1]],header = F, row.names=1,col.names=c("Gene",strsplit(filenames[[1]],"/")[[1]][6])), stringsAsFactors = F)

for (i in filenames[s==F]){
  am <- as.sparse(read.table(i,header= F, row.names=1,col.names=c("Gene",strsplit(i,"/")[[1]][6])))
  if (identical(rownames(am),rownames(sm))){
    sm<-cbind(sm,am)
  }
}

write.table(sm,"astro_adultctmatrix.txt",quote=F,row.names=T,col.names=T)
library(org.Mm.eg.db)
ids=sapply(rownames(sm),function(x)strsplit(x,"[.]")[[1]][1])
fromKey="ENSEMBL"
toKey="SYMBOL"
db=org.Mm.eg.db
selRes<-AnnotationDbi::select(db,keys=ids,keytype=fromKey,columns=c(fromKey,toKey))

##this is really used to get particular keys of the keytype "ensembl" (or any other type)

x=selRes[match(ids,selRes[,1]),1:2]
good_names<-rownames(sm)[!is.na(x[,2])]
rm<-sm[rownames(sm)%in%good_names,]

##am<-am[match(rownames(am),rownames(sm)),]

##l<-merge(d,c,by="row.names",all.x=TRUE)
ids=sapply(rownames(rm),function(x)strsplit(x,"[.]")[[1]][1])
rownames(rm)<-ids
fromKey="ENSEMBL"
toKey="SYMBOL"
db=org.Mm.eg.db
selRes<-AnnotationDbi::select(db,keys=ids,keytype=fromKey,columns=c(fromKey,toKey))
x=selRes[match(ids,selRes[,1]),1:2]


##Need to combine rows that have the same x[,2]
sum(isUnique(x[,2]))
keep<-x[isUnique(x[,2]),]
repeats<-x[!isUnique(x[,2]),]
nm<-rm[rownames(rm)%in%keep$ENSEMBL,]
names(rownames(nm))<-NULL
keys=character()
for(i in 1:nrow(repeats)){
  y=repeats$ENSEMBL[i]
  if (y %in% keys){
    next
  }
  a=repeats$SYMBOL[i]
  newkeys=(repeats[repeats$SYMBOL==a,1])
  keys=c(keys,newkeys)
  total=Matrix::colSums(rm[newkeys,])
  nm<-rbind(nm,total)
  rownames(nm)[nrow(nm)]<-newkeys[1]
}

dim(nm)
all(isUnique(rownames(nm)))
sum(isUnique(rownames(nm)))

ids=rownames(nm)
fromKey="ENSEMBL"
toKey="SYMBOL"
db=org.Mm.eg.db
selRes<-AnnotationDbi::select(db,keys=ids,keytype=fromKey,columns=c(fromKey,toKey))
x=selRes[match(ids,selRes[,1]),1:2]

identical(x[,1],rownames(nm))
rownames(nm)<-x$SYMBOL
all(isUnique(x$SYMBOL))
nm<-nm[,-1]
which(nm[,1]==45775)
which(nm[,5]==max(nm[,5]))
astroadult <- CreateSeuratObject(counts = nm)

##to add a column with gene names to the final spreadsheet



##match finds positions of the first match between the rownames
##(ensemblIDs from DESeq table) and the first column of selRes, which is ensembl
##ID's found in the gene names table from org.Mm.eg.db. 
##so x returns selRes with only the first instance of a given 
##ensemblID and each symbol to go with it. 

astroadult[["percent.mt"]] <- PercentageFeatureSet(object = astroadult, pattern = "^mt-")
VlnPlot(object = astroadult, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot2 <- FeatureScatter(object = astroadult, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

astroadult <- NormalizeData(
object = astroadult,
normalization.method = "LogNormalize",
scale.factor = 1000)


astroadult <- FindVariableFeatures(object = astroadult, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = astroadult), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = astroadult)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(x = astroadult)
astroadult <- ScaleData(object = astroadult, features = all.genes)
astroadult <- RunPCA(object = astroadult, features = VariableFeatures(object = astroadult))
print(x = astroadult[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(object = astroadult, reduction = "pca")
ElbowPlot(object = astroadult)

astroadult <- FindNeighbors(object = astroadult, dims = 1:10)
astroadult <- FindClusters(object = astroadult, resolution = 0.5)
head(x = Idents(object = astroadult), 5)
reticulate::py_install(packages ='umap-learn')

save(astroadult,list=ls(),file="astroadult.RData")
load("~/Desktop/astroadult.RData")
astroadult <- RunUMAP(object = astroadult, dims = 1:10)
DimPlot(object = astroadult, reduction = "umap",label=T)

VlnPlot(object = astroadult, features = c("Aldh1l1", "Cx3cr1")) #0,1,2 are astro, #7 is micro
VlnPlot(object = astroadult, features = c("Nrxn1", "Pdgfra")) #3 maybe opcs, 4/6 neurons
VlnPlot(object = astroadult, features = c("Pecam1", "Olig2")) #5 Endothelial, 8 oligo

VlnPlot(object = astroadult, features = c("Aldh1l1", "Cx3cr1"), slot = "counts", log = TRUE)

FeaturePlot(object = astroadult, features = c("Aldh1l1","Cx3cr1","Nrxn1","Pdgfra","Olig2","Pecam1"))
save(astroadult,list=ls(),file="astroadult_full.RData")

#astrocytes only
astro<-WhichCells(object = astroadult, idents = c("0","1","2"))

aa<-subset(x=astroadult,idents=c("0","1","2"))
aa <- FindVariableFeatures(object = aa, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = aa), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = aa)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(x = aa)
aa <- ScaleData(object = aa, features = all.genes)
aa <- RunPCA(object = aa, features = VariableFeatures(object = aa))
print(x = aa[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(object = aa, reduction = "pca")
ElbowPlot(object = aa)

aa <- FindNeighbors(object = aa, dims = 1:10)
aa <- FindClusters(object = aa, resolution = 0.5)
head(x = Idents(object = aa), 5)


sums<-Matrix::rowSums(aa@assays$RNA@counts)
keep<-names(sums[sums!=0])

FilterGenes <-
function (object, min.value=1, min.cells = 0, genes = NULL) {
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("FilterGenes"))]
    object <- Seurat:::SetCalcParams(object = object, calculation = "FilterGenes", ... = parameters.to.store)
    genes.use <- rownames(object@data)
    
    if (!is.null(genes)) {
        genes.use <- intersect(genes.use, genes)
        object@data <- object@data[genes.use, ]
        return(object)
    } else if (min.cells > 0) {
        num.cells <- Matrix::rowSums(object@data > min.value)
        genes.use <- names(num.cells[which(num.cells >= min.cells)])
        object@data <- object@data[genes.use, ]
        return(object)
    } else {
        return(object)
    }
}

#aa <- FilterGenes(object = aa,min.value = 1, min.cells = 30)
#astro
##allgenes<-list of gene names you care about
#using scaledata which was scaled using all genes expressed in astrocytes
counts<-aa@assays$RNA@scale.data #could use data instead for unscaled counts, or raw.data
genes<-read.table("~/Desktop/UCSF/Mike_Astrocyte analysis/genelist.txt")
genes<-as.character(genes[,1])
newgenes<-genes[genes%in%keep]
genes<-newgenes
counts<-counts[rownames(counts)%in%genes,]

#scale the object for just astrocyte genes



require(Matrix)
lcorr2<-function(x,y=1){
    x<-t(x)
    j=dim(x)[2]
    if(length(y)<2) {
        y=as(Matrix(x[,j]),"dgCMatrix")
    } else {y=y}
    allcor<-corSparse(x,y)
    allcor<-as.data.frame(allcor)
    rownames(allcor)<-colnames(x)
    allcor
}


allcorfun<-function(counts,expression_threshold=-Inf){
  df<-data.frame(1:nrow(counts))
  rownames(df)<-rownames(counts)
  for (i in 1:nrow(counts)){
    genename=rownames(counts)[i]
    expressors=counts[i,]
    expressors=expressors>expression_threshold   #only cells that express that gene will be "True"
    y=as(Matrix(counts[i,expressors]),"dgCMatrix")
    s<-lcorr2(counts[,expressors],y)
    df[,i]<-s
    colnames(df)[i]<-genename
  }
  df
}


library(RColorBrewer)
l<-allcorfun(counts,expression_threshold=-Inf)
l_mat<-data.matrix(l)
#cc <- rainbow(ncol(l_mat), start = 0, end = .3)
#l_heatmap <- heatmap(l_mat, col = cm.colors(256), scale="column", margins=c(5,10))
#jBuPuFun <- colorRampPalette(brewer.pal(n = 9, "BuPu"))
#paletteSize <- 100
#jBuPuPalette <- jBuPuFun(paletteSize)

#col<-colorRampPalette(c("blue","black","red"))(20)
#heatmap(x=l_mat,col=col,symm=T)
#ggsave(filename="cormat_onlyposcells.png",plot=last_plot(),device="png",path="~/Desktop/Mike")

library(corrplot)

l_mat[which(l_mat>0.9)]<-0
corrplot(l_mat,type = "upper", order = "hclust", tl.col = "black", tl.srt = 45,title="Pairwise astrocyte gene correlation_adult",diag=T,tl.cex=.7,cl.lim=c(-0.3,0.3),col=col3(10))


ggsave(filename="cormat.png",plot=last_plot(),device="png",path="~/Desktop/Mike")


#Now make a scatter plot for two genes at a time: (should I be using raw counts, or scaled counts?)







#genes=c("Vamp1","Kcnn2")
#newcounts<-Matrix::colSums(counts[rownames(counts)%in%genes,])
#counts<-rbind(counts,"Adra1"=newcounts)




#2=dep

#make a scatterplot showing the counts for each gene on one axis


#add colors for control vs deprived. 
#library(ggplot2)




#cols.use=c("grey","blue","red","yellow"))

#Note: from raw data, I first normalized the data, then filtered genes, then scaled while regressing out nUMI cause that's very different between my samples. that way cells with more transcripts overall will look similar

#allgenes<-c("Vamp2","Vamp3","Vamp4","Kcna2","Kcna6","Kcnb1","Kcnc4","Kcnd2","Kcnd3","Kcnj10","Kcnj16","Kcnk10","Kcnk2","Kcnn2","Kcnn3","Gria2","Grid1","Grid2","Grik5","Grin1","Grin2b","Grin2c","Grin3a","Grina","Gabbr1","Gabbr2","Gabra2","Gabra4","Gabrb1","Adora1","Adora2a","Adora2b","Adra1a","Adra1b","Adra2a","Adrb1","Adrbk1","Adrbk2","Slc1a2","Slc1a3")
#allgenescombined<-c("Grid","Nmda","Gabbr","Gabra","Adra1")
#genes[order(match(genes,colnames(m)))] ##how to order genes by the colnames of m
allgenes<-genes
makeplots<-function(cmtx=counts,s=aa,expression_threshold=-Inf,cormat=l_mat,genes=allgenes){
  m<-cmtx[rownames(cmtx)%in%c(genes),]
  m<-t(m)
  m<-as.data.frame(m)
  for (i in 1:(ncol(m)-1)){
    g1=colnames(m)[i]
    for (j in (i+1):(ncol(m)-1)){
      g2=colnames(m)[j]
      if (abs(l_mat[g1,g2])<0.1){next}
      partialcounts<-m[m[i]>expression_threshold&m[j]>expression_threshold,c(i,j,ncol(m))]
      scplot<-ggplot(partialcounts, aes_string(x=g1,y=g2)) +
        geom_point()
      ggsave(filename=paste0(g1,"vs",g2,".png"),plot=last_plot(),device="png",path="~/Desktop/Astro/Final_02")
    }
  }
}


makeplots()
  
##Requires that you only include your filtered genes in counts. Filter before running makeplots.

quantile(l_mat,probs=seq(0,1,0.05))

##For cleaning up astrocyte names
orig<-names(astroadult$orig.ident)
s1 <- sapply(strsplit(orig, split='_', fixed=TRUE), function(x) (x[1]))
head(s1)
astroadult$cellnames<-s1
meta<-read.table("Astrocyte_Adult_scRNAseq_SmartSeq2/GSE114000_series_matrix_metadata.txt",header=F, sep="\t",fill=T,stringsAsFactors = F)
meta<-meta[1:14,]
meta[1,1:5]
meta<-meta[c(2,11),]
meta<-meta[,2:ncol(meta)]
rownames(meta)<-c("Cellname","Region")
s1<-as.character(meta[2,])
s1 <- sapply(strsplit(s1, split=': ', fixed=TRUE), function(x) (x[2]))
meta[2,]<-s1
region<-meta[2,]
names(region)<-meta[1,]
region<-as.factor(region)

cellnames<-astroadult$cellnames
names(cellnames)<-NULL
region<-region[match(cellnames,names(region))]
all(cellnames%in%names(region))

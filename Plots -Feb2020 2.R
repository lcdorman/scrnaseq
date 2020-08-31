library(Seurat)

library(ggplot2); theme_set(theme_classic())

#Load MG only (LD) file called "mgAVM02"
#load("/Users/whippoorwill/Dropbox (Anna Molofsky Lab)/2020-Leah-barrelCortex/Manuscript data/Sequencing/LD_AVM02/Data/mgAVM02_filtered_Feb2020.RData")


#Write a function that will print out a graph given a metadata or feature name, seurat object, and split.by/group.by parameters
#sobject = seurat object name
#graphtype = "violin","feature","dim"
#feature = gene name or metadata column to feature
#group = only for violin; metadata to group by
#split = metadata column to split on
#cellnames = vector: metadata column, subset to include (any length); required for violin plot
#namecard = token to include in the plot name (usually the name of the sobject)
setwd("~/Desktop/plots")

PrintSeuratGraph = function(namecard = "a",sobject,graphtype = "feature",feature = NULL,group = NULL,split=NULL,cellnames=NULL){
  if (!is.null(cellnames)){
    Idents(sobject) = cellnames[1]
    cells = colnames(sobject)[Idents(sobject) %in% cellnames[2:length(cellnames)]]} 
  else {cells = cellnames}
  if (graphtype == "feature"){
    graph = FeaturePlot(sobject,features = feature,split.by = split, cells = cells)
  }
  if (graphtype == "violin"){
    graph = VlnPlot(sobject,features = feature, pt.size = 0.1, idents = cellnames[2:length(cellnames)],group.by = group, split.by = split)
  }
  if (graphtype == "dim"){
    graph = DimPlot(sobject,cells = cells, group.by = group, split.by = split)
    
  }
  name = paste0(feature,"_",graphtype,namecard,".eps")
  graph
  setEPS()
  postscript(name)
  print(graph)
  dev.off()
}

#should print a graph called "Ptprc_featureboth.RData" showing expression of ptprc split by id (SvL) using every cell

features = c("nCount_RNA","percent.mito","Mertk","Axl","Tyrobp","P2ry12","P2ry13","Tmem119","Trem2","Cx3cr1","Spp1","Hexb","Fcrls","C3","Ptprc","Itgam","Csf1r","Cd68","Ifitm3","Ctsb","Ctsd")
all(features %in% rownames(mgAVM02))
for(feature in features){
  PrintSeuratGraph(namecard = "MGonly",sobject=mgAVM02,graphtype = "feature",feature = feature)
}

#for(feature in features){
  #PrintSeuratGraph(namecard = "newthresholds_allcells",sobject=MG,graphtype = "feature",feature = feature)
#}
for(feature in features){
  PrintSeuratGraph(namecard = "CP5",sobject=mgAVM02,graphtype = "feature",feature = feature,cellnames = c("sample_description","Control_P5"))
}

for(feature in features){
  PrintSeuratGraph(namecard = "DP5",sobject=mgAVM02,graphtype = "feature",feature = feature,cellnames = c("sample_description","Deprived_P5"))
}
for(feature in features){
  PrintSeuratGraph(namecard = "CP7",sobject=mgAVM02,graphtype = "feature",feature = feature,cellnames = c("sample_description","Control_P7"))
}
for(feature in features){
  PrintSeuratGraph(namecard = "DP7",sobject=mgAVM02,graphtype = "feature",feature = feature,cellnames = c("sample_description","Deprived_P7"))
}

groups = c("age","sex","seurat_clusters")

for(group in groups){
  PrintSeuratGraph(namecard = "CP5",sobject=mgAVM02,graphtype = "dim",group = group, feature = group,cellnames = c("sample_description","Control_P5"))
}

for(group in groups){
  PrintSeuratGraph(namecard = "DP5",sobject=mgAVM02,graphtype = "dim",group = group, feature = group,cellnames = c("sample_description","Deprived_P5"))
}
for(group in groups){
  PrintSeuratGraph(namecard = "CP7",sobject=mgAVM02,graphtype = "dim",group = group, feature = group,cellnames = c("sample_description","Control_P7"))
}
for(group in groups){
  PrintSeuratGraph(namecard = "DP7",sobject=mgAVM02,graphtype = "dim",group = group, feature = group,cellnames = c("sample_description","Deprived_P7"))
}

for(feature in features){
  PrintSeuratGraph(namecard = "MGonlyvlnP5",sobject=mgAVM02,graphtype = "violin",group = "seurat_clusters",feature = feature,split = "sample_description",cellnames = c("age","P5"))
}
for(feature in features){
  PrintSeuratGraph(namecard = "MGonlyvlnP7",sobject=mgAVM02,graphtype = "violin",group = "seurat_clusters",feature = feature,split = "sample_description",cellnames = c("age","P7"))
}
for(feature in features){
  PrintSeuratGraph(namecard = "MGonlyvlnCtrl",sobject=mgAVM02,graphtype = "violin",group = "seurat_clusters",feature = feature,split = "sample_description",cellnames = c("condition","Control"))
}
for(feature in features){
  PrintSeuratGraph(namecard = "MGonlyvlnDep",sobject=mgAVM02,graphtype = "violin",group = "seurat_clusters",feature = feature,split = "sample_description",cellnames = c("condition","Deprived"))
}



r = table(mgAVM02$seurat_clusters,mgAVM02$sample_description)
Control_P5 = r[,1]*100/colSums(r)[1]
Control_P7 = r[,2]*100/colSums(r)[2]
Deprived_P5 = r[,3]*100/colSums(r)[3]
Deprived_P7 = r[,4]*100/colSums(r)[4]
relconpct = cbind(Control_P5,Deprived_P5,Control_P7,Deprived_P7)
relconpct

#notes: P5 - 4,5 P7 - 1,3 P5 Dep = 8 (at the expense of 0), P5 Ctrl - 5,4 
setEPS()
postscript("clustercomppctall.eps")
barplot(relconpct, main="Cluster composition by percent of sample",
        xlab="Cluster", ylab = "Number of Cells", ylim = c(0,100), col=c("darkblue","lightblue","red","orange","purple","lightgrey","darkgrey","magenta","yellow","black","red3","cornflowerblue"),axisnames = T,
        width = .2,xlim = c(0,2),legend = rownames(relconpct), space = 0.6,cex.names = 0.7,axis.lty = 1)
dev.off()

#Find out the relative proportion of male/female cells per cluster
ratio = table(mgAVM02$sex,mgAVM02$seurat_clusters)
Female = ratio[1,]/rowSums(ratio)[1]
Male = ratio[2,]/rowSums(ratio)[2] #normalize each row to 1000 cells per sex
ratio = rbind("Female" = Female*1000,"Male" = Male*1000)
rowSums(ratio)
#find out what percent of each cluster is made by each sex
Female = ratio[1,]*100/colSums(ratio)
Male = ratio[2,]*100/colSums(ratio)
ratio2 = rbind(Female,Male)
setEPS()
postscript("sexpct.eps")
barplot(ratio2, main="Cluster composition by sex",
        xlab="% of cluster", ylab = "Cluster", ylim = c(0,4), col=c("darkblue","lightblue"),axisnames = T,
        width = .2,xlim = c(0,100),legend = rownames(ratio2), space = 0.6,cex.names = 0.7,axis.lty = 1,horiz = T)
dev.off()

Idents(mgAVM02) = "seurat_clusters"
mgAVM02= BuildClusterTree(mgAVM02)
tree = mgAVM02@tools$BuildClusterTree
setEPS()
postscript("newtreeMGonly.eps")
plot.phylo(tree, use.edge.length = T, direction = "rightwards")
dev.off()

FeaturePlot(mgAVM02,features = "percent.mito")
FeaturePlot(mgAVM02,features = "nCount_RNA")

VlnPlot(mgAVM02,features = c("nCount_RNA","percent.mito"),group.by = "sample_description")

Idents(mgAVM02) = "age"
setEPS()
postscript("~/Desktop/mgavm02_condition_P5only.eps")
DimPlot(mgAVM02,group.by = "condition",cells = rownames(mgAVM02@meta.data[mgAVM02$age == "P5",]))
dev.off()

setEPS()
postscript("~/Desktop/mgavm02_condition_P7only.eps")
DimPlot(mgAVM02,group.by = "condition",cells = rownames(mgAVM02@meta.data[mgAVM02$age == "P7",]))
dev.off()


cellnames = mgAVM02$sample_description
cp7 = sample(cellnames[cellnames == "Control_P7"],2955,replace = FALSE)
cp5 = cellnames[cellnames == "Control_P5"]
setEPS()
postscript("~/Desktop/mgavm02_age_controlonly.eps")
DimPlot(mgAVM02,group.by = "age",cells = c(names(cp5),names(cp7)),cols = c("lightskyblue","darkblue"))
dev.off()



setEPS()
postscript("~/Desktop/mgavm02_age_deprivedonly.eps")
DimPlot(mgAVM02,group.by = "age",cells = rownames(mgAVM02@meta.data[mgAVM02$condition == "Deprived",]))
dev.off()

setEPS()
postscript("~/Desktop/mgAVM02_clusters.eps")
DimPlot(mgAVM02,group.by = "seurat_clusters",label = T)
dev.off()
#VolcanoPlot

#fc = read.csv("/Users/whippoorwill/Desktop/plots/allmarkers_vargenesMG_0215.csv",stringsAsFactors = F)
#colnames(fc)[8] = "Gene"

fc = read.csv("/Users/whippoorwill/Dropbox (Anna Molofsky Lab)/2020-Leah-barrelCortex/Manuscript data/Glia RNASeq/P5_P7 MG rnaseq/Spreadsheets/allmarkers_vargenesMG_0501.csv",stringsAsFactors = F)
highlightcluster3 = read.csv("/Users/whippoorwill/Desktop/cluster3highlight.csv",stringsAsFactors = F)
highlightcluster3 = highlightcluster3[1:34,ncol(highlightcluster3)]

highlightcluster4 = read.csv("/Users/whippoorwill/Desktop/cluster4highlight.csv",stringsAsFactors = F,header = F)
highlightcluster4 = highlightcluster4[,ncol(highlightcluster4)-1]
highlightcluster8 = read.csv("/Users/whippoorwill/Desktop/cluster8highlight.csv",stringsAsFactors = F)
highlightcluster8 = highlightcluster8[1:60,ncol(highlightcluster8)]
highlightcluster8 = highlightcluster8[highlightcluster8 != ""]


#fc = markers_all
colnames(fc)[8] = "Gene"
newlist = list()
#Split by cluster
for (i in c(0:9)){
  newlist[[i+1]] = fc[fc$cluster == i,]
}
cluster = 4

#select a single cluster
setwd("~/Desktop/plots")
fc = newlist[[cluster+1]]
#fc = fc[fc$p_val_adj<0.001,]
#fc = fc[order(fc$avg_logFC,decreasing = T),]
#write.csv(fc,file = paste0(cluster,"new_331.csv"))

fc = fc[!is.na(fc$avg_logFC),]
colorkeysdown = fc$Gene[fc$avg_logFC < -log2(1.15) & fc$p_val_adj < 10e-25]
colorkeysup = fc$Gene[fc$avg_logFC > log2(1.15) & fc$p_val_adj < 10e-25]
  

allcolors = rep("darkgrey",length(fc$Gene))
names(allcolors) = fc$Gene
a = highlightcluster3

allcolors[names(allcolors) %in% colorkeysdown] = "blue"
allcolors[names(allcolors) %in% colorkeysup]= "red"
allcolors[names(allcolors)%in% a] = "brown"
names(allcolors)[allcolors == "brown"] = "labeled"

names(allcolors)[allcolors == "red"] = "u"
names(allcolors)[allcolors == "darkgrey"] = "-"
names(allcolors)[allcolors == "blue"] = "d"



  
  
setEPS()
postscript(paste0("~/Desktop/volcano_0505",cluster,"nolabel.eps"))
EnhancedVolcano(fc,
                lab = fc$Gene,
                x = 'avg_logFC',
                y = 'p_val_adj',
                xlim = c(-3, 3),
                title = paste0("LDAVM02_",cluster),
                subtitle = "",
                drawConnectors = F,
                legendPosition = 'right',
                legendVisible = F,
                pCutoff = 10e-25,
                FCcutoff = log2(1.15),
                selectLab = a,
                transcriptPointSize = 1.5,
                transcriptLabSize = 2.0,
                col=c('black', 'black', 'black', 'red3'),
                colCustom = allcolors,
                gridlines.major = F,
                gridlines.minor = F,
                colAlpha = 1)
dev.off()
  
cluster = cluster+1
  
fc = read.csv("/Users/whippoorwill/Desktop/allmarkers_vargenesMG_0501.csv",stringsAsFactors = F)
#fc = markers_all
colnames(fc)[8] = "Gene"

newlist = list()
#Split by cluster
for (i in c(1:9)){
  newlist[[i+1]] = fc[fc$cluster == i,]
}
fc = fc[!is.na(fc$avg_logFC),]
fc = fc[fc$p_val_adj<10e-25,]
fc = fc[abs(fc$avg_logFC)> 0.2,]

library(org.Mm.eg.db)
ids=fc$Gene
fromKey="SYMBOL"
toKey=c("GENENAME")
db=org.Mm.eg.db
selRes<-AnnotationDbi::select(db,keys=ids,keytype=fromKey,columns=c(fromKey,toKey))

##this is really used to get particular keys of the keytype "ensembl" (or any other type)

x=selRes[match(ids,selRes[,1]),1:2]
identical(x$SYMBOL,fc$Gene)
fc$GeneName = x$GENENAME

write.csv(fc,file = "~/Desktop/AllmarkersCutoffe25lfc0_2.csv")

#Save umap embeddings
umap = mgAVM02@reductions$umap@cell.embeddings
write.csv(umap,file = "~/Desktop/Jupyter/mgAVM02/umap_mgAVM02_may2020.csv",row.names = T)
umap2 = read.csv("~/Desktop/Jupyter/mgAVM02/umap_mgAVM02_may2020.csv",header = T,row.names = 1,stringsAsFactors = F)

#save pcs
pc = mgAVM02@reductions$pca@cell.embeddings
write.csv(pc,file = "~/Desktop/Jupyter/mgAVM02/pca_mgAVM02_may2020.csv",row.names = T)

#Find out which clusters are specifically male/female: 4 is male, 3 is female
t = table(mgAVM02$fulldescription,mgAVM02$finalclusters)
t = t[,c(2:7,9,10)]
tP7 = t[c(3,4,7,8),]
totalP7 = rowSums(tP7)
CP7F = tP7[1,]/totalP7[1]
CP7M = tP7[2,]/totalP7[2]
DP7F = tP7[3,]/totalP7[3]
DP7M = tP7[4,]/totalP7[4]
tP7 = rbind(CP7F,CP7M,DP7F,DP7M)
tP7 = round (tP7*1000)
tP7


#LD1 contingency

all = table(LD1$SCT_snn_res.0.3,LD1$sample_description)
total = colSums(all)
clusters = rownames(all)
test = list()
chi = chisq.test(all)


for (cluster in clusters){
  x = rbind(all[cluster,],total-all[cluster,])
  rownames(x) = c(cluster,"total")
  chi = chisq.test(x)
  f = sqrt(chi$statistic/sum(x))
  test[[cluster]] = chi
  test[[cluster]]$f = f
}

#for opcs_2: 0.19
v = sqrt(test$`2`$statistic/sum(all["2",]))
#for microglia_0: 0.11
v = sqrt(test$`0`$statistic/sum(all["0",]))

#for astro_1 = 0.5
v = sqrt(test$astrocytes_1$statistic/sum(all["astrocytes_1",]))
v
#for astro_5 = 0.8
v = sqrt(test$astrocytes_5$statistic/sum(all["astrocytes_5",]))
v

#for opcs_9 = 0.3
v = sqrt(test$opcs_9$statistic/sum(all["opcs_9",]))
v

#for microglia_8 = 0.27
v = sqrt(test$microglia_8$statistic/sum(all["microglia_8",]))
v

#for microglia_6 = 1.13
v = sqrt(test$microglia_6$statistic/sum(all["microglia_6",]))
v

#for microglia_0 = 0.19
v = sqrt(test$microglia_0$statistic/sum(all["microglia_0",]))
v

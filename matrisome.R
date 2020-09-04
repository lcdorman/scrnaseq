#Find out how many of your DE genes are part of the matrisome (http://matrisomeproject.mit.edu/)

library(readr)
library(Seurat)
library(dplyr)

#You will need to add your own data here (note to Leah - put this on GEO so it can be loaded)
seurat_path = file.path("~/Desktop/Sequencing/LD_AVM02/Data/Seurat/mgAVM02_may26.RData")
load(path)

#rename whatever seurat object you have as "sobject"
sobject = mgAVM02

#set thresholds for DE genes per cluster
alpha = 1e-25
foldchange = 0.2

#if public
#de_path = "https://github.com/lcdorman/scrnaseq/blob/master/allmarkers_vargenesMG_0501.csv"
#m_file = "https://github.com/lcdorman/scrnaseq/blob/master/matrisome_mm_masterlist.csv"

#if private
m_file = "matrisome_mm_masterlist.csv"
de_file = "allmarkers_vargenesMG_0501.csv"

#Read in matrisome file
matrisome <- read.csv(m_file,stringsAsFactors = F,fill = T)
matrisome$Category = as.factor(matrisome$Category)
matrisome$Division = as.factor(matrisome$Division)

#show matrisome summary
table(matrisome$Division,matrisome$Category)

#Read in DE genes file
de = read.csv(de_file,stringsAsFactors = F,fill = T)
de$cluster = as.factor(de$cluster)

#Select only the genes that have been used for clustering in sobject
genes = rownames(GetAssayData(sobject,slot = "scale.data"))
all(de$gene %in% genes)

#how many of the dataset genes are in the matrisome, and which subset are they
m_subset = matrisome[matrisome$Gene.Symbol %in% genes,]
dim(m_subset)
table(m_subset$Division) 
table(m_subset$Category) 

#how many differentially expressed genes in each category are actually in this dataset
#select only the genes that pass a threshold
de = de[de$p_val_adj<alpha,]
de = de[de$avg_logFC > foldchange | de$avg_logFC < -foldchange,]

#make a separate data table with all up/down regulated genes
deup = de[de$avg_logFC > foldchange,]
dedown = de[de$avg_logFC < -foldchange,]

#for each cluster, calculate # of matrisome genes in up- and down-regulated gene sets
clusters = levels(de$cluster)
ecmlist = {}
for (cluster in clusters){
  gup = deup[deup$cluster == cluster,"gene"]
  gupecm = mall[mall$Gene.Symbol %in% gup,]
  gdown = dedown[dedown$cluster == cluster,"gene"]
  gdownecm = mall[mall$Gene.Symbol %in% gdown,]
  ecmlist[[cluster]][["up"]] = gupecm
  ecmlist[[cluster]][["down"]] = gdownecm
}
 
str(ecmlist,max.level = 2) 
table(deup$cluster)
table(dedown$cluster)

#Make a graph showing ecm content for each cluster, up and down
up = table(deup$cluster)
down = table(dedown$cluster)
ecmup = rep(0,length(up))
ecmdown = rep(0,length(down))


t = rbind(up,ecmup,"rup" = ecmup,down,ecmdown,"rdown" = ecmdown)

for (cluster in clusters){
  ecmup = nrow(ecmlist[[cluster]]$up)
  ecmdown = nrow(ecmlist[[cluster]]$down)
  t["ecmup",cluster] = ecmup
  t["ecmdown",cluster] = ecmdown
}

t["rup",] = round(t["ecmup",]/t["up",]*100,1)
t["rdown",] = round(t["ecmdown",]/t["down",]*100,1)
  
t = cbind(t,"average" = rowMeans(t,na.rm = T))

#actually make a bargraph of t$rup and t$rdown

#exactly what you want to plot
p = t[c("rdown","rup"),]

setEPS()
postscript("barplot_matrisome_updown.eps")
barplot(p, horiz = F, las=1, xlim = c(0,30),xlab = "cluster", ylab = '% matrisome associated',
        beside=T, col=c('blue','red'),ylim = c(0,30),legend = c("downregulated","upregulated"),axis.lty = 1)
dev.off()

#print out important gene lists (or replace "m" with clusterlevels)
imp = c()

#To select specific clusters
clusters = c('3','4','5','8','9')

#To select all clusters
clusters = levels(de$cluster)

#Make a table that could be printed (from a list)
for (cluster in clusters){
  lu = ecmlist[[cluster]][["up"]][,1:3]
  ld = ecmlist[[cluster]][["down"]][,1:3]
  if (nrow(lu) != 0){
    lu[,"direction"] = "up"
  }
  if (nrow(ld) != 0){
    ld[,"direction"] = "down"
  }
  l = rbind(lu,ld)
  l[,"cluster"] = cluster
  imp = rbind(imp,l)
}

write.csv(imp,file = "ecmlist.csv")




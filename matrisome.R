#matrisome

sobject = mgAVM02
alpha = 1e-25
foldchange = 0.2

matrisome = read.csv("~/Desktop/matrisome_mm_masterlist.csv",fill = T,stringsAsFactors = F)
matrisome$Category = as.factor(matrisome$Category)
matrisome$Division = as.factor(matrisome$Division)

table(matrisome$Division,matrisome$Category)


genes = VariableFeatures(sobject)
de = read.csv("/Users/whippoorwill/Dropbox (Anna Molofsky Lab)/2020-Leah-barrelCortex/Manuscript data/Glia RNASeq/P5_P7 MG rnaseq/Spreadsheets/allmarkers_vargenesMG_0501.csv",stringsAsFactors = F)
all(de$gene %in% genes)
de$cluster = as.factor(de$cluster)

#how many of the variable genes are in the matrisome dataset
mall = matrisome[matrisome$Gene.Symbol %in% genes,]
dim(mall)
table(mall$Division) #123 core, 261 associated
table(mall$Category) #27 collagen, 82 glycoproteins, 74 regulators, 64 ECM-affiliated
                      #14 proteoglycans, 123 secreted factors out of 6000

#how many differentially expressed genes in each category are actually in this dataset
#select only the genes that pass a threshold

de = de[de$p_val_adj<alpha,]
de = de[de$avg_logFC > foldchange | de$avg_logFC < -foldchange,]

#make a separate sheet with all up/down regulated genes
deup = de[de$avg_logFC > foldchange,]
dedown = de[de$avg_logFC < -foldchange,]

#for each cluster, calculate # of ecm genes in up- and down-regulated gene sets
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
#cluster 1: 0,18 out of 454, 160 = 0%, 11.25%
#cluster 2: 0,5 out of 261, 87 = 0%, 1.92%
#cluster 3: 11,1 out of 315, 6 = 3.49%, 17%
#cluster 4: 24,3 out of 230, 420 = 10.4%, 0.7%
#cluster 5(0): 1,3 out of 10, 265 = 10%, 1.13%
#cluster 6: 0,10 out of 131, 88 = 0%, 11.4%
#cluster 8: 5, 0 out of 111, 3 = 4.5%, 0%
#cluster 9: 40, 6 out of 298, 60 = 13.4%, 10%

#allgenes: 384 out of 6000 = 6.4%

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
p = t[c("rup","rdown"),]
setEPS()
postscript("~/Desktop/barplot_matrisome_updown.eps")
barplot(p, horiz = F, las=1, xlim = c(0,30),xlab = "cluster", ylab = '% matrisome associated',
        beside=T, col=c('red','blue'),ylim = c(0,30),legend = c("upregulated","downregulated"),axis.lty = 1)
dev.off()

#print out important gene lists
imp = c()
m = c('3','4','5','8','9')
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

write.csv(imp,file = "~/Desktop/ecmlist.csv")




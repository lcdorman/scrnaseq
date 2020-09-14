#First you will need to ensure all the following programs are installed. Use "BiocLite()" to install new programs. 
source("https://bioconductor.org/biocLite.R")
library(ensembldb)
library(EnsDb.Mmusculus.v79)
library(AnnotationFilter)
library(tximport)
library(readr)
library(DESeq2) #The following code uses the DESEq2 vignette as a base
library(rhdf5)

#Set "dir" to the folder in which your samples can be found. 
dir<-"~/Desktop/RNAseq practice"
contrast=c("condition","1_Pos","2_Neg")
alpha = 0.05
minlfc = 1
name = "practice"

##tximport
##"Samples.txt" must
##be a table with (at least) the sample name and the name of the file containing 
##the corresponding quantification folders. If necessary, create "samples.txt" and write it (below, hashed). 
#samplenames<-c("Control","Deprived","Extra")
#samplefiles<-c("a","b","c")
#condition<-c("treated","control","treated")
#animal<-c("one","two","three")
##add any columns you want; you can use them to design your experiment for DESeq
#samples<-data.frame(cbind(samplenames,samplefiles,condition,animal))

#write.table(samples,file.path(dir,"samples.txt", sep="\t",col.names=T,row.names=F)

##If the file already exists, read it in: 
samples<-read.table(file.path(dir,"samples.txt"),header=TRUE,stringsAsFactors = FALSE)

##creating the tx2gene matrix from ensembleIDs to convert tx_id given by your data to gene id's that you can read
txdf<-transcripts(EnsDb.Mmusculus.v79, return.type="DataFrame")
tx2gene<-as.data.frame(txdf[,c("tx_id","gene_id")])
write.table(tx2gene,"tx2gene.txt",sep="\t", quote=F,col.names=T,row.names=F)

##create a named vector pointing to the quantification files. 
files<-file.path(dir,samples$samplefiles,"abundance.h5") #This opens nested folders dir, then "a.txt", then "abundance.h5". check your folders for naming
names(files)<-samples$samplenames
all(file.exists(files)) #check that this is true
write.table(files,file.path(dir,"files.txt"),sep="\t") #save files so you don't have to recreate the above code

#import kallisto alignment files. If you used STAR or another alignment tool, change "type"
txi<-tximport(files,type='kallisto',ignoreTxVersion=TRUE,tx2gene=tx2gene)
write.table(txi$counts,file.path(dir,"counts.txt"),sep="\t")
write.table(txi,file.path(dir,"txi.txt"),sep="\t")

##change zero length entries to 1 so they don't get dropped
txi$length[txi$length==0]<-1
##http://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
##returns a list with 4 objects: abundance, counts, length, countsfromabundance


#create sampleTable as a framework for DESeq2
samples
##set rows of samples to column names from txicounts (same as file vector) or change to samples$samplenames if you'd rather use descriptions
rownames(samples)<-samples$samplefiles
colnames(txi$counts)<-rownames(samples)

##create object that can be analyzed by DESeq2. The last element in "design" is the experimental variable. See vignette for more complicated designs
dds<-DESeqDataSetFromTximport(txi,colData=samples,design=~animal+condition)
counts<-counts(dds)

##pre-filter for counts higher than 1
dds<-dds[rowSums(counts(dds))>1,]
##set factor levels - the first level is the reference (or use code below for "resDE")
#dds$condition<-factor(dds$condition,levels=c("2_Neg","1_Pos"))

#Differential Expression
ddsDE<-DESeq(dds)

resDE<-results(ddsDE,cooksCutoff=FALSE,contrast=contrast)
##resDE<-results(ddsDE,cooksCutoff=FALSE,addMLE=TRUE,contrast=c("condition","treated","control"))
##cooksCutoff= T means you are removing p-values from extreme count outliers (default)

##to order the results by adjusted p-value
resOrdered<-resDE[order(resDE$padj),]
##number of results that are significant

sig_sum<-sum(resDE$padj<alpha,na.rm=TRUE)
#remove all "na" values and replace with 1
remove_na<-function(dataset){
        for (i in 1:length(dataset)){
                if (is.na(dataset[i])) dataset[i]=1
        }
        dataset
}
resDE$padj<-remove_na(res2DE$padj)

sig<-subset(resDE,log2FoldChange>minlfc)
sig<-subset(sig,padj<alpha)

##to add a column with gene names to the final spreadsheet
library(org.Mm.eg.db)
ids=rownames(sig)
fromKey="ENSEMBL"
toKey="SYMBOL"
db=org.Mm.eg.db
selRes<-AnnotationDbi::select(db,keys=ids,keytype=fromKey,columns=c(fromKey,toKey))

##this is really used to get particular keys of the keytype "ensembl" (or any other type)

x=selRes[match(ids,selRes[,1]),2]
##match finds positions of the first match between the rownames
##(ensemblIDs from DESeq table) and the first column of selRes, which is ensembl
##ID's found in the gene names table from org.Mm.eg.db. 
##so x returns selRes with only the first instance of a given 
##ensemblID and each symbol to go with it. 

DE=data.frame(Ensembl_ID=rownames(sig),Gene_ID=x,sig)
allDE=data.frame(Ensembl_ID=rownames(resDE),Gene_ID=x,resDE)
##write a table
write.csv(DE,file.path(dir,paste0(name,"_DE.csv")),row.names=T,col.names=T)
write.csv(allDE,file.path(dir,paste0(name,"_allDE.csv")),row.names=T,col.names=T)
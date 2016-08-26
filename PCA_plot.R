# PCA plot

# Read data
rm(list=ls(all=TRUE))
setwd("/Data/Computation/microArray/01PCA")

library(lattice)
library(ggplot2)
library(gplots)

foxc1_p47KO_fc=read.table("Foxc1_fc_P47KO_by_gene.txt",header=T)
names(foxc1_p47KO_fc)=c("gene","Foxc1_P47KO")

lin_greco_fc=read.table("Lineage_Greco_HG_Bulge_fc_KO-vs-WT_gene.txt",header=T)
names(lin_greco_fc)[1]="gene"

lin_lien_fc=read.table("Lineage_Lien_aHFSC+qHFSC+TAC_fc_KO-vs-WT_gene.txt",header=T)

lin_zhang=read.table("Lineage_Zhang_P18+P25_bulge_division_fc_KO-vs-WT_gene.txt",header=T)
lin_zhang_fc=data.frame(gene=lin_zhang$gene,
                        p18_bu_div0=apply(lin_zhang[,2:3],1,mean,na.rm=T),
                        p18_bu_div2=apply(lin_zhang[,4:6],1,mean,na.rm=T),
                        p25_bu_div0=apply(lin_zhang[,7:9],1,mean,na.rm=T),
                        p25_bu_div2=apply(lin_zhang[,10:12],1,mean,na.rm=T))


#### Merge all fold change table

files=ls()
fc_files=files[grep("_fc",files,perl=T)]
fc_sublist=fc_files

rm(merged)
merged=lin_lien_fc
for( i in 1:length(fc_sublist)){
  var=as.symbol(fc_sublist[i])
  
  merged=merge(merged,eval(var),all.x=T,all.y=T)
}


### Visualize data 

head(merged)
flatten = stack(merged[,2:ncol(merged)])
celltype=c("Intermediate","Intermediate","Quiescent","Differentiated",
           "Differentiated","Unknown","Activated","Quiescent","Activated",
           "Quiescent","Activated","Quiescent","Intermediate",
           "Quiescent","Quiescent")
names(celltype)= names(merged)[-1]
flatten$celltype=apply(flatten,1,function(x) celltype[x[2]])

### Clean data
#### Exclude genes that contains missing values
merged2=na.omit(merged)


#### Select genes that have more than 2-fold change in at least one of the observations
rowmax=apply(merged2[,2:ncol(merged2)],1,function(x) max(abs(x),na.rm=T))
merged_rowmax=cbind(merged2,rowmax)
union2fold=subset(merged_rowmax,rowmax>1)
row.names(union2fold)=union2fold$gene
union2fold$rowmax=NULL
union2fold=na.omit(union2fold)
gene=union2fold$gene
union2fold$gene=NULL


#### Centering, which returns a transposed matrix
colmean=apply(union2fold,2,mean)
union2fold_std=apply(union2fold, 1, function(x) x-colmean) 

### Principle Component Analysis
#### Compute principle components. Skipping the scaling step gives a better result
res.pca=prcomp(t(union2fold)) 
scores=as.data.frame(res.pca$x)
scores$celltype=as.factor(celltype[rownames(scores)])
scores$celltype=factor(scores$celltype,
                       levels=c("Quiescent","Intermediate","Activated","Differentiated","Unknown" ))

#### Variance explained
par1=round(res.pca$sdev[1]^2/sum(res.pca$sdev^2)*100,2)
par2=round(res.pca$sdev[2]^2/sum(res.pca$sdev^2)*100,2)

####  label colors
color1=c("orange","orange","deepskyblue4","darkgreen","darkgreen","black",
         rep(c("red","deepskyblue4"),3),"orange","deepskyblue4","deepskyblue4")


####  plot PCA in 2D
ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores))) +
  geom_point(aes(colour=celltype),size=4)+
  geom_text(alpha = 1, size = 6,aes(colour=factor(celltype)),show_guide=F,vjust=1) +
  ggtitle("PCA plot HFSC transcriptome")+
  xlab(paste("PC1 ",par1,"%"))+ylab(paste("PC2 ",par2,"%"))+
  theme(legend.title=element_text(size=16),
        legend.text=element_text( size = 12, face = "bold"))

###  Kmeans Clustering
set.seed(120410)
kmeans_out=kmeans(scores[,1:2],4,iter.max = 100)

Cluster=paste("Cluster",kmeans_out$cluster,sep="")
ggplot(data=scores, aes(x=PC1,y=PC2, label=rownames(scores))) +
  geom_text(data=scores, alpha = 1, size = 6,aes(colour=Cluster),show_guide=F) +
  ggtitle("PCA plot HFSC transcriptome")+
  xlab(paste("PC1 ",par1,"%"))+ylab(paste("PC2 ",par2,"%"))


###  Validation 2: Linear Discriminative Analysis (LDA)
library(MASS)
x=scores[rownames(scores)!="Foxc1_P47KO",]
x$celltype=factor(x$celltype)
lda1=lda(celltype~.,data=x[,c(1:14,16)])
y_pred=predict(lda1,x[,c(1:14,16)])
table(y_pred$class,x[,16])
predict(lda1,scores["Foxc1_P47KO",])
#### Prediction is consistent with hypothesis
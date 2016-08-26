print("fold_change_by_gene.R myScript.R [input directory]")
rm(list=ls(all=TRUE))

args <- commandArgs(trailingOnly = TRUE)

setwd(args[1])
#Read data file (with .CEL extention)
library(affy)
library(mouse4302.db)
library(annotate)
library(miscTools)

data=ReadAffy()

processed = expresso(data, bgcorrect.method="rma",
                     normalize.method="quantiles",
                     pmcorrect.method="pmonly",
                     summary.method="medianpolish")
expressions=exprs(processed)

# DO P/M/A calls
calls<-mas5calls(data)
PMcalls<-exprs(calls)
colnames(PMcalls)<-paste(colnames(PMcalls),"PM",sep=".")

# match probes to gene symbols
probe2gene=read.table("/Data/Computation/microArray/00pipeline/Mouse430_2.na33.annot.geneSymbol_multi.txt",header=T)
probe2gene=subset(probe2gene,GeneSymbol!="---")
Probe_ID<-rownames(expressions)
expression_PM=cbind(Probe_ID=Probe_ID,data.frame(expressions),data.frame(PMcalls))

merged=merge(probe2gene,expression_PM)


#output data
output=merged[order(merged[,2],merged[,1]),]

# Read sample table
sample_table=read.table("sample_table.txt",header=T)

# Calculate fold change with P/M/A filter. 
# Compute fold change only when both value in the pair are "P".

pairs=levels(sample_table[,"pair"])

fc_table=data.frame(matrix(ncol = length(pairs), nrow = nrow(output)))
names(fc_table)=pairs
for(i in pairs){
  sub_table=subset(sample_table,pair==i)
  sub_table=sub_table[match(c("KO","WT"),sub_table$group),]
  target=c(paste(sub_table[,1],"CEL",sep="."),paste(sub_table[,1],"CEL","PM",sep="."))
  sub_output=output[,target]
  
  # Compute fold change, the raw value is log 2 transformed
  
  #fc=ifelse(sub_output[,3]=="P" & sub_output[,4]=="P",sub_output[,1]-sub_output[,2],NA)
  fc=sub_output[,1]-sub_output[,2]
  fc_table[,i]=fc
  
}
#####################
# Input pairs_name
#####################
#pairs_name=c("HG_set1","Bulge_set1","HG_set2","Bulge_set2")

if (exists("pairs_name")){
  names(fc_table)=pairs_name
}

fc_table=cbind(Probe_ID=output$Probe_ID, GeneSymbol=output$GeneSymbol,fc_table)
fc_table=fc_table[order(fc_table[,2],fc_table[,1]),]

#write.table(fc_table,file="fc_KO-vs-WT.txt",sep="\t",col.names=TRUE,quote=F,row.names=F)

# Apply median to fold changes of multi-probe genes
gene_list=unique(fc_table$GeneSymbol)

fc_by_gene=data.frame(matrix(ncol = length(pairs), nrow = length(gene_list)))
if (exists("pairs_name")){
  names(fc_by_gene)=pairs_name
}

#group by gene
split=split(fc_table,fc_table$GeneSymbol)
fc_table$GeneSymbol=factor(fc_table$GeneSymbol)
fc_table$Probe_ID=factor(fc_table$Probe_ID)
bygene=sapply(split(fc_table,fc_table$GeneSymbol),function(x) colMedians(x[,3:ncol(fc_table)]))
fc_by_gene=t(bygene)


write.table(fc_by_gene,file="fc_KO-vs-WT_gene071915.txt",sep="\t",col.names=TRUE,quote=F)


library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggbeeswarm)

## Panel=7A

#-------------------------------
#- Select dataset (bulk or NART)
#-------------------------------

#- uncomment for Bulk
#dataset<-"Bulk"
#dataset2<-"Trm-dys"
#file<-"JR_Bulk1TumorTRTvsnTRT.txt"

#- uncomment for nart
dataset<-"NART"
dataset2<-"Neo-CD8"
file<-"L011_Unfiltered_TumNARTneg_vs_TumNARTpos.txt"

#----------------------------
#- Read list of genes to plot
#----------------------------
genes_bulkfinal<-read.table("080419Volcano_annotation_bulk.txt")
genes_bulkfinal<-genes_bulkfinal[,1]
genes_nartfinal<-read.table("080419Volcano_annotation_nart.txt")
genes_nartfinal<-genes_nartfinal[,1]


#----------------------------
#- Read data
#----------------------------


fcthreshold=1
pointsizefactor=1

ha<-read.table(file, sep="\t",header=T)
if(dataset == "NART"){
	ha$logFC <- -ha$logFC	
	names(ha)[3] <- "HGNC" 			#- change from "ext_gene"
	names(ha)[8] <- "log2FoldChange" 	#- change from "logFC"
	names(ha)[10] <- "pvalue"	 	#- change from "PValue"
	names(ha)[11] <- "padj" 		#- change from "FDR"
	fcthreshold=1
	pointsizefactor=5
}
ha <- ha %>% filter(!is.na(pvalue)) 


#-------------------------------------------------------------
#- Add labels/colours for genes according to their FC/P values
#-------------------------------------------------------------
ha = within(ha, {Col="Other"})
ha[abs(ha$log2FoldChange) > fcthreshold, "Col"] <- paste0("|LogFC|>",fcthreshold)
ha[ha$padj <.05, "Col"] <- "FDR<0.05"
ha[ha$padj <.05 & abs(ha$log2FoldChange)> fcthreshold, "Col"] <- paste0("FDR<0.05 & |LogFC|>",fcthreshold)

ha$Col<-factor(ha$Col, 
		levels=c("FDR<0.05",
			paste0("|LogFC|>",fcthreshold),
			paste0("FDR<0.05 & |LogFC|>",fcthreshold),
			"Other"), 
		labels=c("FDR<0.05",
			paste0("|LogFC|>",fcthreshold),
			paste0("FDR<0.05 & |LogFC|>",fcthreshold),
			"Other"))

#------------
#- Plot data
#------------

cols <- c(a = "indianred1", 
	b = "gold2",
	c = "cornflowerblue",
	d="gray47")
names(cols)[1]<-"FDR<0.05"
names(cols)[2]<-paste0("|LogFC|>",fcthreshold)
names(cols)[3]<-paste0("FDR<0.05 & |LogFC|>",fcthreshold)
names(cols)[4]<-"Other"

p<-ggplot(ha, aes(log2FoldChange, -log10(pvalue))) +
	geom_point(aes(col=Col), size=abs(ha$log2FoldChange)/pointsizefactor, alpha=0.5) + 
	ylim(0,10) +
	ggtitle(dataset2) + 
	scale_alpha(guide = "none") +
	xlab(bquote(~Log[2]~ "Fold Change")) + ylab(bquote(~-Log[10]~italic(P))) +
	geom_vline(xintercept = fcthreshold, linetype = 2, alpha = 0.5) +
	geom_vline(xintercept = -fcthreshold, linetype = 2, alpha = 0.5) +
	scale_color_manual(values=cols) + 
	theme_classic() +
	theme(legend.title=element_blank(), legend.position="bottom",legend.text=element_text(size=18),text = element_text(size=20),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) 

#-----------------
#- Add gene labels
#-----------------


if(dataset == "Bulk"){
	ha2 <- ha %>%
	   filter(HGNC %in% genes_bulkfinal) %>% 
	  select(HGNC,log2FoldChange,pvalue,padj,Col)
}else if(dataset =="NART"){
	ha2 <- ha %>%
	   filter(HGNC %in% genes_nartfinal) %>% 
	  select(HGNC,log2FoldChange,pvalue,padj,Col)
}


p<-p+geom_text_repel(data=ha2,aes(label= HGNC), box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"), force=1, size=5, segment.size=0.3,segment.alpha=0.5)
pdf(paste0("volcano_", dataset, ".pdf"), width = 10, height=12)
p
dev.off()


## Panel=7B
comb_dat<-read.table("./combined_dat.txt",header=T,sep="\t",stringsAsFactors=F)
ggplot(comb_dat, aes(x=Log2_FC_CPI,y=Log2_FC_NART)) +geom_point(size=4,stroke=1,shape=23,fill="darkblue",color="grey") + geom_text_repel(aes(x=Log2_FC_CPI,y=Log2_FC_NART, label = gene), size = 3)

## Panel=7C
merge_dat<-read.table("./cxcl13_ccr5_tpm.txt",sep="\t",header=T,stringsAsFactors=F)
cbPalette2 <- c("darkblue","darkred")
ggplot(data=merge_dat,aes(x=response, y=CXCL13,color=response))+geom_violin() +geom_quasirandom(size = 3, shape=18)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+scale_fill_manual(values=cbPalette2)+scale_color_manual(values = cbPalette2)+ theme(plot.margin = unit(c(1,1,1,1), "lines"))+xlab("")+scale_y_log10()
ggplot(data=merge_dat,aes(x=response, y=CCR5,color=response))+geom_violin() +geom_quasirandom(size = 3, shape=18)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+scale_fill_manual(values=cbPalette2)+scale_color_manual(values = cbPalette2)+ theme(plot.margin = unit(c(1,1,1,1), "lines"))+xlab("")+scale_y_log10()





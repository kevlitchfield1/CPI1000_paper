library(pheatmap)
library(CombinePValue)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(nlme)
library(plyr)
require(reshape2)
library(gridExtra)
library(cowplot)
library(survival)
library('survminer')
library(scales)
library('dplyr')
library(ggbeeswarm)
library(meta)
library(ggridges)




rm(list=ls(all=TRUE))

## Define functions

run_meta_analysis_pancan <- function(metrics_output,studies,cols_for_test,num) {

	# Create output matrices for logistic regression results
	ORs_out<-data.frame(matrix(0, ncol = length(studies), nrow = length(cols_for_test)))
	Pvals_out<-data.frame(matrix(1, ncol = length(studies), nrow = length(cols_for_test)))
	Estimates_out<-data.frame(matrix(0, ncol = length(studies), nrow = length(cols_for_test)))
	SEs_out<-data.frame(matrix(0, ncol = length(studies), nrow = length(cols_for_test)))
	names(Pvals_out)<-studies;names(Estimates_out)<-studies;names(SEs_out)<-studies
	row_names<-colnames(metrics_output)[cols_for_test]
	rownames(Estimates_out)<-row_names;rownames(Pvals_out)<-row_names;rownames(SEs_out)<-row_names

	# Run logistic regression analysis for each biomarker one by one
	for(i in 1:length(cols_for_test)){
		this_tab<-metrics_output[,c(1,cols_for_test[i],ncol(metrics_output),ncol(metrics_output)-1,ncol(metrics_output)-2)]
		names(this_tab)[2]<-"metric"
		for(j in 1:length(studies)){
			this_hist<-this_tab[this_tab$histology==as.character(data.frame(strsplit(studies[j],"-"))[1,1]),]
			this_study<-this_hist[this_hist$study==as.character(data.frame(strsplit(studies[j],"-"))[2,1]),]
			this_study<-this_study[complete.cases(this_study),]

			# Only analyse if n=>10
			if(nrow(this_study)>=10 && sum(this_study$metric)){
				this_study<-this_study[!(is.na(this_study$metric)),]
				this_study$metric_test<-scale(this_study$metric)
				# Use logistic regression
				this_study$Response_binary<-1
				this_study$Response_binary[this_study$response=="no_response"]<-0
				b<-summary(glm(Response_binary~metric_test,data=this_study,family=binomial))
				Estimates_out[(i),(j)]<-b$coefficients[2,1]
				SEs_out[(i),(j)]<-b$coefficients[2,2]
				if(class(this_study$metric)=="logical"){p_FET<-fisher.test(table(this_study$metric,this_study$response));Pvals_out[(i),(j)]<-p_FET$p.value}
				if(class(this_study$metric)!="logical"){p_wilcox<-wilcox.test(metric~response,data=this_study);Pvals_out[(i),(j)]<-p_wilcox$p.value}
						}	  
					}
				}

	# Execute meta analysis and get ORs/p's for forest plot 
	meta_set<-data.frame(matrix(1, ncol = 4, nrow = nrow(Estimates_out)))
	rownames(meta_set)<-row_names
	colnames(meta_set)<-c("OR","lower_OR","upper_OR","Meta_Pval")
	for(k in 1:nrow(Estimates_out)){
		meta_dat<-data.frame(t(Estimates_out[(k),]),t(SEs_out[(k),]),studies)
		names(meta_dat)[1]<-"Effect"
		names(meta_dat)[2]<-"SE"
		a<-metagen(Effect,SE, studlab = studies,sm = "OR",data = meta_dat)
		meta_set[k,1]<-exp(a$TE.random)
		meta_set[k,2]<-exp(a$lower.random)
		meta_set[k,3]<-exp(a$upper.random)
		meta_set[k,4]<-a$pval.random
		}
	#Estimates_out_set<-Estimates_out
	#SEs_out_set<-SEs_out
	Estimates_ORs_out_set<-exp(Estimates_out)
	assign(paste0("meta_set",num),meta_set,envir = .GlobalEnv)
	assign(paste0("Estimates_out_set",num),Estimates_out,envir = .GlobalEnv)
	assign(paste0("SEs_out_set",num),SEs_out,envir = .GlobalEnv)
	assign(paste0("Estimates_ORs_out_set",num),Estimates_ORs_out_set,envir = .GlobalEnv)
}


## Figure=4A

# read in signature data per sample:
mutation_sigs<-read.table("./SIGs_metrics_final.txt",header=T,sep="\t",stringsAsFactors=F)
mutation_sigs$signature_apobec<-mutation_sigs$Signature.2+mutation_sigs$Signature.13;mutation_sigs$Signature.13<-NULL;mutation_sigs$Signature.2<-NULL;mutation_sigs$Signature.14<-NULL

# define studies
studies<-c("MELANOMA-SNYDER_NEJM_2014","MELANOMA-VANALLEN_SCIENCE_2015","MELANOMA-HUGO_CELL_2016","MELANOMA-RIAZ_CELL_2017","MELANOMA-CRISTESCU_SCIENCE_2018","LUNG-RIZVI_SCIENCE_2015","LUNG-HELLMAN_ALL_2018","BLADDER-SNYDER_PLOSMED_2017","BLADDER-MARIATHASAN_NATURE_2018","BLADDER-CRISTESCU_SCIENCE_2018","RENAL-MCDERMOT_NMED_2018","CRC-Diaz_MSI_CRC","HEAD AND NECK-CRISTESCU_SCIENCE_2018","BREAST-CRISTESCU_SCIENCE_2018","OTHER-CRISTESCU_SCIENCE_2018")
# Define signature columns for signature analysis
cols_for_test<-c(1:20)
metrics_output<-mutation_sigs[,c(1:19,31,25,26,29,28,30)]
names(metrics_output)[23]<-"histology";names(metrics_output)[24]<-"study";names(metrics_output)[25]<-"response"

# run meta-analysis for mutation signatures
run_meta_analysis_pancan(metrics_output,studies,cols_for_test,1)

# forest plot
meta_results<-meta_set1[1:20,]
meta_results$biomarker<-rownames(meta_results)
meta_results$upper_OR[meta_results$upper_OR>3] <- 3
meta_results$biomarker<-factor(meta_results$biomarker,levels=meta_results$biomarker)
fp <- ggplot(data=meta_results,aes(x=biomarker, y=OR, ymin=lower_OR, ymax=upper_OR))+geom_pointrange(colour="black",size=0.6,shape=22,fill="black") + geom_hline(yintercept=1, lty=2) +xlab("") + ylab("Standardized Odds Ratio (95% CI)")+ theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, hjust = 1))+scale_x_discrete(limits = levels(meta_results$biomarker));fp

# Adjust for TMB
metrics_output$Response_binary<-1
metrics_output$Response_binary[metrics_output$response=="no_response"]<-0
summary(glm(Response_binary~Signature.1A+TMB+histology,data=metrics_output,family=binomial))
summary(glm(Response_binary~Signature.4+TMB+histology,data=metrics_output,family=binomial))
summary(glm(Response_binary~Signature.7+TMB+histology,data=metrics_output,family=binomial))
summary(glm(Response_binary~signature_apobec+TMB+histology,data=metrics_output,family=binomial))


## Figure=4B
ggplot(data=metrics_output,aes(x=histology, y=Signature.7)) + geom_boxplot(outlier.shape=NA,aes(fill = histology))+geom_jitter(shape=16,size=1,position=position_jitter(0.1),color = "grey",fill = "black")+scale_fill_brewer(palette="Blues")+ theme(legend.position = "none",axis.text.x = element_text(angle = 30, hjust = 1,size=8),plot.title = element_text(size = 8),axis.text.y = element_text(size=8),axis.title=element_text(size=8,face="bold"))
ggplot(data=metrics_output,aes(x=histology, y=DNVs)) + geom_boxplot(outlier.shape=NA,aes(fill = histology))+geom_jitter(shape=16,size=1,position=position_jitter(0.1),color = "grey",fill = "black")+scale_fill_brewer(palette="Blues")+ theme(legend.position = "none",axis.text.x = element_text(angle = 30, hjust = 1,size=8),plot.title = element_text(size = 8),axis.text.y = element_text(size=8),axis.title=element_text(size=8,face="bold"))+ scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
ggplot(metrics_output, aes(x=Signature.7,y=DNVs))+geom_point(size=1,stroke=0.3,shape=23,fill="darkblue",color="grey")+theme(axis.text.x=element_text(angle=20, hjust=1,size=8),plot.title = element_text(size = 8),axis.text.y = element_text(size=8),axis.title=element_text(size=8,face="bold"),strip.text = element_text(size = 8))+stat_cor(method="spearman",size=3,label.x.npc=0.4,label.y.npc=0.05)+ scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))+geom_smooth()


## Figure=4C
snv_aa_dat<-read.table("./AA_hm_data_SNVs_data.txt",header=T,sep="\t",stringsAsFactors=F)
paletteLength <- 400;myColor <- colorRampPalette(c("white", "darkred"))(paletteLength)
row.names(snv_aa_dat)<-snv_aa_dat$rowname;snv_aa_dat$rowname<-NULL
pheatmap(log(snv_aa_dat+1),color=myColor, border_color = "grey",cluster_cols = F, cluster_rows = F,fontsize_row=6,fontsize_col=6,)

dnv_aa_dat<-read.table("./AA_hm_data_DNVs_data.txt",header=T,sep="\t",stringsAsFactors=F)
paletteLength <- 400;myColor <- colorRampPalette(c("white", "darkred"))(paletteLength)
row.names(dnv_aa_dat)<-dnv_aa_dat$rowname;dnv_aa_dat$rowname<-NULL
pheatmap(log(dnv_aa_dat+1),color=myColor, border_color = "grey",cluster_cols = F, cluster_rows = F,fontsize_row=6,fontsize_col=6,)

barplot1_dat<-read.table("./bar_plot_dat_number_AAs.txt",header=T,sep="\t",stringsAsFactors=F)
ggplot(barplot1_dat, aes(x=Type,y=Number,fill=Type,colour=Type)) +geom_bar(stat="identity",position="dodge")+scale_fill_manual(values = c("#E7B800", "#FC4E07"))+scale_color_manual(values = c("#E7B800", "#FC4E07"))
barplot2_dat<-read.table("./bar_plot_dat_radical_change.txt",header=T,sep="\t",stringsAsFactors=F)
ggplot(barplot2_dat, aes(x=Type,y=Mean,fill=Type,colour=Type)) +geom_bar(stat="identity",position="dodge")+scale_fill_manual(values = c("#E7B800", "#FC4E07"))+scale_color_manual(values = c("#E7B800", "#FC4E07"))+geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean+Sd), width=.2,position=position_dodge(.9),color="black")


## Figure=4D
dat<-read.table("./hydrophobicity_score.txt",header=T,sep="\t",stringsAsFactors=F)
ggplot(dat, aes(y=Gran_dist,x=Type,fill=Type))+geom_boxplot(outlier.shape=NA)+scale_fill_manual(values = c("#FC4E07","#E7B800"))+scale_color_manual(values = c("#FC4E07","#E7B800"))
ggplot(dat, aes(x=hydrophobicity.change,y=Type,fill=Type))+geom_density_ridges(bandwidth=1)+scale_fill_manual(values = c("#FC4E07","#E7B800"))+scale_color_manual(values = c("#FC4E07","#E7B800"))


## Figure=4E
dat<-read.table("./immunogenic_hydro_scores.txt",header=T,sep="\t",stringsAsFactors=F)
ggplot(dat, aes(y=hydro_score,x=cd8_epitope,fill=cd8_epitope))+geom_boxplot(outlier.shape=NA)+stat_compare_means()+scale_fill_manual(values = c("#E7B800", "#FC4E07"))+scale_color_manual(values = c("#E7B800", "#FC4E07"))


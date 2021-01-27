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
library("DescTools")
library("dplyr")

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
	# Error message "Zero values in seTE replaced by NAs" is expected, this just allows studies/biomarkers with SE exactly zero to be replaced with NA 
	# (which is what we want) as the SE=0 studies/biomarker combinations actually did not have observations n=>10 to allow analysis, and so hence meta package
	# replaces with NA and ignores these studies
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


run_meta_analysis_histology <- function(meta_results,Estimates_out_set1,Estimates_out_set2,studies,this_subgroup) {

	meta_results_by_subgroup<-meta_results[1,]
	meta_results_by_subgroup$Group<-"Group"
	meta_results_by_subgroup<-meta_results_by_subgroup[0,]

	for(k in 1:nrow(Estimates_out_set1)){
		meta_dat<-data.frame(t(Estimates_out_set1[(k),studies]),t(SEs_out_set1[(k),studies]),studies)
		names(meta_dat)[1]<-"Effect";names(meta_dat)[2]<-"SE"
		a<-metagen(Effect,SE, studlab = studies,sm = "OR",data = meta_dat)
		this_result<-data.frame(exp(a$TE.fixed),exp(a$lower.fixed),exp(a$upper.fixed),a$pval.fixed,row.names(Estimates_out_set1)[k],this_subgroup)
		names(this_result)<-names(meta_results_by_subgroup)
		tmp<-meta_results_by_subgroup
		meta_results_by_subgroup<-rbind(tmp,this_result)
		}
	for(k in 1:nrow(Estimates_out_set2)){
		meta_dat<-data.frame(t(Estimates_out_set2[(k),studies]),t(SEs_out_set2[(k),studies]),studies)
		names(meta_dat)[1]<-"Effect";names(meta_dat)[2]<-"SE"
		a<-metagen(Effect,SE, studlab = studies,sm = "OR",data = meta_dat)
		this_result<-data.frame(exp(a$TE.fixed),exp(a$lower.fixed),exp(a$upper.fixed),a$pval.fixed,row.names(Estimates_out_set2)[k],this_subgroup)
		names(this_result)<-names(meta_results_by_subgroup)
		tmp<-meta_results_by_subgroup
		meta_results_by_subgroup<-rbind(tmp,this_result)
		}
	assign(paste0("meta_results_by_subgroup_",this_subgroup),meta_results_by_subgroup,envir = .GlobalEnv)
}

calculate_variance <- function(studies,column) {

	this_tab<-metrics_output
	this_tab$hist_study<-paste0(this_tab$histology,"-",this_tab$study)
	this_study<-this_tab[this_tab$hist_study %in% studies,]
	this_study$Response_binary<-1
	this_study$Response_binary[this_study$response=="no_response"]<-0
	#sources of antigen
	this_input<-this_study[,c(2:14,ncol(this_study))]
	this_mod1<-glm(Response_binary~.,data=this_input,family=binomial) 
	variances_out[1,column]<-PseudoR2(this_mod1)
	#sources of antigen+immune evasasion
	this_input<-this_study[,c(2:24,ncol(this_study))]
	this_mod2<-glm(Response_binary~.,data=this_input,family=binomial) 
	variances_out[2,column]<-(PseudoR2(this_mod2)-PseudoR2(this_mod1))
	#sources of antigen+immune evasasion+host
	this_input<-this_study[,c(2:31,ncol(this_study))]
	this_mod3<-glm(Response_binary~.,data=this_input,family=binomial) 
	variances_out[3,column]<-(PseudoR2(this_mod3)-PseudoR2(this_mod2)) 
	#sources of antigen+immune evasasion+host+immune infiltration (NB need to z-score normalise transcriptome datasets within each study, before combining in the R2 analysis, as some are nanostring, some RNAseq)
	this_input<-cbind(this_study[,c(2:31)],ave(this_study[32],this_study$study,FUN=scale),ave(this_study[33],this_study$study,FUN=scale),ave(this_study[34],this_study$study,FUN=scale),ave(this_study[35],this_study$study,FUN=scale),this_study[,ncol(this_study)]);names(this_input)[35]<-"Response_binary"
	this_mod4<-glm(Response_binary~.,data=this_input,family=binomial) 
	variances_out[4,column]<-(PseudoR2(this_mod4)-PseudoR2(this_mod3))
	assign(paste0("variances_out",column),variances_out,envir = .GlobalEnv)
}


## Panel=Fig2A

# read in biomarker data per sample:
metrics_output<-read.table("./meta_analysis_input_data.txt",sep="\t",header=T,stringsAsFactors=F)

# Part 1) First analyse biomarkers which can be measured using exome genomic data 
# Define studies with exome data
studies<-c("MELANOMA-SNYDER_NEJM_2014","MELANOMA-VANALLEN_SCIENCE_2015","MELANOMA-HUGO_CELL_2016","MELANOMA-RIAZ_CELL_2017","MELANOMA-CRISTESCU_SCIENCE_2018","LUNG-RIZVI_SCIENCE_2015","LUNG-HELLMAN_ALL_2018","BLADDER-SNYDER_PLOSMED_2017","BLADDER-MARIATHASAN_NATURE_2018","BLADDER-CRISTESCU_SCIENCE_2018","RENAL-MCDERMOT_NMED_2018","CRC-Diaz_MSI_CRC","HEAD AND NECK-CRISTESCU_SCIENCE_2018","BREAST-CRISTESCU_SCIENCE_2018")
# Define biomarkers from exome data
cols_for_test<-c(2:31)

# run meta-analysis
run_meta_analysis_pancan(metrics_output,studies,cols_for_test,1)

# Part 2) Analyse biomarkers which can be measured using transcriptome data 
# Define studies with transcriptome data
studies<-c("MELANOMA-SNYDER_NEJM_2014","MELANOMA-VANALLEN_SCIENCE_2015","MELANOMA-HUGO_CELL_2016","MELANOMA-RIAZ_CELL_2017","MELANOMA-CRISTESCU_SCIENCE_2018","BLADDER-SNYDER_PLOSMED_2017","BLADDER-MARIATHASAN_NATURE_2018","BLADDER-CRISTESCU_SCIENCE_2018","RENAL-MCDERMOT_NMED_2018","HEAD AND NECK-CRISTESCU_SCIENCE_2018","BREAST-CRISTESCU_SCIENCE_2018")
# Define biomarkers measured with transcriptome data
cols_for_test<-c(32:35)

# run meta-analysis
run_meta_analysis_pancan(metrics_output,studies,cols_for_test,2)

# Make merged matrix of meta-analysis results, with blank columns/rows between tumour types/biomarker groups 
# All these outputs are returned from the function above
Estimates_ORs_out_set1$blank1<-1;Estimates_ORs_out_set1$blank2<-1;Estimates_ORs_out_set1$blank3<-1;Estimates_ORs_out_set1$blank4<-1;Estimates_ORs_out_set1$blank5<-1;Estimates_ORs_out_set1$blank6<-1;Estimates_ORs_out_set1$blank7<-1
Estimates_ORs_out_set1<-Estimates_ORs_out_set1[,c(1,2,3,4,5,15,8,9,10,16,11,17,13,18,14,19,6,7,17,12)]
Estimates_ORs_out_set2$blank1<-1;Estimates_ORs_out_set2$blank2<-1;Estimates_ORs_out_set2$blank3<-1;Estimates_ORs_out_set2$blank4<-1;Estimates_ORs_out_set2$blank5<-1;Estimates_ORs_out_set2$blank6<-1;Estimates_ORs_out_set2$blank7<-1;Estimates_ORs_out_set2$blank8<-1;Estimates_ORs_out_set2$blank9<-1
Estimates_ORs_out_set2<-Estimates_ORs_out_set2[,c(1,2,3,4,5,12,6,7,8,13,9,14,10,15,11,16,17,18,19,20)]
blank_row1<-Estimates_ORs_out_set1[1,];blank_row1[,]<-1;rownames(blank_row1)<-"blank1"
blank_row2<-Estimates_ORs_out_set1[1,];blank_row2[,]<-1;rownames(blank_row2)<-"blank2"
blank_row3<-Estimates_ORs_out_set1[1,];blank_row3[,]<-1;rownames(blank_row3)<-"blank3"
names(Estimates_ORs_out_set2)<-names(Estimates_ORs_out_set1)
matrix<-rbind(Estimates_ORs_out_set1[1:13,],blank_row1,Estimates_ORs_out_set1[14:23,],blank_row2,Estimates_ORs_out_set1[24:nrow(Estimates_ORs_out_set1),],blank_row3,Estimates_ORs_out_set2)

# Make heatmap plot
# For clearer plotting, ORs of non-significant extreme Odds Ratios values above 10 at capped OR=10, and those below 0.1 capped at OR=0.1 
# none of these OR values are significant, they are unstable point estimates resulting from small sample sizes. This approach is clearly explained in the legend.   
matrix[matrix>=10] <- 10; matrix[matrix<=0.1] <- 0.1
paletteLength <- 400
myColor <- colorRampPalette(c("darkred", "white", "darkblue"))(paletteLength)
myBreaks <- c(seq(min(log2(matrix)), 0, length.out=ceiling(paletteLength/2) + 1),seq(max(log2(matrix))/paletteLength, max(log2(matrix)), length.out=floor(paletteLength/2)))
pheatmap(log2(matrix),color=myColor, breaks=myBreaks,border_color = "white",cluster_cols = F, cluster_rows = F,fontsize_row=6,fontsize_col=6,)

# forest plot
# meta_set1 and meta_set2 are returned from the function above
meta_results<-rbind(meta_set1,meta_set2)
meta_results$biomarker<-rownames(meta_results)
meta_results$biomarker<-factor(meta_results$biomarker,levels=meta_results$biomarker)
fp <- ggplot(data=meta_results,aes(x=biomarker, y=OR, ymin=lower_OR, ymax=upper_OR))+geom_pointrange(colour="black",size=0.6,shape=22,fill="black") + geom_hline(yintercept=1, lty=2) +xlab("") + ylab("Odds Ratio (95% CI)")+ theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, hjust = 1))+scale_x_discrete(limits = rev(levels(meta_results$biomarker)));fp+theme(axis.text.x = element_text(angle = 30, hjust = 1))+coord_flip()


## Panel=Fig2B - meta-analysis by histology

# Define studies for first subgroup - Melanoma - anti-PD1/L1
studies<-c("MELANOMA-HUGO_CELL_2016","MELANOMA-RIAZ_CELL_2017","MELANOMA-CRISTESCU_SCIENCE_2018")
this_subgroup<-"Melanoma_anti_PD1_L1"
run_meta_analysis_histology(meta_results,Estimates_out_set1,Estimates_out_set2,studies,this_subgroup)

# Define studies for second subgroup - Melanoma - anti-CTLA-4
studies<-c("MELANOMA-SNYDER_NEJM_2014","MELANOMA-VANALLEN_SCIENCE_2015")
this_subgroup<-"Melanoma_anti_CTLA_4"
run_meta_analysis_histology(meta_results,Estimates_out_set1,Estimates_out_set2,studies,this_subgroup)

# Define studies for third subgroup - Urothelial anti-PD1/L1
studies<-c("BLADDER-SNYDER_PLOSMED_2017","BLADDER-MARIATHASAN_NATURE_2018","BLADDER-CRISTESCU_SCIENCE_2018")
this_subgroup<-"Urothelial_anti_PD1_L1"
run_meta_analysis_histology(meta_results,Estimates_out_set1,Estimates_out_set2,studies,this_subgroup)

# Define studies for fourth subgroup - Lung anti-PD1/L1 (no transcriptome data)
studies<-c("LUNG-RIZVI_SCIENCE_2015","LUNG-HELLMAN_ALL_2018")
this_subgroup<-"Lung_anti_PD1_L1"
meta_results_by_subgroup<-meta_results[1,]
meta_results_by_subgroup$Group<-"Group"
meta_results_by_subgroup_Lung_anti_PD1_L1<-meta_results_by_subgroup[0,]
for(k in 1:nrow(Estimates_out_set1)){
	meta_dat<-data.frame(t(Estimates_out_set1[(k),studies]),t(SEs_out_set1[(k),studies]),studies)
	names(meta_dat)[1]<-"Effect";names(meta_dat)[2]<-"SE"
	a<-metagen(Effect,SE, studlab = studies,sm = "OR",data = meta_dat)
	this_result<-data.frame(exp(a$TE.random),exp(a$lower.random),exp(a$upper.random),a$pval.random,row.names(Estimates_out_set1)[k],this_subgroup)
	names(this_result)<-names(meta_results_by_subgroup)
	tmp<-meta_results_by_subgroup_Lung_anti_PD1_L1
        meta_results_by_subgroup_Lung_anti_PD1_L1<-rbind(tmp,this_result)
	}

# Filter results either significant in pan-cancer analysis (Fig1A), or significant in one subgroup meta-analysis
meta_results_by_subgroup<-rbind(meta_results_by_subgroup_Melanoma_anti_PD1_L1,meta_results_by_subgroup_Melanoma_anti_CTLA_4,meta_results_by_subgroup_Urothelial_anti_PD1_L1,meta_results_by_subgroup_Lung_anti_PD1_L1)
subgroup_significant<-meta_results_by_subgroup[meta_results_by_subgroup$Meta_Pval<0.05,]
pancan_significant<-meta_results_by_subgroup[meta_results_by_subgroup$biomarker %in% rownames(meta_results[meta_results$Meta_Pval<0.05,]),]
merged_significant<-rbind(subgroup_significant,pancan_significant)
merged_significant<-merged_significant[!duplicated(merged_significant),]
merged_significant<-merged_significant[complete.cases(merged_significant$biomarker),]
# Correct missing NA value (no SERPINB3 mutation in lung subgroup) with null result (OR=1) 
merged_significant[is.na(merged_significant)]<-1
merged_significant$Group<-factor(merged_significant$Group,levels = c("Melanoma_anti_PD1_L1","Melanoma_anti_CTLA_4","Urothelial_anti_PD1_L1","Lung_anti_PD1_L1"),ordered = TRUE)
ggplot(merged_significant,aes(x=biomarker, y=OR,fill=Group,color="black"))+geom_point(shape=24,size=3)+facet_grid(~Group,scales = "free_x",space = "free_x")+theme(strip.background = element_rect(fill="lightgrey"),axis.text.x=element_text(angle=45,hjust=1,size=8),plot.title = element_text(size = 8),strip.text = element_text(size = 8),legend.position="none",panel.grid.major.x = element_line(colour = "grey50",size=0.5,linetype = "dashed"))+scale_fill_manual(values=c("#F1BB7B", "#FD6467", "#5B1A18", "black","darkblue"))+ geom_hline(yintercept=1)+ylim(c(0,5))


## Panel=Fig2C - Correlation plot
# use samples with at least TMB and CD8A data
metrics_output2<-metrics_output[!(is.na(metrics_output$TMB)) & !(is.na(metrics_output$CD8A)),]
# Load purity data
coverage<-read.table("./coverage_purity_stats.csv",header=T,sep=",",stringsAsFactors=F)
metrics_output2<-merge(metrics_output2,coverage,by.x="case",by.y="Patient")
metrics_output2$coverage<-NULL

# Select continous variables for correlation analysis
metrics_output2<-select_if(metrics_output2, is.numeric)
corr_out<-data.frame(matrix(0, ncol = ncol(metrics_output2), nrow = ncol(metrics_output2)))

# Calculate correlation coefficients
for(i in 1:ncol(metrics_output2)){
	for(j in 1:ncol(metrics_output2)){
		a<-cor.test(metrics_output2[,i],metrics_output2[,j],method="spearman")
		corr_out[i,j]<-a$estimate
					}			
				}
row_names<-colnames(metrics_output2)
rownames(corr_out)<-row_names
names(corr_out)<-row_names

# Plot correlation heatmap
paletteLength <- 400
myColor <- colorRampPalette(c("darkblue", "white", "darkred"))(paletteLength)
myBreaks <- c(seq(-1, 0, length.out=ceiling(paletteLength/2) + 1),seq(max(corr_out)/paletteLength, max(corr_out), length.out=floor(paletteLength/2)))
pheatmap(corr_out,color=myColor, breaks=myBreaks,border_color = "white",cluster_cols = F, cluster_rows = F,fontsize_row=6,fontsize_col=6,)


## Panel=Fig2D
# proportion of variance explained plot - create output matrix
variances_out<-data.frame(matrix(0, ncol = 4, nrow = 4))
names(variances_out)<-c("Melanoma_anti_PD1","Melanoma_anti_CTLA4","Urothelial_anti_PDL1_PD1","Lung_anti_PD1")
rownames(variances_out)<-c("sources_of_antigen","immune_evasion","host_factors","infiltration")

# Calculate variance for first subgroup - Melanoma - anti-PD-1
studies<-c("MELANOMA-HUGO_CELL_2016","MELANOMA-RIAZ_CELL_2017","MELANOMA-CRISTESCU_SCIENCE_2018")
calculate_variance(studies,1)

# Calculate variance for second subgroup - Melanoma - anti-CTLA-4
studies<-c("MELANOMA-SNYDER_NEJM_2014","MELANOMA-VANALLEN_SCIENCE_2015")
calculate_variance(studies,2)

# Calculate variance for third subgroup
studies<-c("BLADDER-SNYDER_PLOSMED_2017","BLADDER-MARIATHASAN_NATURE_2018","BLADDER-CRISTESCU_SCIENCE_2018")
calculate_variance(studies,3)

# Calculate variance for fourth subgroup 
studies<-c("LUNG-RIZVI_SCIENCE_2015","LUNG-HELLMAN_ALL_2018","LUNG-CRISTESCU_SCIENCE_2018")
this_tab<-metrics_output
this_tab$hist_study<-paste0(this_tab$histology,"-",this_tab$study)
this_study<-this_tab[this_tab$hist_study %in% studies,]
this_study$Response_binary<-1
this_study$Response_binary[this_study$response=="no_response"]<-0
#sources of antigen
this_input<-this_study[,c(2:14,ncol(this_study))]
this_mod1<-glm(Response_binary~.,data=this_input,family=binomial) 
variances_out[1,4]<-PseudoR2(this_mod1)
#sources of antigen+immune evasasion
this_input<-this_study[,c(2:24,ncol(this_study))]
this_mod2<-glm(Response_binary~.,data=this_input,family=binomial) 
variances_out[2,4]<-(PseudoR2(this_mod2)-PseudoR2(this_mod1))
#sources of antigen+immune evasasion+host
this_input<-this_study[,c(2:31,ncol(this_study))]
this_mod3<-glm(Response_binary~.,data=this_input,family=binomial) 
variances_out[3,4]<-(PseudoR2(this_mod3)-PseudoR2(this_mod2))

#plot
variances_out[1]<-variances_out1[1];variances_out[2]<-variances_out2[2];variances_out[3]<-variances_out3[3]
variances_out$metric<-row.names(variances_out)
variances_out_long<-melt(variances_out,id.vars=c("metric"))
myColor <- colorRampPalette(c("lightgreen", "darkgreen"))(4)
ggplot(variances_out_long, aes(factor(variable), value, fill = factor(metric))) +geom_bar(stat = "identity") +theme(axis.text.x = element_text(angle = 45,hjust = 1))+ scale_fill_manual(values =myColor)

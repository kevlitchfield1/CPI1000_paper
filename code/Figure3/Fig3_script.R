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
library("glmnet")
library("ROCR")
library("xgboost")
library("DiagrammeR")
library("pROC")
rm(list=ls(all=TRUE))

#### Panel=Fig3A

# read in biomarker data per sample:
metrics_output<-read.table("./meta_analysis_input_data.txt",sep="\t",header=T,stringsAsFactors=F)
# define and extract biomarkers which are significant in last column of Fig 2A:
cols_uni_sig <- c("TMB","CXCL9","Clonal_TMB","CD8A","CD274","NMD_escape_TMB","Indel_TMB","Signature.4","Signature.7","signature_apobec","T_inflam_GEP","Is_Male")
# all biomarkers for the multivariable model to be scaled, TMB is left as absolute, for benchmarking purposes
cols_to_scale<-c(2:12)
mva_ind_set<-metrics_output[,colnames(metrics_output) %in% cols_uni_sig]
mva_ind_set<-cbind(mva_ind_set,metrics_output$case,metrics_output$histology,metrics_output$study,metrics_output$response)
names(mva_ind_set)[ncol(mva_ind_set)-3]<-"case";names(mva_ind_set)[ncol(mva_ind_set)-2]<-"histology";names(mva_ind_set)[ncol(mva_ind_set)-1]<-"study";names(mva_ind_set)[ncol(mva_ind_set)]<-"response"
mva_ind_set$hist_study<-paste0(mva_ind_set$histology,"-",mva_ind_set$study)
mva_set_scaled<-mva_ind_set
mva_set_scaled[cols_to_scale]<-ave(mva_set_scaled[cols_to_scale],mva_set_scaled$hist_study,FUN=scale)
mva_set_scaled$Response_binary<-1;mva_set_scaled$Response_binary[mva_set_scaled$response=="no_response"]<-0

# Show distribution of feature importance scores per study using Monte Carlo sampling, in the biggest cohort with DNA+transcriptmoe for the main 4 tumour types (mela,bladd,H&N,renal) with n>25
studies<-c("BLADDER-MARIATHASAN_NATURE_2018","HEAD AND NECK-CRISTESCU_SCIENCE_2018","MELANOMA-CRISTESCU_SCIENCE_2018","RENAL-MCDERMOT_NMED_2018")
for(j in 1:length(studies)){
	opt_tmb_auc<-c(); tmb_auc<-c()	
	coefs_template<-data.frame(matrix(0, ncol = 1, nrow = length(cols_uni_sig)))
	coefs_template$Feature<-names(mva_ind_set)[1:length(cols_uni_sig)]
	coefs_out<-data.frame(matrix(0, ncol = 1, nrow = length(cols_uni_sig)))
	coefs_out$Feature<-names(mva_ind_set)[1:length(cols_uni_sig)]
	coefs_out <- coefs_out[order(coefs_out$Feature),]
	coefs_out<-coefs_out[,c(2,1)]
	names(coefs_out)[2]<-"Gain"
	
for (i in 1:1000){
	this_study<-mva_ind_set[mva_ind_set$hist_study==as.character(studies[j]),]
	this_study$Response_binary<-1;this_study$Response_binary[this_study$response=="no_response"]<-0	
	#scale variables within each study
	for(k in 1:length(cols_to_scale)){this_study[cols_to_scale[k]]<-as.numeric(scale(this_study[cols_to_scale[k]]))}	
	train_ind <- sample(seq_len(nrow(this_study)), size = (nrow(this_study)*0.75))
	train<-this_study[train_ind,];test<-this_study[-train_ind,]
	while(sum(test$Response_binary)==0){train_ind <- sample(seq_len(nrow(this_study)), size = (nrow(this_study)*0.75));train<-this_study[train_ind,];test<-this_study[-train_ind,]}
	mva_mod <- xgboost(data = as.matrix(train[,2:length(cols_uni_sig)]), label = train$Response_binary, nthread = 2,max.depth = 2,nrounds = 15,eta=0.2, objective = "binary:logistic",verbose=F)
	importance_matrix <- xgb.importance(model = mva_mod)
	these_coefs<-merge(coefs_template,importance_matrix,all.x=T)
	these_coefs[is.na(these_coefs)] <- 0
	these_coefs <- these_coefs[order(these_coefs$Feature),]
	coefs_out_tmp<-coefs_out; coefs_out<-rbind(coefs_out_tmp,these_coefs[,c(1,3)])
	}
	coefs_out<-coefs_out[-(1:length(cols_uni_sig)),]
	coefs_out<-coefs_out[coefs_out$Feature!="TMB",]
	plot_coe<-ggplot(data=coefs_out,aes(x = reorder(Feature, -Gain), y=Gain)) + geom_boxplot(fill = "darkblue",outlier.shape=NA)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+ theme(plot.margin = unit(c(1,1,1,1), "lines"),plot.title = element_text(size = 7, face = "bold"),axis.text.x = element_text(angle = 30,hjust = 1,size=7))+ ggtitle(paste0("Cohort: ",studies[j]))+ylab("")+xlab("")+scale_y_continuous()
	assign(paste("c",j, sep = "_"), plot_coe)		
			}
iteration_names<-c("tmp")
for(k in 1:length(studies)){iteration_names<-c(iteration_names,rep(as.character(studies[k]),i))}
plot_grid(c_1,c_2,c_3,c_4,nrow=2)
ggsave("./coeff_plots_cell_aug2020.pdf.pdf",dpi=600)

#### Panel=Fig3B/C
# Train final multivariable model, and identical model with just TMB for comparison purposes 
mva_mod <- xgboost(data = as.matrix(mva_set_scaled[,2:length(cols_uni_sig)]), label = mva_set_scaled$Response_binary,objective = "binary:logistic",verbose=F,max.depth = 2,nrounds = 15,eta=0.2)
TMB_mod <- xgboost(data = as.matrix(mva_set_scaled[,1]), label = mva_set_scaled$Response_binary,objective = "binary:logistic",verbose=F,max.depth = 2, nrounds = 15,eta=0.2)
importance_matrix <- xgb.importance(model = mva_mod)
xgb.plot.importance(importance_matrix = importance_matrix)

#### Panel=Fig3D

# Read in and analyse first test cohort
mva_test_set1<-read.table("./keynote_other_tumour_type_test_set.txt",sep="\t",header=T,stringsAsFactors=F)
mva_test_set1_scaled<-mva_test_set1
mva_test_set1_scaled[cols_to_scale]<-scale(mva_test_set1_scaled[cols_to_scale])
mva_test_set1_scaled$Response_binary<-1;mva_test_set1_scaled$Response_binary[mva_test_set1_scaled$response=="no_response"]<-0
#test mva model
pred <- predict(mva_mod, as.matrix(mva_test_set1_scaled[,2:length(cols_uni_sig)]))
pr <- prediction(as.numeric(pred),mva_test_set1_scaled$Response_binary)
auc <- performance(pr, measure = "auc")
auc_val<-auc@y.values[[1]]
#test tmb model
pred2 <- predict(TMB_mod, as.matrix(mva_test_set1_scaled[,1]))
pr2 <- prediction(as.numeric(pred2),mva_test_set1_scaled$Response_binary)
auc2 <- performance(pr2, measure = "auc")
auc_val2<-auc2@y.values[[1]]
#ROC curves 
auc <- performance(pr,"tpr","fpr");auc2 <- performance(pr2,"tpr","fpr")
plot(auc,col="darkblue",lwd=3,xaxt="n")+axis(2, at = seq(0.0,1.0, by = 0.2), las=2)+mtext(paste0("Multivariate AUC: ",round(auc_val,2)),col="darkblue",line=-8,cex=0.8)+ abline(coef = c(0,1))
plot(auc2,add=T,col="darkred",lwd=3)+mtext(paste0("TMB AUC: ",round(auc_val2,2)),col="darkred",line=-7,cex=0.8)
# compare curves
roc1 <- roc(mva_test_set1_scaled$Response_binary, as.numeric(pred))
roc2 <- roc(mva_test_set1_scaled$Response_binary,as.numeric(pred2))
print(roc.test(roc1, roc2))

# Read in and analyse second test cohort
mva_test_set2<-read.table("./mva_data_EVA.txt",header=T,stringsAsFactors=F,sep="\t")
mva_test_set2_scaled<-mva_test_set2[,colnames(mva_test_set2) %in% cols_uni_sig]
mva_test_set2_scaled<-cbind(mva_test_set2_scaled,mva_test_set2$case,mva_test_set2$histology,mva_test_set2$study,mva_test_set2$response)
names(mva_test_set2_scaled)[ncol(mva_test_set2_scaled)-3]<-"case";names(mva_test_set2_scaled)[ncol(mva_test_set2_scaled)-2]<-"histology";names(mva_test_set2_scaled)[ncol(mva_test_set2_scaled)-1]<-"study";names(mva_test_set2_scaled)[ncol(mva_test_set2_scaled)]<-"response"
mva_test_set2_scaled$hist_study<-paste0(mva_test_set2_scaled$histology,"-",mva_test_set2_scaled$study)
mva_test_set2_scaled$Response_binary<-1;mva_test_set2_scaled$Response_binary[mva_test_set2_scaled$response=="no_response"]<-0
for(k in 1:length(cols_to_scale)){mva_test_set2_scaled[cols_to_scale[k]]<-as.numeric(scale(mva_test_set2_scaled[cols_to_scale[k]]))}	
pred <- predict(mva_mod, as.matrix(mva_test_set2_scaled[,2:length(cols_uni_sig)]))
pr <- prediction(as.numeric(pred),mva_test_set2_scaled$Response_binary)
auc <- performance(pr, measure = "auc")
auc_val3<-auc@y.values[[1]]
pred2 <- predict(TMB_mod, as.matrix(mva_test_set2_scaled[,1]))
pr2 <- prediction(as.numeric(pred2),mva_test_set2_scaled$Response_binary)
auc2 <- performance(pr2, measure = "auc")
auc_val4<-auc2@y.values[[1]]
#ROC curves for test cohort 2
auc <- performance(pr,"tpr","fpr");auc2 <- performance(pr2,"tpr","fpr")
plot(auc,col="darkblue",lwd=3,xaxt="n")+axis(2, at = seq(0.0,1.0, by = 0.2), las=2)+mtext(paste0("Multivariate AUC: ",round(auc_val3,2)),col="darkblue",line=-8,cex=0.8)+ abline(coef = c(0,1))
plot(auc2,add=T,col="darkred",lwd=3)+mtext(paste0("TMB AUC: ",round(auc_val4,2)),col="darkred",line=-7,cex=0.8)
# compare curves
roc1 <- roc(mva_test_set2_scaled$Response_binary, as.numeric(pred))
roc2 <- roc(mva_test_set2_scaled$Response_binary,as.numeric(pred2))
print(roc.test(roc1, roc2))

# Read in and analyse third test cohort (which only has TMB/PDL1/Smoking data/Sex)
mva_test_set3<-read.table("./mva_data_AO.txt",header=T,stringsAsFactors=F,sep="\t")
cols_test_cohort3 <- c("TMB","CD274","Signature.4","Is_Male")
cols_to_scale<-c(1:4)
mva_test_set3_scaled<-mva_test_set3[,colnames(mva_test_set3) %in% cols_test_cohort3]
mva_test_set3_scaled<-cbind(mva_test_set3_scaled,mva_test_set3$case,mva_test_set3$histology,mva_test_set3$study,mva_test_set3$response)
names(mva_test_set3_scaled)[ncol(mva_test_set3_scaled)-3]<-"case";names(mva_test_set3_scaled)[ncol(mva_test_set3_scaled)-2]<-"histology";names(mva_test_set3_scaled)[ncol(mva_test_set3_scaled)-1]<-"study";names(mva_test_set3_scaled)[ncol(mva_test_set3_scaled)]<-"response"
mva_test_set3_scaled$hist_study<-paste0(mva_test_set3_scaled$histology,"-",mva_test_set3_scaled$study)
mva_test_set3_scaled$Response_binary<-1
mva_test_set3_scaled$Response_binary[mva_test_set3_scaled$response=="no_response"]<-0
for(k in 1:length(cols_to_scale)){mva_test_set3_scaled[cols_to_scale[k]]<-as.numeric(scale(mva_test_set3_scaled[cols_to_scale[k]]))}	
#Build new model with the four variables where we have data
names(mva_set_scaled)[names(mva_set_scaled) == "TMB"] <- "TMB_unscaled"
mva_set_scaled$TMB<-ave(mva_set_scaled$TMB_unscaled,mva_set_scaled$hist_study,FUN=scale)
mva_mod <- xgboost(data = as.matrix(mva_set_scaled[,cols_test_cohort3]), label = mva_set_scaled$Response_binary,objective = "binary:logistic",nthread = 2,verbose=F,max.depth = 2,nrounds = 15,eta=0.2)
importance_matrix <- xgb.importance(model = mva_mod)
xgb.plot.importance(importance_matrix = importance_matrix)
#Predict for third test cohort
pred <- predict(mva_mod, as.matrix(mva_test_set3_scaled[,cols_test_cohort3]))
pr <- prediction(as.numeric(pred),mva_test_set3_scaled$Response_binary)
auc <- performance(pr, measure = "auc")
auc_val3<-auc@y.values[[1]]
pred2 <- predict(TMB_mod, as.matrix(mva_test_set3[,2]))
pr2 <- prediction(as.numeric(pred2),mva_test_set3_scaled$Response_binary)
auc2 <- performance(pr2, measure = "auc")
auc_val4<-auc2@y.values[[1]]

#ROC curves
auc <- performance(pr,"tpr","fpr");auc2 <- performance(pr2,"tpr","fpr")
plot(auc,col="darkblue",lwd=3,xaxt="n")+axis(2, at = seq(0.0,1.0, by = 0.2), las=2)+mtext(paste0("Multivariate AUC: ",round(auc_val3,2)),col="darkblue",line=-8,cex=0.8)+ abline(coef = c(0,1))
plot(auc2,add=T,col="darkred",lwd=3)+mtext(paste0("TMB AUC: ",round(auc_val4,2)),col="darkred",line=-7,cex=0.8)
# compare curves
roc1 <- roc(mva_test_set3_scaled$Response_binary, as.numeric(pred))
roc2 <- roc(mva_test_set3_scaled$Response_binary,as.numeric(pred2))
print(roc.test(roc1, roc2))

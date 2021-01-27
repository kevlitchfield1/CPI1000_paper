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

## Fig=6A

dat<-read.table("./amps_to_plot_data.txt",header=T,sep="\t",stringsAsFactors=F)
dat_m<-melt(dat,measure.vars = c("RR_mut","RR_wt"))
ggplot(data=dat_m, aes(x=reorder(gene,-Diff), y=value,color=variable,fill=variable)) +geom_bar(stat="identity", position=position_dodge())+ theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values=c("black", "lightgrey"))+scale_color_manual(values=c("black", "lightgrey"))

dat<-read.table("./deep_dels_to_plot_data.txt",header=T,sep="\t",stringsAsFactors=F)
dat_m<-melt(dat,measure.vars = c("RR_mut","RR_wt"))
ggplot(data=dat_m, aes(x=reorder(gene,-Diff), y=value,color=variable,fill=variable)) +geom_bar(stat="identity", position=position_dodge())+ theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values=c("black", "lightgrey"))+scale_color_manual(values=c("black", "lightgrey"))


## Fig=6B
dat<-read.table("./CCND1_histology_counts.txt",header=T,stringsAsFactors=F,sep="\t")
ggplot(data=dat, aes(x=reorder(Histology,Count), y=Count)) +geom_bar(stat="identity", position=position_dodge(),color="darkgreen",fill="darkgreen")+ theme(axis.text.x = element_text(angle = 45,hjust = 1))+coord_flip()


## Fig=6C
dat<-read.table("./CCND1.txt",header=T,stringsAsFactors=F,sep="\t")
ggplot(data=dat,aes(x=response, y=CCND1_TPM,color=response))+geom_boxplot() +geom_quasirandom(size = 3, shape=18)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means()+ theme(plot.margin = unit(c(1,1,1,1), "lines"))+xlab("")+scale_y_log10()


## Fig=6D
splots <- list()
dat<-read.table("./genie_IO_bladder_CCND1.txt",header=T,stringsAsFactors=F,sep="\t")
summary(coxph(formula=Surv(SURVIVAL_MONTHS,SURVIVAL_EVENT)~CCND1_amp, dat))
timestrata.surv1 <- survfit(Surv(SURVIVAL_MONTHS,SURVIVAL_EVENT)~ CCND1_amp, dat)
splots[[1]]<-ggsurvplot(timestrata.surv1,risk.table = TRUE, risk.table.height=0.25,risk.table.title="",font.x=9,font.y=9,font.tickslab=10,legend = "none",xlab="Time (months)",ylab="",tables.theme = theme_cleantable(),palette = c("#D55E00","#009E73"))


## Fig=6E
dat<-read.table("./genie_nonIO_bladder_CCND1.txt",header=T,stringsAsFactors=F,sep="\t")
summary(coxph(formula=Surv(SURVIVAL_MONTHS,SURVIVAL_EVENT)~CCND1_amp, dat))
timestrata.surv1 <- survfit(Surv(SURVIVAL_MONTHS,SURVIVAL_EVENT)~ CCND1_amp, dat)
splots[[2]]<-ggsurvplot(timestrata.surv1,risk.table = TRUE, risk.table.height=0.25,risk.table.title="",font.x=9,font.y=9,font.tickslab=10,legend = "none",xlab="Time (months)",ylab="",tables.theme = theme_cleantable(),palette = c("#D55E00","#009E73"))


## Fig=6F
dat<-read.table("./genie_IO_pancan_CCND1.txt",header=T,stringsAsFactors=F,sep="\t")
summary(coxph(formula=Surv(SURVIVAL_MONTHS,SURVIVAL_EVENT)~CCND1_amp, dat))
timestrata.surv1 <- survfit(Surv(SURVIVAL_MONTHS,SURVIVAL_EVENT)~ CCND1_amp, dat)
splots[[3]]<-ggsurvplot(timestrata.surv1,risk.table = TRUE, risk.table.height=0.25,risk.table.title="",font.x=9,font.y=9,font.tickslab=10,legend = "none",xlab="Time (months)",ylab="",tables.theme = theme_cleantable(),palette = c("#D55E00","#009E73"))


## Fig=6G
dat<-read.table("./genie_nonIO_pancan_CCND1.txt",header=T,stringsAsFactors=F,sep="\t")
summary(coxph(formula=Surv(SURVIVAL_MONTHS,SURVIVAL_EVENT)~CCND1_amp, dat))
timestrata.surv1 <- survfit(Surv(SURVIVAL_MONTHS,SURVIVAL_EVENT)~ CCND1_amp, dat)
splots[[4]]<-ggsurvplot(timestrata.surv1,risk.table = TRUE, risk.table.height=0.25,risk.table.title="",font.x=9,font.y=9,font.tickslab=10,legend = "none",xlab="Time (months)",ylab="",tables.theme = theme_cleantable(),palette = c("#D55E00","#009E73"))

arrange_ggsurvplots(splots, print = TRUE,ncol = 4, nrow = 1, risk.table.height = 0.4)

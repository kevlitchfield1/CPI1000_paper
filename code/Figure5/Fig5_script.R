library("ggplot2")
library("cowplot")

rm(list=ls(all=TRUE))

## Fig=5A
#plotted from raw data on HPC using R copynumber package


## Fig=5B
# read in data
case_list_for_scna<-read.table("./case_scna_analysis.txt",header=T,sep="\t",stringsAsFactors=F)
counts<-table(case_list_for_scna$Response)
count_non_resp<-counts[1]
count_resp<-counts[2]
cyto_dat<-read.table("./cytoband_counts.txt",header=T,sep="\t",stringsAsFactors=F)

# do Fishers Exact test per cytoband
cyto_dat$del_pval<-1
cyto_dat$amp_pval<-1
for(i in 1:nrow(cyto_dat)){
	res_loss <- fisher.test(matrix(c(cyto_dat[i,2],count_resp-cyto_dat[i,2],cyto_dat[i,3],count_non_resp-cyto_dat[i,3]),nrow = 2))
	cyto_dat[i,6]<-res_loss$p.value
	res_gain <- fisher.test(matrix(c(cyto_dat[i,4],count_resp-cyto_dat[i,4],cyto_dat[i,5],count_non_resp-cyto_dat[i,5]),nrow = 2))
	cyto_dat[i,7]<-res_gain$p.value
			}
cyto_dat$del_qval<-p.adjust(cyto_dat$del_pval, "fdr")
cyto_dat$amp_qval<-p.adjust(cyto_dat$amp_pval, "fdr")

del_hits<-head(cyto_dat[order(cyto_dat$del_pval),],10)
amp_hits<-head(cyto_dat[order(cyto_dat$amp_pval),],10)
del_hits$freq_diff<-((del_hits$loss_count_responders/count_resp)-(del_hits$loss_count_non_responders/count_non_resp))*100
amp_hits$freq_diff<-((amp_hits$gain_count_responders/count_resp)-(amp_hits$gain_count_non_responders/count_non_resp))*100

# plot results
p1<-ggplot(data=del_hits, aes(x=reorder(cytoband,freq_diff), y=freq_diff)) +geom_bar(stat="identity", position=position_dodge(),fill="darkred")+ theme_minimal()+coord_flip()
p2<-ggplot(data=amp_hits, aes(x=reorder(cytoband,freq_diff), y=freq_diff)) +geom_bar(stat="identity", position=position_dodge(),fill="darkblue")+ theme_minimal()+coord_flip()
plot_grid(p1,p2)

#ggsave("./panel_5b.pdf",dpi=600)


## Panel=Fig5C

#make 9q34 plot
out_del_to_plot<-read.table("./9q34_deletion_freqs.txt",header=T,sep="\t",stringsAsFactors=F)
out_del_to_plot_m<-melt(out_del_to_plot,measure.vars = c("resp_del_freq","non_resp_del_freq"))
ggplot(out_del_to_plot_m) + geom_line(aes(y = -freq_diff, x = chr9_Mb),color="black",size=0.6,data = out_del_to_plot, stat="identity")+scale_y_continuous()+theme_bw()

## Panel=Fig5D

#make freq plot by histology
hist_freq<-read.table("./TRAF2_hist_freq.txt",header=T,stringsAsFactors=F)
hist_freq_m<-melt(hist_freq,measure.vars = c("X._del_resp","X._del_non_resp"))
hist_freq_m$Tumour_Type<-ordered(hist_freq_m$Tumour_Type,levels=c("Other","Melanoma","Bladder","Pan-cancer"))
hist_freq_m$variable<-ordered(hist_freq_m$variable,levels=c("X._del_non_resp","X._del_resp"))
ggplot(data=hist_freq_m, aes(x=Tumour_Type, y=value,color=variable,fill=variable)) +geom_bar(stat="identity", position=position_dodge())+ theme_classic()+coord_flip()+scale_fill_manual(values=c("darkred","darkblue"))+scale_color_manual(values=c("darkred","darkblue"))


## Panel=Fig5E
dat<-read.table("./haplo_dat.txt",header=T,sep="\t",stringsAsFactors=F)
dat$order<-rownames(dat)
dat$order<-as.numeric(rownames(dat))
ggplot(dat, aes(order,prob_haploinsufficient,color = prob_haploinsufficient))+geom_point(shape=1, size=0.75)+theme_minimal()+scale_color_gradient(low = "pink", high = "darkred")


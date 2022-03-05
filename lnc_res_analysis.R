library(ggplot2)
library(Biostrings)
setwd("D:\\sorf_models\\model\\tests")
data<-read.csv("lnc_res.txt",sep="\t",header = T,stringsAsFactors = F,quote = "")[,-1]
data<-data[data$lncRNA!="HSALNT0245581",]
overall_cutoff<-quantile(data$overall_foldchange,0.95)
self_cutoff<-quantile(data$self_foldchange,0.95)

print(paste0("overall_foldchange mean: ",mean(data$overall_foldchange)))
print(paste0("self_foldchange mean: ",mean(data$self_foldchange)))

lncs<-readDNAStringSet("all_lnc_new_plane.fa")
my_lnc<-read.csv("expression_profiles_HPA/max10955lnc.txt",header = F,stringsAsFactors = F)[,1]
lncs<-lncs[names(lncs) %in% my_lnc]
lncs_s<-lncs[width(lncs)<250]
lncs_l<-lncs[width(lncs)>=250]
lncs_l<-subseq(lncs_l,start=1,width=250)
lncs1<-c(lncs_s,lncs_l)

GCs<-(as.integer(letterFrequency(lncs1,letters = "G"))+as.integer(letterFrequency(lncs1,letters = "C")))/width(lncs)
GCs<-data.frame(lncRNA=names(lncs1),GC=GCs)
data1<-merge(data,GCs)
data1<-data1[order(data1$overall_foldchange,decreasing = T),]
R=cor(x = data1$overall_foldchange,y=data1$GC)

graphics.off()
pdf("lnc_res.pdf",width = 6,height = 5)
ggplot()+
  geom_rect(aes(xmin=overall_cutoff,xmax=Inf,ymin=1,ymax=Inf),fill="#ffe36c",alpha=0.25)+
  geom_rect(aes(xmin=1,xmax=Inf,ymin=self_cutoff,ymax=Inf),fill="#69c2ff",alpha=0.25)+
  geom_point(data=data1,aes(x=overall_foldchange,y=self_foldchange),size=0.5,alpha=0.5,shape=16)+
  geom_vline(xintercept = overall_cutoff,size=0.1)+
  geom_hline(yintercept = self_cutoff,size=0.1)+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  labs(x="mTIS foldchange",y="Self foldchange")+
  coord_cartesian(xlim = c(1,1.5))
dev.off()


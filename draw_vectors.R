library(Biostrings)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(ggpmisc)
setwd("D:\\sorf_models\\model\\tests\\final_true_vs_predicts")
fas<-readDNAStringSet("../gencode.v38.annotation.gff3_gffreadto_simple_plane_new.fa")

if (TRUE) {
  numpy<-read.csv("inspect.txt",header = T)
  t=strsplit(names(numpy),split='X..')[[1]][2]
  numpy<-numpy[,1]
  
  plotdata<-data.frame(x=1:length(numpy),y=numpy)
  graphics.off()
  pdf(paste0("draw_vector_",t,".pdf"),width = 50,height = 1)
  print(ggplot(data=plotdata[1:min(length(numpy),2000),])+
    geom_bar(aes(x=x,y=y),stat = "identity")+
    scale_x_continuous(breaks=1:length(numpy),labels=strsplit(as.character(fas[t]),split="")[[1]])+
    coord_cartesian(expand = c(0,0))+
    theme(line = element_blank(),
          title = element_blank(),
          axis.text.y = element_text(size=4),
          axis.text.x = element_text(size=0.7,vjust=17)))
  dev.off()
}


if (TRUE) {
  numpy<-read.csv("inspect.txt",header = T)
  output<-read.csv("inspect_output.txt",header = T)
  pos<-read.csv("inspect_pos.txt",header = T)
  t=strsplit(names(numpy),split='X..')[[1]][2]
  numpy<-numpy[,1]
  output<-output[,1]
  pos<-pos[,1]
  
  plotdata0<-data.frame(x=1:length(numpy),y=numpy,class=rep("true",length(numpy)))
  plotdata0=plotdata0[1:min(length(numpy),2000),]
  plotdata1<-data.frame(x=pos+1,y=-output*max(numpy)/max(output),class=rep("predicted",length(pos)))
  plotdata<-rbind(plotdata0,plotdata1)
  graphics.off()
  pdf(paste0("draw_vector_",t,"_true_vs_predicted_finalrun_5_ratio_best.pdf"),width = 10,height = 1)
  print(ggplot(data=plotdata)+
          geom_bar(aes(x=x,y=y,group=class,fill=class),stat = "identity",position = 'identity',alpha=1)+
          scale_x_continuous(breaks=1:length(numpy),labels=strsplit(as.character(fas[t]),split="")[[1]])+
          coord_cartesian(expand = c(0,0))+
          theme(line = element_blank(),
                title = element_blank(),
                axis.text.y = element_text(size=4),
                axis.text.x = element_text(size=0.7,vjust=17)))
  dev.off()
}

sigs<-c()
iters<-c()
outputs<-list()
extra<-c()
anno<-read.csv("../mTIS_identity/all_tran_CDS_positions.txt",sep = "\t",header = F,stringsAsFactors = F,quote = "")[,1:2]
names(anno)<-c("transcripts","annotated_positions")
ii<-1

if (TRUE) {
  seq<-read.csv("FXII_seq.txt",sep = "\t",stringsAsFactors = F,header = F)[ii,]
  iter<-seq[,1]
  start<-seq[,2]
  seq<-seq[,3]
  numpy<-read.csv("inspect.txt",header = T)
  output<-read.csv("inspect_output.txt",header = T,sep="\t")
  for (i in 1:ncol(output)) {
    output[,i]<-output[,i]/sum(abs(output[,i]))
  }
  pos<-read.csv("inspect_pos.txt",header = T)
  t=strsplit(names(numpy),split='X..')[[1]][2]
  #start<-anno$annotated_positions[anno$transcripts==t]
  #start<-82
  #seq<-"ATTAGAGTCTGTGCTTCACTTCCGTTCCAGCCTCAGCGGCAGCTGGATCGCTCGACGGAGTGCCTCTGGTAGTTGGCCAAGACGCCGAATATCAAAATCTTCAGCGGCAGCTCCCACCAGGACTTATCCCAGAAAATTGCTGACCGCCTGGGCCTGGAGCTAGGCAAGGTGGTGACTAAGAAATTCAGCAACCAGGAGACCTGCGTGGAAATTGATGAGAGTGTGCGTGGAGAGGATGTCTACATCGTTCAGAGTGGTTGTGGCGAAATCAACGACAGTCTAATGGAGCTTTTGATCATGATTAATGCCTGCAAGATTGCTTCAGCTAGCCGAGTTACTGCAGTCATCCCATGCTTCCCTTATGCCCGACAGGATAAGAAGGATAAGAGCCGGTCCCCAATCTCTGCCAAGCTTGTTGCAAATATG"
  #seq<-as.character(fas[t])
  numpy<-numpy[,1]
  pos<-pos[,1]
  
  plotdata1<-data.frame(x=pos+1,y=rowMeans(output),fill="#636363")
  plotdata1[which.max(plotdata1$y),"fill"]<-"red"
  mymax<-function(x,thresh){max(x,thresh)}
  mymax1<-Vectorize(mymax,vectorize.args = "x")
  cutoff<-5*mean(mymax1(plotdata1$y[order(plotdata1$y,decreasing = F)][1:(length(plotdata1$y)/1.5)],0.0001))
  sTIS_index<-which(plotdata1$y>cutoff)
  sTIS<-plotdata1$y[sTIS_index]
  sTIS_index<-sTIS_index[order(sTIS,decreasing = T)][2:min(8,length(sTIS_index))]
  plotdata1[sTIS_index,"fill"]<-"orange"
  
  plotdata2<-cbind(x=rep(pos+1,ncol(output)),y=melt(output,id.vars = NULL))
  plotdata3<-data.frame(x=rep(start,2),y=c(-0.5,-1))
  plotdata4<-data.frame(x=start,y=0)
  
  custom<-ggplotGrob(ggplot()+geom_line(data=plotdata3,aes(x=x-0.5,y=y),color="dark green",lineend = "butt", linejoin = "mitre",
                                        arrow = arrow(length=unit(1,"mm"), ends="first"))+
                       geom_text(data=plotdata4,aes(label=x,x=x-0.5,y=y,vjust=1,hjust=0.5),color="black",size=2)+
                       scale_x_continuous(limits = c(0,max(plotdata1$x)+2.5),expand = c(0,0))+
                       theme(panel.background = element_blank(),
                             panel.grid = element_blank(),
                             axis.title = element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             plot.margin = unit(c(0,0,0,0),"cm"),
                             axis.line = element_blank(),
                             panel.border = element_blank()))
  
  
  graphics.off()
  scale_factor<-0.1929636/max(plotdata2$y.value)
  pdf(paste0("draw_vector_Calvo_",t,"_predicted_finalrun_all10_ratio_best.pdf"),width = 20,height = 2)
  print(ggplot()+
          geom_bar(data=plotdata1,aes(x=x,y=y,fill=fill),stat = "identity",position = 'identity',alpha=1)+
          geom_point(data=plotdata2,aes(x=x,y=y.value,color=y.variable))+
          geom_line(data=plotdata2,aes(x=x,y=y.value,color=y.variable),alpha=0.5)+
          geom_bar(data=plotdata1,aes(x=x,y=y,fill=fill),stat = "identity",position = 'identity',alpha=0.5)+
          geom_abline(intercept = cutoff,slope = 0,size=0.2,color="#034e7b")+
          annotation_custom(custom,xmin = 0,xmax = max(plotdata1$x)+2.5,ymin = -0.04/scale_factor,ymax = -0.01/scale_factor)+
          scale_x_continuous(limits = c(0,max(plotdata1$x)+2.5),breaks=1:length(numpy),labels=strsplit(as.character(seq),split="")[[1]],expand = c(0,0))+
          scale_fill_manual(breaks = unique(plotdata1$fill),values=unique(plotdata1$fill))+
          coord_cartesian(ylim=c(0,max(plotdata2$y.value)*1.1),expand = c(0,0))+
          theme_bw()+
          theme(line = element_blank(),
                title = element_blank(),
                axis.text.y = element_text(size=8),
                axis.text.x = element_text(size=4,vjust=5),
                legend.position = "none",
                plot.margin = unit(c(0.2,0.2,0.5,0.2),"cm")))
  dev.off()
  
}
ii<-ii+1

if (start %in% plotdata1$x) {
  if (iter!="F12_mut3") {
    sigs<-c(sigs,plotdata1$y[plotdata1$x==start])
  }else{
    sigs<-c(sigs,plotdata1$y[plotdata1$x==start]+plotdata1$y[plotdata1$x==(start-9)])
  }
  iters<-c(iters,t)
}else{
  extra<-c(extra,iter)
}
outputs<-c(outputs,plotdata1)
#sigs2<-c(0.14429542,0.07132571,0.16368895,0.12500914,0.17635284,0.16075970,0.10729650,0.10880718)
#iters2<-c("allele_C_2","allele_T_2","mut1_2","mut2_2","mut3_2","uORF1_2","uORF2_2","uORF3_2")


#sigs<-c(0.14429542,0.07132571,0.16368895,0.12500914,0.17635284,0.16075970,0.10729650,0.10880718)
#iters<-read.csv("FXII_seq.txt",sep = "\t",stringsAsFactors = F,header = F)[,1]
#reals<-c(1,0.58,1.11,0.87,1.31,0.35,0.39,0.29)
res<-data.frame(iters,sigs)
#save(res,file = "res.RData")
#save(outputs,file="Calvo_outputs.RData")
load("res.RData")
load("Calvo_outputs.RData")
res1<-res[1:8,]
res2<-res[9:50,]
res1$sigs<-res1$sigs/res1$sigs[1]
res1<-res1[-1,]
res2_0<-res2[seq(1,nrow(res2),2),]
res2<-res2[seq(2,nrow(res2),2),]
res2$sigs<-res2$sigs/res2_0$sigs
res3<-rbind(res1,res2)
reals<-c(0.58,1.11,0.87,1.31,0.35,0.39,0.29,0.85,0.62,0.59,0.69,0.59,0.44,0.4,0.6,0.39,0.35,0.23,0.33,0.2,0.2,0.59,0.56,0.44,0.27,0.05,0.01,0)
#reals<-c(0.5,1.1,0.87,1.24,0.45,0.44,0.32,0.65,0.62,0.6,0.58,0.51,0.5,0.45,0.42,0.39,0.33,0.26,0.27,0.21,0.2,0.63,0.5,0.42,0.3,0.05,0,0)
res4<-cbind(res3,reals)

mytruncate<-function(x,trun="."){
  return(unlist(lapply(x, function(m){strsplit(m,split=paste0("\\",trun))[[1]][1]})))
}
res4$iters<-mytruncate(res4$iters,"_n|_u$")

R<-cor(x=res4$sigs,y=res4$reals)

graphics.off()
pdf(file = "Calvo_res3.pdf",width = 7,height = 7)
ggplot(data=res4,aes(x=sigs,y=reals))+theme_bw()+
  geom_smooth(formula = "y ~ x",method = lm, se=F,color="blue")+
  geom_smooth(data=res5,formula = "y ~ x",method = lm, se=F,color="red")+
  geom_point()+
  geom_text_repel(aes(label=iters))+
  geom_text(aes(x=1.08,y=0.735,label=paste0("R=",round(R,2))),size=5,hjust=-0.05,vjust=0.5,color="blue")+
  geom_text(aes(x=0.95,y=0.82,label=paste0("R=",round(R3,2))),size=5,hjust=-0.05,vjust=0.5,color="red")+
  coord_fixed(xlim=c(0,max(res4$sigs)),ylim=c(0,max(res4$reals)))+
  labs(x="ILTP Signal Foldchange",y="Normalized Luciferase Activity Foldchange")
dev.off()











library(ggplot2)
library(dplyr)
setwd("D:\\sorf_models\\model\\tests\\expression_profiles_HPA")
data<-read.csv("HPATranscriptTPM.tsv",sep="\t",header = T,stringsAsFactors = F,quote = "")
colon<-rowMeans(data[,names(data)[as.integer(gregexpr(pattern = "colon",names(data)))>=1]])
mytruncate<-function(x,trun="."){
  return(unlist(lapply(x, function(m){strsplit(m,split=paste0("\\",trun))[[1]][1]})))
}
new_data<-data[,1,drop=F]
for (i in unique(mytruncate(names(data),trun="_"))[-1]) {
  new_data<-cbind(new_data,new=rowMeans(data[,names(data)[as.integer(gregexpr(pattern = i,names(data)))>=1]]))
  names(new_data)[ncol(new_data)]<-i
}
data<-new_data
rm(new_data)

meta<-read.csv("Lnc_meta_info.txt",sep = "\t",stringsAsFactors = F,header = F,quote = "")[,c(2,4)]
names(meta)<-c("ids","geneID")
data<-data[data$transcriptid %in% meta$ids,]
maxexp<-apply(data[,2:ncol(data)],1, max)
meanexp<-rowMeans(data[,2:ncol(data)])

if (FALSE) {
  graphics.off()
  pdf(file="lncexp.pdf",width = 10,height = 10)
  ggplot(data=data.frame(maxexp,meanexp))+geom_point(aes(x=log10(maxexp+1),y=log10(meanexp+1)),size=0.1,alpha=0.2)+theme_bw()
  dev.off()
}

#取maxexp前10000的lncRNA转录本：
sum(maxexp>4)
used<-merge(meta,data.frame(ids=data$transcriptid[maxexp>4],maxexp=maxexp[maxexp>4]))
used<-used %>% group_by(geneID) %>% summarise(ids=ids[which.max(maxexp)],maxexp=max(maxexp))

write.table(used$ids,file="max10955lnc.txt",col.names = F,row.names = F,quote = F,sep = "\t")

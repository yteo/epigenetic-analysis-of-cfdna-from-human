# adjust %s for sample name accordingly
multi<-read.table(sprintf("%s\\%s_L1HS_multi_frac.txt",sample,sample,sep="\t"))
count<-read.table(sprintf("%s_L1HS_consensuscoord_count.txt",sample,sample,sep="\t"))
num<-read.table(sprintf("%s_L1HS_consensuscoord_num.txt",sample,sample,sep="\t"))
colnames(count)<-c("Pos","Count")
colnames(num)<-c("Pos","Num")
colnames(multi)<-c("repeat","Pos","multi_frac")


# count = unique count
#num=num of L1HS mapped at the consensus
# multi_frac = multi count fraction
count_num<-merge(count,num,by="Pos")
count_num<-merge(count_num,multi,by="Pos")

# combined mulfi fraction and unique counts
count_num$multi_unique<-count_num$Count+count_num$multi_frac


i<-i+1
lib<-libsize[grep(paste(sample),libsize$V1),]$V2
count_num$norm<-count_num$multi_unique/lib
count_num2<-merge(count_num,counter,all=T,by="Pos")
count_num2[is.na(count_num2)] <- 0
count_num2$FFT = filterFFT(count_num2$norm, pcKeepComp=0.01)


lines(count_num2$Pos,count_num2$FFT,col=paste("deepskyblue",i,sep=""),lty=2,lwd=0.1)

dat<-cbind(dat,count_num2$FFT)

dat$ave<-rowMeans(dat[,2:4])
lines(dat$pos,dat$ave,lty=1,col=paste("deepskyblue4"),lwd=2)
print(cor.test(L1HS_map$V4,dat$ave))
## Figure 2
# compartment values
library(FSA)
library(dunn.test)
library(rcompanion)
library(reshape2)
sample<-"HealthyCen1_full_compartment.bed"
cnt2 <- read.table(sprintf(sample, header=F))
cnt2<-cnt2[,-1]

cnt2<-cnt2[!grepl ("B4",cnt2$V5),]

samples<-c("HealthyCen1_full_compartment.bed","HealthyCen2_full_compartment.bed","HealthyCen3_full_compartment.bed","Old1_full_compartment.bed","Old2_full_compartment.bed","Old3_full_compartment.bed","UnhealthyCen1_full_compartment.bed","UnhealthyCen2_full_compartment.bed","UnhealthyCen3_full_compartment.bed","Young1_full_compartment.bed","Young2_compartment.bed","Young3_full_compartment.bed")

for (sample in samples){
  cnt<-read.table(sprintf(sample, header=F))
  cnt<-cnt[!grepl ("B4",cnt$V5),]
  colnames(cnt)<-c(paste(sample),"chr","start","end","compartment")
  cnt2<-merge(cnt2,cnt,by.x=c("V2","V3","V4","V5"),by.y=c("chr","start","end","compartment"))
  
}


colnames(cnt2)[1:4]<-c("chrom", "start","stop","compartment")
cnt2<-cnt2[complete.cases(cnt2), ]

cnt2<-melt(cnt2, id.vars=c("chrom", "start","stop","compartment"))
cnt2$age<-gsub("_.*","",cnt2$variable)

#  Healthy centenarian
comp_medHC<-cnt2[grep("HealthyCen",cnt2$age),]

comp_medHC<-comp_medHC[complete.cases(comp_medHC),]
comp_medHC_lm = lm(value ~ compartment, data = comp_medHC)
# diagnostic
HC.mod = data.frame(Fitted = fitted(comp_medHC_lm),
                    Residuals = resid(comp_medHC_lm), compartment = comp_medHC$compartment)
ggplot(HC.mod, aes(Fitted, Residuals, colour = compartment)) + geom_point()
# there are a few points with high residuals. try removing them
remrow <- function(x, rows) x[-rows,, drop = FALSE]
HC.mod1<-remrow(HC.mod,(which.max(HC.mod$Residuals)))
HC.mod2<-remrow(HC.mod1,(which.max(HC.mod1$Residuals)))
HC.mod3<-remrow(HC.mod2,(which.max(HC.mod2$Residuals)))
HC.mod4<-remrow(HC.mod3,(which.max(HC.mod3$Residuals)))
HC.mod5<-remrow(HC.mod4,(which.max(HC.mod4$Residuals)))
HC.mod6<-remrow(HC.mod5,(which.max(HC.mod5$Residuals)))
HC.mod7<-remrow(HC.mod6,(which.max(HC.mod6$Residuals)))
HC.mod8<-remrow(HC.mod7,(which.max(HC.mod7$Residuals)))
HC.mod9<-remrow(HC.mod8,(which.max(HC.mod8$Residuals)))
ggplot(HC.mod9, aes(Fitted, Residuals, colour = compartment)) + geom_point()

comp_medHC1<-remrow(comp_medHC,(which.max(HC.mod$Residuals)))
comp_medHC2<-remrow(comp_medHC1,(which.max(HC.mod1$Residuals)))
comp_medHC3<-remrow(comp_medHC2,(which.max(HC.mod2$Residuals)))
comp_medHC4<-remrow(comp_medHC3,(which.max(HC.mod3$Residuals)))
comp_medHC5<-remrow(comp_medHC4,(which.max(HC.mod4$Residuals)))
comp_medHC6<-remrow(comp_medHC5,(which.max(HC.mod5$Residuals)))
comp_medHC7<-remrow(comp_medHC6,(which.max(HC.mod6$Residuals)))
comp_medHC8<-remrow(comp_medHC7,(which.max(HC.mod7$Residuals)))
comp_medHC9<-remrow(comp_medHC8,(which.max(HC.mod8$Residuals)))


# Kruskal wallis test

### Dunn test


kruskal.test(value ~ compartment, data=comp_medHC9)
dunn_HC = dunnTest(value  ~ compartment,
                   data=comp_medHC9,
                   method="bh")


# Unhealthy centenarian
##
comp_medUHC<-cnt2[grep("UnhealthyCen",cnt2$age),]

comp_medUHC<-comp_medUHC[complete.cases(comp_medUHC),]
comp_medUHC_lm = lm(value ~ compartment, data = comp_medUHC)
# diagnostic
UHC.mod = data.frame(Fitted = fitted(comp_medUHC_lm),
                     Residuals = resid(comp_medUHC_lm), compartment = comp_medUHC$compartment)
ggplot(UHC.mod, aes(Fitted, Residuals, colour = compartment)) + geom_point()
# there are a few points with high residuals. try removing them
remrow <- function(x, rows) x[-rows,, drop = FALSE]
UHC.mod1<-remrow(UHC.mod,(which.max(UHC.mod$Residuals)))
UHC.mod2<-remrow(UHC.mod1,(which.max(UHC.mod1$Residuals)))
UHC.mod3<-remrow(UHC.mod2,(which.max(UHC.mod2$Residuals)))
UHC.mod4<-remrow(UHC.mod3,(which.max(UHC.mod3$Residuals)))
UHC.mod5<-remrow(UHC.mod4,(which.max(UHC.mod4$Residuals)))
UHC.mod6<-remrow(UHC.mod5,(which.max(UHC.mod5$Residuals)))
UHC.mod7<-remrow(UHC.mod6,(which.max(UHC.mod6$Residuals)))
UHC.mod8<-remrow(UHC.mod7,(which.max(UHC.mod7$Residuals)))
UHC.mod9<-remrow(UHC.mod8,(which.max(UHC.mod8$Residuals)))
ggplot(UHC.mod9, aes(Fitted, Residuals, colour = compartment)) + geom_point()

comp_medUHC1<-remrow(comp_medUHC,(which.max(UHC.mod$Residuals)))
comp_medUHC2<-remrow(comp_medUHC1,(which.max(UHC.mod1$Residuals)))
comp_medUHC3<-remrow(comp_medUHC2,(which.max(UHC.mod2$Residuals)))
comp_medUHC4<-remrow(comp_medUHC3,(which.max(UHC.mod3$Residuals)))
comp_medUHC5<-remrow(comp_medUHC4,(which.max(UHC.mod4$Residuals)))
comp_medUHC6<-remrow(comp_medUHC5,(which.max(UHC.mod5$Residuals)))
comp_medUHC7<-remrow(comp_medUHC6,(which.max(UHC.mod6$Residuals)))
comp_medUHC8<-remrow(comp_medUHC7,(which.max(UHC.mod7$Residuals)))
comp_medUHC9<-remrow(comp_medUHC8,(which.max(UHC.mod8$Residuals)))

kruskal.test(value ~ compartment, data=comp_medUHC9)
dunn_UHC = dunnTest(value  ~ compartment,
                   data=comp_medUHC9,
                   method="bh")





# Young
###
comp_medYoung<-cnt2[grep("Young",cnt2$age),]

comp_medYoung<-comp_medYoung[complete.cases(comp_medYoung),]
comp_medYoung_lm = lm(value ~ compartment, data = comp_medYoung)
# diagnostic
Young.mod = data.frame(Fitted = fitted(comp_medYoung_lm),
                       Residuals = resid(comp_medYoung_lm), compartment = comp_medYoung$compartment)
ggplot(Young.mod, aes(Fitted, Residuals, colour = compartment)) + geom_point()
# there are a few points with high residuals. try removing them
remrow <- function(x, rows) x[-rows,, drop = FALSE]
Young.mod1<-remrow(Young.mod,(which.max(Young.mod$Residuals)))
Young.mod2<-remrow(Young.mod1,(which.max(Young.mod1$Residuals)))
Young.mod3<-remrow(Young.mod2,(which.max(Young.mod2$Residuals)))
Young.mod4<-remrow(Young.mod3,(which.max(Young.mod3$Residuals)))
Young.mod5<-remrow(Young.mod4,(which.max(Young.mod4$Residuals)))
Young.mod6<-remrow(Young.mod5,(which.max(Young.mod5$Residuals)))
Young.mod7<-remrow(Young.mod6,(which.max(Young.mod6$Residuals)))
Young.mod8<-remrow(Young.mod7,(which.max(Young.mod7$Residuals)))
Young.mod9<-remrow(Young.mod8,(which.max(Young.mod8$Residuals)))
ggplot(Young.mod9, aes(Fitted, Residuals, colour = compartment)) + geom_point()

comp_medYoung1<-remrow(comp_medYoung,(which.max(Young.mod$Residuals)))
comp_medYoung2<-remrow(comp_medYoung1,(which.max(Young.mod1$Residuals)))
comp_medYoung3<-remrow(comp_medYoung2,(which.max(Young.mod2$Residuals)))
comp_medYoung4<-remrow(comp_medYoung3,(which.max(Young.mod3$Residuals)))
comp_medYoung5<-remrow(comp_medYoung4,(which.max(Young.mod4$Residuals)))
comp_medYoung6<-remrow(comp_medYoung5,(which.max(Young.mod5$Residuals)))
comp_medYoung7<-remrow(comp_medYoung6,(which.max(Young.mod6$Residuals)))
comp_medYoung8<-remrow(comp_medYoung7,(which.max(Young.mod7$Residuals)))
comp_medYoung9<-remrow(comp_medYoung8,(which.max(Young.mod8$Residuals)))


kruskal.test(value ~ compartment, data=comp_medYoung9)
dunn_Young = dunnTest(value  ~ compartment,
                    data=comp_medYoung9,
                    method="bh")

# Old
###

comp_medOld<-cnt2[grep("Old",cnt2$age),]

comp_medOld<-comp_medOld[complete.cases(comp_medOld),]
comp_medOld_lm = lm(value ~ compartment, data = comp_medOld)
# diagnostic
Old.mod = data.frame(Fitted = fitted(comp_medOld_lm),
                     Residuals = resid(comp_medOld_lm), compartment = comp_medOld$compartment)
ggplot(Old.mod, aes(Fitted, Residuals, colour = compartment)) + geom_point()
# there are a few points with high residuals. try removing them
remrow <- function(x, rows) x[-rows,, drop = FALSE]
Old.mod1<-remrow(Old.mod,(which.max(Old.mod$Residuals)))
Old.mod2<-remrow(Old.mod1,(which.max(Old.mod1$Residuals)))
Old.mod3<-remrow(Old.mod2,(which.max(Old.mod2$Residuals)))
Old.mod4<-remrow(Old.mod3,(which.max(Old.mod3$Residuals)))
Old.mod5<-remrow(Old.mod4,(which.max(Old.mod4$Residuals)))
Old.mod6<-remrow(Old.mod5,(which.max(Old.mod5$Residuals)))
Old.mod7<-remrow(Old.mod6,(which.max(Old.mod6$Residuals)))
Old.mod8<-remrow(Old.mod7,(which.max(Old.mod7$Residuals)))
Old.mod9<-remrow(Old.mod8,(which.max(Old.mod8$Residuals)))
ggplot(Old.mod9, aes(Fitted, Residuals, colour = compartment)) + geom_point()

comp_medOld1<-remrow(comp_medOld,(which.max(Old.mod$Residuals)))
comp_medOld2<-remrow(comp_medOld1,(which.max(Old.mod1$Residuals)))
comp_medOld3<-remrow(comp_medOld2,(which.max(Old.mod2$Residuals)))
comp_medOld4<-remrow(comp_medOld3,(which.max(Old.mod3$Residuals)))
comp_medOld5<-remrow(comp_medOld4,(which.max(Old.mod4$Residuals)))
comp_medOld6<-remrow(comp_medOld5,(which.max(Old.mod5$Residuals)))
comp_medOld7<-remrow(comp_medOld6,(which.max(Old.mod6$Residuals)))
comp_medOld8<-remrow(comp_medOld7,(which.max(Old.mod7$Residuals)))
comp_medOld9<-remrow(comp_medOld8,(which.max(Old.mod8$Residuals)))


kruskal.test(value ~ compartment, data=comp_medOld9)
dunn_Old = dunnTest(value  ~ compartment,
                      data=comp_medOld9,
                      method="bh")




allsignal_removedout<-data.frame(compartment=comp_medYoung9$compartment, Young=comp_medYoung9$value,Old=comp_medOld9$value,HealthyCen=comp_medHC9$value,UnhealthyCen=comp_medUHC9$value)
allsignal_removedout<-melt(allsignal_removedout,id.vars = c("compartment"))

ggplot()+facet_wrap(~variable,ncol = 1)+geom_boxplot(data=allsignal_removedout,aes(x=compartment, y=value),outlier.shape=NA,notch=T)+theme_bw()+
  scale_y_continuous(limits = c(1000,3000))+ylab("Average cfDNA signals in 100kb spanning regions")+xlab("Subcompartment")




# variance of all compartment signals

allsignal_removedout<-data.frame(compartment=comp_medYoung9$compartment, Young=comp_medYoung9$value,Old=comp_medOld9$value,HealthyCen=comp_medHC9$value,UnhealthyCen=comp_medUHC9$value)
allsignal_removedout<-allsignal_removedout[complete.cases(allsignal_removedout), ]
variance<-apply(allsignal_removedout[,-1], 2, var)
allsignal_removedout<-melt(allsignal_removedout,id.vars = c("compartment"))



library(PairedData)
library(car)
# levene's test for variance 
allsignal_removedout2<-allsignal_removedout[,-1]
leveneTest(value ~ variable, data = allsignal_removedout)

levene.Var.test(allsignal_removedout2$value[allsignal_removedout2$variable=="Young"], allsignal_removedout2$value[allsignal_removedout2$variable=="Old"],location=c("median")) 
levene.Var.test(allsignal_removedout2$value[allsignal_removedout2$variable=="Young"], allsignal_removedout2$value[allsignal_removedout2$variable=="HealthyCen"],location=c("median")) 
levene.Var.test(allsignal_removedout2$value[allsignal_removedout2$variable=="Young"], allsignal_removedout2$value[allsignal_removedout2$variable=="UnhealthyCen"],location=c("median")) 


levene.Var.test(allsignal_removedout2$value[allsignal_removedout2$variable=="HealthyCen"], allsignal_removedout2$value[allsignal_removedout2$variable=="UnhealthyCen"],location=c("median")) 
levene.Var.test(allsignal_removedout2$value[allsignal_removedout2$variable=="HealthyCen"], allsignal_removedout2$value[allsignal_removedout2$variable=="Old"],location=c("median")) 
levene.Var.test(allsignal_removedout2$value[allsignal_removedout2$variable=="UnhealthyCen"], allsignal_removedout2$value[allsignal_removedout2$variable=="Old"],location=c("median")) 

p.adjust(c("2.2e-16","2.2e-16","2.2e-16","2.2e-16","2.2e-16","1.061e-05"))
#######################

# chr 11 ideogram



# variance of all compartment signals

allsignal_removedout<-data.frame(chrom=comp_medYoung9$chrom,start=comp_medYoung9$start,stop=comp_medYoung9$stop,compartment=comp_medYoung9$compartment, Young=comp_medYoung9$value,Old=comp_medOld9$value,HealthyCen=comp_medHC9$value,UnhealthyCen=comp_medUHC9$value)

chr11<-allsignal_removedout[grep("11",allsignal_removedout$chrom),]
chr11<-melt(chr11,id.vars = c("chrom","start","stop","compartment"))
library(dplyr)
chr11<-data.frame(chr11 %>% group_by(chrom,start,stop,compartment,variable) %>%
  summarise(mean=mean(value)))


ggplot(chr11, aes(start, mean,color=compartment, alpha=0.5)) + 
  geom_point(size=0.5) + xlab("") + ylab("") +facet_wrap(~variable,ncol=1)+
  
  scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black")) +theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ scale_x_continuous(labels = comma)+ylab("Average cfDNA signals")+
  scale_y_continuous(limit=c(1200,2500))




# 3D
chr11_wide <- spread(chr11, variable, mean)


pos<-read.table("Data\\position_chr11_100kb.txt")
dat<-read.table("Data\\chr11_100kb_3Dcoor.xyz",sep="\t")
dat<-dat[-1317,]
dat<-dat[,-1]
dat$pos<-pos$V1
chr11_wide<-merge(chr11_wide,dat,by.x="start",by.y="pos")


library("plot3D")
meansignalster3D(chr11_wide$V2,chr11_wide$V3,chr11_wide$V4,lwd=10,colvar = as.numeric(as.character(chr11_wide$Young)),type="l",main = "Young",colkey = list(side = 1),bty="n",clim=c(600,2400))
meansignalster3D(chr11_wide$V2,chr11_wide$V3,chr11_wide$V4,lwd=10,colvar = as.numeric(as.character(chr11_wide$Old)),type="l",main = "Old",colkey = list(side = 1),bty="n",clim=c(640,2400))
meansignalster3D(chr11_wide$V2,chr11_wide$V3,chr11_wide$V4,lwd=10,colvar = as.numeric(as.character(chr11_wide$HealthyCen)),type="l",main = "Healthy centenarian",colkey = list(side = 1),bty="n",clim=c(640,2400))
meansignalster3D(chr11_wide$V2,chr11_wide$V3,chr11_wide$V4,lwd=10,colvar = as.numeric(as.character(chr11_wide$UnhealthyCen)),type="l",main = "Unhealthy centenarian",colkey = list(side = 1),bty="n",clim=c(640,2400))





allsignal_removedout<-data.frame(chrom=comp_medYoung9$chrom,start=comp_medYoung9$start,stop=comp_medYoung9$stop,compartment=comp_medYoung9$compartment, Young=comp_medYoung9$value,Old=comp_medOld9$value,HealthyCen=comp_medHC9$value,UnhealthyCen=comp_medUHC9$value)
meansignals<-data.frame(allsignal_removedout %>% group_by(chrom,start,stop,compartment) %>%
                    summarise(meanYoung=mean(Young),meanOld=mean(Old),meanHC=mean(HealthyCen),meanUHC=mean(UnhealthyCen)))
## meansignalsterplot of average signals

smoothScatter(meansignals$meanOld,meansignals$meanYoung,xlim=c(640,2500),ylim=c(640,2500),xlab=("cfDNA signals in 100kb regions in old"),ylab=("cfDNA signals in 100kb regions in young"))
smoothScatter(meansignals$meanUHC,meansignals$meanYoung,xlim=c(640,2500),ylim=c(640,2500),xlab=("cfDNA signals in 100kb regions in unhealthy centenarian"),ylab=("cfDNA signals in 100kb regions in young"))
smoothScatter(meansignals$meanHC,meansignals$meanYoung,xlim=c(640,2500),ylim=c(640,2500),xlab=("cfDNA signals in 100kb regions in healthy centenarian"),ylab=("cfDNA signals in 100kb regions in young"))
smoothScatter(meansignals$meanHC,meansignals$meanOld,xlim=c(640,2500),ylim=c(640,2500),xlab=("cfDNA signals in 100kb regions in healthy centenarian"),ylab=("cfDNA signals in 100kb regions in old"))
smoothScatter(meansignals$meanUHC,meansignals$meanOld,xlim=c(640,2500),ylim=c(640,2500),xlab=("cfDNA signals in 100kb regions in unhealthy centenarian"),ylab=("cfDNA signals in 100kb regions in old"))
smoothScatter(meansignals$meanUHC,meansignals$meanHC,xlim=c(640,2500),ylim=c(640,2500),xlab=("cfDNA signals in 100kb regions in unhealthy centenarian"),ylab=("cfDNA signals in 100kb regions in healthy centenarian"))


cor.test(meansignals$meanOld,meansignals$meanYoung)
cor.test(meansignals$meanUHC,meansignals$meanYoung)
cor.test(meansignals$meanHC,meansignals$meanYoung)
cor.test(meansignals$meanHC,meansignals$meanOld)
cor.test(meansignals$meanUHC,meansignals$meanOld)
cor.test(meansignals$meanUHC,meansignals$meanHC)


## Figure 3
##################################################################################
# Fold change between conditions
library("DESeq2")
library("plot3D")

sample<-"HealthyCen1_full_compartment.bed"
cnt2 <- read.table(sprintf(sample, header=F))
cnt2<-cnt2[,-1]

  
samples<-c("HealthyCen1_full_compartment.bed","HealthyCen2_full_compartment.bed","HealthyCen3_full_compartment.bed","Old1_full_compartment.bed","Old2_full_compartment.bed","Old3_full_compartment.bed","UnhealthyCen1_full_compartment.bed","UnhealthyCen2_full_compartment.bed","UnhealthyCen3_full_compartment.bed","Young1_full_compartment.bed","Young2_full_compartment.bed","Young3_full_compartment.bed")
  
for (sample in samples){
  cnt<-read.table(sprintf(sample, header=F))
  colnames(cnt)<-c(paste(sample),"chr","start","end","compartment")
  cnt2<-merge(cnt2,cnt,by.x=c("V2","V3","V4","V5"),by.y=c("chr","start","end","compartment"))
  
}



rownames(cnt2)<-paste(cnt2$V5,cnt2$V2,cnt2$V3,sep="_")



cnt2<-cnt2[,-1]
cnt2<-cnt2[,-1]
cnt2<-cnt2[,-1]
cnt2<-cnt2[,-1]
cnt2<-cnt2[apply(cnt2, 1, function(x) !all(x==0)),]



condition=gsub("_.*","",colnames(cnt2))

cnt2<-cnt2[!grepl ("B4",rownames(cnt2)),]

(coldata <- data.frame(row.names=colnames(cnt2), condition=gsub("_.*","",colnames(cnt2))))
  countdata<-cnt2
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds
# pre-filtering
dds <- dds[ rowSums(counts(dds)) > 1, ]
# Run the DESeq pipeline
dds <- DESeq(dds, test="LRT", reduced= ~ 1)

# Plot dispersions

plotDispEsts(dds, main="Dispersion plot")


########## To test whether something is changed due to condition in general
#Likelihood ratio test ~ One-way anova test
ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
resLRT <- results(ddsLRT)
resLRT
resLRT1<-data.frame(resLRT)
sigresLRT<-resLRT1[which(resLRT1$padj<0.05),]
########################################

########################################
# individual comparison
dds<-DESeq(dds)
res <- results(dds)
sizeFactors(dds)
res

res_healthycen.young<-data.frame(results(ddsLRT, contrast=c("condition","HealthyCen","Young"),test="Wald"))
res_unhealthycen.young<-data.frame(results(ddsLRT, contrast=c("condition","UnhealthyCen","Young"),test="Wald"))
res_old.young<-data.frame(results(ddsLRT, contrast=c("condition","Old","Young"),test="Wald"))

res_unhealthycen.healthycen<-data.frame(results(ddsLRT, contrast=c("condition","UnhealthyCen","HealthyCen"),test="Wald"))
res_unhealthycen.old<-data.frame(results(ddsLRT, contrast=c("condition","UnhealthyCen","Old"),test="Wald"))
res_old.healthycen<-data.frame(results(ddsLRT, contrast=c("condition","Old","HealthyCen"),test="Wald"))



write.table(cbind(rownames(res_healthycen.young),res_healthycen.young$log2FoldChange),file="logFChealthycen.young_DESEQ2.rnk",sep="\t",quote=F,row.names = F,col.names = F)
write.table(cbind(rownames(res_unhealthycen.young),res_unhealthycen.young$log2FoldChange),file="logFCunhealthycen.young_DESEQ2.rnk",sep="\t",quote=F,row.names = F,col.names = F)
write.table(cbind(rownames(res_old.young),res_old.young$log2FoldChange),file="logFCold.young_DESEQ2.rnk",sep="\t",quote=F,row.names = F,col.names = F)
write.table(cbind(rownames(res_old.healthycen),res_old.healthycen$log2FoldChange),file="logFCold.healthycen_DESEQ2.rnk",sep="\t",quote=F,row.names = F,col.names = F)
write.table(cbind(rownames(res_unhealthycen.old),res_unhealthycen.old$log2FoldChange),file="logFCunhealthycen.old_DESEQ2.rnk",sep="\t",quote=F,row.names = F,col.names = F)
write.table(cbind(rownames(res_unhealthycen.healthycen),res_unhealthycen.healthycen$log2FoldChange),file="logFCunhealthycen.healthycen_DESEQ2.rnk",sep="\t",quote=F,row.names = F,col.names = F)





##############################################################################
# 3D
# stitching differential signals together to color 3D structure

peak.meds <- read.table("logFChealthycen.young_DESEQ2.rnk", header=F,row.names=1)
peak.meds1<-data.frame(compartment=gsub("_.*","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(peak.meds))),chrom=gsub(".*_","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(peak.meds))),start=gsub(".*_","\\1", rownames(peak.meds)),HealthyCen_Young=peak.meds$V2)
peak.meds1<-peak.meds1[peak.meds1$chrom=="11",]






peak.meds <- read.table("logFCold.young_DESEQ2.rnk", header=F,row.names=1)
peak.meds2<-data.frame(compartment=gsub("_.*","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(peak.meds))),chrom=gsub(".*_","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(peak.meds))),start=gsub(".*_","\\1", rownames(peak.meds)),Old_Young=peak.meds$V2)
peak.meds2<-peak.meds2[peak.meds2$chrom=="11",]



peak.meds <- read.table("logFCunhealthycen.young_DESEQ2.rnk", header=F,row.names=1)
peak.meds3<-data.frame(compartment=gsub("_.*","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(peak.meds))),chrom=gsub(".*_","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(peak.meds))),start=gsub(".*_","\\1", rownames(peak.meds)),UnhealthyCen_Young=peak.meds$V2)
peak.meds3<-peak.meds3[peak.meds3$chrom=="11",]





peak.meds <- read.table("logFCold.healthycen_DESEQ2.rnk", header=F,row.names=1)
peak.meds4<-data.frame(compartment=gsub("_.*","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(peak.meds))),chrom=gsub(".*_","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(peak.meds))),start=gsub(".*_","\\1", rownames(peak.meds)),HealthyCen_Old=peak.meds$V2)
peak.meds4<-peak.meds4[peak.meds4$chrom=="11",]





peak.meds <- read.table("logFCunhealthycen.old_DESEQ2.rnk", header=F,row.names=1)
peak.meds5<-data.frame(compartment=gsub("_.*","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(peak.meds))),chrom=gsub(".*_","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(peak.meds))),start=gsub(".*_","\\1", rownames(peak.meds)),UnhealthyCen_Old=peak.meds$V2)
peak.meds5<-peak.meds5[peak.meds5$chrom=="11",]




peak.meds <- read.table("logFCunhealthycen.healthycen_DESEQ2.rnk", header=F,row.names=1)
peak.meds6<-data.frame(compartment=gsub("_.*","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(peak.meds))),chrom=gsub(".*_","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(peak.meds))),start=gsub(".*_","\\1", rownames(peak.meds)),UnhealthyCen_HealthyCen=peak.meds$V2)
peak.meds6<-peak.meds6[peak.meds6$chrom=="11",]




all<-cbind(as.character(peak.meds1$chrom),as.character(peak.meds1$start),as.numeric(as.character(peak.meds1$stop)),as.character(peak.meds1$compartment),as.numeric(peak.meds1$HealthyCen_Young),as.numeric(peak.meds2$Old_Young),as.numeric(peak.meds3$UnhealthyCen_Young),as.numeric(peak.meds4$HealthyCen_Old),as.numeric(peak.meds5$UnhealthyCen_Old),as.numeric(peak.meds6$UnhealthyCen_HealthyCen))

colnames(all)<-c("chrom","start","compartment","HealthyCen_Young","Old_Young","UnhealthyCen_Young","Old_HealthyCen","UnhealthyCen_Old","UnhealthyCen_HealthyCen")
all<-data.frame(all)
#write.table(all,file="diffcfDNA_100kb.bed",quote=F,sep="\t",row.names=F)
all<-all[order(all$start),]
all$start<-as.numeric(as.character(all$start))
pos<-read.table("Z:\\data\\yteo\\BGI_cfDNA_cleanrun\\HiC\\scripts\\position_chr11_100kb.txt")

all2<-merge(pos,all,by.x="V1",by.y="start")

colnames(all2)[1]<-"start"
#write.table(all2,file="diffcfDNA_100kb.bed",quote=F,sep="\t",row.names=F)



dat<-read.table("Data\\chr11_100kb_3Dcoor.xyz",sep="\t")
## caution~ our rank file is missing chr 11 54600000 --it is at linek514 of dat file
dat<-dat[,-1]
dat<-dat[-1317,]
dat<-dat[-514,]

colnames(dat)<-c("X","Y","Z")



scatter3D(dat$X,dat$Y,dat$Z,lwd=5,colvar = as.numeric(as.character(all2$HealthyCen_Young)),type="l",main = "Healthy centenarian - Young",colkey = list(side = 1),bty="n",clim=c(-0.3,0.3),col =ramp.col (col = c("yellow", "blue"), n = 100, alpha = 1))
scatter3D(dat$X,dat$Y,dat$Z,lwd=5,colvar = as.numeric(as.character(all2$UnhealthyCen_Young)),type="l",main = "Unhealthy centenarian - Young",colkey = list(side = 1),bty="n",clim=c(-0.3,0.3),col =ramp.col (col = c("yellow", "blue"), n = 100, alpha = 1))
scatter3D(dat$X,dat$Y,dat$Z,lwd=5,colvar = as.numeric(as.character(all2$Old_Young)),type="l",main = "Old - Young",colkey = list(side = 1),bty="n",clim=c(-0.3,0.3),col =ramp.col (col = c("yellow", "blue"), n = 100, alpha = 1))

scatter3D(dat$X,dat$Y,dat$Z,lwd=5,colvar = as.integer(all2$compartment),type="l",col = c("darkred", "red", "mediumblue", "royalblue1", "skyblue3"), colkey = list(at = c(1,2, 3, 4,5),side = 1, labels = c("A1", "A2", "B1","B2","B3")),bty="n")



scatter3D(dat$X,dat$Y,dat$Z,lwd=5,colvar = as.numeric(as.character(all2$Old_HealthyCen)),type="l",main = "Old - Healthy centenarian",colkey = list(side = 1),bty="n",clim=c(-0.2,0.2),col =ramp.col (col = c("yellow", "blue"), n = 100, alpha = 1))
scatter3D(dat$X,dat$Y,dat$Z,lwd=5,colvar = as.numeric(as.character(all2$UnhealthyCen_HealthyCen)),type="l",main = "Unhealthy centenarian - Healthy centenarian",colkey = list(side = 1),bty="n",clim=c(-0.2,0.2),col =ramp.col (col = c("yellow", "blue"), n = 100, alpha = 1))
scatter3D(dat$X,dat$Y,dat$Z,lwd=5,colvar = as.numeric(as.character(all2$UnhealthyCen_Old)),type="l",main = "Unhealthy centenarian - Old",colkey = list(side = 1),bty="n",clim=c(-0.2,0.2),col =ramp.col (col = c("yellow", "blue"), n = 100, alpha = 1))



######## plotting ideograms 

library(ggplot2)
library(scales)
peak.meds <- read.table("logFCold.healthycen_DESEQ2.rnk", header=F)
peak.meds<-data.frame(compartment=gsub("_.*","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(x))),chrom=gsub(".*_","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(x))),start=gsub(".*_","\\1", rownames(x)),med=peak.meds$V2)



mychrom<-"11"
peak.meds.sub <- subset(peak.meds, peak.meds$chrom == mychrom)

d<-ggplot(peak.meds.sub ,aes(as.numeric(as.character(peak.meds.sub$start)),med,color=compartment,alpha=0.5))+geom_point(size=0.5)+ scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  ggtitle("Old - Healthy centenarian",paste("Chr",mychrom,sep=""))+theme_bw()+ scale_x_continuous(labels = scales::comma)+ ylab("Log fold change of cfDNA signals") + xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=0,linetype="dashed")+scale_y_continuous(limits=c(-0.15,0.15))

############

peak.meds <- read.table("logFCunhealthycen.healthycen_DESEQ2.rnk", header=F)
peak.meds<-data.frame(compartment=gsub("_.*","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(x))),chrom=gsub(".*_","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(x))),start=gsub(".*_","\\1", rownames(x)),med=peak.meds$V2)



mychrom<-"11"
peak.meds.sub <- subset(peak.meds, peak.meds$chrom == mychrom)

e<-ggplot(peak.meds.sub ,aes(as.numeric(as.character(peak.meds.sub$start)),med,color=compartment,alpha=0.5))+geom_point(size=0.5)+ scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  ggtitle("Unhealthy centenarian -  Healthy centenarian",paste("Chr",mychrom,sep=""))+theme_bw()+ scale_x_continuous(labels = scales::comma)+ ylab("Log fold change of cfDNA signals") + xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=0,linetype="dashed")+scale_y_continuous(limits=c(-0.15,0.15))

############
peak.meds <- read.table("logFCunhealthycen.old_DESEQ2.rnk", header=F)
peak.meds<-data.frame(compartment=gsub("_.*","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(x))),chrom=gsub(".*_","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(x))),start=gsub(".*_","\\1", rownames(x)),med=peak.meds$V2)



mychrom<-"11"
peak.meds.sub <- subset(peak.meds, peak.meds$chrom == mychrom)

f<-ggplot(peak.meds.sub ,aes(as.numeric(as.character(peak.meds.sub$start)),med,color=compartment,alpha=0.5))+geom_point(size=0.5)+ scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  ggtitle("Unhealthy centenarian -  Old ",paste("Chr",mychrom,sep=""))+theme_bw()+ scale_x_continuous(labels = scales::comma)+ ylab("Log fold change of cfDNA signals") + xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=0,linetype="dashed")+scale_y_continuous(limits=c(-0.15,0.15))

#######
peak.meds <- read.table("logFCold.young_DESEQ2.rnk", header=F)
peak.meds<-data.frame(compartment=gsub("_.*","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(x))),chrom=gsub(".*_","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(x))),start=gsub(".*_","\\1", rownames(x)),med=peak.meds$V2)



mychrom<-"11"
peak.meds.sub <- subset(peak.meds, peak.meds$chrom == mychrom)

a<-ggplot(peak.meds.sub ,aes(as.numeric(as.character(peak.meds.sub$start)),med,color=compartment,alpha=0.5))+geom_point(size=0.5)+ scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  ggtitle("Old - Young",paste("Chr",mychrom,sep=""))+theme_bw()+ scale_x_continuous(labels = scales::comma)+ ylab("Log fold change of cfDNA signals") + xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=0,linetype="dashed")+scale_y_continuous(limits=c(-0.15,0.15))

##########
peak.meds <- read.table("logFCunhealthycen.young_DESEQ2.rnk", header=F)
peak.meds<-data.frame(compartment=gsub("_.*","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(x))),chrom=gsub(".*_","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(x))),start=gsub(".*_","\\1", rownames(x)),med=peak.meds$V2)



mychrom<-"11"
peak.meds.sub <- subset(peak.meds, peak.meds$chrom == mychrom)

b<-ggplot(peak.meds.sub ,aes(as.numeric(as.character(peak.meds.sub$start)),med,color=compartment,alpha=0.5))+geom_point(size=0.5)+ scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  ggtitle("Unhealthy centenarian - Young",paste("Chr",mychrom,sep=""))+theme_bw()+ scale_x_continuous(labels = scales::comma)+ ylab("Log fold change of cfDNA signals") + xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=0,linetype="dashed")+scale_y_continuous(limits=c(-0.15,0.15))

#########
peak.meds <- read.table("logFChealthycen.young_DESEQ2.rnk", header=F)
peak.meds<-data.frame(compartment=gsub("_.*","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(x))),chrom=gsub(".*_","\\1",gsub(".*_[a-z]|_\\d+$", "", rownames(x))),start=gsub(".*_","\\1", rownames(x)),med=peak.meds$V2)



mychrom<-"11"
peak.meds.sub <- subset(peak.meds, peak.meds$chrom == mychrom)

c<-ggplot(peak.meds.sub ,aes(as.numeric(as.character(peak.meds.sub$start)),med,color=compartment,alpha=0.5))+geom_point(size=0.5)+ scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  ggtitle("Healthy centenarian - Young",paste("Chr",mychrom,sep=""))+theme_bw()+ scale_x_continuous(labels = scales::comma)+ ylab("Log fold change of cfDNA signals") + xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=0,linetype="dashed")+scale_y_continuous(limits=c(-0.15,0.15))


multiplot(a,b,c,ncol=1)
multiplot(d,e,f,ncol=1)



# Distribution plot all chromosomes

y<-read.table("logFChealthycen.young_DESEQ2.rnk",sep="\t",header=F,row.names=1)
y$compartment<-gsub("_.*","",rownames(y))
y<-y[(y$compartment!="NA"),]
b<-ggplot(data=y,aes(x=V2,color=compartment)) + geom_density()+theme_bw()+scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  scale_x_continuous(limits = c(-0.35,0.25))+xlab("Log fold change")+ylab("Density")+ggtitle("Healthy centenarian - Young")

y<-read.table("logFCold.young_DESEQ2.rnk",sep="\t",header=F,row.names=1)
y$compartment<-gsub("_.*","",rownames(y))
y<-y[(y$compartment!="NA"),]
a<-ggplot(data=y,aes(x=V2,color=compartment)) + geom_density()+theme_bw()+scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  scale_x_continuous(limits = c(-0.35,0.25))+xlab("Log fold change")+ylab("Density")+ggtitle("Old - Young")


y<-read.table("logFCunhealthycen.young_DESEQ2.rnk",sep="\t",header=F,row.names=1)
y$compartment<-gsub("_.*","",rownames(y))
y<-y[(y$compartment!="NA"),]
c<-ggplot(data=y,aes(x=V2,color=compartment)) + geom_density()+theme_bw()+scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  scale_x_continuous(limits = c(-0.35,0.25))+xlab("Log fold change")+ylab("Density")+ggtitle("Unhealthy centenarian - Young")

multiplot(a,b,c,ncol=1)





y<-read.table("logFCold.healthycen_DESEQ2.rnk",sep="\t",header=F,row.names=1)
y$compartment<-gsub("_.*","",rownames(y))
y<-y[(y$compartment!="NA"),]
d<-ggplot(data=y,aes(x=V2,color=compartment)) + geom_density()+theme_bw()+scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  scale_x_continuous(limits = c(-0.35,0.25))+xlab("Log fold change")+ylab("Density")+ggtitle("Old - Healthy centenarian")

y<-read.table("logFCunhealthycen.healthycen_DESEQ2.rnk",sep="\t",header=F,row.names=1)
y$compartment<-gsub("_.*","",rownames(y))
y<-y[(y$compartment!="NA"),]
e<-ggplot(data=y,aes(x=V2,color=compartment)) + geom_density()+theme_bw()+scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  scale_x_continuous(limits = c(-0.35,0.25))+xlab("Log fold change")+ylab("Density")+ggtitle("Unhealthy centenarian - Healthy centenarian")


y<-read.table("logFCunhealthycen.old_DESEQ2.rnk",sep="\t",header=F,row.names=1)
y$compartment<-gsub("_.*","",rownames(y))
y<-y[(y$compartment!="NA"),]
f<-ggplot(data=y,aes(x=V2,color=compartment)) + geom_density()+theme_bw()+scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  scale_x_continuous(limits = c(-0.35,0.25))+xlab("Log fold change")+ylab("Density")+ggtitle("Unhealthy centenarian - Old")

multiplot(d,e,f,ncol=1)





# Distribution plot chromosome 11

y<-read.table("logFChealthycen.young_DESEQ2.rnk",sep="\t",header=F,row.names=1)
y<-y[grep("_11_",rownames(y)),,drop=F]
y$compartment<-gsub("_.*","",rownames(y))
y<-y[(y$compartment!="NA"),]
b<-ggplot(data=y,aes(x=V2,color=compartment)) + geom_density()+theme_bw()+scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  scale_x_continuous(limits = c(-0.35,0.25))+xlab("Log fold change")+ylab("Density")+ggtitle("Healthy centenarian - Young")

y<-read.table("logFCold.young_DESEQ2.rnk",sep="\t",header=F,row.names=1)
y<-y[grep("_11_",rownames(y)),,drop=F]
y$compartment<-gsub("_.*","",rownames(y))
y<-y[(y$compartment!="NA"),]
a<-ggplot(data=y,aes(x=V2,color=compartment)) + geom_density()+theme_bw()+scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  scale_x_continuous(limits = c(-0.35,0.25))+xlab("Log fold change")+ylab("Density")+ggtitle("Old - Young")


y<-read.table("logFCunhealthycen.young_DESEQ2.rnk",sep="\t",header=F,row.names=1)
y<-y[grep("_11_",rownames(y)),,drop=F]
y$compartment<-gsub("_.*","",rownames(y))
y<-y[(y$compartment!="NA"),]
c<-ggplot(data=y,aes(x=V2,color=compartment)) + geom_density()+theme_bw()+scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  scale_x_continuous(limits = c(-0.35,0.25))+xlab("Log fold change")+ylab("Density")+ggtitle("Unhealthy centenarian - Young")

multiplot(a,b,c,ncol=1)





y<-read.table("logFCold.healthycen_DESEQ2.rnk",sep="\t",header=F,row.names=1)
y<-y[grep("_11_",rownames(y)),,drop=F]
y$compartment<-gsub("_.*","",rownames(y))
y<-y[(y$compartment!="NA"),]
d<-ggplot(data=y,aes(x=V2,color=compartment)) + geom_density()+theme_bw()+scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  scale_x_continuous(limits = c(-0.35,0.25))+xlab("Log fold change")+ylab("Density")+ggtitle("Old - Healthy centenarian")

y<-read.table("logFCunhealthycen.healthycen_DESEQ2.rnk",sep="\t",header=F,row.names=1)
y<-y[grep("_11_",rownames(y)),,drop=F]
y$compartment<-gsub("_.*","",rownames(y))
y<-y[(y$compartment!="NA"),]
e<-ggplot(data=y,aes(x=V2,color=compartment)) + geom_density()+theme_bw()+scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  scale_x_continuous(limits = c(-0.35,0.25))+xlab("Log fold change")+ylab("Density")+ggtitle("Unhealthy centenarian - Healthy centenarian")


y<-read.table("logFCunhealthycen.old_DESEQ2.rnk",sep="\t",header=F,row.names=1)
y<-y[grep("_11_",rownames(y)),,drop=F]
y$compartment<-gsub("_.*","",rownames(y))
y<-y[(y$compartment!="NA"),]
f<-ggplot(data=y,aes(x=V2,color=compartment)) + geom_density()+theme_bw()+scale_color_manual(values=c("darkred", "red", "mediumblue", "royalblue1", "skyblue3","black"))+
  scale_x_continuous(limits = c(-0.35,0.25))+xlab("Log fold change")+ylab("Density")+ggtitle("Unhealthy centenarian - Old")

multiplot(d,e,f,ncol=1)


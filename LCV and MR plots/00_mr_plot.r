R
library(data.table)
library(plotrix)
library(scales)

caro_orange <- rgb(249,164,24,maxColorValue=255)
beauty_red <- rgb(237,33,36,maxColorValue=255)
adam_blue <- rgb(41,153,208,maxColorValue=255)


d1 <- fread('PGC.PTSD_MRlap_results_adam.csv',data.table=F)
d1$color <-'black'
d1[d1$Type=="MR",]$color <- adam_blue
d1[d1$Type=="LAP",]$color <- caro_orange




d1$mr <- d1$IVW_beta
d1$LCImr <- d1$IVW_beta - 1.96*d1$IVW_se
d1$UCImr <- d1$IVW_beta + 1.96*d1$IVW_se


cis3 <- subset(d1,Exposure_Name == "PTSD") #axis gets flipped if changing to horizontal..

pdf('PTSD_MR.pdf',7,5)
par(mar=c(5, 14, 4, 6) + 0.5)
plotCI(y=cis3$Placement,x=cis3$mr,li=cis3$LCImr,ui=cis3$UCImr,err="x",lwd=2,pch=19,cex.axis=1.25,xlab="MR effect estimate",main="",cex.lab=1.45,col=cis3$color,scol=alpha(cis3$color,.4),sfrac=0,xaxt='n',yaxt='n',cex=1.8,ylab="",xlim=c(-0.5,0.5))

abline(v=0,col='black',lwd=2,lty=2)

axis(1,at=c(-1,-0.5,0,0.5,1),cex.axis=1.25)
axis(2,at=c(1:dim(cis3)[1]), labels=cis3$Outcome_name,cex.axis=1,las=1)

dev.off()



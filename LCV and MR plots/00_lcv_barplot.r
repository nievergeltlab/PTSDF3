R
library(data.table)
library(plotrix)
library(scales)


#this data isa lready filtered to bonferroni significant. I am missing 3 due to not having all from renato:
#Adult education class
#Duration to first press of snap-button in each round
#Volume level set by participant (right)

res2 <- fread('rg_lcv_pan.ukb.adam_v2.csv',data.table=F)

#815 traist in the database... nee
res2 <- subset(res2, h2_z >= 6 | h2_z <=-6)

#Where are highest correlations
res2 <- res2[order(res2$rg_w_ptsd),]

bon1 <- 818

#818 have significant rg, so we do LCV on all tehse

#Count number of phenotypes with rg significant
res2 <- subset(res2,pval.gcpzero.2tailed <= 0.05/818)


res2$OR <- res2$gcp.pm
res2$LCI <- res2$gcp.pm - 1.96*res2$gcp.pse
res2$UCI <- res2$gcp.pm + 1.96*res2$gcp.pse

res2$rg <- res2$rg_w_ptsd
res2$LCIrg <- res2$rg_w_ptsd - 1.96*res2$rg_se
res2$UCIrg <- res2$rg_w_ptsd  + 1.96*res2$rg_se

res2$color <- NA
caro_orange <- rgb(249,164,24,maxColorValue=255)
beauty_red <- rgb(237,33,36,maxColorValue=255)
adam_blue <- rgb(41,153,208,maxColorValue=255)


#this will be fore the second data

#Nagelkerke is 0.0096 on the observed scale - 0.0157 assuming the sample prevalence matches the population prevalence.

cis2 <- cis2[nrow(cis2):1,]

pdf('PTSD_phecoded2.pdf',7,7)
par(mar=c(5, 12, 4, 6) + 0.5)
plotCI(y=1:dim(cis2)[1],x=cis2$OR,li=cis2$LCI,ui=cis2$UCI,err="x",lwd=2,pch=19,cex.axis=1.25,xlab="Genetic causality proportion (GCP)",main="",ylab='',cex.lab=1.45,col=adam_blue,scol=alpha(adam_blue,.4),sfrac=0,xaxt='n',yaxt='n',cex=1.8)

#legend("topleft",col=c("white",beauty_red,caro_orange,adam_blue),legend=c("Training -> Target", "MVP -> UKB", "UKB -> MVP", "Meta-analysis -> MVP R4"),bty="n",pch=19,cex=1.5)



abline(v=0.0,lty=2,col="black",lwd=2)
#axis(2, 1:nrow(cis2), cis2$category_adam_arbitrary, line = 5,las=1)
axis(2,at=c(1:dim(cis2)[1]), labels=cis2$Name_adam,cex.axis=1,las=2)

axis(1,at=c(-1,-0.5,0,0.5,1),cex.axis=1.25)

dev.off()


#Plot rgs with all of these phenotypes




#all rgs range from 0.2-0.4, shitty plot
cis3 <- cis2 #axis gets flipped if changing to horizontal..
pdf('PTSD_rg2.pdf',7,7)
par(mar=c(5, 12, 4, 6) + 0.5)
plotCI(y=1:dim(cis3)[1],x=cis3$rg,li=cis3$LCIrg,ui=cis3$UCIrg,err="x",lwd=2,pch=19,cex.axis=1.25,xlab="Genetic Correlation",main="",cex.lab=1.45,col=beauty_red,scol=alpha(beauty_red,.4),sfrac=0,xaxt='n',yaxt='n',cex=1.8,ylab="",xlim=c(-1,1))
abline(v=0,col='black',lwd=2,lty=2)
#legend("topleft",col=c("white",beauty_red,caro_orange,adam_blue),legend=c("Training -> Target", "MVP -> UKB", "UKB -> MVP", "Meta-analysis -> MVP R4"),bty="n",pch=19,cex=1.5)

axis(1,at=c(-1,-0.5,0,0.5,1),cex.axis=1.25)
axis(2,at=c(1:dim(cis3)[1]), labels=cis3$Name_adam,cex.axis=0.6,las=1)
#axis(2, 1:nrow(cis3), cis3$category_adam_arbitrary, line = 5,las=2)
dev.off()


#MRlap

d1 <- fread('PGC.PTSD_MRlap_results (002).csv',data.table=F)
d1$color <- caro_orange

d1$mr <- d1$IVW
d1$LCImr <- d1$IVW - 1.96*d1$SE
d1$UCImr <- d1$IVW  + 1.96*d1$SE


cis3 <- d1 #axis gets flipped if changing to horizontal..
pdf('PTSD_MR.pdf',7,7)
par(mar=c(5, 12, 4, 6) + 0.5)
plotCI(y=1:dim(cis3)[1],x=cis3$IVW,li=cis3$LCImr,ui=cis3$UCImr,err="x",lwd=2,pch=19,cex.axis=1.25,xlab="MR effect estimate",main="",cex.lab=1.45,col=cis3$color,scol=alpha(cis3$color,.4),sfrac=0,xaxt='n',yaxt='n',cex=1.8,ylab="",xlim=c(-0.5,0.5))
abline(v=0,col='black',lwd=2,lty=2)
#legend("topleft",col=c("white",beauty_red,caro_orange,adam_blue),legend=c("Training -> Tamret", "MVP -> UKB", "UKB -> MVP", "Meta-analysis -> MVP R4"),bty="n",pch=19,cex=1.5)

axis(1,at=c(-1,-0.5,0,0.5,1),cex.axis=1.25)
axis(2,at=c(1:dim(cis3)[1]), labels=cis3$Outcome,cex.axis=1,las=1)

dev.off()


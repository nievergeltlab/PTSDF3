
R
library(data.table)
library(Hmisc)
library(fmsb)
library(data.table)
library(plotrix)
library(lmtest)

aqua_ppt <- rgb(0,176,240,max=255)
beauty_red <- rgb(237,33,36,maxColorValue=255)
orange_ppt <- rgb(255,192,0,max=255)
purple=rgb(179,162,199,max=255)


res2 <- fread('prs_f3_to_rep.csv',data.table=F)
res2$color <- "black"

res2$OR <- exp(res2$Beta)
res2$LCI <- exp(res2$Beta - 1.96*res2[,"SE"])
res2$UCI <- exp(res2$Beta + 1.96*res2[,"SE"])

res2[res2$Ancestry=="European",]$Decile <- res2[res2$Ancestry=="European",]$Decile  
res2[res2$Ancestry=="African",]$Decile <- res2[res2$Ancestry=="African",]$Decile  + 0.075
res2[res2$Ancestry=="HNA",]$Decile <- res2[res2$Ancestry=="HNA",]$Decile + 0.15

res2[res2$Ancestry=="European",]$color <- aqua_ppt
res2[res2$Ancestry=="African",]$color <- beauty_red
res2[res2$Ancestry=="HNA",]$color <- purple


res2A <- subset(res2,select=c(Decile,OR,LCI,UCI,color))

replicationres <- as.data.frame(matrix(nrow=5,ncol=5))
names(replicationres) <- c("Decile","OR","LCI","UCI","color")
replicationres$Decile <-c (1:5) - 0.075


#F2 

replicationres$beta <- c(0,0.13,0.23,0.31,0.41)
replicationres$OR <- exp(replicationres$beta)
replicationres$LCI <- exp( log(replicationres$OR) - 1.96*0.035)
replicationres$UCI<- exp( log(replicationres$OR) + 1.96*0.035)


#F2.5
# replicationres$beta <- c(0,0.16,0.33,0.43,0.61)
# replicationres$OR <- exp(replicationres$beta)
# replicationres$LCI <- exp( log(replicationres$OR) - 1.96*0.035)
# replicationres$UCI<- exp( log(replicationres$OR) + 1.96*0.035)

#F2
#replicationres$OR <- c(1, 1.07,1.13,1.23,1.38)
#replicationres$LCI <- exp( log(replicationres$OR) - 1.96*0.06)
#replicationres$UCI<- exp( log(replicationres$OR) + 1.96*0.06)

replicationres[1,]$LCI <- 1
replicationres[1,]$UCI <- 1

replicationres$color <- orange_ppt 

replicationres2A <- subset(replicationres,select=c(Decile,OR,LCI,UCI,color))


cis2 <- rbind(res2A,replicationres2A)

# layout(matrix(c(1,1,2),nrow=1),widths=c(.65,.75,1))
# #Put space between 
# plotCI(x=cis3$Decile,y=cis3$OR,li=cis3$LCI,ui=cis3$UCI,lwd=2,ylim=c(.8,1.8),pch=19,cex.axis=1.25,xlab="PRS Quintile",ylab="Quintile Odds Ratio (95% CI)",main="",cex.lab=1.45,col=cis3$color,scol=alpha(cis3$color,.4),sfrac=0,xaxt='n')
# legend("topleft",col=c("white",colorlist[c(2,3,1,5)]),legend=c("Training -> Target", "UKBB Women -> UKBB Men", "UKBB Men -> UKBB Women", "UKBB -> PGC 1.5"),bty="n",pch=19,cex=1.5)
# axis(1,at=c(1:5),labels=c(1,2,3,4,5),cex.axis=1.25)

pdf('prs_decile_f3_to_mvpr4_v2.pdf',7,7)
par(mar=c(5, 5, 5, 5) + 0.5,cex.axis=5)
plotCI(x=cis2$Decile,y=cis2$OR,li=cis2$LCI,ui=cis2$UCI,lwd=2,ylim=c(1,3),pch=19,cex.axis=1.25,xlab="PRS Quintile",ylab="Quintile Odds Ratio (95% CI)",main="",cex.lab=1.45,col=cis2$color,scol=alpha(cis2$color,.4),sfrac=0,xaxt='n',yaxt='n',cex=1.8)
legend("topleft",col=c("white",orange_ppt,aqua_ppt,beauty_red,purple),legend=c("Training -> Target","PGC2 -> EA", "PGC3 -> EA",  "PGC3 -> AA", "PGC3 -> LAT"),bty="n",pch=19,cex=1.5)
axis(1,at=c(1:5),cex.axis=1.25)
axis(2,at=c(1,1.25,1.5,1.75,2,2.25,2.5,2.75,3), labels=c("1","1.25","1.5","1.75","2","2.25","2.5","2.75", "3.0"),cex.axis=1.25)

dev.off()

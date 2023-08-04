R
library(data.table)
library(plotrix)
library(scales)

#this needs to be updated for what is in hte mixer plots, i think bipoalr may be old!
res2 <- fread('ptsd_psych_rgs.csv',data.table=F) 


res2$OR <- res2$rg
res2$LCI <- res2$rg - 1.96*res2$se
res2$UCI <- res2$rg + 1.96*res2$se

res2$color <- NA
caro_orange <- rgb(249,164,24,maxColorValue=255)
beauty_red <- rgb(237,33,36,maxColorValue=255)
adam_blue <- rgb(41,153,208,maxColorValue=255)


res2$color <- beauty_red
res2[which(res2$h2_obs /res2$h2_obs_se <=4),]$color <- "grey"

res2 <- res2[nrow(res2):1,]

# res2[which(res2$p2 == "Bipolar (????)"),]$color <- caro_orange
# res2[which(res2$p2 == "MDD (Howard)"),]$color <- adam_blue
# res2[which(res2$p2 == "SCZ 3"),]$color <- beauty_red



#this will be fore the second data


cis2 <- rbind(res2)

pdf('PTSD_rgs.pdf',11,9)
#jpeg('PTSD_rgs.jpeg',width=900,heigh=700,quality=100)
par(mar=c(5,35,3,3) + 0.5)
plotCI(y=1:dim(cis2)[1],xlim=c(0,1),x=cis2$OR,li=cis2$LCI,ui=cis2$UCI,err="x",lwd=2,pch=19,cex.axis=1.25,xlab="Genetic correlation",main="",cex.lab=1.45,col=res2$color,scol=alpha(res2$color,.4),sfrac=0,xaxt='n',yaxt='n',cex=1.8,ylab="")

#legend("topleft",col=c("white",beauty_red,caro_orange,adam_blue),legend=c("Training -> Target", "MVP -> UKB", "UKB -> MVP", "Meta-analysis -> MVP R4"),bty="n",pch=19,cex=1.5)

axis(1,at=c(-1,-0.5,0,0.5,1),cex.axis=1.25)
axis(2,at=c(1:dim(cis2)[1]), labels=cis2$p2,cex.axis=2,las=2)

#axis(2, 1:nrow(cis2), cis2$category_adam_arbitrary, line = 2.5,las=2)

dev.off()



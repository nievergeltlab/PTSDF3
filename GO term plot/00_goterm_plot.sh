R
library(data.table)
library(plotrix)
library(scales)


caro_orange <- rgb(249,164,24,maxColorValue=255)
beauty_red <- rgb(237,33,36,maxColorValue=255)
adam_blue <- rgb(41,153,208,maxColorValue=255)
adam_blue2 <-  rgb(60,160,255,maxColorValue=255)

aqua_ppt <- rgb(0,176,240,max=255)
red_ptt <- red
orange_ppt <- rgb(255,192,0,max=255)


goterm <- fread('eur_f3_goterms.csv',data.table=F)
goterm$color <- adam_blue2

goterm[which(goterm$Term == "Cellular Component"),]$color <- beauty_red
goterm[which(goterm$Term == "Molecular Function"),]$color <- aqua_ppt
goterm[which(goterm$Term == "Biological Process"),]$color <-  orange_ppt # rgb(127,253,166,maxColorValue=255) 


dm <- goterm

dm$logp <- -log10(dm$P)

dm <- dm[nrow(dm):1,]

#pdf('PTSD_goterm.pdf',7,8)
jpeg('PTSD_goterm.jpeg',2000,2000,quality=100,antialias='none')
par(mar=c(15, 5, 5, 5) + 0.5,lwd=5,cex.axis=5)
myplot <- barplot(height=dm$logp,names.arg=dm$FULL_NAME,col=dm$color,las=1,horiz=TRUE,cex.axis=2,cex.lab=5,xaxt='n')
abline(v=-log10(0.05/15483),col='black',lwd=7,lty=2)
axis(1,c(0,2,4,6,8),cex.axis=5,cex.lab=5,lwd=7,lwd.ticks=7,mgp=c(0,5,0))


dev.off()
,xlab="-log10 p-value",
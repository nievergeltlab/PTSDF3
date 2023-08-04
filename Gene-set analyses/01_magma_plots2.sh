R
library(data.table)
library(plotrix)
library(scales)


caro_orange <- rgb(249,164,24,maxColorValue=255)
beauty_red <- rgb(237,33,36,maxColorValue=255)
adam_blue <- rgb(41,153,208,maxColorValue=255)
adam_blue2 <- rgb(60 ,160,254,maxColorValue=255)


magma_whole <- fread('magma_eur_subtissues.csv',data.table=F)

dm <- magma_whole
dm$logp <- -log10(dm$P)
dm$color <- beauty_red
dm <- dm[order(dm$P),]
dm[which(dm$P >= 0.05/54),]$color <- "grey"

pdf('PTSD_tissues.pdf',15,8)
#jpeg('PTSD_tissues.jpg',1500,800)
par(mar=c(16, 5, 5, 5) + 0.5)


plots <- barplot(height=dm$logp,names.arg=dm$FULL_NAME,col=dm$color,horiz=FALSE,ylab="-log10 p-value",srt=75,las=2,xaxt='n',ylim=c(0,12),cex.axis=1.2)
abline(h=-log10(0.05/54),col='black',lwd=2,lty=2)

text(plots,rep(-0.3 ,nrow(dm)),labels=dm$FULL_NAME,srt=70,las=2,xpd=TRUE,adj=1)


dev.off()

https://www.tenderisthebyte.com/blog/2019/04/25/rotating-axis-labels-in-r/

linn <- fread('magma_celltype_Linnarsson_GSE76381_Human_Midbrain.gsa.out',data.table=F)
linn$color <- adam_blue
dm <- rbind(linn) # ,karo
dm <- dm[order(dm$P),]

dm[which(dm$P >= 0.05/268),]$color <- "grey"

dm$logp <- -log10(dm$P)



pdf('PTSD_gsa_linn.pdf',10,8)
#jpeg('PTSD_gsa_linn.jpeg')

par(mar=c(15, 5, 5, 5) + 0.5)
plots <- barplot(height=dm$logp,names.arg=dm$VARIABLE,col=dm$color,las=2,horiz=FALSE,ylab="-log10 p-value",ylim=c(0,8),xaxt='n',srt=75,las=2,cex.axis=1.2)

abline(h=-log10(0.05/268),col='black',lwd=2,lty=2)
text(plots,rep(-0.3 ,nrow(dm)),labels=dm$VARIABLE,srt=70,las=2,xpd=TRUE,adj=1)


dev.off()

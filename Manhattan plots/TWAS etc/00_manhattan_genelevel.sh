library(data.table)

#Read the main GWAS results
 variantgwas <- fread('../freeze3_manhattan/eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.mh',data.table=F)
 
#Read gene-based results
 genegwas=fread('freeze3_eur_genelevel.csv',data.table=F)
 
#Read sizes of each chromosome
 chrfile <- fread('chromosome_lengths.txt',data.table=F)
 

#Get the cumulative position of a chromosome to give the correct plotting
 offsets <- c(0, cumsum(as.numeric(chrfile$bp)))
 
 chrfile$cumpos <-   offsets[-length(offsets)]



 # need to subtract out the initial length to get the start
 chrfile$CHR <- chrfile$Chromosome

#Get genome length, to use as an offset factor
 genome_length=max(chrfile$cumpos)

#Apply offset to map new plot positions for both sets of data (matrix expansion method - works because chromosome number matches row number)
 variantgwas$xpos <- chrfile[variantgwas$CHR,]$cumpos + variantgwas$BP

 genegwas$xpos_start <- chrfile[genegwas$CHR,]$cumpos + genegwas$START
 genegwas$xpos_stop <- chrfile[genegwas$CHR,]$cumpos + genegwas$STOP
 genegwas$color
 
 #Also use matrix expansion to provide colors

 
 
#consider truncating the loci positions liek a circos plot...

#Design an empty plot
pdf('loci_plot.pdf',14,7)
plot(x=c(1,max(offsets )),y=c(0,10),col="white",xaxt="no",yaxt="no",xlab='',ylab='')
axis(1,at=chrfile$cumpos,labels=c(chrfile$Chromosome),cex=1.25) #tick marks for each chroomosome. Right now it starts at 1, which is kind of inconveninted

#For each gene, make it on the plot
 rect(xleft=genegwas$xpos_start,xright=genegwas$xpos_stop ,ybottom=-log10(genegwas$P),ytop=-log10(genegwas$P),col="blue",lwd=1) 
dev.off()

#this is a cool way to do a manhattan plot but may as well just use the other plot at this point...

 genegwas=fread('../freeze3_manhattan_genelevel/freeze3_eur_genelevel.csv',data.table=F)
 
dsub <- genegwas
dsub$BP <- (dsub$START + dsub$STOP) /2
dsub$SNP <- dsub$SYMBOL

dsub2 <- subset(dsub, select=c(CHR,BP,P,SNP))
write.table(dsub2,file="eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.GENEBASED.mh",quote=F,row.names=F)


#14,935 genes in TWAS
twas <- fread('TWAS_Brain.csv',data.table=F)
twas$CHR <- twas$Chr
twas$BP <- (twas$bp_s +twas$bp_e) /2
twas$SNP <- twas$Name
twas$P <- twas$Pval


dsub2 <- subset(twas, select=c(CHR,BP,P,SNP))
write.table(dsub2,file="eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.TWAS.mh",quote=F,row.names=F)

#9003 in SMR
esmr <- fread('eSMR_Brain.csv',data.table=F)
esmr$CHR <- esmr$topSNP_chr
esmr$BP <- esmr$topSNP_bp
esmr$SNP <- esmr$Gene
esmr$P <- esmr$p_SMR # how to color by p heidi?
esmr2 <- subset(esmr,p_HEIDI > 0.05)

dsub2 <- subset(esmr2, select=c(CHR,BP,P,SNP))
write.table(dsub2,file="eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.esmr.mh",quote=F,row.names=F)

#1,209 in pQTL
pqtl <- fread('ptsdf3_pQTL.csv',data.table=F)
pqtl$CHR <- pqtl$topSNP_chr
pqtl$BP <- pqtl$topSNP_bp_GRCh37
pqtl$SNP <- pqtl$Gene
pqtl$P <- apply(pqtl[,c("p_SMR","p_SMR_multi")],1,min)
pqtl2 <- subset(pqtl,p_HEIDI > 0.05)

dsub2 <- subset(pqtl2, select=c(CHR,BP,P,SNP))
write.table(dsub2,file="eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.pqtl.mh",quote=F,row.names=F)



#Annotations table for PLI, etc...
genemaps <- fread('eur_genemapped_proteincoding_sig.csv',data.table=F)
genemaps$SNP <- genemaps$symbol
genemaps$CHR <- genemaps$chr
genemaps$BP <- (genemaps$start + genemaps$end) /2

pli <- genemaps
pli$P <- 10^-pli$pLI

write.table(subset(pli,select=c(CHR,BP,P,SNP)),file="eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.pli.mh",quote=F,row.names=F)

cadd <- genemaps
cadd$P <- 10^-genemaps$posMapMaxCADD 
write.table(subset(cadd,select=c(CHR,BP,P,SNP)),file="eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.cadd.mh",quote=F,row.names=F)

#load finemapping posterior probabilities


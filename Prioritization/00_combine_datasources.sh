
library(data.table)
library(plyr)

ucsc <- fread('zcat ucscgenes.tsv.gz',data.table=F)

#Master list of symbols and ensembl ids
ensemble <- fread('zcat ensemble_genes.tsv.gz',data.table=F)
ensemble2 <- fread('zcat ensemble_genes2.tsv.gz',data.table=F)

names(ensemble) <- c("name","SYMBOL")
ensemble2$EnsembleID <- ensemble2$name2

ensemble3 <- merge(ensemble,ensemble2,by='name',all=TRUE)

ensemble4 <- subset(ensemble3,select=c("EnsembleID","SYMBOL"))
ensemble5 <- ensemble4[!duplicated(ensemble4$SYMBOL),]






 
#gene-based 

 genebased <- read.csv('f3_eur_magma.csv',stringsAsFactors=F,header=T)
 genebased$EnsembleID <- genebased$GENE
 
 genebased$gene_evidence <- 0
 genebased[which(genebased$P <= 2.5e-6),]$gene_evidence <- 2
 
 

#chiayen
 pqtl <- fread('PGCPTSDv3_UKBPPP1500.msmr.summary.txt',data.table=F)
  pqtl$pqtl_yes <- 0
  pqtl[which((pqtl$p_SMR_multi <= 0.05/1209 | pqtl$p_SMR <= 0.05/1209) & pqtl$p_HEIDI > 0.05 ),]$pqtl_yes  <- 1
  
  #consider taking minimum p-value, or convert to wide format, or wide-esque variable that lists protein types that are significant.
  
 pqtl2 <- ddply(pqtl, ~ Gene, colwise(sum,'pqtl_yes'))
 
 pqtl2$pqtl_evidence <- 0
 pqtl2[which(pqtl2$pqtl_yes >= 1),]$pqtl_evidence <- 3
 
 pqtl2$SYMBOL <- pqtl2$Gene
 
 
#Genemapping
 genemap <- fread('f3_eur_genemapping.csv',data.table=F)
 
 
 #Add ensemble ID name
  genemap$EnsembleID <- genemap$ensg

 #Positional mapping
  genemap$pos_evidence <- 0
  genemap[which(genemap$posMapSNPs >0),]$pos_evidence <- 2

 #eQTL mapping 
  genemap$eqtlmap_evidence <- 0
  genemap[which(genemap$eqtlMapSNPs >0),]$eqtlmap_evidence <- 2

 #CI map
  genemap$cimap_evidence <- 0
  genemap[which(genemap$ciMap =="Yes"),]$cimap_evidence <- 2
  
 #CADD score
  genemap$cadd_evidence <- 0
  genemap[which(genemap$posMapMaxCADD >= 12.37 ),]$cadd_evidence <- 15
  
 #PLI score
  genemap$pli_evidence <- 0
  genemap[which(genemap$pLI >= 0.9 ),]$pli_evidence <- 1
  
 #If a gene is annotated due to being implicated by multiple loci, score highly
  
  genemap$locus_evidence <- 0
  genemap[grep(":",genemap$GenomicLocus),]$locus_evidence <- 2
  
#Exonic
  snpstxt <- fread('snps.txt',data.table=F)
  snpstxt <- subset(snpstxt,gwasP <= 5e-8 & dist==0)
  
  snpstxt$exon_yes <- 0
  snpstxt[which(snpstxt$func == "exonic" ),]$exon_yes <- 1
  snpstxt$SYMBOL <- snpstxt$nearestGene
  
  
  snpstxt$rdb_yes <- 0
  snpstxt[which(snpstxt$RDB %in% c("1a"  , "1b" ,  "1c" ,  "1d",  "1f")),]$rdb_yes <- 1
  snpstxt[which(snpstxt$RDB %in% c("2a"  , "2b" ,  "2c")),]$rdb_yes <- 1
  
  snps2a <- ddply(snpstxt, ~ SYMBOL, colwise(sum,'exon_yes'))
  snps2b <- ddply(snpstxt, ~ SYMBOL, colwise(sum,'rdb_yes'))
   
  snps3 <- merge(snps2a,snps2b,by="SYMBOL")
 
  snps3$exon_evidence <- 0
  snps3[which(snps3$exon_yes >= 1),]$exon_evidence <- 1

  snps3$rdb_evidence <- 0
  snps3[which(snps3$rdb_yes >= 1),]$rdb_evidence <- 1

  
  
 #PIPs
  pips <- fread('finemap_pips.txt',data.table=F)
  pips2 <- merge(pips,snpstxt,by="uniqID",suffixes=c("","_snps"))
  pips3 <- ddply(pips2, ~ SYMBOL, colwise(max,'PIP'))
  pips3$pip_evidence <- 1
  pips3[which(pips3$PIP >= 0.7),]$pip_evidence <- 2
  # pips3[which(pips3$PIP >= 0.9),] $pip_evidence <- 1
   
  
 #5 contributions from this file

 #CI map across multiple loci?
 


#TWAS 

 twas <- fread('TWAS_Brain.csv',data.table=F)
 twas$twas_yes <- 0
 twas[which(twas$Pval <= 0.05/14935 ),]$twas_yes <- 1

 
 twas2 <- ddply(twas, ~ EnsembleID, colwise(sum,'twas_yes'))
 
 twas2$twas_evidence <- 0
 twas2[which(twas2$twas_yes >= 1),]$twas_evidence <- 3

#Some of these do not have ensemble IDs...

#Collapse this, needs to be done so that we only count a gene once, 


#ESMR

 smr <- fread('eSMR_Brain.csv',data.table=F)
 smr$EnsembleID <- smr$probeID
 smr$smr_yes <- 0
 smr[which(smr$p_SMR < 0.05/9003 & smr$p_HEIDI > 0.05),]$smr_yes <- 1
 
  smr2 <- ddply(smr, ~ EnsembleID , colwise(sum,'smr_yes'))
 
 smr2$smr_evidence <- 0
 smr2[which(smr2$smr_yes >= 1),]$smr_evidence <- 3
 
 
 #There is an empty ensemble ID somewhere, we don't care about it necessarily..
 #Merge left on ensemble ID
 d1a <- merge(ensemble5,genebased,by="EnsembleID",all=TRUE,suffixes=c("","_genebased"))
 d1b <- merge(d1a, twas2,by="EnsembleID" ,all=TRUE)
 d1c <- merge(d1b,smr2,by="EnsembleID" ,all=TRUE) 
 d1d <- merge(d1c,genemap,by="EnsembleID" ,all=TRUE) #not sure yet if any of the duplicates are ones i care about anyway
 d1e <- merge(d1d,pqtl2,by="SYMBOL" ,all=TRUE,suffixes=c("","_gpwas")) 

 d1f <- merge(d1e,snps3,by="SYMBOL" ,all=TRUE)
 d1g <- merge(d1f,pips3,by="SYMBOL" ,all=TRUE)
  
 #This needs to be condensed because some things are in multiple times. I'll score first, then take the maximum score of each gene..
 
 
 grep("evidence",names(d1),value=TRUE)
 
 #twas2[-which(twas2$EnsembleID %in% ensemble5$EnsembleID),] #sometimes there is evidence, but some are lincrnas and stuff, not just protein coding genes..
# snpstxt[-which(snpstxt$SYMBOL %in% ensemble5$SYMBOL),] #sometimes there is evidence, but some are lincrnas and stuff, not just protein coding genes..
 
 
#A few symbols in twice, due to different transcripts. e.g. d1a[d1a$SYMBOL== "SLC25A10",]
 
 d1 <- d1g
 #Sum evidence columns..
 d1$sumscore <- apply(d1[,grep("evidence",names(d1))],1,sum,na.rm=T)
 d2 <- subset(d1, !is.na(GenomicLocus)) #Only take genes mapped to loci in the main -- the others can pretty much be ignored..
 
 
 # #Only a few dont have scores..
 #Over each locus, get the quantiles of gene scores
 #Tier 1: We want a >20% relative difference over the best and second best
 

 bestgenes <- ddply(d2,~GenomicLocus, .fun=colwise(max,'sumscore'))
 names(bestgenes)[2] <- "best_score"
 secondbest <- function(x)
 {
  sort(na.omit(x),decreasing=TRUE)[2]
 }
 
 secondbestgenes <- ddply(d2,~GenomicLocus, .fun=colwise(secondbest,'sumscore'))
 
 names(secondbestgenes)[2] <- "second_score"

 genecount <- ddply(d2,~GenomicLocus, .fun=colwise(length,'sumscore'))
 names(genecount)[2] <-"Ngenes"
 
  bestcom <- merge(bestgenes,secondbestgenes,by="GenomicLocus")
 bestcom$ratio <- (bestcom$best_score - bestcom$second_score)/bestcom$second_score
 bestcom2 <- merge(bestcom,genecount,by="GenomicLocus")
 bestcom2$gt4 <- 0
 bestcom2[which(bestcom2$best_score >=4),]$gt4 <- 1
 
 bestcom2$tier <- 3
 
 bestcom2[which(bestcom2$best_score >= 4 & (bestcom2$ratio < 0.2 )),]$tier <- 2
  
 bestcom2[which(bestcom2$best_score >= 4 & (bestcom2$ratio >=0.2 | is.na(bestcom2$ratio))),]$tier <- 1
 

 #bestcom2[which(bestcom2$tier != 1 & bestcom2$best_score >= 2 (bestcom2$ratio >=0.2 & bestcom2$ratio >=0.5 )),]$tier <- 2

 #if the leader is not more than 50% better, its good enough to be tier 2
 
 Furthermore, other genes in a locus harboring a tier 1 gene were classified as tier 2 prioritized genes if the relative score difference versus the highest ranked (tier 1)
 gene was between 20% and 50%. Lastly, when the relative score difference between the highest ranked gene and other genes in the same locus was <20%,
 then both the highest ranked gene and all genes with a score difference <20% were classified as tier 2 prioritized genes in the investigated locus; based on the current evidence, it is difficult to prioritize two or more similarly scored genes.
 
 tier1s <- subset(bestcom2,tier==1)
 
 d2_tier1s <- merge(tier1s,d2,by="GenomicLocus")
 d2_tier1s <- d2_tier1s[order(d2_tier1s$sumscore,decreasing=TRUE),]
 
 topcol <- function(x)
 {
  x[1,]
 }
 d2_tier1s_named <- ddply(d2_tier1s,~GenomicLocus,.fun=topcol)
 write.table(d2_tier1s_named,file='tier1.txt',row.names=F)
 
 d3_tier1s <- merge(bestcom2,d2,by="GenomicLocus")
 rowiterate <- function(x){
  x$row <- 1:nrow(x)
  x}
  d3_tier1s <-  d3_tier1s[order( d3_tier1s$GenomicLocus, d3_tier1s$sumscore,decreasing=TRUE),]
  d3_tier1s2 <- ddply(d3_tier1s, ~ GenomicLocus, .fun=rowiterate)

 d3_tier1s2$ratio2 <- (d3_tier1s2$best_score - d3_tier1s2$sumscore)/d3_tier1s2$sumscore
 

  d3_tier1s2[which(d3_tier1s2$ratio2 >= 0.2 & d3_tier1s2$ratio2 <0.5),]$tier <- 2
  d3_tier1s2[which(d3_tier1s2$ratio2 >0.5),]$tier <- 3
  
 d3_tier1s2_fi <- subset(d3_tier1s2, tier <= 2)
 
  write.table(d3_tier1s2,file='tierall.txt',row.names=F)
 
 write.table(d3_tier1s2_fi,file='tiert1t2.txt',row.names=F)
 
 
 #For tier 2 genes, note anything that is close in ratio to tier 1
 
 tier2s <- subset(bestcom2,tier==2)
  d2_tier2s <- merge(tier2s,d2,by="GenomicLocus")
 d2_tier2s <- d2_tier2s[order(d2_tier2s$sumscore,decreasing=TRUE),]
  d2_tier2s_named <- ddply(d2_tier2s,~GenomicLocus,.fun=topcol)
 write.table(d2_tier2s,file='tier2leads.txt')
 
 
 
 #need to expand loci that overlapped, or combine them.
 
 
 
 d3 <- subset(d2,sumscore>=1) # 

 #I need the raw positions of each gene to assign them to loci

 table(d1$sumscore)
 
  head(d1[order(d1$sumscore,decreasing=TRUE),],10)
 
 
 head(d1[order(d1$sumscore,decreasing=TRUE),grep("evidence",names(d1))])
 
 d2 <- subset(d1,sumscore>5)
 
 d2 <- d2[order(d2$sumscore,decreasing=FALSE),]
 
 d3 <- subset(d1,sumscore>=7 )

 d3$GenomicLocus2 <- d3$GenomicLocus
 d3[which(d3$GenomicLocus=="40:41"),]$GenomicLocus2 <-40
  d3[which(d3$GenomicLocus=="27:28"),]$GenomicLocus2 <-27

  
 d3 <- d3[order(d3$chr,d3$START,decreasing=c(F,F)),]
 d3[,c(2,grep('evidence',names(d3)))]
 write.table(d3,file="genes_evidnece_feb17_2023.txt",row.names=F)
 
 d4 <- subset(d1,sumscore==6 )
 d4 <- d4[order(d4$chr,d4$START,decreasing=c(F,F)),]
 #order by genomic locus, but some are 40:41 - need to fix that..
 
 
 #Why are some eQTLs yet NOT in Nikos SMR analysis?
 

 d3$index <- 1:dim(d3)[1]
 d3$pli_color <- "grey"
 d3[which(d3$pli_evidence ==1),]$pli_color <- "red"
 
 d3$cadd_color <- "grey"
 d3[which(d3$cadd_evidence ==1),]$cadd_color <- "red"
 
#Color code this based on where mapping is. maybe have a big brain in the corner??

 d3$cimap_color <- "grey"
d3[which(d3$cimap_evidence ==1),]$cimap_color <- "red"


 #taking hte log for plotting will cause negatives that dont get plotted - adjust


 d3$gene_color <- "grey"
 d3[which(d3$gene_evidence == 1),]$gene_color <- "red"

 d3$eqtlmap_color <- "grey"
 d3[which(d3$eqtlmap_evidence ==1),]$eqtlmap_color <- "red"
 
  d3$finemap_color <- "grey"
 d3[which(d3$finemap_evidence ==1),]$finemap_color <- "red"
 
  d3$twas_color <- "grey"
 d3[which(d3$twas_evidence ==1),]$twas_color <- "red"
 
   d3$smr_color <- "grey"
 d3[which(d3$smr_evidence ==1),]$smr_color <- "red"
 
   d3$pqtl_color <- "grey"
 d3[which(d3$pqtl_evidence ==1),]$pqtl_color <- "red"
 
 
 #not fair to filter to only finemapped!
 
 plot(y=c(1,dim(d3)[1]),x=c(1,10),col="white",yaxt='n',ylab='')
 axis(2,at=1:dim(d3)[1],labels=d3$SYMBOL,las=1)
 #Consider making small points slightly larger
 
 #Finemapping gene
  points(rep(1,dim(d3)[1]), d3$index,cex=1.5,col=d3$finemap_color,pch=15) #x position will be hte slot number, y will  be the row, size/color shape refer to the actual thing


 #Gene based
 points(rep(2,dim(d3)[1]), d3$index,cex=1.5,col=d3$gene_color,pch=15) #x position will be hte slot number, y will  be the row, size/color shape refer to the actual thing



 #PLI scores
 points(rep(3,dim(d3)[1]), d3$index,cex=d3$pLI*1.5,col=d3$pli_color,pch=19) #x position will be hte slot number, y will  be the row, size/color shape refer to the actual thing

 #CADD scores
 points(rep(3.5,dim(d3)[1]), d3$index,cex=abs(log(d3$posMapMaxCADD/9))*2,col=d3$cadd_color,pch=19) #x position will be hte slot number, y will  be the row, size/color shape refer to the actual thing

 #pQTL
  points(rep(5,dim(d3)[1]), d3$index,cex=1.5,col=d3$pqtl_color,pch=15) #x position will be hte slot number, y will  be the row, size/color shape refer to the actual thing

 #TWAS
  points(rep(5.5,dim(d3)[1]), d3$index,cex=1.5,col=d3$twas_color,pch=15) #x position will be hte slot number, y will  be the row, size/color shape refer to the actual thing
  
 #SMR
  points(rep(6,dim(d3)[1]), d3$index,cex=1.5,col=d3$smr_color,pch=15) #x position will be hte slot number, y will  be the row, size/color shape refer to the actual thing
 
 #eQTL
  points(rep(7,dim(d3)[1]), d3$index,cex=1.5,col=d3$eqtlmap_color,pch=15) #x position will be hte slot number, y will  be the row, size/color shape refer to the actual thing

 #CI maps
  points(rep(8,dim(d3)[1]), d3$index,cex=1.5,col=d3$cimap_color,pch=15) #x position will be hte slot number, y will  be the row, size/color shape refer to the actual thing

 #add positional maps back


 dev.off()
 
library(VennDiagram)

ds1  <- subset(d1,sumscore>=6,select=grep("evidence",names(d1),value=TRUE))


venn.diagram(ds1, filename = "venn-4-dimensions.png")

library(ggvenn)
ggvenn(
  ds1, 
  stroke_size = 0.5, set_name_size = 4
  )
  
 #TWAS
  #Consider arrows for up or down regulation (if uniform across sig tissues)

 #pSMR
 
 #eSMR
 
 #can stratify over brain region

 #points(rep(4,dim(d3)[1]), d3$index,cex=abs(d3$ZSTAT)/6,col=d3$smr_color,pch=19) #x position will be hte slot number, y will  be the row, size/color shape refer to the actual thing

 #see if you can get a summary column made somehow - steal from AD
 
 
#I shouldnt have more than 20k genes... yet here I am?
 
 #TRAF3 and CD40..
 #protein analysis??
 
 #can I make a protein bubble plot?? 
 
#add finemapping results
 #1 if gene in finemapping region, 0 if not. Maybe this adds weight?
 
library(data.table)
library(janitor)
library(plyr)
 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }
  
 tdata="pts1"


 ##harmonize to 0-1 based on theoretical ranges
 
 #split up teic kaufman etc based on substudy
 substd <- fread('data/pts1_substudies.txt',data.table=F)
 names(substd)[3] <- "study_abbrev"
 
 fam <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc1.fam',sep=''),data.table=F)
 
 names(fam) <- c("FID","IID","M","F","G","P")
 for (ancgroup in c("eur","aam"))
 {
  pcs <-  read.table(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/covariates/pts_',tdata,'_mix_am-qc1-',ancgroup,'_pca.menv.mds_cov',sep=''),header=T,stringsAsFactors=F)

  dma <- merge(pcs,fam,by=c("FID","IID"),suffixes=c("_pc","")) 
  #dm1x <- merge(dm1,dma,by=c("FID","IID"),suffixes=c("_fam","")) 
  dm1x <- merge(dma,substd,by=c("FID","IID"))
  dm1x$Case <- dm1x$P # Current_PTSD_Dx
  dm1_exp <- subset(dm1x,!is.na(Case),select=c(FID,IID,Case))
 
  print(table(is.na(dm1_exp$Case)))
  
  print(table(is.na(dm1_exp$Case)))
    print(table(dm1_exp$Case))
     counts <- table(dm1_exp$Case)
  4/((1/counts[1]) + (1/counts[2]))


  table(is.na(dm1x$Case),dm1x$G,useNA='always')
  table(is.na(dm1x$Case),dm1x$study_abbrev,useNA='always')

  dm1_exp2 <- subset(dm1x,!is.na(Case),select=c(FID,IID,study_abbrev,Case))
 
  table(dm1_exp2$Case,dm1_exp2$study_abbrev,useNA='always')

  
  write.table(dm1_exp,file=paste("pheno/p2_",tdata,"_pts1","_",ancgroup,".pheno",sep=""),quote=F,row.names=F)
   
  #for the cnv analysis, will export an additional var with study_abbrev
   dm1_exp2 <- subset(dm1x,!is.na(Case),select=c(FID,IID,Case,study_abbrev))
   write.table(dm1_exp2,file=paste("pheno/p2_",tdata,"_pts1","_",ancgroup,".pheno_cnv",sep=""),quote=F,row.names=F)
      
  dm_cov <- dm1x  

  covexp4 <- subset(dm_cov,select=c("FID","IID","C1","C2","C3","C4","C5")) #mrs doesnt need sex!
  
  write.table(covexp4,file=paste("pheno/p2_",tdata,"_",ancgroup,"_pts1","_pcs.cov",sep=""),quote=F,row.names=F)
 }  
  
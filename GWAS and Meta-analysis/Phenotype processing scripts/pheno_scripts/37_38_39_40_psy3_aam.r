library(data.table)
library(janitor)
 
 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }
  
tdata="psy3"

 #dm1 <- fread(paste('data/p2_bry2','.csv',sep=''),data.table=F,na.strings=c(NA,"-9"))
 
 #dm1 <- remove_empty(dm1, which = c("cols"))
 studs <- fread('data/psy3_studies.csv',data.table=F)
 
 fam <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc.fam',sep=''),data.table=F)
 names(fam) <- c("FID","IID","M","F","G","P")
 for (ancgroup in c("aam"))# ,"lat","aam"))
 {
  pcs <-  read.table(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/covariates/pts_',tdata,'_mix_am-qc-',ancgroup,'_pca.menv.mds_cov',sep=''),header=T,stringsAsFactors=F)
  
  dm1x <- merge(pcs,fam,by=c("FID","IID"),suffixes=c("_pc","")) 
    dm1x$Case <- dm1x$P
  
  
  dm1_exp <- subset(dm1x,!is.na(Case),select=c(FID,IID,Case))
 
  table(is.na(dm1_exp$Case))
  
  dm1xa <- merge(dm1x,studs,by="IID")
  dm1_exp2 <- subset(dm1xa,!is.na(P) & P != -9, select=c(FID,IID,Cohort,P))
 
  table(dm1_exp2$Cohort,dm1_exp2$P)


  write.table(dm1_exp,file=paste("pheno/p2_",tdata,"_psy3","_",ancgroup,".pheno",sep=""),quote=F,row.names=F)
    
  dm_cov <- dm1x  

  covexp4 <- subset(dm_cov,select=c("FID","IID","C1","C2","C3","C4","C5")) #mrs doesnt need sex!
  
  write.table(covexp4,file=paste("pheno/p2_",tdata,"_",ancgroup,"_psy3","_pcs.cov",sep=""),quote=F,row.names=F)
 }  
  
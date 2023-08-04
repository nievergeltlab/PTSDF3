library(data.table)
 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }
  
tdata="mamh"

 dm1 <- fread(paste('data/p3_MAMH_Benjet_final.csv',sep=''),data.table=F,na.strings=c(NA,"-9"))
 
 
 wm2 <- function(x,...)
 {
  ff <- which.max(x)
  if(length(ff) > 0)
  {
   return(ff)
  } 
   if(length(ff) == 0)
  {
   return(NA)
  }  
 }
  
 dm1$Current_PTSD_Continuous_whichmax <- apply(dm1[,c("PTSsym_count05",	"PTSsym_count13")],1,wm2)

 takecol <- function(x)
 {
  column_to_pick <- as.numeric(x[length(x)])
  if(!is.na(column_to_pick))
  { 
   return(x[column_to_pick])
  } 
   if(is.na(column_to_pick))
  { 
   return(NA)
  } 
   
 }
 
 dm1$Current_PTSD_Continuous <- apply(dm1[,c("PTSsym_count05",	"PTSsym_count13","Current_PTSD_Continuous_whichmax")],1,takecol)


 fam <- fread(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/qc/pts_',tdata,'_mix_am-qc1.fam',sep=''),data.table=F)
 names(fam) <- c("FID","IID","M","F","G","P")
 
 for (ancgroup in c("lat")) # ,"lat","aam"))
 {
  pcs <-  read.table(paste('/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/',tdata,'/covariates/pts_',tdata,'_mix_am-qc1-natlat_pca.menv.mds_cov',sep=''),header=T,stringsAsFactors=F)
  
  dma <- merge(pcs,fam,by=c("FID","IID"),suffixes=c("_pc","")) 
  
  #Only have IID in this file..
  dm1x <- merge(dm1,dma,by=c("IID"),suffixes=c("_fam","")) 
  
  
  dm1_exp <- subset(dm1x,!is.na(Current_PTSD_Continuous),select=c(FID,IID,Current_PTSD_Continuous))
 
  print(table(is.na(dm1_exp$Current_PTSD_Continuous)))
  
  
  write.table(dm1_exp,file=paste("pheno/p2_",tdata,"_mamh","_",ancgroup,".pheno",sep=""),quote=F,row.names=F)
    
  dm_cov <- dm1x  

  covexp4 <- subset(dm_cov,select=c("FID","IID","C1","C2","C3","C4","C5")) #mrs doesnt need sex!
  
  write.table(covexp4,file=paste("pheno/p2_",tdata,"_",ancgroup,"_mamh","_pcs.cov",sep=""),quote=F,row.names=F)
 }  
   
  
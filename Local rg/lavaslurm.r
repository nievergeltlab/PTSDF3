args <- commandArgs(trailingOnly = TRUE)
 jobnum <- args[1]
 infile <- args[2]
 
 print(jobnum)
library(LAVA)
library(data.table)
#mdd (reformatted) 


# input = process.input(input.info.file="input_ptsd_mdd.txt",           # input info file
                      # sample.overlap.file="sample_overlap.txt",   # sample overlap file (can be set to NULL if there is no overlap)
                      # ref.prefix="g1000_eur",                    # reference genotype data prefix
                      # phenos=c("ptsd","mdd"))       # subset of phenotypes listed in the input info file that we want to process

load(infile)

loci = read.loci(paste('temp/',jobnum,sep=''))
resmat <- as.data.frame(matrix(nrow=dim(loci)[1],ncol=15))
names(resmat) <-c("SNP","CHR","START","STOP","h2_pts","p_pts","h2_mdd","p_mdd","rho","rho.lower","rho.upper","r2","r2.lower","r2.upper","p")

for (i in 1:dim(loci)[1])
{
 testo <- try({
  locusX = process.locus(loci[i,], input)
  uvr <- run.univ(locusX)
  bvr <- run.bivar(locusX)
  resmat[i,1:4] <- loci[i,]
  resmat[i,5] <- uvr[1,2]
  resmat[i,6] <- uvr[1,3]
  resmat[i,7] <- uvr[2,2]
  resmat[i,8] <- uvr[2,3]
  resmat[i,9:15] <- bvr[3:length(bvr)]
    #rm(locusX,uvr,bvr)
  },silent=T)

}
write.table(resmat,file=paste("output_aitd/ptsd_aitd_rg_gw_",jobnum,".txt",sep=""),row.names=F)

cat PGC_UKB_depression_genome-wide.txt_2  | awk '{N="N";Z="Z"; if (NR==1) {$1="SNP";$4="A1";$5="A2";N="N";Z="Z";$9="P"}; if (NR>1) { N="807553"  ; Z=$7/$8; }; print $1,toupper($4),toupper($5),N,Z,$9}' > PGC_UKB_depression_genome-wide.txt_2.lava


library(LAVA)
library(data.table)
#mdd (reformatted) 



library(LAVA)
library(data.table)
#mdd (reformatted) 


input = process.input(input.info.file="input_ptsd_mdd2019.txt",           # input info file
                      sample.overlap.file="sample_overlap2.txt",   # sample overlap file (can be set to NULL if there is no overlap)
                      ref.prefix="g1000_eur",                    # reference genotype data prefix
                      phenos=c("ptsd","mdd"))       # subset of phenotypes listed in the input info file that we want to process

# inspect the processed input data
ls(input)                     # this is actually an environment; hence ls() rather than str()
ls(input$sum.stats)   # processed summary statistics
head(input$sum.stats$bmi)
input$info            # processed input info file, with additional variables computed by process.input()
input$sample.overlap  # sample overlap file


head(input$ref$bim)   # bim file from reference data
  

loci = read.loci("f3_autosome.loci")

resmat <- as.data.frame(matrix(nrow=dim(loci)[1],ncol=15))
names(resmat) <-c("SNP","CHR","START","STOP","h2_pts","p_pts","h2_mdd","p_mdd","rho","rho.lower","rho.upper","r2","r2.lower","r2.upper","p")

for (i in 1:dim(loci)[1])
{
 print (i)
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
write.table(resmat,file="ptsd_mdd2019_rg_5eneg8.txt",row.names=F)



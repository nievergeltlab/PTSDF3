
https://github.com/Kyoko-wtnb/FUMA_scRNA_data

#Combine datasets together



##Linnarsson Midbrain data
/mnt/ukbb/adam/tinnitus_gwas/magma --gene-results ../magma.genes.raw \
 --gene-covar FUMA_scRNA_data-master/processed_data/Linnarsson_GSE76381_Human_Midbrain.txt  \
 --out results/linnarsson_midbrain

#Conditional on average (like website)
/mnt/ukbb/adam/tinnitus_gwas/magma --gene-results ../magma.genes.raw \
 --gene-covar FUMA_scRNA_data-master/processed_data/Linnarsson_GSE76381_Human_Midbrain.txt \
 --model condition-hide=Average   \
 --out  results/linnarsson_midbrain_condavg #human



for tissue in $(cat allKI.expressionlog2.txt | head -n1 | cut -d " " -f 2-)
do

/mnt/ukbb/adam/tinnitus_gwas/magma --gene-results magma.genes.raw \
--gene-covar processed_data/Linnarsson_GSE76381_Human_Midbrain.txt.gz \
--model condition-residualize=$tissue   \
--out  results/ptsdf3_allKI_magma_condition_$tissue #human
done


### Gtex main tissues conditional


#Get genome-wide expression average for other tissues form gtex v8
#cat ../magma.genes.raw | cut -d \  -f1 | tail -n+3 > genes_examined.txt

library(data.table)
#magma <- fread('genes_examined.txt',data.table=F,header=F)

 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )
	return(t(toret))
  }
  
d1 <- fread('/mnt/ukbb/adam/tinnitus_gwas/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',data.table=F,skip=2)
#Remove isoforms
d1$Name <- t(sapply(d1$Name,unlist_split,"[.]"))[,1]

d1 <- d1[,-2] #cut the description column

#Average first
d1$avgcount <- apply(d1[,-c(1)],1,mean)
d1$brainavg <- apply(d1[,grep("Brain",names(d1))],1,mean) #if working with log2, need a different average.
d1$othavg <- apply(d1[,-grep("Brain",names(d1))],1,mean)


#get log2+1 expression
d1[,-1] <- log2(d1[,-1]+1)

#Rename the columns
colnames(d1) <- gsub(" ", ".", colnames(d1))
colnames(d1) <- gsub("-", ".", colnames(d1))

#curiosity - what is the cor between cerebellum and cerebellar hemisphere?
cor(d1$Brain...Cerebellar.Hemisphere,d1$Brain...Cerebellum)

cor(d1$Brain...Cerebellar.Hemisphere,d1$Brain...Hypothalamus)
#conditional correlation is still perfect
summary(lm(d1$Brain...Cerebellar.Hemisphere~d1$Brain...Cerebellum +d1$avgcount))

summary(lm(d1$Brain...Cerebellar.Hemisphere~d1$Brain...Cerebellum +d1$brainavg))

write.table(d1,file="gtexv8_avg.txt",quote=F,row.names=F,sep="\t")

/mnt/ukbb/adam/tinnitus_gwas/magma --gene-results ../magma.genes.raw \
 --gene-covar gtexv8_avg.txt  \
 --model direction=pos condition-hide=avgcount   \
 --out  gtexv8_avg_test #human
 

/mnt/ukbb/adam/tinnitus_gwas/magma --gene-results ../magma.genes.raw \
 --gene-covar gtexv8_avg.txt  \
 --model direction=pos condition-hide=avgcount,brainavg   \
 --out  gtexv8_avg_test2 #human
 
 
#cant do both frontal cortex and avg, will not get a result..
/mnt/ukbb/adam/tinnitus_gwas/magma --gene-results ../magma.genes.raw \
 --gene-covar gtexv8_avg.txt  \
 --model direction=pos condition-hide=avgcount,Brain...Frontal.Cortex.\(BA9\)   \
 --out  gtexv8_avg_test3 #human
 
/mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/eur_ptsd_pcs_v4_aug3_2021.fuma.gz
#zcat /mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/eur_ptsd_pcs_v4_aug3_2021.fuma.gz | awk 'NR>1{print $1,$3,0,$2,toupper($4),toupper($5)}' > eur_ptsd_pcs_v4_aug3_2021.fuma.gz.bim


#With newest magma (to check differences between versions, supposed bias in earlier versions)


library(data.table)

d1 <- fread('/mnt/ukbb/adam/tinnitus_gwas/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',data.table=F,skip=2)
#Remove isoforms
d1$Name <- d1[,2]
d1 <- d1[,-2] #cut the description column

#Average first
d1$avgcount <- apply(d1[,-c(1)],1,mean)
d1$brainavg <- apply(d1[,grep("Brain",names(d1))],1,mean) #if working with log2, need a different average.
#d1$othavg <- apply(d1[,-grep("Brain",names(d1))],1,mean)


#get log2+1 expression
d1[,-1] <- log2(d1[,-1]+1)

#Rename the columns
colnames(d1) <- gsub(" ", ".", colnames(d1))
colnames(d1) <- gsub("-", ".", colnames(d1))

d1 <- d1[!duplicated(d1$Name),]

write.table(d1,file="gtexv8_avg_hugo.txt",quote=F,row.names=F,sep="\t")


 awk '{print $6,$2,$3,$4,$5,$1}' NCBI37.3.gene.loc >  NCBI37.3.gene.loc2
 ./magma --annotate --snp-loc  g1000_eur.bim --gene-loc NCBI37.3.gene.loc2 --out hg19annot
 
 ./magma  --bfile g1000_eur --gene-annot  hg19annot.genes.annot \
 --pval eur_ptsd_pcs_v4_aug3_2021.fuma use=MarkerName,P-value ncol=Weight \
 --out f3_genes
 
  
 ./magma --gene-results f3_genes.genes.raw \
 --gene-covar gtexv8_avg_hugo.txt  \
 --model direction=pos condition-hide=avgcount   \
 --out  gtexv8_avg_test1_magma110 #human
 
 
 ./magma --gene-results f3_genes.genes.raw \
 --gene-covar gtexv8_avg_hugo.txt  \
 --model direction=pos condition-hide=avgcount,brainavg   \
 --out  gtexv8_avg_test2_magma110 #human
 
 
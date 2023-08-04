library(data.table)
d1 <- fread('replications_stein.csv',data.table=F)

d2 <- fread('../lava/f3_autosome.loci',data.table=F)
d2 <- fread('trans_loci.txt',data.table=F)

#Scan to see if the prior hits are within any risk locus

d1$replicated <- NA
#All 3 must be true
for (i in 1:dim(d1)[1])
{
 d1[i,]$replicated <- any(d1[i,]$chr == d2$CHR & d1[i,]$pos >= d2$START & d1[i,]$pos <= d2$STOP)
 }
 
 d2$prior <- NA
 
 for (j in 1:dim(d2)[1])
 {
 d2[j,]$prior <- any(d2[j,]$CHR == d1$chr & d2[j,]$START <= d1$pos & d2[j,]$STOP >= d1$pos)
}
 
 write.table(d2, file="prior_found_trans.txt",quote=F,row.names=F)
 
  write.table(d1, file="replicated_trans.txt",quote=F,row.names=F)
 
 #Scan to see if trans loci are overlapping regular loci
 

d1 <- fread('../lava/f3_autosome.loci',data.table=F)
d2 <- fread('trans_loci.txt',data.table=F)

#Scan to see if the prior hits are within any risk locus

d1$overlap <- NA
#All 3 must be true
for (i in 1:dim(d1)[1])
{
 d1[i,]$overlap <- any(d1[i,]$CHR == d2$CHR & d1[i,]$pos >= d2$START & d1[i,]$pos <= d2$STOP)
 }
 
 d2$prior <- NA
 
 for (j in 1:dim(d2)[1])
 {
  for (k in 1:dim(d1)[1])
  {
   if( d2[j,]$CHR == d1[k,]$CHR & any(d2[j,]$START:d2[j,]$STOP %in% d1[k,]$START:d1[k,]$STOP) )
    {
   d2[j,]$prior <- TRUE
  }
  }
}

 
 write.table(d1, file="trans_locus_overlap.txt",quote=F,row.names=F)
 
#  write.table(d1, file="overlap_trans.txt",quote=F,row.names=F)
 
 
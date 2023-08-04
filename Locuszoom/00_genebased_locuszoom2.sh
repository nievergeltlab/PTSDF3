
#Whole GWAS

     zcat /mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz | \
   awk -v STARTPOS=$bp1 -v STOPPOS=$bp2 -v chr=$chr 'BEGIN{OFS="\t"}{$1=$1;if (NR == 1 ) print "CHR","POS", "A1","A2","MAF","P"; else  print $1,$2,$4,$5,$6,$9}' | sed 's/X/23/g' | sort -g -k1,1 -k2,2 \
   > eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.lz.website.txt
  
  
  
#for genes
 IFS=$'\n'
for gene in $(cat f3_geneset_manytargets.loci | awk '{print $1}' | sort -u  ) 
do 

 locus=$( grep $gene f3_geneset_manytargets.loci | head -n1)
 
 locusnum=$(echo $locus | awk '{print $1}')
 chr=$(echo $locus | awk '{print $3}')
 bp1=$(echo $locus | awk '{print $4}')
 bp2=$(echo $locus | awk '{print $5}')

  #  zcat /mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz | \
 #  awk -v STARTPOS=$bp1 -v STOPPOS=$bp2 -v chr=$chr 'BEGIN{OFS="\t"}{if (NR == 1 ) print "SNP", "P"; else if ($1 == chr && $2 >= STARTPOS - 5000000 && $2 <= STOPPOS + 5000000) print $3,$9}' \
 #  > "$locusnum".lz
   window=0 #for plotting, set window width
  for refsnp in $(grep $gene  f3_geneset_manytargets.loci | awk '{print $2}')
  do

  python2 /mnt/ukbb/locuszoom/bin/locuszoom  \
  --build hg19 --metal "$locusnum".lz --markercol SNP --pvalcol P --ignore-vcf-filter --markercol SNP --refsnp $refsnp --chr $chr --start $(($bp1 - $window))  --end $(($bp2 + $window)) --no-cleanup    \
  --ld-vcf ../../1000_unrel_vcf/UKB_tinnitus_eur_unrelated_1000_"$chr".vcf.gz  --plotonly --prefix "$locusnum"_"$refsnp"


 done
 
 done

 


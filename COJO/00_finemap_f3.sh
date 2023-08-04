#need pairwise LD matrix for all markers. 

#Reformat PTSD GWAS.
#THe program complains if MAF > 0.5. Have to apply some flipping to deal with this..


 #Finemap format
 zcat ../eur_ptsd_pcs_v4_aug3_2021.fuma.gz | awk 'BEGIN{OFS="\t"}{if (NR==1) {$1="chromosome"; $2="position"; $3="rsid";A1="allele1";A2="allele2";maf="maf";BETA="beta";SE="se";}; if(NR>1) {BETA=($8)*3.09/sqrt(2*$7*$6*(1-$6)); SE=BETA/($8+0.000001); A1=toupper($4); A2=toupper($5);maf=$6}; print $3,$1,$2,A1,A2,maf,BETA,SE} ' > eur_ptsd_pcs_v4_aug3_2021.fuma.gz.finemap
 
 #mtcojo format
  zcat ../eur_ptsd_pcs_v4_aug3_2021.fuma.gz | awk 'BEGIN{OFS="\t"}{if (NR==1) {$1="chromosome"; $2="position"; $3="snp";A1="A1";A2="A2";maf="freq";BETA="b";SE="se"; p="p";N="N"}; if(NR>1) {BETA=($8)*3.09/sqrt(2*$7*$6*(1-$6)); SE=BETA/($8+0.000001); A1=toupper($4); A2=toupper($5);maf=$6; p=$9; N=$7}; print $3,A1,A2,maf,BETA,SE,p,N} ' > eur_ptsd_pcs_v4_aug3_2021.fuma.gz.cojo
 
 
 ptsd_gwas=eur_ptsd_pcs_v4_aug3_2021.fuma.gz.cojo
 
#Get list of relevant loci (GenomicRiskLOci.txt)
 cp ../lava/f3_autosome.loci .

#SNP list from GWAS, to find mutual intersection with LD reference
 cat "$ptsd_gwas" | awk '{print $1}' > gwas_snps.txt

#Sort GWAS results by SNP for later join step
 cat "$ptsd_gwas" | LC_ALL=C sort -k1b,1 > "$ptsd_gwas".ordered 
 
#So merged SNP results have a header (rsid is the SNP name in the GWAS output)
 head -n1 "$ptsd_gwas" | awk '{OFS=" "; $1=$1}{print}' > gwas_header.txt
 
#Todo: If there are duplicate snps in the PLINK data, may cause problems!


#Analysis steps
IFS=$'\n'
for locus in $(cat f3_geneset.loci | grep DRD2 ) #f3_autosome.loci for GWAS hits. f3_geneset for non HLA region gene-set results
do 

 #Parse locus file so data can be subset
 locusnum=$(echo $locus | awk '{print $1}')
 chr=$(echo $locus | awk '{print $2}')
 bp1=$(echo $locus | awk '{print $3}')
 bp2=$(echo $locus | awk '{print $4}')
 
  
 #Subset the genotype data to the overlapping SNPs
  ../../plink --bfile /mnt/ukbb/adam/ptsd/UKBB/bgeneur/ukb_unrel_"$chr" --chr $chr --from-bp $bp1 --to-bp $bp2 --extract gwas_snps.txt --make-bed --out ukb_ref/locus_"$locusnum"
 
 #Establish an ordering for GWAS results to corresopnd to the pairwise LD matrix.
  awk '{print NR, $2}' ukb_ref/locus_"$locusnum".bim |  LC_ALL=C sort -k2b,2 > ukb_ref/locus_"$locusnum".bim.joinable
 
 #Subset GWAS results to overlapping SNPs in the correct order. This is sent to the finemap software
  LC_ALL=C join -1 2 -2 1 ukb_ref/locus_"$locusnum".bim.joinable "$ptsd_gwas".ordered  | LC_ALL=C sort -g -k2 | cut -d " " -f1,3- | cat gwas_header.txt - > sumstats_reformat/"$ptsd_gwas".locus"$locusnum".z
  


 # #Pairwise LD. This is sent to the finemap software
  # ../plink --bfile ukb_ref/locus_"$locusnum" --r square --out corrmat/"$ptsd_gwas".locus"$locusnum"
  
  # #Have to convert format from tab delim..
  # awk 'OFS=" "{$1=$1; print}' corrmat/"$ptsd_gwas".locus"$locusnum".ld > corrmat/"$ptsd_gwas".locus"$locusnum".ld2 
  # mv corrmat/"$ptsd_gwas".locus"$locusnum".ld2  corrmat/"$ptsd_gwas".locus"$locusnum".ld
  
 # #error checking: 
 # wc -l ukb_ref/locus_"$locusnum".bim #should be n
 # wc -l "$ptsd_gwas".locus"$locusnum".z #Should be n+1
 
 # echo "z;ld;snp;config;cred;log;n_samples" > master.head
 # echo "sumstats_reformat/"$ptsd_gwas".locus"$locusnum".z;corrmat/"$ptsd_gwas".locus"$locusnum".ld;output/"$ptsd_gwas".locus"$locusnum".snp;output/"$ptsd_gwas".locus"$locusnum".config;dataset1.cred;output/"$ptsd_gwas".locus"$locusnum".log;641533" | cat master.head - > finemapmaster/locus"$locusnum".master

 # #run finemap
 # /mnt/ukbb/adam/ptsd/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64 --cond --in-files finemapmaster/locus"$locusnum".master
 
 # #Sort results by conditional p-value. Take top result. See if this is significant... there are multiple z-resid columns potentially. needs update!
 # awk 'NR>1{print sqrt($3*$3)}' output/eur_ptsd_pcs_v4_aug3_2021.fuma.gz.finemap.locus"$locusnum".snp | sort -gr | head -n1 | awk -v locus=$locusnum '{print locus,$1}' >> condmapping.txt


 #run mtcojo 
  #set threshold to suggestive SNPs for gene-based analysis...
  pthresh=1e-3
   gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --cojo-file sumstats_reformat/"$ptsd_gwas".locus"$locusnum".z --cojo-slct --cojo-p $pthresh --bfile ukb_ref/locus_"$locusnum" --out cojotest_genes/locus"$locusnum"_$pthresh

done
 
#Analysis steps

#For loci
IFS=$'\n'
for locus in $(cat f3_autosome.loci | tail -n+2 )
do 
 locusnum=$(echo $locus | awk '{print $1}')
 chr=$(echo $locus | awk '{print $2}')
 bp1=$(echo $locus | awk '{print $3}')
 bp2=$(echo $locus | awk '{print $4}')
 
 awk -v locus=$locusnum '{print locus,$0}'  cojotest/locus"$locusnum".jma.cojo  >> cojo.results
 done
 
 
 #for genes
 IFS=$'\n'
for locus in $(cat f3_geneset.loci  )
do 
 locusnum=$(echo $locus | awk '{print $1}')
 chr=$(echo $locus | awk '{print $2}')
 bp1=$(echo $locus | awk '{print $3}')
 bp2=$(echo $locus | awk '{print $4}')
 
 awk -v locus=$locusnum '{print locus,$0}'  cojotest_genes/locus"$locusnum".jma.cojo  >> cojo_genes.results
 done
 
  cat cojo_genes.results | sort -g -k 2 | awk '{if (NR==1 || $2!="Chr")print}' >cojo_genes.results.1e-3.txt
  
 
#No genes selected (no p < 1e-5) for the following :
#awk: fatal: cannot open file `cojotest_genes/locusGRM7.jma.cojo' for reading (No such file or directory)
#awk: fatal: cannot open file `cojotest_genes/locusKCNMB2.jma.cojo' for reading (No such file or directory)
#awk: fatal: cannot open file `cojotest_genes/locusCACNA2D3.jma.cojo' for reading (No such file or directory)

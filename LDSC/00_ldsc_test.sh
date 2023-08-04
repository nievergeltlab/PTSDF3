 conda activate ldsc
 
 LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' w_hm3.snplist | LC_ALL=C sort -k1b,1 ) <(zcat eur_ehr_PTSD_Broad_sswv2_pgcptsdf251.tbl.fuma.gz | awk '{if (NR == 1) $3="SNP"; print}'   | LC_ALL=C sort -k3b,3 ) | gzip > eur_ehr_PTSD_Broad_sswv2_pgcptsdf251.tbl.fuma2.gz
 python /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --sumstats eur_ehr_PTSD_Broad_sswv2_pgcptsdf251.tbl.fuma2.gz --signed-sumstats Zscore,0 --snp MarkerName --out eur_ehr_PTSD_Broad_sswv2_pgcptsdf251.tbl.fuma.munge.gz #add --N-col OBS_CT for the sample size
 
 
  /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/ldsc.py \
 --h2 eur_ehr_PTSD_Broad_sswv2_pgcptsdf251.tbl.fuma.munge.gz.sumstats.gz \
 --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
 --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
 --out eur_ehr_PTSD_Broad_sswv2_pgcptsdf251.tbl.fuma.munge.gz.ldsc
 
 
 
  
 LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' w_hm3.snplist | LC_ALL=C sort -k1b,1 ) <(zcat eur_ehr_PTSD_Broad_sswv2_pgcptsdf25_mvp.tbl.fuma.gz | awk '{if (NR == 1) $3="SNP"; print}'   | LC_ALL=C sort -k3b,3 ) | gzip > eur_ehr_PTSD_Broad_sswv2_pgcptsdf25_mvp.tbl.fuma2.gz
 python /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --sumstats eur_ehr_PTSD_Broad_sswv2_pgcptsdf25_mvp.tbl.fuma2.gz --signed-sumstats Zscore,0 --snp MarkerName --out eur_ehr_PTSD_Broad_sswv2_pgcptsdf25_mvp.tbl.fuma.munge.gz #add --N-col OBS_CT for the sample size
 
 
  /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/ldsc.py \
 --h2 eur_ehr_PTSD_Broad_sswv2_pgcptsdf25_mvp.tbl.fuma.munge.gz.sumstats.gz \
 --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
 --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
 --out eur_ehr_PTSD_Broad_sswv2_pgcptsdf25_mvp.tbl.fuma.munge.gz.ldsc
 
 
 
  LC_ALL=C join -1 1 -2 1 <(awk '{print $1}' w_hm3.snplist | LC_ALL=C sort -k1b,1 ) <(zcat  eur_dec28_2017_maf01_info6.results_nefff.gz | awk '{if (NR == 1) $1="SNP"; print}'   | LC_ALL=C sort -k1b,1 ) | gzip > eur_dec28_2017_maf01_info6.results_nefff.fuma2.gz
 python /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --sumstats eur_dec28_2017_maf01_info6.results_nefff.fuma2.gz   --signed-sumstats Effect,0 --N-col Neff --out eur_dec28_2017_maf01_info6.results_nefff.fuma.munge.gz #add --N-col OBS_CT for the sample size
 
   LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' w_hm3.snplist | LC_ALL=C sort -k1b,1 ) <(zcat eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz | awk '{if (NR == 1) $3="SNP"; print}'   | LC_ALL=C sort -k3b,3 ) | gzip > eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma2.gz
 python /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --sumstats eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma2.gz --signed-sumstats Zscore,0 --snp MarkerName --out eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.munge.gz #add --N-col OBS_CT for the sample size
 
  /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/ldsc.py \
 --rg eur_dec28_2017_maf01_info6.results_nefff.fuma.munge.gz.sumstats.gz,eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.munge.gz.sumstats.gz \
 --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
 --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
 --out freeze2_and_eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.munge.gz.ldsc
 
 
 
 

 
 
 
 #Format for LDSC
LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat eur_broad_mar172020_allchr.broad_tinnitus.maf01.resultsa.txt.fuma.gz   |  LC_ALL=C sort -u -k2b,2 ) > eur_broad_mar172020_allchr.broad_tinnitus.maf01.resultsa.txt.fuma.gz.premunge
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz    | LC_ALL=C sort -u -k1b,1 ) > f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.premunge
LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat metal_results/ukbbcoding3_related_nocov_mvp_eur_broad1.tbl.fuma.gz | awk '{if (NR == 1) $3="SNP"; print}'  | LC_ALL=C sort -u -k3b,3 ) > metal_results/ukbbcoding3_related_nocov_mvp_eur_broad1.tbl.fuma.gz.premunge


conda activate ldsc

awk '{if (NR ==1) {$1="SNP";$5="A1";$6="A2"} ; print $1,$5,$6}' eur_broad_mar172020_allchr.broad_tinnitus.maf01.resultsa.txt.fuma.gz.premunge  > eur_broad_mar172020_allchr.broad_tinnitus.maf01.resultsa.txt.fuma.gz.premunge.alleles
awk '{if (NR ==1) {$1="SNP";$4="A1";$5="A2"} ; print $1,$4,$5}' f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.premunge  > f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.premunge.alleles


 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --sumstats f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.premunge --merge-alleles eur_broad_mar172020_allchr.broad_tinnitus.maf01.resultsa.txt.fuma.gz.premunge.alleles --N 175995 --a1 ALLELE1 --a2 ALLELE0 --out f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.munge.gz
 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --sumstats eur_broad_mar172020_allchr.broad_tinnitus.maf01.resultsa.txt.fuma.gz.premunge --N 281398 --a2 AX --merge-alleles f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.premunge.alleles --out eur_broad_mar172020_allchr.broad_tinnitus.maf01.resultsa.txt.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size
 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --sumstats metal_results/ukbbcoding3_related_nocov_mvp_eur_broad1.tbl.fuma.gz.premunge  --out metal_results/ukbbcoding3_related_nocov_mvp_eur_broad1.tbl.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size


#TBI

 LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' /mnt/ukbb/adam/ptsd/w_hm3.snplist | LC_ALL=C sort -k1b,1 ) <(zcat eur_tbistrict_dec172020_allchr.tbistrict.maf01.resultsa.fuma.gz | awk '{if (NR == 1) $1="SNP"; print}'   | LC_ALL=C sort -k1b,1 ) | gzip > eur_tbistrict_dec172020_allchr.tbistrict.maf01.resultsa.ldsc1.gz
 python /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --sumstats  eur_tbistrict_dec172020_allchr.tbistrict.maf01.resultsa.ldsc1.gz --a1 A1 --a2 AX --frq A1_FREQ  --N-col OBS_CT  --out eur_tbistrict_dec172020_allchr.tbistrict.maf01.resultsa.munge.gz #add --N-col OBS_CT for the sample size
 
 zcat eur_tbistrict_dec172020_allchr.tbistrict.maf01.resultsa.ldsc1.gz | awk '{if (NR==1){ $2="CHR"; $3="BP";$5="A2";$6="FRQ";$7="INFO";$8="N"} print $1,$2,$3,$4,$5,$6,$7,$8,$9,$12}' > eur_tbistrict_dec172020.ldhub
 
 
 zip -1 eur_tbistrict_dec172020.sumstats.zip eur_tbistrict_dec172020.ldhub
 
  /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/ldsc.py \
 --h2 eur_tbistrict_dec172020_allchr.tbistrict.maf01.resultsa.munge.gz.sumstats.gz \
 --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
 --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
 --out eur_tbistrict_dec172020_allchr.tbistrict.maf01.resultsa.ldsclog
 
 
 
   /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/ldsc.py \
 --rg eur_tbistrict_dec172020_allchr.tbistrict.maf01.resultsa.munge.gz.sumstats.gz,/mnt/ukbb/adam/ptsd/eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.gz1.tbl.munge.gz.sumstats.gz \
 --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
 --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
 --out tbi_ptsd
 
 
 #Freeze 2 v 2.5
 
 
 #F3
 
 
  /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/ldsc.py \
 --h2 eur_ptsd_pcs_v4_aug3_2021.fuma.gz.munge.gz.sumstats.gz \
 --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
 --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
 --out eur_ptsd_pcs_v4_aug3_2021.fuma.gz.ldsc
 
 Total Observed scale h2: 0.0532 (0.002)
Lambda GC: 1.5511
Mean Chi^2: 1.7194
Intercept: 1.0524 (0.0097)
Ratio: 0.0729 (0.0134)
Analysis finished at Mon Aug 22 12:34:26 2022
Total time elapsed: 8.75s

 
 
 
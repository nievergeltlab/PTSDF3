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

 
 #Stratified ldsc with the 1.2 baseline annotations from finucane et al.
 
  python2 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/ldsc.py \
 --h2 eur_ptsd_pcs_v4_aug3_2021.fuma.gz.munge.gz.sumstats.gz \
 --ref-ld-chr baseline_v1.2/baseline. \
 --w-ld-chr weights_hm3_no_hla/weights. \
 --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
 --overlap-annot \
 --out f3_ldsc/eur_ptsd_pcs_v4_aug3_2021.fuma.gz.ldsc.stratified
 


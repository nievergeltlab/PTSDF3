  /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/ldsc.py \
 --rg  eur_ptsd_pcs_v4_females_oct11_2021.tbl.munge.gz.sumstats.gz,eur_ptsd_pcs_v4_males_oct11_2021.tbl.munge.gz.sumstats.gz \
 --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
 --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
 --out ptsdf3_males_femalesrg
 
 
 
 

 
 
 LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat eur_ptsd_pcs_v4_aug3_2021_noukbb.fuma.gz  | sed 's/MarkerName/SNP/g'  | LC_ALL=C sort -u -k3b,3 ) > eur_ptsd_pcs_v4_aug3_2021_noukbb.fuma.gz.premunge
python2 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats eur_ptsd_pcs_v4_aug3_2021_noukbb.fuma.gz.premunge  --out eur_ptsd_pcs_v4_aug3_2021_noukbb.fuma.gz.munge.gz


  python2 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/ldsc.py \
 --h2 eur_ptsd_pcs_v4_aug3_2021_noukbb.fuma.gz.munge.gz.sumstats.gz \
 --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
 --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
 --out ptsdf3_nokbb_h2
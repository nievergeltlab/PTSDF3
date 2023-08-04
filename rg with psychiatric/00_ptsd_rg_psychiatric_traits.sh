#Perform genetic correlations between a trait and all existing PGC psychiatric traits (as of july 2021)

#need to make a file where the possible alleles are given - will crash otherwise, because sometimes (rare) the same rs has different coded alleles 
#Basically just make a file consisting of the SNP, A1, A2 from your main dataset where you want to see what traits are genetically correlated to it.

LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/eur_ptsd_pcs_v4_aug3_2021.fuma.gz | awk '{if(NR==1) {$4="A1"; $5="A2"}; print $3,toupper($4),toupper($5)}' | LC_ALL=C sort -k1b,1 | sed 's/MarkerName/SNP/g')  > eur_ptsd_pcs_v4_aug3_2021.fuma.alleles

#Reformat main summary statistics
LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/eur_ptsd_pcs_v4_aug3_2021.fuma.gz  | sed 's/MarkerName/SNP/g'  | LC_ALL=C sort -u -k3b,3 ) > eur_ptsd_pcs_v4_aug3_2021.fuma.gz.premunge
 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats eur_ptsd_pcs_v4_aug3_2021.fuma.gz.premunge  --out eur_ptsd_pcs_v4_aug3_2021.fuma.gz.munge.gz



#Reformat all other summary statistics for ldsc

LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/adhd_eur_jun2017.gz   |  LC_ALL=C sort -u -k2b,2 ) > rgpsych/adhd_eur_jun2017.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/adhd_eur_jun2017.gz.premunge --N 53293 --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/adhd_eur_jun2017.gz.munge.gz
 
LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/pgc_alcdep.eur_discovery.aug2018_release_FIXED.txt.gz  |  LC_ALL=C sort -u -k2b,2 | cut -d " " -f1-8 ) > rgpsych/pgc_alcdep.eur_discovery.aug2018_release.txt.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/pgc_alcdep.eur_discovery.aug2018_release.txt.gz.premunge --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/pgc_alcdep.eur_discovery.aug2018_release.txt.gz.munge.gz
 
LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/pgcAN2.2019-07.vcf.tsv.gz | grep -v "#" | awk '{if (NR==1) $3="SNP"; print}'  |  LC_ALL=C sort -u -k3b,3 ) > rgpsych/pgcAN2.2019-07.vcf.tsv.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/pgcAN2.2019-07.vcf.tsv.gz.premunge --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --a1 REF --a2 ALT --N 72517 --out rgpsych/pgcAN2.2019-07.vcf.tsv.gz.munge.gz
 
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/anxiety.meta.full.cc.tbl.gz | awk '{if (NR==1) {$1="SNP"; $10="N"}; print}'  |  LC_ALL=C sort -u -k1b,1 ) > rgpsych/anxiety.meta.full.cc.tbl.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/anxiety.meta.full.cc.tbl.gz.premunge --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/anxiety.meta.full.cc.tbl.gz.munge.gz
 
LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/iPSYCH-PGC_ASD_Nov2017.gz   |  LC_ALL=C sort -u -k2b,2 ) > rgpsych/iPSYCH-PGC_ASD_Nov2017.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/iPSYCH-PGC_ASD_Nov2017.gz.premunge --N 46351 --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/iPSYCH-PGC_ASD_Nov2017.gz.munge.gz
 
LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/daner_PGC_BIP32b_mds7a_0416a.gz   |  LC_ALL=C sort -u -k2b,2 ) > rgpsych/daner_PGC_BIP32b_mds7a_0416a.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/daner_PGC_BIP32b_mds7a_0416a.gz.premunge --N 51710 --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/daner_PGC_BIP32b_mds7a_0416a.gz.munge.gz
 
 LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/ptsd/freeze3/freeze3_psychiatric_rgs/bip2021.gz  | tail -n+73  | sed 's/#//g' | sed 's/ID/SNP/' |  LC_ALL=C sort -u -k3b,3 ) > rgpsych/bip2021.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/bip2021.gz.premunge --N 413466 --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/bip2021.gz.munge.gz
 
 
 
LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/Cannabis_ICC_UKB_het.txt.gz   |  LC_ALL=C sort -u -k2b,2 ) > rgpsych/Cannabis_ICC_UKB_het.txt.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/Cannabis_ICC_UKB_het.txt.gz.premunge --ignore Z --N 162082 --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/Cannabis_ICC_UKB_het.txt.gz.munge.gz
 
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(cat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/PGC_UKB_depression_genome-wide.txt  | awk '{if (NR==1) $1="SNP"; print}'  |  LC_ALL=C sort -u -k1b,1 ) > rgpsych/PGC_UKB_depression_genome-wide.txt.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/PGC_UKB_depression_genome-wide.txt.premunge --signed-sumstats LogOR,0 --N 500199 --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/PGC_UKB_depression_genome-wide.txt.munge.gz
 
LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/ocd_aug2017.gz   |  LC_ALL=C sort -u -k2b,2 ) > rgpsych/ocd_aug2017.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/ocd_aug2017.gz.premunge --N 9725 --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/ocd_aug2017.gz.munge.gz
 
LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz | awk '{if (NR==1) $3="SNP"; print}'  |  LC_ALL=C sort -u -k3b,3 ) > rgpsych/eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz.premunge --N 55374 --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz.munge.gz
 
LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/PGC3_SCZ_wave3_public.v2.tsv.gz   |  LC_ALL=C sort -u -k2b,2 ) > rgpsych/PGC3_SCZ_wave3_public.v2.tsv.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/PGC3_SCZ_wave3_public.v2.tsv.gz.premunge --N 161405 --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/PGC3_SCZ_wave3_public.v2.tsv.gz.munge.gz
 
 LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/SA_in_MDD_BIP_SCZ_2019.gz   |  LC_ALL=C sort -u -k2b,2 ) > rgpsych/SA_in_MDD_BIP_SCZ_2019.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/SA_in_MDD_BIP_SCZ_2019.gz.premunge --N 23801 --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/SA_in_MDD_BIP_SCZ_2019.gz.munge.gz

LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/TS_Oct2018.gz  |  LC_ALL=C sort -u -k1b,1 ) > rgpsych/TS_Oct2018.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/TS_Oct2018.gz.premunge --N 14307 --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/TS_Oct2018.gz.munge.gz
 
 LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' /mnt/ukbb/software/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/ptsd/freeze3/freeze3_psychiatric_rgs/ADHD2022_iPSYCH_deCODE_PGC.meta.gz  |  LC_ALL=C sort -u -k2b,2 ) > rgpsych/ADHD2022_iPSYCH_deCODE_PGC.meta.gz.premunge
/mnt/ukbb/software/ldsc/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/ADHD2022_iPSYCH_deCODE_PGC.meta.gz.premunge --N 225534 --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/ADHD2022_iPSYCH_deCODE_PGC.meta.gz.munge.gz
 
 
#Hearing

LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(cat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/GERA-EUR-ARHI.mtcojo  |  LC_ALL=C sort -u -k1b,1 ) > rgpsych/GERA-EUR-ARHI.mtcojo.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/GERA-EUR-ARHI.mtcojo.premunge --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/GERA-EUR-ARHI.mtcojo.munge.gz
 

#Insomnia
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' /mnt/ukbb/adam/tinnitus_gwas/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/Insomnia_sumstats_Jansenetal.txt.gz  | awk '{if (NR == 1 || (length($5) == 1 && length($6) == 1)) print}' |  LC_ALL=C sort -u -k1b,1 ) > rgpsych/Insomnia_sumstats_Jansenetal.txt.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgpsych/Insomnia_sumstats_Jansenetal.txt.gz.premunge --merge-alleles eur_ptsd_pcs_v4_aug3_2021.fuma.alleles --out rgpsych/Insomnia_sumstats_Jansenetal.txt.gz.munge.gz
 

  /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/ldsc.py \
 --h2  eur_ptsd_pcs_v4_aug3_2021.fuma.gz.munge.gz.sumstats.gz \
 --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
 --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
 --out rgpsych/ptsobsh2
 
 --samp-prev 0.3427468 \
 --pop-prev 0.3427468 \

#Get genetic correlations
for analysis in ADHD2022_iPSYCH_deCODE_PGC.meta.gz.munge.gz.sumstats.gz # PGC_UKB_depression_genome-wide.txt.munge.gz.sumstats.gz #GERA-EUR-ARHI.mtcojo.munge.gz.sumstats.gz   Insomnia_sumstats_Jansenetal.txt.gz.munge.gz.sumstats.gz SA_in_MDD_BIP_SCZ_2019.gz.munge.gz.sumstats.gz  adhd_eur_jun2017.gz.munge.gz.sumstats.gz pgc_alcdep.eur_discovery.aug2018_release.txt.gz.munge.gz.sumstats.gz pgcAN2.2019-07.vcf.tsv.gz.munge.gz.sumstats.gz anxiety.meta.full.cc.tbl.gz.munge.gz.sumstats.gz iPSYCH-PGC_ASD_Nov2017.gz.munge.gz.sumstats.gz bip2021.gz.munge.gz.sumstats.gz Cannabis_ICC_UKB_het.txt.gz.munge.gz.sumstats.gz  ocd_aug2017.gz.munge.gz.sumstats.gz eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz.munge.gz.sumstats.gz PGC3_SCZ_wave3_public.v2.tsv.gz.munge.gz.sumstats.gz  TS_Oct2018.gz.munge.gz.sumstats.gz 
do
 outname=$(echo $analysis | cut -d "." -f1)
 echo $outname
  /mnt/ukbb/software/ldsc/ldsc.py \
 --rg  eur_ptsd_pcs_v4_aug3_2021.fuma.gz.munge.gz.sumstats.gz,rgpsych/$analysis \
 --ref-ld-chr /mnt/ukbb/software/ldsc//eur_w_ld_chr/ \
 --w-ld-chr  /mnt/ukbb/software/ldsc//eur_w_ld_chr/ \
 --out rgs/ptsd_$outname
 done
  
#Send rg outputs to a single file
for analysis in Insomnia_sumstats_Jansenetal.txt.gz.munge.gz.sumstats.gz GERA-EUR-ARHI.mtcojo.munge.gz.sumstats.gz  SA_in_MDD_BIP_SCZ_2019.gz.munge.gz.sumstats.gz  adhd_eur_jun2017.gz.munge.gz.sumstats.gz pgc_alcdep.eur_discovery.aug2018_release.txt.gz.munge.gz.sumstats.gz pgcAN2.2019-07.vcf.tsv.gz.munge.gz.sumstats.gz anxiety.meta.full.cc.tbl.gz.munge.gz.sumstats.gz iPSYCH-PGC_ASD_Nov2017.gz.munge.gz.sumstats.gz bip2021.gz.munge.gz.sumstats.gz Cannabis_ICC_UKB_het.txt.gz.munge.gz.sumstats.gz PGC_UKB_depression_genome-wide.txt.munge.gz.sumstats.gz ocd_aug2017.gz.munge.gz.sumstats.gz eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz.munge.gz.sumstats.gz PGC3_SCZ_wave3_public.v2.tsv.gz.munge.gz.sumstats.gz  TS_Oct2018.gz.munge.gz.sumstats.gz
do
outname=$(echo $analysis | cut -d "." -f1)
 echo $outname
  grep "Summary of Genetic Correlation Results" -A2  rgpsych/rgs/tbi_$outname.log | tail -n1 >> rgpsych/mvptbi_rgs.txt
   done
  
  
  
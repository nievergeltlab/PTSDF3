#Now including flashPCA PCs in eur (see 00_pca_ukbb_freeze3.sh)
#the f.22000.0.0 being -9 is fixed

    
    
  BOLT-LMM_v2.3.5/bolt \
    --bed=pca/UKB_ptsd_eur_related_{1..22}.bed \
    --bim=pca/UKB_ptsd_eur_related_{1..22}.bim \
    --fam=pca/UKB_ptsd_eur_related_1.fam \
    --LDscoresFile=BOLT-LMM_v2.3.5/tables/LDSCORE.1000G_EUR.tab.gz \
    --geneticMapFile=BOLT-LMM_v2.3.5/tables/genetic_map_hg19_withX.txt.gz \
    --numThreads=6 \
    --bgenFile=/mnt/ukbb/adam/bgen/ukb_imp_chr{1..22}_v3.bgen \
    --bgenMinMAF=1e-2 \
    --bgenMinINFO=0.8 \
    --sampleFile=/mnt/ukbb/adam/ptsd/preimputation_genotypes/ukb41209_imp_chr1_v3_s487320.sample \
    --verboseStats \
    --lmmForceNonInf \
    --statsFile=pclc_ukbb_may13_2021_unrelated.stats.gz \
    --statsFileBgenSnps=pts_ukbb_may13_2021_unrelated.bgen.stats.gz \
    --phenoFile=UKB_ptsd_eur_unrelated_m0_may13_2021.pheno \
    --phenoCol=PCL_sum_imp  \
    --covarFile=UKB_ptsd_eur_unrelated_m0_may13_2021.pheno \
    --covarCol=f.22000.0.0  --covarCol=f.54.0.0  --qCovarCol=PC{1:10}   \
    --covarMaxLevels=107 \
    --remove bolt.in_plink_but_not_imputed.FID_IID.969.txt \
    --noBgenIDcheck
    
    
    #if this doesnt run, change map file .
    # Number of samples in BGEN header does not match sample file
    #god fucking damn it. I may need to just line up and convert.
  BOLT-LMM_v2.3.5/bolt \
    --bed=pca/UKB_ptsd_eur_related_{1..22}.bed \
    --bim=pca/UKB_ptsd_eur_related_{1..22}.bim \
    --fam=pca/UKB_ptsd_eur_related_1.fam \
    --LDscoresFile=BOLT-LMM_v2.3.5/tables/LDSCORE.1000G_EUR.tab.gz \
    --geneticMapFile=BOLT-LMM_v2.3.5/tables/genetic_map_hg19_withX.txt.gz \
    --numThreads=6 \
    --bgenFile=/mnt/ukbb/adam/bgen/ukb_imp_chr23_v3.bgen \
    --bgenMinMAF=1e-2 \
    --bgenMinINFO=0.8 \
    --sampleFile=/mnt/ukbb/adam/ptsd/preimputation_genotypes/ukb41209_imp_chrX_v3_s486561.sample \
    --verboseStats \
    --lmmForceNonInf \
    --statsFile=pclc_ukbb_may13_2021_unrelated.stats.gz \
    --statsFileBgenSnps=pts_ukbb_may13_2021_unrelated_chrX.bgen.stats.gz \
    --phenoFile=UKB_ptsd_eur_unrelated_m0_may13_2021.pheno \
    --phenoCol=PCL_sum_imp  \
    --covarFile=UKB_ptsd_eur_unrelated_m0_may13_2021.pheno \
    --covarCol=f.22000.0.0  --covarCol=f.54.0.0  --qCovarCol=PC{1:10}   \
    --covarMaxLevels=107 \
    --remove bolt.in_plink_but_not_imputed.FID_IID.120.txt \
    --noBgenIDcheck
    
    
    
  BOLT-LMM_v2.3.5/bolt \
    --bed=pca/UKB_ptsd_eur_related_{1..22}.bed \
    --bim=pca/UKB_ptsd_eur_related_{1..22}.bim \
    --fam=pca/UKB_ptsd_eur_related_1.fam \
    --LDscoresFile=BOLT-LMM_v2.3.5/tables/LDSCORE.1000G_EUR.tab.gz \
    --geneticMapFile=BOLT-LMM_v2.3.5/tables/genetic_map_hg19_withX.txt.gz \
    --numThreads=6 \
    --bgenFile=/mnt/ukbb/adam/bgen/ukb_imp_chr{1..22}_v3.bgen \
    --bgenMinMAF=1e-2 \
    --bgenMinINFO=0.8 \
    --sampleFile=/mnt/ukbb/adam/ptsd/preimputation_genotypes/ukb41209_imp_chr1_v3_s487320.sample \
    --verboseStats \
    --lmmForceNonInf \
    --statsFile=pts_ukbb_may13_2021_unrelated_males.stats.gz \
    --statsFileBgenSnps=pts_ukbb_may13_2021_unrelated_males.bgen.stats.gz \
    --phenoFile=UKB_ptsd_eur_unrelated_m0_may13_2021_males.pheno \
    --phenoCol=PCL_sum_imp  \
    --covarFile=UKB_ptsd_eur_unrelated_m0_may13_2021_males.pheno \
    --covarCol=f.22000.0.0  --covarCol=f.54.0.0  --qCovarCol=PC{1:10}   \
    --covarMaxLevels=107 \
    --remove bolt.in_plink_but_not_imputed.FID_IID.969.txt \
    --noBgenIDcheck
    
awk '{if (NR ==1 || $20 == "1") print }' UKB_ptsd_eur_unrelated_m0_may13_2021.pheno > UKB_ptsd_eur_unrelated_m0_may13_2021_males.pheno
awk '{if (NR ==1 || $20 == "0") print }' UKB_ptsd_eur_unrelated_m0_may13_2021.pheno > UKB_ptsd_eur_unrelated_m0_may13_2021_females.pheno


  BOLT-LMM_v2.3.5/bolt \
    --bed=pca/UKB_ptsd_eur_related_{1..22}.bed \
    --bim=pca/UKB_ptsd_eur_related_{1..22}.bim \
    --fam=pca/UKB_ptsd_eur_related_1.fam \
    --LDscoresFile=BOLT-LMM_v2.3.5/tables/LDSCORE.1000G_EUR.tab.gz \
    --geneticMapFile=BOLT-LMM_v2.3.5/tables/genetic_map_hg19_withX.txt.gz \
    --numThreads=6 \
    --bgenFile=/mnt/ukbb/adam/bgen/ukb_imp_chr{10..22}_v3.bgen \
    --bgenMinMAF=1e-2 \
    --bgenMinINFO=0.8 \
    --sampleFile=/mnt/ukbb/adam/ptsd/preimputation_genotypes/ukb41209_imp_chr1_v3_s487320.sample \
    --verboseStats \
    --lmmForceNonInf \
    --statsFile=pts_ukbb_may13_2021_unrelated_females.stats.gz \
    --statsFileBgenSnps=pts_ukbb_may13_2021_unrelated_females2.bgen.stats.gz \
    --phenoFile=UKB_ptsd_eur_unrelated_m0_may13_2021_females.pheno \
    --phenoCol=PCL_sum_imp  \
    --covarFile=UKB_ptsd_eur_unrelated_m0_may13_2021_females.pheno \
    --covarCol=f.22000.0.0  --covarCol=f.54.0.0  --qCovarCol=PC{1:10}   \
    --covarMaxLevels=107 \
    --remove bolt.in_plink_but_not_imputed.FID_IID.969.txt \
    --noBgenIDcheck
    
    cat <(zcat pts_ukbb_may13_2021_unrelated_females.bgen.stats.gz | awk '{if (NR==1 || $2 <= 9) print}' ) <(zcat pts_ukbb_may13_2021_unrelated_females2.bgen.stats.gz | tail -n+2) | gzip > pts_ukbb_may13_2021_unrelated_females.bgen.stats.gz.complete.gz


    #as a case/control trait
  BOLT-LMM_v2.3.5/bolt \
    --bed=pca/UKB_ptsd_eur_related_{1..22}.bed \
    --bim=pca/UKB_ptsd_eur_related_{1..22}.bim \
    --fam=pca/UKB_ptsd_eur_related_1.fam \
    --LDscoresFile=BOLT-LMM_v2.3.5/tables/LDSCORE.1000G_EUR.tab.gz \
    --geneticMapFile=BOLT-LMM_v2.3.5/tables/genetic_map_hg19_withX.txt.gz \
    --numThreads=6 \
    --bgenFile=/mnt/ukbb/adam/bgen/ukb_imp_chr{1..22}_v3.bgen \
    --bgenMinMAF=1e-2 \
    --bgenMinINFO=0.8 \
    --sampleFile=/mnt/ukbb/adam/ptsd/preimputation_genotypes/ukb41209_imp_chr1_v3_s487320.sample \
    --verboseStats \
    --lmmForceNonInf \
    --statsFile=pcld_ukbb_may13_2021_unrelated.stats.gz \
    --statsFileBgenSnps=pcldx_ukbb_may13_2021_unrelated.bgen.stats.gz \
    --phenoFile=UKB_ptsd_eur_unrelated_m0_may13_2021.pheno \
    --phenoCol=PCL_dx  \
    --covarFile=UKB_ptsd_eur_unrelated_m0_may13_2021.pheno \
    --covarCol=f.22000.0.0  --covarCol=f.54.0.0  --qCovarCol=PC{1:10}   \
    --covarMaxLevels=107 \
    --remove bolt.in_plink_but_not_imputed.FID_IID.969.txt \
    --noBgenIDcheck
    
    dm <- fread('UKB_ptsd_eur_unrelated_m0_may13_2021.pheno',data.table=F)
    table(dm$PCL_dx)
    


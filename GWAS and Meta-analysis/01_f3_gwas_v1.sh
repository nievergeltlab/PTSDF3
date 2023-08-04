
### 1)Study Level Analysis steps: 

##Group 1 and 2 GWAS only: GWAS

cov=pcs
sex=all # females males
working_dir=$(pwd)

#may need to do gtpc as well
#Make sure time codes are correct
IFS=$'\n'
for line in $(cat dosage_locations_f3_aam_cont.csv  | grep saf2)
do
 study_1=$(echo $line | awk 'BEGIN{FS=","}  {print $1}')
 study_2=$(echo $line | awk 'BEGIN{FS=","}  {print $2}')

 ancgroup=$(echo $line | awk 'BEGIN{FS=","} {print $3}')
 timecode=$(echo $line | awk 'BEGIN{FS=","} {print $4}')
 exclude=$(echo $line | awk 'BEGIN{FS=","}  {print $5}')

 
#analysis of phenotypes coded into .pheno files

 #if [ $exclude == "0" ] 
 #then
  echo gwas for $study_1 $study_2 $timecode 
  sbatch --time=$timecode --error errandout/"$study_1"_"$study_2"_"$cov"_"$ancgroup"_"$sex".e --output errandout/"$study_1"_"$study_2"_"$cov"_"$ancgroup"_"$sex".o  \
   --export=ALL,study="$study_1",cov="$cov",study_2="$study_2",ancgroup="$ancgroup",sex="$sex" run_trauma_gwas_v2_freeze3.slurm -D $working_dir 
# fi


#Case/Control analysis of built in phenotypes (For Grotzinger) - beware, phenotypes in files may not be right  or at least not optimal - such as for army starrs data!
 # if [ $exclude == "1" ] 
 # then
  # echo gwas for $study_1 $study_2 $timecode 
  # sbatch --time=$timecode --error errandout/"$study_1"_"$study_2"_"$cov"_"$ancgroup"_"$sex".e --output errandout/"$study_1"_"$study_2"_"$cov"_"$ancgroup"_"$sex".o  \
   # --export=ALL,study="$study_1",cov="$cov",study_2="$study_2",ancgroup="$ancgroup",sex="$sex" run_trauma_gwas_v2_freeze3_case_control.slurm -D $working_dir 
 # fi

done


##Group 3 and 4 and 5 only: Summary stat conversion from genome-wide to split by chromosome
 #Caution: these data may need to be reformatted first for meta-analysis. Do this PRIOR to running this script
 
 #To reformat by chromosome, supply a list .gz files. Each row should be a file name followed by the column number for the chromosome
 IFS=$'\n'
 for line in $(cat sumstat_reformat.csv  )
 do
  infile=$(echo $line | awk ' {print $1}')
  colno=$(echo $line | awk '{print $2}')
  freqcol=$(echo $line | awk ' {print $3}')
  echo $infile $colno $freqcol
  sbatch --time=00:25:00 --error errandout/"$infile"_reformat.e --output errandout/"$infile"_reformat.o  \
   --export=ALL,colno="$colno",freqcol="$freqcol",infile="$infile" reformat.slurm -D $working_dir 
 done
 
 #Some X chromosome data needs to be converted, from showing 23 to showing X
 for files in ptsd_qt_vetsa_may12_2021_related_filtered.imputed.stats.gz_23.gz pts_ukbb_may13_2021_unrelated_chrX.bgen.stats.gz # daner_iPSYCH2015_PTSDbroad_chrX_HRC_MAF01.gz PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_With_X_WG.txt_Broad.regenie.gz_23.gz PTSD_broad_EstBB_GWAS_saige_chr23.txt.noinfo.gz
 do
 # zcat sumstats/bychr/$files | sed 's/^23[[:blank:]]/X /g' | gzip > sumstats/bychr/"$files".23.gz
    zcat sumstats/bychr/$files | sed 's/[[:blank:]]23[[:blank:]]/ X /g' | gzip > sumstats/bychr/"$files".23.gz
  
  zcat sumstats/bychr/"$files".23.gz | head
 done

### 2)Meta analysis

#User: Make a meta-analysis script. Follow format of eur_ptsd_m0.mi
#User: Studies should be weighted by effective sample size. Weights are set in the meta-script
#00_weight_checkr.r contains code to take input phenotypes and estimate a weighting factor


#Adjust meta-analysis script for each chromosome (replaces "XXX" with a chromosome number)
 for chr in {1..22} #X
 do 
  #sed s/XXX/ f3_gwas_may13_2021_aam_nodnhs.mi
   # sed s/XXX/$chr.gz/g f3_gwas_may13_2021_aam_nodnhs.mi  >   metal_scripts/f3_gwas_may13_2021_aam_nodnhs.mi_$chr

  #Basic analysis of all data

  sed s/XXX/$chr.gz/g f3_gwas_may13_2021_no_mdd_overlap.mi  >   metal_scripts/f3_gwas_may13_2021_no_mdd_overlap.mi_$chr 
 #   sed s/XXX/$chr.gz/g f3_gwas_may13_2021_ehr_only.mi  >   metal_scripts/f3_gwas_may13_2021_ehr_only.mi_$chr
#  sed s/XXX/$chr.gz/g f3_gwas_may13_2021_pgc_only.mi  >   metal_scripts/f3_gwas_may13_2021_pgc_only.mi_$chr
  # sed s/XXX/$chr.gz/g f3_gwas_may13_2021_nomvp.mi  >   metal_scripts/f3_gwas_may13_2021_nomvp.mi_$chr
e



#List all chromosome meta-analysis files into metafilelist.txt file
 ls  metal_scripts/Rf3_gwas_may13_2021.mi_* > metafilelist.txt
 
#User: Give a name for the error log files
 dataset=PTSD_25

 sbatch -t 01:50:00  --error errandout/"$dataset".e --output errandout/"$dataset".o   --export=ALL,metafile=metafilelist.txt -D /home/maihofer/freeze3_gwas run_meta_v2_loo_v2.slurm
 
#Concatenate meta-analysis results


#Format meta results for fuma #Note: COLUMNS WILL CHANGE if ANALYZE HET IS ON!

#User: find out what the max N is make sure only decently covered loci are returned (by default, 90% of total Neff)
#TBD: X chromosome needs its own cutoffs

 #Main analysis
  #european #added R to not delete original file, i'm getting out the het values for plots
   cat metal_results/eur_ptsd_pcs_v4_aug3_2021_[0-9].gz1.tbl metal_results/eur_ptsd_pcs_v4_aug3_2021_[0-9][0-9].gz1.tbl metal_results/eur_ptsd_pcs_v4_aug3_2021_chrX.gz1.tbl | awk '{if (NR == 1 || ($3 != "MarkerName" )) print}' > metal_results/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl
   totalN=641553
   totalNX=424500.00
   percentN=0.8
   cat metal_results/eur_ptsd_pcs_v4_aug3_2021_*.gz1.tbl   |  awk  -v totalN=$totalN -v totalNX=$totalNX -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && (($1 != "X" && $10 >= percentN*totalN )||( $1 == "X" && $10 >= percentN*totalNX ))) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/eur_ptsd_pcs_v4_aug3_2021.fuma.gz

#meta analyis with het (positions off. pull positions in from other file by merge)

   cat metal_results/Reur_ptsd_pcs_v4_aug3_2021_*.gz1.tbl   |  awk  '{ if (NR ==1 || ($4 >= 0.01 && $4 <= 0.99 && $1 != "MarkerName" )) print}'   | grep -v : | sort -g -k 10 | gzip > results_filtered/Reur_ptsd_pcs_v4_aug3_2021.allchr.het.fuma.gz


#Test: low cutoff
   totalN=641553
   percentN=0.25
  cat metal_results/eur_ptsd_pcs_v4_aug3_2021_[0-9].gz1.tbl metal_results/eur_ptsd_pcs_v4_aug3_2021_[0-9][0-9].gz1.tbl   |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/eur_ptsd_pcs_v4_aug3_2021_LOWCUTOFF.fuma.gz

  #Add the X (if starting from scratch, beware the previous *, which may have already added the X. In that case, you have to code an if statement
   totalN=424500.00
   percentN=0.8
   cat metal_results/eur_ptsd_pcs_v4_aug3_2021_chrX.gz1.tbl   |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/eur_ptsd_pcs_v4_aug3_2021.X.fuma.gz

   zcat results_filtered/eur_ptsd_pcs_v4_aug3_2021.fuma.gz results_filtered/eur_ptsd_pcs_v4_aug3_2021.X.fuma.gz | sort -g -k 9 | awk '{if (NR == 1 || $1 != "Chromosome") print}' | gzip > results_filtered/eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz

  #AAM
   cat metal_results/aam_ptsd_pcs_v5_jan4_2021_*.gz1.tbl | awk '{if (NR == 1 || ($3 != "MarkerName" )) print}' > metal_results/aam_ptsd_pcs_v5_jan4_2022.gz1.tbl
   totalN=42804
   percentN=0.8
   cat metal_results/aam_ptsd_pcs_v5_jan4_2022.gz1.tbl  |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/aam_ptsd_pcs_v5_jan4_2022.fuma.gz
  
   totalN=13246.00
   percentN=0.8
   cat metal_results/aam_ptsd_pcs_v5_jan4_2021_X.gz1.tbl  |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/aam_ptsd_pcs_v5_jan4_2021_X.fuma.gz

   zcat results_filtered/aam_ptsd_pcs_v5_jan4_2022.fuma.gz results_filtered/aam_ptsd_pcs_v5_jan4_2021_X.fuma.gz | sort -g -k 9 | awk '{if (NR == 1 || $1 != "Chromosome") print}' | gzip > results_filtered/aam_ptsd_pcs_v5_jan4_2022.allchr.fuma.gz
  
   #LAT
   cat metal_results/hna_ptsd_pcs_v4_aug3_2021_*.gz1.tbl | awk '{if (NR == 1 || ($3 != "MarkerName" )) print}' > metal_results/hna_ptsd_pcs_v4_aug3_2021.gz1.tbl
   totalN=6530
   percentN=0.8
   cat metal_results/hna_ptsd_pcs_v4_aug3_2021.gz1.tbl  |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/hna_ptsd_pcs_v4_aug3_2021.fuma.gz
  
   totalN=6530
   percentN=0.8
   cat metal_results/hna_ptsd_pcs_v4_aug3_2021_X.gz1.tbl  |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/hna_ptsd_pcs_v4_aug3_2021_X.fuma.gz
  
   zcat results_filtered/hna_ptsd_pcs_v4_aug3_2021.fuma.gz results_filtered/hna_ptsd_pcs_v4_aug3_2021_X.fuma.gz | sort -g -k 9 | awk '{if (NR == 1 || $1 != "Chromosome") print}' | gzip > results_filtered/hna_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz
    
   
   
   #TRANS
   totalN=$((641553+42804+6530)
   percentN=0.8
   cat metal_results/trans_ptsd_pcs_v5_jan4_20221.tbl  |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/trans_ptsd_pcs_v5_jan4_2022.fuma.gz
   
    
   
  
###LDSC analysis
 conda activate ldsc

 #Filter data down to just LDSC SNPs then munge
 LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' ~/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(cat metal_results/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl | awk '{if (NR == 1) $3="SNP"; print}'   | LC_ALL=C sort -k3b,3 ) | gzip > results_filtered/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.premunge.gz 
 python2 ~/trauma_gwas/ldsc-master/munge_sumstats.py --sumstats results_filtered/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.premunge.gz  --N-col Weight --out results_filtered/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.gz #add --N-col OBS_CT for the sample size
 
 #Get format for LDhub
 zcat results_filtered/eur_ptsd_pcs_v4_aug3_2021.fuma.gz | awk '{if (NR==1) print "CHR","BP","SNP","A1","A2","FRQ","N","Z","P"; else print $0}' > results_filtered/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.txt
 cd results_filtered
 zip eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.txt.zip eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.txt
 
 #Run LDSC
  python2 /home/maihofer/trauma_gwas/ldsc-master/ldsc.py \
 --h2 results_filtered/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.gz.sumstats.gz \
 --ref-ld-chr /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --w-ld-chr  /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --out results_filtered/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.gz.tbl.ldsc


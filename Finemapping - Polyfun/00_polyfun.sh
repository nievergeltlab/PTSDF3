
#install polyfun
 git clone https://github.com/omerwe/polyfun
 cd polyfun
 sudo /usr/local/bin/anaconda3/condabin/conda env create -f polyfun.yml
 python test_polyfun.py #seems to run

#Get existing functional annotations
 #wget https://data.broadinstitute.org/alkesgroup/LDSCORE/baselineLF_v2.2.UKB.polyfun.tar.gz
 
#download an LD matrix
mkdir  LD_cache
mkdir  output
mkdir -p LD_temp
cd LD_temp
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr1_46000001_49000001.npz
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr1_46000001_49000001.gz
cd ..

#Start polyfun
conda activate polyfun
 
 
#Munge Sumstats

#Premunge - capitalize A1 and A2, name MAF and N columns, change X to 23.
#zcat eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz | awk '{if(NR==1) { $3="SNP"; $4="A1";$5="A2";$6="MAF";$7="N"; $8="Z";$9="P"}; if (NR>1) {$4=toupper($4); $5=toupper($5)}; gsub("X",23); print}' > eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.premunge

#The X chromosome is not present in polyfun.. just remove it.
zcat eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz | awk '{if(NR==1) { $3="SNP"; $4="A1";$5="A2";$6="MAF";$7="N"; $8="Z";$9="P"}; if (NR>1) {$4=toupper($4); $5=toupper($5)};  print}' | grep -v X > eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.premunge

python /home/genetics/polyfun/munge_polyfun_sumstats.py --sumstats eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.premunge --out eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.polyfun

#Calculate priors - using approach 1 of pre-computed prior stuff
 python /home/genetics/polyfun/extract_snpvar.py --sumstats eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.polyfun --out eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.polyfun_prior
 
#Get the LD files

cd /home/storage/polyfun_ld
for files in $(cat /home/adam/ptsdf3_polyfun/ldfiles_require.txt)
do
 if [ ! -f "$files".gz ]
 then
  wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/"$files".gz
  wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/"$files".npz
 fi
done



#We should set max-num causal to 1, per author recommendations as we don't ahve good LD, but this generates NULL results! 
dos2unix finemaplist_withrefs.csv

#Need SNPVAR column
IFS=$'\n'
for snpset in $(cat finemaplist_withrefs.csv )
do

snp=$(echo $snpset | awk 'BEGIN{FS=","}{print $1}')
chr=$(echo $snpset | awk 'BEGIN{FS=","}{print $2}')
start=$(echo $snpset | awk 'BEGIN{FS=","}{print $3}')
stop=$(echo $snpset | awk 'BEGIN{FS=","}{print $4}')
fileld=$(echo $snpset | awk 'BEGIN{FS=","}{print $5}')


  wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/"$fileld".gz
  wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/"$fileld".npz
  
    python /home/genetics/polyfun/finemapper.py \
    --ld $fileld \
    --sumstats eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.polyfun_prior \
    --chr $chr \
    --n 641533 \
    --start  $start 	  \
    --end $stop   \
    --method susie \
    --max-num-causal 2 \
    --out output/finemapLD.inform.EUR."$snp"."$chr"."$start"."$stop".gz
      rm "$fileld".gz "$fileld".npz
done


IFS=$'\n'
for snpset in $(cat finemaplist_withrefs.csv )
do

snp=$(echo $snpset | awk 'BEGIN{FS=","}{print $1}')
chr=$(echo $snpset | awk 'BEGIN{FS=","}{print $2}')
start=$(echo $snpset | awk 'BEGIN{FS=","}{print $3}')
stop=$(echo $snpset | awk 'BEGIN{FS=","}{print $4}')
fileld=$(echo $snpset | awk 'BEGIN{FS=","}{print $5}')

#Get markers in credible set

zcat output/finemapLD.inform.EUR."$snp"."$chr"."$start"."$stop".gz | awk -v leadsnp=$snp '{if (NR==1 || $15!=0) print leadsnp, $0}' > output_filtered/finemapLD.inform.EUR."$snp"."$chr"."$start"."$stop".credible

zcat output/finemapLD.inform.EUR."$snp"."$chr"."$start"."$stop".gz | awk -v leadsnp=$snp '{print leadsnp, $0}' > output_filtered/finemapLD.inform.EUR."$snp"."$chr"."$start"."$stop".allsnps

#calculate credible set widths
done




#
cat output_filtered/*.credible  | awk '{if(NR==1){$1="Locus"}; if(NR==1 || $2 != "CHR") print}' > finemapLD.inform.EUR.complete.txt

cat output_filtered/*.allsnps  | awk '{if(NR==1){$1="Locus"}; if(NR==1 || $2 != "CHR") print}' > finemapLD.inform.EUR.complete.allsnps.txt

    python /home/genetics/polyfun/finemapper.py \
    --sumstats eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.polyfun_prior \
    --chr $chr \
    --n 641533 \
    --start  $start 	  \
    --end $stop   \
    --method susie \
    --max-num-causal 1 \
    --out output/finemapNOLD.inform.EUR."$snp"."$chr"."$start"."$stop".gz
    
  
    for files in $(ls )
    do 
    newname=$(echo $files | sed 's/\r//g')
    mv $files $newname
    done
    
    for files in $(ls | grep gz$)
    do
     zcat $files | awk '{if(NR==1 || $14==1)print}' 
     done
     
    
    zcat output/finemap.inform.EUR.rs78201023.1.35664657.36375226.gz| awk '{if(NR==1 || $14==1)print}' 
    
    #Tet locus
    1	 38198744 	 38459210 


    #Without the reference panel we find NOTHING
    
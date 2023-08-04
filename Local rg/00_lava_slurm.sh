#!/bin/bash

# module load 2021
# module load R/4.1.0-foss-2021a
#split -a3 --numeric-suffixes=100 -l 10 blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile

#echo "LOC CHR START STOP" > lava_loci_header.txt
cd  /home/maihofer/freeze3_gwas/lava

cp ptsdf3_mddf2_lava.R  "$TMPDIR"/.

i=1
IFS=$'\n'
for lineno in $(cat jobfiles/x${SLURM_ARRAY_TASK_ID} )
do
 echo $lineno | cat lava_loci_header.txt - > temp/lava_js_x${SLURM_ARRAY_TASK_ID}_$i.txt
 
  /sara/eb/AVX2/Debian10/EB_production/2021/software/R/4.1.0-foss-2021a/bin/Rscript lavaslurm.r lava_js_x${SLURM_ARRAY_TASK_ID}_$i.txt "$TMPDIR"/ptsd_mddf2_lava.r  &
 
i=$((i+1))
 
done

wait

#cat lava_loci_header.txt  jobfiles/x${SLURM_ARRAY_TASK_ID} > temp/lava_js_x${SLURM_ARRAY_TASK_ID}.txt


#  /sara/eb/AVX2/Debian10/EB_production/2021/software/R/4.1.0-foss-2021a/bin/Rscript lavaslurm.r lava_js_x${SLURM_ARRAY_TASK_ID}.txt "$TMPDIR"/ptsdf3_mddf3_lava.R  
 
#  sbatch --array=101-349 --time=0:17:00 --error errandout/ptsd_mdd_lava_%a.e --output errandout/ptsd_mdd_lava_%a.o --export=ALL 00_lava_slurm.sh

  

cat output/ptsd_mdd2019_rg_gw_lava_js_*  | sed 's/\"//g' | sort -g -k 15 | awk '{if (NR==1 || $15 != "p") print}' > ptsd_mdd2019.allrg


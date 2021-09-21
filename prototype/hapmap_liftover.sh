#!/bin/bash

module load python/2.7
module load plink


plink --bfile /data/LNG/hampton_temp/liftover_files/hapmap/HAPMAP_hg19_new --recode --out /data/LNG/hampton_temp/liftover_files/hapmap/hapmap_recode

python /data/LNG/hampton_temp/liftover_files/liftOverPlink.py --map /data/LNG/hampton_temp/liftover_files/hapmap/hapmap_recode.map --out hapmap_lifted --chain /data/LNG/hampton_temp/liftover_files/hg19ToHg38.over.chain.gz --bin /data/LNG/hampton_temp/liftover_files/liftOver

python /data/LNG/hampton_temp/liftover_files/rmBadLifts.py --map /data/LNG/hampton_temp/liftover_files/hapmap/hapmap_lifted.map --out /data/LNG/hampton_temp/liftover_files/hapmap/good_lifted.map --log /data/LNG/hampton_temp/liftover_files/hapmap/bad_lifted.dat

module load python/3.7

#this script needs to be where your good_lifted.map file is
python /data/LNG/hampton_temp/liftover_files/hapmap/good_lifted_snp_list.py

plink --bfile /data/LNG/hampton_temp/liftover_files/hapmap/HAPMAP_hg19_new --extract /data/LNG/hampton_temp/liftover_files/hapmap/snplist.txt --recode --out /data/LNG/hampton_temp/liftover_files/hapmap/lifted

plink --ped lifted.ped --map good_lifted.map --make-bed --out hapmap_liftedhg38
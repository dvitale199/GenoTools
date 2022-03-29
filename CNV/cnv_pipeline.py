from QC.utils import shell_do


'''This pipeline will need to run on each unique sentrix barcode and then split to individual samples later'''

    

def idat_to_cnv(idat_path, bpm, bpm_csv, egt, ref_fasta, out_path):
    '''
    current structure of idat storage is such that a directory of each SentrixBarcode_A with all idats for that barcode in it
    for ex.
    1112223334
        --> 1112223334_R01C01_Red.idat
        --> 1112223334_R01C01_Grn.idat
        --> 1112223334_R01C02_Red.idat
        --> 1112223334_R01C02_Grn.idat
        etc.
        
    '''
    
    out_tmp = f'{out_path}/tmp'
    os.mkdir(out_tmp, mode)
    
    
    idat_to_gtc_cmd = f'\
{iaap} gencall \
{bpm} \
{egt} \
{out_path} \
-f {idat_path} \
-g \
-t 8'   
    
    # export path to plugins temporarily for biowulf. will figure this out later
    gtc2vcf_cmd = f'\
export BCFTOOLS_PLUGINS="/data/vitaled2/bin"; bcftools +gtc2vcf \
--no-version -Ob \
--bpm {bpm} \
--csv {bpm_csv} \
--egt {egt} \
--gtcs {out_tmp} \
--fasta-ref {ref_fasta} \
bcftools norm --no-version -Oz -c w -f {grch38_fasta} > {vcf_path}/gp2_snps_{code}.vcf.gz'

    # --extra {vcf_path}/gp2_snps_{code}_metadata.tsv | \


# sort vcf
tmp_path = {vcf_path}/tmp
os.mkdir(tmp_path, mode)
with open(f'{swarm_scripts_dir}/bcftools_sort.swarm', 'w') as f:
    for code in key.SentrixBarcode_A.unique():
        sort_cmd = f'cd {vcf_path}; bcftools sort gp2_snps_{code}.vcf.gz -T ./tmp -Oz -o gp2_snps_sorted_{code}.vcf.gz'
        f.write(f'{sort_cmd}\n')
f.close()


# split indels and snps in vcf
with open(f'{swarm_scripts_dir}/vcftools_remove_indels.swarm', 'w') as f:
    for code in key.SentrixBarcode_A.unique():
        ext_snps_cmd = f'cd {vcf_path}; vcftools --gzvcf gp2_snps_sorted_{code}.vcf.gz --remove-indels --recode --recode-INFO-all --out gp2_snps_sorted_snps_only_{code}'
        f.write(f'{ext_snps_cmd}\n')
f.close()


# split indels and snps in vcf
with open(f'{swarm_scripts_dir}/vcftools_keep_indels.swarm', 'w') as f:
    for code in key.SentrixBarcode_A.unique():
        keep_indels_cmd = f'cd {vcf_path}; vcftools --gzvcf gp2_snps_sorted_{code}.vcf.gz --keep-only-indels --recode --recode-INFO-all --out gp2_snps_sorted_indels_only_{code}'
        f.write(f'{keep_indels_cmd}\n')
f.close()


# get snp info from each vcf
with open(f'{swarm_scripts_dir}/get_logr_baf.swarm', 'w') as f:
    for code in key.SentrixBarcode_A.unique():
        get_logr_baf = f'python3 process_vcf_snps.py --vcf {vcf_path}/gp2_snps_sorted_snps_only_{code}.recode.vcf --gene_ref /data/CARD/PD/GP2/ref_panel/glist-hg38 --out {vcf_path}/gp2_snp_metrics_{code}.txt'
        f.write(f'{get_logr_baf}\n')
f.close()


# get snp info from each vcf
with open(f'{swarm_scripts_dir}/clean_snp_metrics.swarm', 'w') as f:
    for code in key.SentrixBarcode_A.unique():
        clean_snp_metrics = f'python3 clean_snp_metrics.py --infile {vcf_path}/gp2_snp_metrics_{code}.txt --outpath {vcf_path}/gp2_snp_metrics'
        f.write(f'{clean_snp_metrics}\n')
f.close()
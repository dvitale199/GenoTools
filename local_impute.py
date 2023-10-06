import pandas as pd
from genotools.impute import *
from genotools.utils import shell_do


project_path = '/data/GP2/projects/gp2_imputation'
ref_path = f'{project_path}/refs'
eagle_map = f'{ref_path}/genetic_map_hg38_withX.txt.gz'

minimac_path = '/data/vitaled2/GP2_data_processing/bin/minimac4'
harmonizer_path = '/data/vitaled2/GP2_data_processing/GenotypeHarmonizer-1.4.25-SNAPSHOT/GenotypeHarmonizer.jar'
geno_in = f'{project_path}/raw_genotypes/EUR/EUR_maf_hwe_release5'
eagle_path = '/data/vitaled2/GP2_data_processing/bin/Eagle_v2.4.1/eagle'


chroms = [str(x) for x in range(1, 23)] + ["X.PAR1","X.nonPAR","X.PAR2"]

# will loop through these after test
chrom = '1'
label = 'EUR'

# create paths for harmonized and phased outputs
harm_dir = f'{project_path}/harmonized/{label}'
harm_geno_in = f"{project_path}/raw_genotypes/EUR/{label}_maf_hwe_release5_chr{chrom}"
harm_geno_out = f"{harm_dir}/{label}_chr{chrom}_harmonized"
chunk_geno_out = f"{harm_geno_out}_chunk"


phased_dir = f'{project_path}/phased/{label}'
eagle_geno_out = f'{phased_dir}/{label}_chr{chrom}_phased'

harmonizer_ref = f'{ref_path}/vcfs/ALL.chr{chrom}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'
eagle_ref = f'{ref_path}/bcfs/1kg_GRCh38_chr{chrom}.bcf' 
eagle_map = f'{ref_path}/genetic_map_hg38_withX.txt.gz'
impute_ref = f'{ref_path}/mvcfs/1kg_GRCh38_chr{chrom}.msav'

# run harmonizer
print('running harmonizer')
# harmonize(f'{harm_geno_in}', harm_geno_out, harmonizer_ref, harmonizer_path)
print('harmonizer complete')
print()

print('chunking genotypes')
chunk_output = chunk_genotypes(harm_geno_out, chunk_geno_out, chunk_size=20000000, chrom=chrom)
print('chunking complete')
print()

# eventually run this in parallel- either launch sbatch job for each chunk or use multiprocessing
print('phasing + imputing genotypes')
for chunk_start, chunk_end in chunk_output.keys():

    eagle_input = chunk_output[(chunk_start, chunk_end)]

    run_eagle(
        eagle_input, 
        eagle_geno_out, 
        eagle_ref, 
        eagle_map, 
        chrom, 
        chunk_start, 
        chunk_end, 
        overlap=5000000, 
        eagle_path=eagle_path
        )

    run_minimac4(
        geno_in=eagle_geno_out, 
        geno_out=f'{eagle_geno_out}_imputed', 
        ref=impute_ref, 
        window=500000, 
        start=chunk_start, 
        end=chunk_end, 
        min_ratio=0.00001, 
        cpus=16, 
        out_format='GT,DS,HDS,GP,SD', 
        minimac_path=minimac_path
        )
        
print('phasing + imputation complete')
print()

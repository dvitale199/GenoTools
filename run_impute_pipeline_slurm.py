from genotools.impute import *
from genotools.impute import run_phasing_and_imputation
import argparse
from concurrent.futures import ProcessPoolExecutor
import os
import shutil

def main(args):
    chrom = args.chrom
    label = args.label
    memory = args.memory
    threads = args.threads
    project_path = args.project_path
    # ref_path = arg.ref_path

    # leave as-is for testing
    # project_path = '/data/jamesml/genotool_trials'
    ref_path = f'{project_path}/refs'
    #minimac_path = '/data/vitaled2/GP2_data_processing/bin/minimac4'
    #load minimac module
    harmonizer_path = '/data/jamesml/GP2_data_processing/GenotypeHarmonizer-1.4.25-SNAPSHOT/GenotypeHarmonizer.jar'
    eagle_path = '/data/jamesml/GP2_data_processing/bin/Eagle_v2.4.1/eagle'
   
   

    harm_dir = f'{project_path}/harmonized/{label}'
    if os.path.exists(harm_dir):
        shutil.rmtree(harm_dir)
    os.makedirs(harm_dir)
    harm_geno_in = f"{project_path}/raw_genotypes/{label}/{label}_maf_hwe_release5_chr{chrom}"
    harm_geno_out = f"{harm_dir}/{label}_chr{chrom}_harmonized"
    chunk_geno_out = f"{harm_geno_out}_chunk"


    phased_dir = f'{project_path}/phased/{label}'
    if os.path.exists(phased_dir):
        shutil.rmtree(phased_dir)
    os.makedirs(phased_dir)
    eagle_geno_out = f'{phased_dir}/{label}_chr{chrom}_phased'

    impute_dir = f'{project_path}/imputed/{label}'
    if os.path.exists(impute_dir):
        shutil.rmtree(impute_dir)
    os.makedirs(impute_dir)
    impute_out = f'{impute_dir}/{label}_chr{chrom}_imputed'

    harmonizer_ref = f'{ref_path}/vcfs/ALL.chr{chrom}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'
    eagle_ref = f'{ref_path}/bcfs/1kg_GRCh38_chr{chrom}.bcf' 
    eagle_map = f'{ref_path}/genetic_map_hg38_withX.txt.gz'
    impute_ref = f'{ref_path}/mvcfs/1kg_GRCh38_chr{chrom}.msav'


    # run harmonizer
    print('running harmonizer')
    harmonize(f'{harm_geno_in}', harm_geno_out, harmonizer_ref, harmonizer_path, memory=memory)
    print('harmonizer complete')
    print()

    print('chunking genotypes')
    #chunk_outs = chunk_genotypes(harm_geno_out, chunk_geno_out, chunk_size=20000000, chrom=chrom)
    print('chunking complete')
    print()


    # # Initialize ProcessPoolExecutor
    # with ProcessPoolExecutor() as executor:
    #     # Create a list of futures
    #     futures = []
    #     for chunk_start, chunk_end in chunk_outs.keys():
    #         eagle_input = chunk_outs[(chunk_start, chunk_end)]
    #         eagle_output = f'{eagle_geno_out}_{chunk_start}_{chunk_end}'
    #         impute_output = f'{impute_out}_{chunk_start}_{chunk_end}'
    #         future = executor.submit(
    #             run_phasing_and_imputation,
    #             geno_in=eagle_input,
    #             phase_out=eagle_output,
    #             impute_out=impute_output,
    #             phase_ref=eagle_ref,
    #             impute_ref=impute_ref,
    #             eagle_map=eagle_map,
    #             chrom=chrom,
    #             chunk_start=chunk_start,
    #             chunk_end=chunk_end,
    #             threads=threads,
    #             eagle_path=eagle_path,
    #             minimac_path=minimac_path
    #         )
    #         futures.append(future)

    #     # Wait for all futures to complete and collect the results
    #     for future in futures:
    #         future.result()  # This will raise an exception if the function encountered any errors


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Harmonizer, Phasing, and Imputation')
    
    parser.add_argument('--chrom', required=True, help='Chromosome number')
    parser.add_argument('--label', required=True, help='Population label')
    parser.add_argument('--memory', required=True, help='Memory to allocate to each job')
    parser.add_argument('--threads', required=True, help='Number of threads to use')
    parser.add_argument('--project_path', required=True, help='String file path to working directory')
    # parser.add_argument('--ref_path', required=True, help='')

    args = parser.parse_args()
    main(args)

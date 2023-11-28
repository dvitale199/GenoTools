import pandas as pd
import subprocess
import shutil
import os
import time
from concurrent.futures import ThreadPoolExecutor
from genotools.utils import shell_do
from genotools.dependencies import check_plink, check_plink2

plink_exec = check_plink()
plink2_exec = check_plink2()

def split_chroms(input, output, chrom, filetype='pgen'):
    """
    Splits genotype data by chromosome using PLINK2.

    Parameters:
    -----------
    input : str
        Path to the input genotype file (either pgen or bed format).
    output : str
        Prefix for the output file.
    chrom : int
        Chromosome number. Note: 23 refers to the X chromosome.
    filetype : str, optional
        The format of the genotype file. Either 'pgen' (default) or 'bed'.
    """

    # Base command for both pgen and bed
    cmd_base = "plink2 --max-alleles 2 --set-missing-var-ids @:#$1:$2 --allow-extra-chr --chr {} --make-bed --out {}_chr{}"

    # Determine file type
    if filetype == 'pgen':
        cmd_base = cmd_base.replace("plink2", f"plink2 --pfile {input}")
    elif filetype == 'bed':
        cmd_base = cmd_base.replace("plink2", f"plink2 --bfile {input}")

    # Determine chromosome
    if chrom == 23:
        cmd = cmd_base.format(f"{chrom},X", input, chrom)
    else:
        cmd = cmd_base.format(chrom, output, chrom)
    
    
    shell_do(cmd)


def add_chr(geno_in, geno_out):
    """
    Add 'chr' prefix to the chromosome column of the .bim file and copy the associated .bed and .fam files.

    Parameters:
    - geno_in: Path to the input .bim file without the '.bim' extension.
    - geno_out: Path to the output .bim file without the '.bim' extension.
    """
    # Read the .bim file
    geno = pd.read_csv(f'{geno_in}.bim', sep='\t', header=None, names=['chr', 'snp', 'cm', 'bp', 'a1', 'a2'])
    
    # Add 'chr' prefix to the chromosome column
    geno['chr'] = geno['chr'].apply(lambda x: f'chr{x}')
    
    # Save the modified .bim file
    geno.to_csv(f'{geno_out}.bim', sep='\t', header=None, index=False)
    
    # Copy the associated .bed and .fam files
    shutil.copy(f'{geno_in}.bed', f'{geno_out}.bed')
    shutil.copy(f'{geno_in}.fam', f'{geno_out}.fam')


def harmonize(geno_in, geno_out, ref_path, harmonizer_path, memory):

    """
    Harmonizes genotype data using genotype harmonizer: https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer.

    Parameters:
    -----------
    geno_in : str
        Path to the input genotype file (PLINK BED format).
    geno_out : str
        Path for the harmonized output genotype file.
    ref_path : str
        Path to the reference data (in VCF format).
    harmonizer_path : str
        Path to the harmonizer tool (JAR file).

    Notes:
    ------
    This function first adds chromosome info to the input genotype data and saves it to a temp file.
    It then runs the harmonizer tool to harmonize the input genotype data against the reference.
    Post harmonization, chromosome info is added again and the output is saved.
    Temporary files created during the process are removed at the end.

    The function assumes 'add_chr' is a predefined function used to append chromosome information.

    The overall run time of the function is printed at the end.
    """

    start_time = time.time()

    tmp1 = f'{geno_in}_tmp1'
    tmp2 = f'{geno_in}_tmp2'

    add_chr(geno_in, tmp1)

    # Run commands
    cmd = f"""\
java -Xmx{memory}g -jar {harmonizer_path} \
--keep \
--input {tmp1} \
--ref {ref_path} \
--inputType PLINK_BED \
--callRateFilter 0.90 \
--refType VCF \
--output {tmp2}
"""

    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    output, error = process.communicate()

    add_chr(tmp2, geno_out)

    # # Remove temporary files
    # for suffix in ['.bim', '.bed', '.fam']:
    #     os.remove(f'{tmp1}{suffix}')
    #     os.remove(f'{tmp2}{suffix}')

    print(f"Overall run time: {time.time() - start_time}")


def threaded_execute(cmds):
    print("Threaded_execute called")
    def run_cmd(cmd):
        subprocess.run(cmd, shell=True, check=True)

    with ThreadPoolExecutor() as executor:
        executor.map(run_cmd, cmds)


def chunk_genotypes(geno_in, geno_out, chrom, chunk_size=20000000):

    bim = pd.read_csv(f'{geno_in}.bim', sep='\s+', header=None, names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'], dtype={'BP': int})
    chr_end = bim.BP.max()

    cmds = []
    chunk_output = {} 

    for i in range(1, chr_end + 1, chunk_size):
        chunk_start = i
        chunk_end = min(i + chunk_size - 1, chr_end)

        recode_vcf_cmd = f"plink2 --bfile {geno_in} --chr {chrom} --from-bp {chunk_start} --to-bp {chunk_end} --export vcf-4.2 --output-chr chrM --set-missing-var-ids @:#\$1:\$2 --out {geno_out}_{chunk_start}_{chunk_end}"
        bcftools_sort_cmd = f'bcftools sort {geno_out}_{chunk_start}_{chunk_end}.vcf -Oz -o {geno_out}_{chunk_start}_{chunk_end}.vcf.gz'
        index_vcf_cmd = f'tabix -f -p vcf {geno_out}_{chunk_start}_{chunk_end}.vcf.gz'

        full_cmd = f'{recode_vcf_cmd} && {bcftools_sort_cmd} && {index_vcf_cmd} && rm {geno_out}_{chunk_start}_{chunk_end}.vcf'
        
        # Add the output filename to the chunk_output dictionary
        chunk_output[(chunk_start, chunk_end)] = f'{geno_out}_{chunk_start}_{chunk_end}.vcf.gz'
        
        cmds.append(full_cmd)

    # Run all commands in parallel
    print("Commands to run:", cmds)
    try:
        threaded_execute(cmds)
    except Exception as e:
        print("An exception occurred:", e)
    
    return chunk_output



# def chunk_genotypes(geno_in, geno_out, chrom, chunk_size=20000000):
#     """
#     Splits genotype data into chunks based on base pair positions and exports them as VCF files.

#     Parameters:
#     -----------
#     geno_in : str
#         Path to the input genotype file in PLINK BED format.
#     geno_out : str
#         Prefix for the output VCF files.
#     chunk_size : int
#         Size of each chunk in base pair units.
#     chrom : int
#         Chromosome number to be considered for chunking.
#     overlap : int, optional
#         Number of base pairs to overlap between consecutive chunks. Default is 5,000,000.

#     Returns:
#     --------
#     chunk_output : dict
#         Dictionary mapping chunk start and end positions to the corresponding output VCF filenames.

#     Notes:
#     ------
#     The function reads base pair positions from the provided genotype data, then divides the data 
#     into chunks of a specified size with an optional overlap. Each chunk is exported as a VCF file 
#     using PLINK2, sorted with bcftools, and indexed using tabix. 

#     The function assumes necessary tools like PLINK2, bcftools, and tabix are available in the 
#     environment where it's executed.
#     """

#     # Read the bim file
#     bim = pd.read_csv(f'{geno_in}.bim', sep='\s+', header=None, names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'], dtype={'BP': int})
#     chr_end = bim.BP.max()

#     cmds = []
#     chunk_output = {} 

#     for i in range(1, chr_end + 1, chunk_size):
#         chunk_start = i
#         chunk_end = min(i + chunk_size - 1, chr_end)

#         recode_vcf_cmd = f"plink2 --bfile {geno_in} --chr {chrom} --from-bp {chunk_start} --to-bp {chunk_end} --export vcf-4.2 --output-chr chrM --set-missing-var-ids @:#\$1:\$2 --out {geno_out}_{chunk_start}_{chunk_end}"
#         bcftools_sort_cmd = f'bcftools sort {geno_out}_{chunk_start}_{chunk_end}.vcf -Oz -o {geno_out}_{chunk_start}_{chunk_end}.vcf.gz'
#         index_vcf_cmd = f'tabix -f -p vcf {geno_out}_{chunk_start}_{chunk_end}.vcf.gz'
        
#         cmds.extend([recode_vcf_cmd, bcftools_sort_cmd, index_vcf_cmd])

#         # Add the output filename to the chunk_output dictionary
#         chunk_output[(chunk_start, chunk_end)] = f'{geno_out}_{chunk_start}_{chunk_end}.vcf.gz'
    
#     for cmd in cmds:
#         shell_do(cmd)
    
#     return chunk_output


def run_eagle(geno_in, geno_out, ref_path, map_path, chrom, chunk_start, chunk_end, overlap=5000000, threads=16, eagle_path=None):
    """
    Generate and run the eagle command based on the provided parameters.

    Parameters:
    - eagle_geno_in: Input path for the eagle command
    - chunk_start: Start of the chunk
    - chunk_end: End of the chunk
    - chrom: Chromosome number

    Returns:
    - eagle_cmd: Generated eagle command
    """
    start_time = time.time()
    if chunk_start == 1:
        bp_start = 1
    else:
        bp_start = chunk_start - overlap

    bp_end = chunk_end + overlap

    if eagle_path:
        eagle_cmd = (f'{eagle_path} --vcfRef {ref_path} '
                 f'--vcfTarget {geno_in}  --geneticMapFile {map_path} '
                 f'--outPrefix {geno_out} --chrom chr{chrom} --bpStart {bp_start} --bpEnd {bp_end} '
                 f'--allowRefAltSwap --Kpbwt=100000 --numThreads={threads}  --vcfOutFormat z')
    else:
        eagle_cmd = (f'eagle --vcfRef {ref_path} '
                 f'--vcfTarget {geno_in}  --geneticMapFile {map_path} '
                 f'--outPrefix {geno_out} --chrom chr{chrom} --bpStart {bp_start} --bpEnd {bp_end} '
                 f'--allowRefAltSwap --Kpbwt=100000 --numThreads={threads} --vcfOutFormat z')
    
    os.system(eagle_cmd)
    

    # Run the bcftools index command
    index_cmd = f'bcftools index {geno_out}.vcf.gz'
    os.system(index_cmd)
    
    print(f"Overall run time: {time.time() - start_time}")


def run_minimac4(geno_in, geno_out, ref, window, region, min_ratio, threads, out_format='GT,DS,HDS,GP,SD', minimac_path=None):

    start_time = time.time()

    if minimac_path:
        minimac = minimac_path
        
    else:
        minimac = 'minimac4'

    impute_cmd = (
        f"{minimac} "
        f"{ref} "
        f"{geno_in} "
        f"--output {geno_out}.sav "
        f"--overlap {window} "
        f"--region {region} "
        f"--min-ratio {min_ratio} "
        f"--format {out_format} "
        f"--threads {threads} "
        f"--all-typed-sites --empirical-output {geno_out}.empirical.sav ")

    os.system(impute_cmd)
    print(f"Overall run time: {time.time() - start_time}")

def run_minimac4(geno_in, geno_out, ref, overlap, region, min_ratio, threads, out_format='bcf', info_format='GT,DS,HDS,GP,SD', minimac_path=None):
    """
    Runs the Minimac4 program for genotype imputation.

    Parameters:
    - geno_in (str): Path to the input genotype file.
    - geno_out (str): Prefix for the output file names.
    - ref (str): Path to the reference panel.
    - overlap (bool): Whether to consider overlapping markers.
    - region (str): Specifies the region for imputation.
    - min_ratio (float): Minimum ratio to consider.
    - threads (int): Number of threads for parallel processing.
    - out_format (str, optional): Desired output format, defaults to 'bcf'.
    - info_format (str, optional): Information format for the output, defaults to various info types.
    - minimac_path (str, optional): Custom path to the Minimac4 executable. If not provided, assumes 'minimac4' is in PATH.

    Outputs:
    - Prints the overall runtime of the imputation.
    """
    
    start_time = time.time()

    minimac = minimac_path if minimac_path else 'minimac4'

    impute_cmd = (
        f"{minimac} "
        f"{ref} "
        f"{geno_in} "
        f"--output {geno_out}.sav "
        f"--overlap {overlap} "
        f"--region {region} "
        f"--min-ratio {min_ratio} "
        f"--format {info_format} "
        f"--output-format {out_format}"
        f"--threads {threads} "
        f"--all-typed-sites --empirical-output {geno_out}_empirical.txt "
    )

    os.system(impute_cmd)
    print(f"Overall run time: {time.time() - start_time}")


def run_phasing_and_imputation(geno_in, phase_out, impute_out, phase_ref, impute_ref, eagle_map, chrom, chunk_start, chunk_end, threads, eagle_path=None, minimac_path=None):

    run_eagle(
        geno_in,
        phase_out,
        phase_ref,
        eagle_map,
        chrom,
        chunk_start,
        chunk_end,
        overlap=5000000,
        threads=threads,
        eagle_path=eagle_path
    )
    
    run_minimac4(
        geno_in=f'{phase_out}.vcf.gz',
        geno_out=f'{impute_out}_imputed',
        ref=impute_ref,
        window=500000,
        region=f'chr{chrom}:{chunk_start}-{chunk_end}',
        min_ratio=0.00001,
        threads=threads,
        info_format='GT,DS,HDS,GP,SD',
        out_format='bcf',
        minimac_path=minimac_path
    )
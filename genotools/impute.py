from genotools.utils import shell_do
import pandas as pd
import subprocess
import shutil
import os
import time


def split_chroms(input, output, chrom, filetype='pgen'):

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
    
    cmd = load_plink_str + cmd
    
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


def harmonize(geno_in, geno_out, ref_path, harmonizer_path):
    start_time = time.time()

    tmp1 = f'{geno_in}_tmp1'
    tmp2 = f'{geno_in}_tmp2'

    add_chr(geno_in, tmp1)

    # Run commands
    cmd = f"""\
java -Xmx16g -jar {harmonizer_path} \
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

    # Remove temporary files
    for suffix in ['.bim', '.bed', '.fam']:
        os.remove(f'{tmp1}{suffix}')
        os.remove(f'{tmp2}{suffix}')


    print(f"Overall run time: {time.time() - start_time}")


def chunk_genotypes(geno_in, geno_out, chunk_size, chrom):
    # Read the bim file
    bim = pd.read_csv(f'{geno_in}.bim', sep='\s+', header=None, names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'], dtype={'BP': int})
    chr_start = bim.BP.min()
    chr_end = bim.BP.max()

    cmds = []
    chunk_output = {}  # Dictionary to store the output filenames indexed by (chunk_start, chunk_end)
    for i in range(chr_start, chr_end + 1, chunk_size):
        chunk_start = i
        chunk_end = min(i + chunk_size - 1, chr_end)
        recode_vcf_cmd = f"plink2 --bfile {geno_in} --chr {chrom} --from-bp {chunk_start} --to-bp {chunk_end} --export vcf-4.2 --output-chr chrM --set-missing-var-ids @:#\$1:\$2 --out {geno_out}_{chunk_start}_{chunk_end}"
        bcftools_sort_cmd = f'bcftools sort {geno_out}_{chunk_start}_{chunk_end}.vcf -Oz -o {geno_out}_{chunk_start}_{chunk_end}.vcf.gz'
        index_vcf_cmd = f'tabix -f -p vcf {geno_out}_{chunk_start}_{chunk_end}.vcf.gz'
        
        cmds.extend([recode_vcf_cmd, bcftools_sort_cmd, index_vcf_cmd])

        # Add the output filename to the chunk_output dictionary
        chunk_output[(chunk_start, chunk_end)] = f'{geno_out}_{chunk_start}_{chunk_end}.vcf.gz'
    
    for cmd in cmds:
        # shell_do(cmd)
        !{cmd}
    
    return chunk_output


def run_eagle(geno_in, geno_out, ref_path, map_path, chrom, chunk_start, chunk_end, overlap=5000000):
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

    eagle_cmd = (f'eagle --vcfRef {ref_path} '
                 f'--vcfTarget {geno_in}  --geneticMapFile {map_path} '
                 f'--outPrefix {geno_out} --chrom chr{chrom} --bpStart {bp_start} --bpEnd {bp_end} '
                 '--allowRefAltSwap --Kpbwt=100000 --numThreads=16  --vcfOutFormat z')
    
    os.system(eagle_cmd)
    

    # Run the bcftools index command
    index_cmd = 'bcftools index {geno_out}'
    os.system(index_cmd)
    
    print(f"Overall run time: {time.time() - start_time}")


def run_minimac4(geno_in, geno_out, ref, window, start, end, minRatio, cpus, out_format='GT,DS,HDS,GP,SD'):
    """
    Runs minimac4 for genotype imputation.
    
    Parameters:
    - refHaps: Path to reference haplotypes
    - haps: Path to the target haplotypes
    - window: Size of the imputation window
    - start: Start position
    - end: End position
    - prefix: Output prefix for the imputed files
    - minRatio: Minimum ratio
    - cpus: Number of CPUs to use
    - format: Output format for the imputed genotypes
    - allTypedSites: Flag for including all typed sites
    - meta: Flag for meta info
    
    Returns:
    - None
    """
    start_time = time.time()
    impute_cmd = (
        f"minimac4 "
        f"--haps {haps} "
        f"--prefix {out} "
        f"--refHaps {refHaps} "
        f"--window {window} "
        f"--start {start} --end {end} "
        f"--minRatio {minRatio} "
        f"--format {out_format} "
        f"--cpus {cpus} "
        f"--allTypedSites --meta"
        f"--noPhoneHome ")

    os.system(impute_cmd)
    print(f"Overall run time: {time.time() - start_time}")

    
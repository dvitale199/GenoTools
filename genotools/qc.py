import pandas as pd
import numpy as np
import os
import shutil
from genotools.utils import shell_do, count_file_lines, concat_logs, bfiles_to_pfiles
from genotools.dependencies import check_plink, check_plink2

plink_exec = check_plink()
plink2_exec = check_plink2()


class SampleQC:

    def __init__(self, geno_path=None, out_path=None):
        self.geno_path = geno_path
        self.out_path = out_path

    def run_callrate_prune(self, mind=0.02):

        """
        Execute call rate pruning on genotype data using PLINK.

        This function performs call rate pruning based on a user-specified threshold for 
        the missing genotype rate. If an individual's missing genotype rate surpasses the 
        given threshold, the individual will be excluded from the dataset.

        Parameters:
        - mind (float, optional): Threshold for permissible missing genotype rate. Default is set to 0.02, which corresponds to 2%.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('callrate_prune').
            * 'metrics': Metrics associated with the pruning, such as the 'outlier_count'.
            * 'output': Dictionary containing paths to the generated output files like 'pruned_samples' and 'plink_out'.
        """
        
        geno_path = self.geno_path
        out_path = self.out_path

        # Check that paths are set
        if geno_path is None or out_path is None:
            raise ValueError("Both geno_path and out_path must be set before calling this method.")

        # Check path validity
        if not os.path.exists(f'{geno_path}.bed'):
            raise FileNotFoundError(f"{geno_path} does not exist.")
        
        # Check type of mind
        if not isinstance(mind, (int, float)):
            raise TypeError("mind should be of type int or float.")
        
        # Check valid range for mind
        if mind < 0 or mind > 1:
            raise ValueError("mind should be between 0 and 1.")

        step = "callrate_prune"
        
        outliers_out = f'{out_path}.outliers'
        
        plink_cmd1 = f"{plink2_exec} --bfile {geno_path} --mind {mind} --make-bed --out {out_path}"

        shell_do(plink_cmd1)

        listOfFiles = [f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)
        
        if os.path.isfile(f'{out_path}.mindrem.id'):
            irem = pd.read_csv(f'{out_path}.mindrem.id', sep='\s+', header=None, names=['FID','IID'])
            irem.to_csv(outliers_out, sep='\t', header=True, index=False)

            outlier_count = sum(1 for line in open(f'{outliers_out}'))
            
        else:
            outlier_count = 0
            
        process_complete = True
        
        outfiles_dict = {
            'pruned_samples': f'{outliers_out}',
            'plink_out': f'{out_path}',
        }

        metrics_dict = {
            'outlier_count': outlier_count
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict



    def run_sex_prune(self, check_sex=[0.25,0.75]):

        """
        Execute sex-based pruning on genotype data using PLINK.
        
        Parameters:
        - check_sex (list of float, optional): Two values indicating the lower and upper bounds for sex discrepancy checks. Default values are [0.25, 0.75].
        
        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('sex_prune').
            * 'metrics': Metrics associated with the pruning, such as 'outlier_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """
        
        geno_path = self.geno_path
        out_path = self.out_path

        step = "sex_prune"

        # create filenames
        sex_tmp1 = f"{out_path}_tmp1"
        sex_tmp2 = f"{out_path}_tmp2"
        sex_fails = f"{out_path}.outliers"

        # check sex 2 methods
        plink_cmd1 = f"{plink_exec} --bfile {geno_path} --check-sex 0.25 0.75 --maf 0.05 --out {sex_tmp1}"
        plink_cmd2 = f"{plink_exec} --bfile {geno_path} --chr 23 --from-bp 2781479 --to-bp 155701383 --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex  0.25 0.75 --out {sex_tmp2}"

        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd)

        listOfFiles = [f'{sex_tmp1}.log', f'{sex_tmp2}.log']
        concat_logs(step, out_path, listOfFiles)

        if os.path.isfile(f'{sex_tmp1}.sexcheck') and os.path.isfile(f'{sex_tmp2}.sexcheck'):

            # grab fails from .sexcheck files
            sex1 = pd.read_csv(f'{sex_tmp1}.sexcheck', sep='\s+')
            sex_fail1 = sex1[sex1.STATUS=='PROBLEM']

            sex2 = pd.read_csv(f'{sex_tmp2}.sexcheck', sep='\s+')
            sex_fail2 = sex2[sex2.STATUS=='PROBLEM']

            # combine and output
            sex_fail_df = pd.concat([sex_fail1, sex_fail2], ignore_index=True)
            sex_fail_ids = sex_fail_df.loc[:,['FID','IID']].drop_duplicates(subset=['FID','IID'])
            sex_fail_count = sex_fail_ids.shape[0]
            sex_fail_ids.to_csv(sex_fails, sep='\t', header=True, index=False)

            # remove sex fail samples from geno
            plink_cmd3 = f"{plink2_exec} --bfile {geno_path} --remove {sex_fails} --make-bed --out {out_path}"
            
            shell_do(plink_cmd3)

            listOfFiles = [f'{out_path}.log']
            concat_logs(step, out_path, listOfFiles)
            
            process_complete = True

            # log outputs
            outfiles_dict = {
                'pruned_samples': sex_fails,
                'plink_out': out_path
            }

            metrics_dict = {
                'outlier_count': sex_fail_count
            }
        
        else:
            print('Sex Prune Failed!')
            print(f'Check {sex_tmp1}.log and {sex_tmp2}.log for more information')

            process_complete = False

            outfiles_dict = {
                'pruned_samples': 'Sex Prune Failed!',
                'plink_out': out_path
            }

            metrics_dict = {
                'outlier_count': 0
            }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict
    

    def run_het_prune(self, het_filter=[-0.25,0.25]):

        """
        Execute heterozygosity-based pruning on genotype data using PLINK.
        
        Parameters:
        - het_filter (list of float, optional): Two values indicating the bounds for heterozygosity outlier checks. Default values are [-0.25, 0.25].
        
        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('het_prune').
            * 'metrics': Metrics associated with the pruning, such as 'outlier_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """

        geno_path = self.geno_path
        out_path = self.out_path

        step = "het_prune"

        het_tmp = f"{out_path}_tmp"
        het_tmp2 = f"{out_path}_tmp2"
        het_tmp3 = f"{out_path}_tmp3"
        outliers_out = f"{out_path}.outliers"
        
        # variant(maf=0.05, geno=0.01, indep_pairwise=[50,5,0.5])

        plink_cmd1 = f"{plink2_exec} --bfile {geno_path} --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out {het_tmp}"
        plink_cmd2 = f"{plink2_exec} --bfile {geno_path} --extract {het_tmp}.prune.in --make-bed --out {het_tmp2}"
        plink_cmd3 = f"{plink2_exec} --bfile {het_tmp2} --het --out {het_tmp3}"

        cmds1 = [plink_cmd1, plink_cmd2, plink_cmd3]

        for cmd in cmds1:
            shell_do(cmd)

        listOfFiles = [f'{het_tmp}.log', f'{het_tmp2}.log', f'{het_tmp3}.log']
        concat_logs(step, out_path, listOfFiles)
        
        # check if het_tmp3 is created. if not, skip this.
        # NOTE: there may be legitimate reasons for this, for example, if only one sample in genotype file (sometimes happens in ancestry split method)
        hetpath = f'{het_tmp3}.het'
        if os.path.isfile(hetpath):
            het = pd.read_csv(hetpath, sep='\s+')
            het_outliers = het[((het.F <= het_filter[0]) | (het.F >= het_filter[1]))]
            outlier_count = het_outliers.shape[0]
            het_outliers.to_csv(f'{outliers_out}', sep='\t', header=True, index=False)
        
            plink_cmd4 = f"{plink2_exec} --bfile {geno_path} --remove {outliers_out} --make-bed --out {out_path}"

            shell_do(plink_cmd4)

            listOfFiles = [f'{out_path}.log']
            concat_logs(step, out_path, listOfFiles)

            if os.path.isfile(f'{out_path}.bed'):
                outfiles_dict = {
                    'pruned_samples': outliers_out,
                    'plink_out': out_path
                }

                metrics_dict = {
                    'outlier_count': outlier_count
                }

                process_complete = True
            
            else:
                print(f'Heterozygosity pruning failed!')
                print(f'Check {out_path}.log for more information')

                outfiles_dict = {
                    'pruned_samples': 'Heterozygosity Pruning Failed!',
                    'plink_out': [het_tmp, het_tmp2, het_tmp3, out_path]
                }

                metrics_dict = {
                    'outlier_count': 0
                }

                process_complete = False

            out_dict = {
                'pass': process_complete,
                'step': step,
                'metrics': metrics_dict,
                'output': outfiles_dict
            }          
        
        else:
            print(f'Heterozygosity pruning failed!')
            print(f'Check {het_tmp}.log, {het_tmp2}.log, or {het_tmp3}.log for more information')
            
            outfiles_dict = {
                'pruned_samples': 'Heterozygosity Pruning Failed!',
                'plink_out': [het_tmp, het_tmp2, het_tmp3]
            }

            metrics_dict = {
                'outlier_count': 0
            }
        
            process_complete = False
        
        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict


    def run_related_prune(self,related_cutoff=0.0884, duplicated_cutoff=0.354, prune_related=True, prune_duplicated=True):

        """
        Execute pruning based on relatedness and duplication checks on genotype data using PLINK and KING.
        
        Parameters:
        - related_cutoff (float, optional): Threshold for relatedness check. Default value is 0.0884.
        - duplicated_cutoff (float, optional): Threshold for duplication check. Default value is 0.354.
        - prune_related (bool, optional): Whether to prune related samples. Default is True.
        - prune_duplicated (bool, optional): Whether to prune duplicated samples. Default is True.
        
        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('related_prune').
            * 'metrics': Metrics associated with the pruning.
            * 'output': Dictionary containing paths to the generated output files.
        """

        geno_path = self.geno_path
        out_path = self.out_path

        step = "related_prune"

        grm1 = f"{out_path}_total_grm"
        grm2 = f"{out_path}_related_grm"
        grm3 = f"{out_path}_duplicated_grm"

        related_out = f"{out_path}_pairs.related"
        related_pruned_out = f"{out_path}.pruned"

        # create pfiles
        king_cmd1 = f'{plink2_exec} --bfile {geno_path} --hwe 0.0001 --mac 2 --make-pgen --out {grm1}'
        # create table of related pairs
        king_cmd2 = f'{plink2_exec} --pfile {grm1} --make-king-table --make-king triangle bin --king-table-filter {related_cutoff} --out {out_path}_pairs'
        # see if any samples are related (includes duplicates)
        king_cmd3 = f'{plink2_exec} --pfile {grm1} --king-cutoff {out_path}_pairs {related_cutoff} --out {grm2}' 
        # see if any samples are duplicated (grm cutoff >= 0.354)
        king_cmd4 = f'{plink2_exec} --pfile {grm1} --king-cutoff {out_path}_pairs {duplicated_cutoff} --out {grm3}' 

        cmds = [king_cmd1, king_cmd2, king_cmd3, king_cmd4]
        for cmd in cmds:
            shell_do(cmd)

        listOfFiles = [f'{grm1}.log', f'{out_path}_pairs.log', f'{grm2}.log', f'{grm3}.log']
        concat_logs(step, out_path, listOfFiles)

        if os.path.isfile(f'{out_path}_pairs.kin0') and os.path.isfile(f'{grm2}.king.cutoff.out.id') and os.path.isfile(f'{grm3}.king.cutoff.out.id'):

            # create .related related pair sample files
            shutil.copy(f'{out_path}_pairs.kin0',f'{out_path}_pairs.related')

            # create .related and .duplicated single sample files
            shutil.copy(f'{grm2}.king.cutoff.out.id',f'{grm2}.related')
            related_count = sum(1 for line in open(f'{grm2}.related'))

            shutil.copy(f'{grm3}.king.cutoff.out.id',f'{grm3}.duplicated')
            duplicated_count = sum(1 for line in open(f'{grm3}.duplicated'))

            related_count = related_count - duplicated_count
            duplicated = pd.read_csv(f'{grm3}.duplicated', sep = '\s+')

            # concat duplicated sample ids to related sample ids, drop_duplicates(keep='last) because all duplicated would also be considered related
            if prune_related and prune_duplicated:
                plink_cmd1 = f'{plink2_exec} --pfile {grm1} --remove {grm2}.king.cutoff.out.id --make-bed --out {out_path}'
                shell_do(plink_cmd1) 

                related = pd.read_csv(f'{grm2}.related', sep = '\s+')
                grm_pruned = pd.concat([related, duplicated], ignore_index=True)

                if '#FID' in grm_pruned:
                    grm_pruned.rename(columns={"#FID": "FID", "#IID": "IID"}, inplace = True)
                else:
                    grm_pruned['FID'] = 0
                    grm_pruned.rename(columns={"#IID": "IID"}, inplace = True)
                grm_pruned.drop_duplicates(subset=['FID','IID'], keep='last', inplace=True)
                grm_pruned.to_csv(related_pruned_out, sep='\t', header=True, index=False)
                process_complete = True
            
            if prune_duplicated and not prune_related:
                plink_cmd1 = f'{plink2_exec} --pfile {grm1} --remove {grm3}.king.cutoff.out.id --make-bed --out {out_path}'
                shell_do(plink_cmd1) 

                grm_pruned = duplicated
                
                if '#FID' in grm_pruned:
                    grm_pruned.rename(columns={"#FID": "FID", "#IID": "IID"}, inplace = True)
                else:
                    grm_pruned['FID'] = 0
                    grm_pruned.rename(columns={"#IID": "IID"}, inplace = True)
                grm_pruned.drop_duplicates(subset=['FID','IID'], keep='last', inplace=True)
                grm_pruned.to_csv(related_pruned_out, sep='\t', header=True, index=False)
                process_complete = True
                
            if not prune_related and not prune_duplicated:
                plink_cmd1 = f'echo prune_related and prune_duplicated set to False. Pruning passed'
                shell_do(plink_cmd1)

                process_complete = True
                related_pruned_out = None
                process_complete = True
                
            if not prune_duplicated and prune_related:
                print('This option is invalid. Cannot prune related without also pruning duplicated')
                process_complete = False


            if os.path.isfile('{out_path}.log'):
                listOfFiles = [f'{out_path}.log']
                concat_logs(step, out_path, listOfFiles)

            # remove intermediate files
            os.remove(f'{grm1}.pgen')
            os.remove(f'{grm1}.pvar')
            os.remove(f'{grm1}.psam')
            os.remove(f'{out_path}_pairs.kin0')

            outfiles_dict = {
                'pruned_samples': related_pruned_out,
                'related_samples': related_out,
                'plink_out': out_path
            }

            metrics_dict = {
                'related_count': related_count,
                'duplicated_count': duplicated_count
            }

        else:
            print(f'Relatedness pruning failed!')
            print(f'Check all_plink_logs.gtlog for more information')
            process_complete = False

            outfiles_dict = {
                'pruned_samples': 'Related Pruning Failed',
                'related_samples': 0, 
                'plink_out': [grm1, grm2, grm3]
            }

            metrics_dict = {
                'related_count': 0,
                'duplicated_count': 0
            }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict
    

class VariantQC:

    def __init__(self, geno_path, out_path):
        self.geno_path = geno_path
        self.out_path = out_path

    def run_geno_prune(self, geno_threshold=0.05):

        """
        Prunes SNPs based on missing genotype data.

        Parameters:
        - geno_threshold (float, optional): Maximum allowable proportion of missing genotyping for a SNP. Default is 0.05 (5%).

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('geno_prune').
            * 'metrics': Metrics associated with the pruning, such as 'geno_removed_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """

        geno_path = self.geno_path
        out_path = self.out_path

        step = "geno_prune"

        # get initial snp count
        initial_snp_count = count_file_lines(f'{geno_path}.bim')

        # variant missingness
        plink_cmd1 = f"{plink2_exec} --bfile {geno_path} --geno {geno_threshold} --make-bed --out {out_path}"
        shell_do(plink_cmd1)

        listOfFiles = [f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)

        # geno pruned count
        geno_snp_count = count_file_lines(f'{out_path}.bim')
        geno_rm_count = initial_snp_count - geno_snp_count

        process_complete = True

        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'geno_removed_count': geno_rm_count
        }
        
        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict
        

    def run_case_control_prune(self, p_threshold=1e-4):

        """
        Prunes SNPs based on missing genotype data, stratified by case-control status.

        Parameters:
        - p_threshold (float, optional): P-value threshold for testing the difference in missing rates between cases and controls. Default is 1e-4.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('case_control_missingness_prune').
            * 'metrics': Metrics associated with the pruning, such as 'mis_removed_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """
        
        geno_path = self.geno_path
        out_path = self.out_path

        step = "case_control_missingness_prune"
        mis_tmp = f'{out_path}_mis_tmp'
        
        # get initial snp count
        initial_snp_count = count_file_lines(f'{geno_path}.bim')

        fam = pd.read_csv(f'{geno_path}.fam', sep='\s+', header=None, usecols=[5], names=['case'])
        # check if this contains both cases and controls
        if  all(x in fam['case'].unique() for x in [1, 2]):
            # missingness by case control (--test-missing), using P > 1E-4
            plink_cmd1 = f"{plink_exec} --bfile {geno_path} --test-missing --out {mis_tmp}"
            shell_do(plink_cmd1)

            listOfFiles = [f'{mis_tmp}.log']
            concat_logs(step, out_path, listOfFiles)

            mis1 = f'{mis_tmp}.missing'
            if os.path.isfile(mis1):

                mis = pd.read_csv(mis1, sep='\s+')
                exclude = mis[mis.P <= p_threshold].loc[:,'SNP']
                exclude.to_csv(f'{mis_tmp}.exclude', sep='\t', header=False, index=False)

                plink_cmd2 = f"{plink2_exec} --bfile {geno_path} --exclude {mis_tmp}.exclude --make-bed --out {out_path}"
                shell_do(plink_cmd2)

                listOfFiles = [f'{out_path}.log']
                concat_logs(step, out_path, listOfFiles)

                # mis pruned count
                mis_snp_count = count_file_lines(f'{out_path}.bim')
                mis_rm_count = initial_snp_count - mis_snp_count

                process_complete = True

            else:
                print(f'Case/Control Missingness pruning failed!')
                print(f'Check {mis_tmp}.log for more information')

                process_complete = False
                mis_snp_count = count_file_lines(f'{geno_path}.bim')
                mis_rm_count = 0
        else:
            print(f'Case/control pruning failed! May be missing Controls!')
            print(f'Check number of Cases/Controls')
            
            process_complete = False

        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'mis_removed_count': mis_rm_count
        }

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict


    def run_haplotype_prune(self, p_threshold=1e-4):
        
        """
        Prunes SNPs based on missing haplotype data.

        Parameters:
        - p_threshold (float, optional): P-value threshold for missingness by haplotype. Default is 1e-4.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('haplotype_prune').
            * 'metrics': Metrics associated with the pruning, such as 'haplotype_removed_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """
            
        geno_path = self.geno_path
        out_path = self.out_path

        step = "haplotype_prune"

        # make tmp names
        # missingness by haplotype
        hap_tmp = f'{out_path}_hap_tmp'

        # get initial snp count
        initial_snp_count = count_file_lines(f'{geno_path}.bim')

        # missingness by haplotype (--test-mishap), using P > 1E-4
        plink_cmd1 = f"{plink_exec} --bfile {geno_path} --test-mishap --out {hap_tmp}"
        shell_do(plink_cmd1)

        listOfFiles = [f'{hap_tmp}.log']
        concat_logs(step, out_path, listOfFiles)

        # read .missing.hap file and grab flanking snps for P <= 0.0001. write flanking snps to file to exclude w bash
        mis_hap = pd.read_csv(f'{hap_tmp}.missing.hap', sep='\s+')
        mis_hap_snps = list(mis_hap[mis_hap.P <= p_threshold].loc[:,'FLANKING'].str.split('|'))
        snp_ls_df = pd.DataFrame({'snp':[rsid for ls in mis_hap_snps for rsid in ls]})
        snp_ls_df['snp'].to_csv(f'{hap_tmp}.exclude',sep='\t', header=False, index=False)

        plink_cmd2 = f"{plink2_exec} --bfile {geno_path} --exclude {hap_tmp}.exclude --make-bed --out {out_path}"
        shell_do(plink_cmd2)

        listOfFiles = [f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)

        # hap pruned count
        hap_snp_count = count_file_lines(f'{out_path}.bim')
        hap_rm_count = initial_snp_count - hap_snp_count

        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'haplotype_removed_count': hap_rm_count
        }

        process_complete = True

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict


    def run_hwe_prune(self, hwe_threshold=1e-4, filter_controls=False):
        
        """
        Prunes SNPs based on Hardy-Weinberg equilibrium (HWE) p-values.

        Parameters:
        - hwe_threshold (float, optional): P-value threshold for SNP filtering based on HWE. Default is 1e-4.
        - filter_controls (bool, optional): If set to True, perform HWE test only on control samples. Default is False.

        Returns:
        - dict: A structured dictionary containing:
            * 'pass': Boolean indicating the successful completion of the process.
            * 'step': The label for this procedure ('hwe_prune').
            * 'metrics': Metrics associated with the pruning, such as 'hwe_removed_count'.
            * 'output': Dictionary containing paths to the generated output files.
        """
        
        geno_path = self.geno_path
        out_path = self.out_path

        step = "hwe_prune"

        hwe_tmp = f'{out_path}_hwe_tmp'

        # get initial snp count
        initial_snp_count = count_file_lines(f'{geno_path}.bim')
        
        if filter_controls:
            # HWE from controls only using P > 1E-4
            plink_cmd1 = f"{plink_exec} --bfile {geno_path} --filter-controls --hwe {hwe_threshold} --write-snplist --out {hwe_tmp}"
        else:
            # HWE using P > 1E-4
            plink_cmd1 = f"{plink_exec} --bfile {geno_path} --hwe {hwe_threshold} --write-snplist --out {hwe_tmp}"

        plink_cmd2 = f"{plink2_exec} --bfile {geno_path} --extract {hwe_tmp}.snplist --make-bed --out {out_path}"

        cmds = [plink_cmd1, plink_cmd2]
        for cmd in cmds:
            shell_do(cmd) 

        listOfFiles = [f'{hwe_tmp}.log', f'{out_path}.log']
        concat_logs(step, out_path, listOfFiles)

        # hwe pruned count
        final_snp_count = count_file_lines(f'{out_path}.bim')
        hwe_rm_count = initial_snp_count - final_snp_count

        outfiles_dict = {
            'plink_out': out_path
        }

        metrics_dict = {
            'hwe_removed_count': hwe_rm_count
        }

        process_complete = True

        out_dict = {
            'pass': process_complete,
            'step': step,
            'metrics': metrics_dict,
            'output': outfiles_dict
        }

        return out_dict

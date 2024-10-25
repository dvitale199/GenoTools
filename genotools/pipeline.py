# Copyright 2023 The GenoTools Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================


import os
import shutil
import argparse
import warnings
import pathlib
import platform
import pandas as pd


def gt_argparse():
    # definte arg parse
    parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')

    # file i/o arguments
    parser.add_argument('--bfile', type=str, nargs='?', default=None, const=None, help='Genotype: String file path to PLINK 1.9 format genotype file (everything before the *.bed/bim/fam)')
    parser.add_argument('--pfile', type=str, nargs='?', default=None, const=None, help='Genotype: String file path to PLINK 2 format genotype file (everything before the *.pgen/pvar/psam)')
    parser.add_argument('--vcf', type=str, nargs='?', default=None, const=None, help='Genotype: String file path to VCF format genotype file')
    parser.add_argument('--out', type=str, nargs='?', default=None, const=None, help='Prefix for output (including path)', required=True)
    parser.add_argument('--full_output', type=str, nargs='?', default='False', const='True', help='Output everything')
    parser.add_argument('--skip_fails', type=str, nargs='?', default='False', const='True', help='Skip up front check for fails')
    parser.add_argument('--warn', type=str, nargs='?', default='True', const='True', help='Warn of error and continue running pipeline')

    # ancerstry arguments
    parser.add_argument('--ancestry', type=str, nargs='?', default='False', const='True', help='Split by ancestry')
    parser.add_argument('--ref_panel', type=str, nargs='?', default=None, const=None, help='Genotype: (string file path). Path to PLINK format reference genotype file, everything before the *.bed/bim/fam.')
    parser.add_argument('--ref_labels', type=str, nargs='?', default=None, const=None, help='tab-separated plink-style IDs with ancestry label (FID  IID label) with no header')
    parser.add_argument('--model', type=str, nargs='?', default=None, const='path', help='Path to pickle file with trained ancestry model for passed reference panel')
    parser.add_argument('--container', type=str, nargs='?', default='False', const='True', help='Run predictions in container')
    parser.add_argument('--singularity', type=str, nargs='?', default='False', const='True', help='Run containerized precitions via singularity')
    parser.add_argument('--cloud', type=str, nargs='?', default='False', const='True', help='Run predictions in GCloud')
    parser.add_argument('--cloud_model', type=str, nargs='?', default='NeuroBooster', help='Model for GCloud predictions')
    parser.add_argument('--subset_ancestry', nargs='*', help='Subset to continue analysis for')
    parser.add_argument('--min_samples', type=int, nargs='?', default=0, const=50, help='Minimum number of samples in an ancestry group required for subsequent analyses to be performed')

    # sample-level qc arguments
    parser.add_argument('--callrate', type=float, nargs='?', default=None, const=0.02, help='Minimum Callrate threshold for QC')
    parser.add_argument('--sex', nargs='*', help='Sex prune with cutoffs')
    parser.add_argument('--related', type=str, nargs='?', default='False', const='True', help='Relatedness prune')
    parser.add_argument('--related_cutoff', type=float, nargs='?', default=0.0884, const=0.0884, help='Relatedness cutoff')
    parser.add_argument('--duplicated_cutoff', type=float, nargs='?', default=0.354, const=0.354, help='Relatedness cutoff')
    parser.add_argument('--prune_related', type=str, nargs='?', default='False', const='True', help='Relatedness prune')
    parser.add_argument('--prune_duplicated', type=str, nargs='?', default='True', const='True', help='Relatedness prune')
    parser.add_argument('--het', nargs='*', help='Het prune with cutoffs')
    parser.add_argument('--amr_het', type=str, nargs='?', default='False', const='True', help='Custom het prune for GP2 AMR samples')
    parser.add_argument('--kinship_check', type=str, nargs='?', default='False', const='True', help='Confirming familial labels')
    parser.add_argument('--all_sample', type=str, nargs='?', default='False', const='True', help='Run all sample-level QC')

    # variant-level qc arguments
    parser.add_argument('--geno', type=float, nargs='?', default=None, const=0.05, help='Minimum Missingness threshold for QC')
    parser.add_argument('--case_control', type=float, nargs='?', default=None, const=1e-4, help='Case control prune')
    parser.add_argument('--haplotype', type=float, nargs='?', default=None, const=1e-4, help='Haplotype prune')
    parser.add_argument('--hwe', type=float, nargs='?', default=None, const=1e-4, help='HWE pruning')
    parser.add_argument('--filter_controls', type=str, nargs='?', default='False', const='True', help='Control filter for HWE prune')
    parser.add_argument('--ld', nargs='*', help='LD prune with window size, step size, r2 threshold')
    parser.add_argument('--all_variant', type=str, nargs='?', default='False', const='True', help='Run all variant-level QC')

    # GWAS and PCA argument
    parser.add_argument('--pca', type=int, nargs='?', default=None, const=10, help='PCA and number of PCs')
    parser.add_argument('--build', type=str, nargs='?', default='hg38', const='hg38', help='Build for PCA')
    parser.add_argument('--gwas', type=str, nargs='?', default='False', const='True', help='Run GWAS')
    parser.add_argument('--covars', type=str, nargs='?', default=None, const=None, help='Path to external covars')
    parser.add_argument('--covar_names', type=str, nargs='?', default=None, const=None, help='Covar names to use from external file')
    parser.add_argument('--maf_lambdas', type=str, nargs='?', default='False', const='True', help='MAF prune before lambda calculations')


    # parse args and turn into dict
    args = parser.parse_args()

    return args


def execute_ancestry_predictions(geno_path, out_path, args, ancestry, tmp_dir):
    if not args['full_output']:
        out_path_pathlib = pathlib.PurePath(out_path)
        out_path_name = out_path_pathlib.name
        out_path = f'{tmp_dir.name}/{out_path_name}'

    ancestry.geno_path = geno_path
    ancestry.out_path = out_path
    ancestry.ref_panel = args['ref_panel']
    ancestry.ref_labels = args['ref_labels']
    ancestry.model_path = args['model']
    ancestry.containerized = args['container']
    ancestry.singularity = args['singularity']
    ancestry.cloud = args['cloud']
    ancestry.cloud_model = args['cloud_model']
    ancestry.subset = args['subset_ancestry']
    ancestry.min_samples = args['min_samples']
    
    ancestry_dict = ancestry.run_ancestry()

    return ancestry_dict


def execute_pipeline(steps, steps_dict, geno_path, out_path, samp_qc, var_qc, assoc, args, tmp_dir):
    # to know which class to call
    samp_steps = ['callrate','sex','related','het','kinship_check']
    var_steps = ['case_control','haplotype','hwe','geno','ld']

    # if full output requested, go to out path
    if args['full_output']:
        step_paths = [out_path]

    # otherwise tmpdir
    else:
        geno_path_pathlib = pathlib.PurePath(geno_path)
        geno_path_name = geno_path_pathlib.name
        geno = f'{tmp_dir.name}/{geno_path_name}'
        out_path_pathlib = pathlib.PurePath(out_path)
        out_path_name = out_path_pathlib.name
        out = f'{tmp_dir.name}/{out_path_name}'
        step_paths = [out]

    # initialize pass/fail and out dicts
    pass_fail = dict()
    out_dict = dict()

    # loop through steps
    for step in steps:
        # if warn is true, find last passed step
        if args['warn']:
            # find last passed step
            last_passed = None
            for completed_step in pass_fail:
                if pass_fail[completed_step]['status']:
                    last_passed = completed_step
            # if the last passed step exists, point to its output
            if last_passed:
                step_input = pass_fail[last_passed]['output']
                step_output = f'{step_input}_{step}'
            # otherwise no steps have passed so go back to geno path
            else:
                step_input = geno_path if ((args['full_output']) or (not args['ancestry'])) else geno
                step_output = f'{out_path}_{step}' if args['full_output'] else f'{out}_{step}'
            
            # last step case
            if step == steps[-1]:
                step_output = f'{out_path}'
        
        # otherwise just go in order
        else:
            if step != steps[0]:
                step_input = f'{step_paths[-1]}'
            elif args['full_output'] or (not args['ancestry']):
                step_input = geno_path
            else:
                step_input = geno

            # step_input = f'{step_paths[-1]}' if step != steps[0] else geno_path
            step_output = f'{step_paths[-1]}_{step}' if step != steps[-1] else out_path
        
        print(f'Running: {step} with input {step_input} and output: {step_output}')
        step_paths.append(step_output)

        # if warn is on and input doesn't exist, all samples or variants were pruned in a previous step
        if args['warn'] and (not os.path.isfile(f'{step_input}.pgen')):
            print(f'Step {step} cannot be run! All samples or variants were pruned in a previous step!')
            pass_fail[step] = {'status':False, 'input':step_input, 'output':step_output}
        
        # otherwise run the qc step
        else:
            # samp qc setup and call
            if step in samp_steps:
                samp_qc.geno_path = step_input
                samp_qc.out_path = step_output

                # related has more than one parameter
                if step == 'related':
                    out_dict[step] = steps_dict[step](related_cutoff=args['related_cutoff'], duplicated_cutoff=args['duplicated_cutoff'],
                                    prune_related=args['prune_related'], prune_duplicated=args['prune_duplicated'])

                elif step == 'kinship_check':
                    # check that OS is not macOS
                    if platform.system() != 'Linux':
                        print('Relatedness Assessment can only run on a linux or windows OS!')
                        out_dict.pop(step, None)
                    elif platform.system() == 'Linux':
                        out_dict[step] = steps_dict[step]()

                else:
                    out_dict[step] = steps_dict[step](args[step])

            # var qc setup and call
            if step in var_steps:
                var_qc.geno_path = step_input
                var_qc.out_path = step_output

                # hwe and ld have extra parameters
                if step == 'hwe':
                    out_dict[step] = steps_dict[step](hwe_threshold=args['hwe'], filter_controls=args['filter_controls'])

                elif step == 'ld':
                    out_dict[step] = steps_dict[step](window_size=args['ld'][0], step_size=args['ld'][1], r2_threshold=args['ld'][2])

                else:
                    out_dict[step] = steps_dict[step](args[step])

            # assoc setup and call
            if step == 'assoc':
                assoc.geno_path = step_input
                assoc.out_path = step_output
                assoc.pca = args['pca']
                assoc.build = args['build']
                assoc.gwas = args['gwas']
                assoc.covar_path = args['covars']
                assoc.covar_names = args['covar_names']
                assoc.maf_lambdas = args['maf_lambdas']
                out_dict[step] = steps_dict[step]()
            
            pass_fail[step] = {'status':out_dict[step]['pass'], 'input':step_input, 'output':step_output}

            # remove old files when appropriate
            if (not args['full_output']) and (not args['warn']) and (step != 'assoc') and (step != 'ancestry') and (step != 'kinship_check'):
                # when warn is True and step fails, don't remove old file
                if args['warn'] and ('pass' in out_dict[step].keys()) and (not out_dict[step]['pass']):
                    remove = False
                else:
                    remove = True
                    # remove_step_index = step_paths.index(step_output) - 1
                    remove_path = pass_fail[step]['input']

                if remove:
                    # remove_path = step_paths[remove_step_index]
                    # make sure we're not removing the output
                    if os.path.isfile(f'{remove_path}.pgen') and (remove_path != out_path) and (remove_path != geno_path):
                        os.remove(f'{remove_path}.pgen')
                        os.remove(f'{remove_path}.psam')
                        os.remove(f'{remove_path}.pvar')
    
    # if the pipeline has more than one step, warn is on, and the last step of the pipeline fails, move the last passed files to output
    if args['warn'] and (len(steps) > 1) and (not pass_fail[steps[-1]]['status']):
        if last_passed:
            if os.path.isfile(f"{pass_fail[last_passed]['output']}.pgen"):
                os.rename(f"{pass_fail[last_passed]['output']}.pgen", f"{out_path}.pgen")
                os.rename(f"{pass_fail[last_passed]['output']}.psam", f"{out_path}.psam")
                os.rename(f"{pass_fail[last_passed]['output']}.pvar", f"{out_path}.pvar")
            # cases when step is passed but all samples or vars get pruned
            else:
                os.rename(f"{pass_fail[last_passed]['input']}.pgen", f"{out_path}.pgen")
                os.rename(f"{pass_fail[last_passed]['input']}.psam", f"{out_path}.psam")
                os.rename(f"{pass_fail[last_passed]['input']}.pvar", f"{out_path}.pvar")
        # cases when no steps are passed
        else:
            move_path = geno_path if ((args['full_output']) or (not args['ancestry'])) else geno
            os.rename(f"{move_path}.pgen", f"{out_path}.pgen")
            os.rename(f"{move_path}.psam", f"{out_path}.psam")
            os.rename(f"{move_path}.pvar", f"{out_path}.pvar")

    out_dict['paths'] = step_paths
    out_dict['pass_fail'] = pass_fail

    return out_dict


def build_metrics_pruned_df(metrics_df, pruned_df, gwas_df, related_df, dictionary, out, ancestry='all'):
    for step in ['callrate', 'sex', 'related', 'het', 'case_control', 'haplotype', 'hwe', 'geno','ld']:
        if step in dictionary.keys():
            qc_step = dictionary[step]['step']
            pf = dictionary[step]['pass']

            if step in ['callrate', 'sex', 'related', 'het']:
                level = 'sample'
                samplefile = dictionary[step]['output']['pruned_samples']
                if (samplefile is not None) and os.path.isfile(samplefile):
                    pruned = pd.read_csv(samplefile, sep='\t')
                    if pruned.shape[0] > 0:
                        pruned.loc[:,'step'] = step
                        pruned_df = pd.concat([pruned_df, pruned[['#FID','IID','step']]], ignore_index=True)

                if step == 'related':
                    relatedfile = dictionary[step]['output']['related_samples']
                    if (relatedfile is not None) and os.path.isfile(relatedfile):
                        related = pd.read_csv(relatedfile, sep=',')
                        if related.shape[0] > 0:
                            if ancestry == 'all':
                                related_out_path = f'{out}.related'
                            else:
                                related_out_path = f'{out}_{ancestry}.related'
                            related.to_csv(related_out_path, index=False)

                            related['ancestry'] = ancestry
                            related_df = pd.concat([related_df, related], ignore_index=True)

            else:
                level = 'variant'

            for metric, value in dictionary[step]['metrics'].items():
                tmp_metrics_df = pd.DataFrame({'step':[qc_step], 'pruned_count':[value], 'metric':[metric], 'ancestry':[ancestry], 'level':[level], 'pass': [pf]})
                metrics_df = pd.concat([metrics_df, tmp_metrics_df], ignore_index=True)

    if ('assoc' in dictionary.keys()) and ('gwas' in dictionary['assoc'].keys()):
        for metric, value in dictionary['assoc']['gwas']['metrics'].items():
            tmp_gwas_df = pd.DataFrame({'value':[value], 'metric':[metric], 'ancestry':[ancestry]})
            gwas_df = pd.concat([gwas_df, tmp_gwas_df], ignore_index=True)

    return metrics_df, pruned_df, gwas_df, related_df
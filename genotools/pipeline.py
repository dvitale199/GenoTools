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


    # parse args and turn into dict
    args = parser.parse_args()
    
    return args


def execute_pipeline(steps, steps_dict, geno_path, out_path, samp_qc, var_qc, ancestry, assoc, args, tmp_dir):
    # to know which class to call
    samp_steps = ['callrate','sex','related','het']
    var_steps = ['case_control','haplotype','hwe','geno','ld']

    # if full output requested, go to out path
    if args['full_output']:
        step_paths = [out_path]

    # otherwise tmpdir
    else:
        out_path_pathlib = pathlib.PurePath(out_path)
        out_path_name = out_path_pathlib.name
        step_paths = [f'{tmp_dir.name}/{out_path_name}']

    # if first step is ancestry, make new steps list to call within-ancestry
    if steps[0] == 'ancestry':
        steps_ancestry = steps[1:]
        steps = [steps[0]]

    out_dict = dict()

    # loop through steps
    for step in steps:
        # use geno_path for first step, out_path for last step
        step_input = f'{step_paths[-1]}' if step != steps[0] else geno_path
        step_output = f'{step_paths[-1]}_{step}' if step != steps[-1] else out_path
        print(f'Running: {step} with input {step_input} and output: {step_output}')

        # ancestry setup and call
        if step == 'ancestry':
            # output goes to out_path if no steps requested after ancestry, otherwise put in at out_path_ancestry (for tmps)
            step_output = f'{step_paths[-1]}_{step}' if len(steps_ancestry) > 0 else out_path
            step_paths.append(step_output)

            ancestry.geno_path = step_input
            ancestry.out_path = step_output
            ancestry.ref_panel = args['ref_panel']
            ancestry.ref_labels = args['ref_labels']
            ancestry.model_path = args['model']
            ancestry.containerized = args['container']
            ancestry.singularity = args['singularity']
            ancestry.cloud = args['cloud']
            ancestry.cloud_model = args['cloud_model']
            ancestry.subset = args['subset_ancestry']
            ancestry.min_samples = args['min_samples']
            out_dict[step] = steps_dict[step]()

            # call ancestry specific steps within each group
            if len(steps_ancestry) > 0:
                for geno, label in zip(out_dict[step]['output']['split_paths'], out_dict[step]['data']['labels_list']):
                    out_dict[label] = execute_pipeline(steps_ancestry, steps_dict, geno, f'{out_path}_{label}', samp_qc, var_qc, ancestry, assoc, args, tmp_dir)

        else:
            # if warn is True and step input doesn't exist print error and reset step input
            if args['warn'] and (not os.path.isfile(f'{step_input}.pgen')) and (len(step_paths) > 1):
                warnings.warn(f'{step_input}.pgen was not created. Continuing to next step...', stacklevel=2)
                step_input = f'{step_paths[-2]}' if step != steps[1] else geno_path
                
                # very rare edge case when multiple steps fail
                if not os.path.isfile(f'{step_input}.pgen'):
                    step_input = f'{step_paths[-3]}'

                step_output = f'{step_paths[-2]}_{step}' if step != steps[-1] else out_path
                step_paths.append(step_output)
            
            # otherwise keep track of paths
            else:
                step_paths.append(step_output)

        # samp qc setup and call
        if step in samp_steps:
            samp_qc.geno_path = step_input
            samp_qc.out_path = step_output

            # related has more than one parameter
            if step == 'related':
                out_dict[step] = steps_dict[step](related_cutoff=args['related_cutoff'], duplicated_cutoff=args['duplicated_cutoff'],
                                 prune_related=args['prune_related'], prune_duplicated=args['prune_duplicated'])
            
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
            out_dict[step] = steps_dict[step]()
        
        # remove old files when appropriate 
        if (not args['full_output']) and (step != 'assoc') and (step != 'ancestry'):
            # when warn is True and step fails, don't remove old file
            if args['warn'] and ('pass' in out_dict[step].keys()) and (not out_dict[step]['pass']):
                remove = False
            else:
                remove = True
                remove_step_index = step_paths.index(step_output) - 1

            if remove:
                remove_path = step_paths[remove_step_index]
                # make sure we're not removing the output
                if os.path.isfile(f'{remove_path}.pgen') and (remove_path != out_path):
                    os.remove(f'{remove_path}.pgen')
                    os.remove(f'{remove_path}.psam')
                    os.remove(f'{remove_path}.pvar')

    out_dict['paths'] = step_paths

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
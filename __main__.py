import subprocess
import sys
import shutil
import argparse
import os
import tempfile
import pandas as pd
import json

from genotools.utils import shell_do, upfront_check
from genotools.qc import SampleQC, VariantQC
from genotools.ancestry import Ancestry
from genotools.gwas import Assoc

def execute_pipeline(steps, steps_dict, geno_path, out_path, samp_qc, var_qc, ancestry, assoc, args, tmp_dir):
    # to know which class to call
    samp_steps = ['callrate','sex','related','het']
    var_steps = ['case_control','haplotype','hwe','geno']

    # if full output requested, go to out path
    if args['full_output']:
        step_paths = [out_path]

    # otherwise tmpdir
    else:
        step_paths = [f'{tmp_dir.name}/out']

    # if first step is ancestry, make new steps list to call within-ancestry
    if steps[0] == 'ancestry':
        steps_ancestry = steps[1:]
        steps = [steps[0]]

        # move common snps file to tmp dir if full output is not requested
        if not args['full_output']:
            common_snps_file = f'{os.path.dirname(out_path)}/ref_common_snps.common_snps'

            if os.path.isfile(common_snps_file):
                shutil.copy(common_snps_file, f'{tmp_dir.name}/ref_common_snps.common_snps')
            else:
                raise FileNotFoundError(f'{common_snps_file} does not exist.')

    out_dict = dict()

    # loop through steps
    for step in steps:
        # use geno_path for first step, out_path for last step
        step_input = f'{step_paths[-1]}' if step != steps[0] else geno_path
        step_output = f'{step_paths[-1]}_{step}' if step != steps[-1] else out_path
        print(f'Running: {step} with input {step_input} and output: {step_output}')

    #     # ancestry setup and call
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
            ancestry.subset = args['subset_ancestry']
            out_dict[step] = steps_dict[step]()

            # call ancestry specific steps within each group
            if len(steps_ancestry) > 0:
                for geno, label in zip(out_dict[step]['output']['split_paths'], out_dict[step]['data']['labels_list']):
                    out_dict[label] = execute_pipeline(steps_ancestry, steps_dict, geno, f'{out_path}_{label}', samp_qc, var_qc, ancestry, assoc, args, tmp_dir)

        # keep track of paths
        else:
            step_paths.append(step_output)

        # samp qc setup and call
        if step in samp_steps:
            samp_qc.geno_path = step_input
            samp_qc.out_path = step_output

            # related has more than one parameter
            if step == 'related':
                out_dict[step] = steps_dict[step](related_cutoff=args['related_cutoff'], duplicate_cutoff=args['duplicate_cutoff'],
                                 prune_related=args['prune_related'], prune_duplicated=args['prune_duplicated'])
            
            else:
                out_dict[step] = steps_dict[step](args[step])
        
        # var qc setup and call
        if step in var_steps:
            var_qc.geno_path = step_input
            var_qc.out_path = step_output

            # hwe has more than one parameter
            if step == 'hwe':
                out_dict[step] = steps_dict[step](hwe_threshold=args['hwe'], filter_controls=args['filter_controls'])
            
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
        if (not args['full_output']):
            remove_step_index = step_paths.index(step_output) - 1
            remove_path = step_paths[remove_step_index]
            if os.path.isfile(f'{remove_path}.pgen'):
                os.remove(f'{remove_path}.pgen')
                os.remove(f'{remove_path}.psam')
                os.remove(f'{remove_path}.pvar')

    out_dict['paths'] = step_paths

    return out_dict


def build_metrics_pruned_df(metrics_df, pruned_df, dictionary, ancestry='all'):
    #TODO: Add association output
    for step in ['callrate', 'sex', 'realted', 'het', 'case_control', 'haplotype', 'hwe', 'geno']:
        if step in dictionary.keys():
            qc_step = dictionary[step]['step']
            pf = dictionary[step]['pass']
            ancestry_label = ancestry

            if step in ['callrate', 'sex', 'realted', 'het']:
                level = 'sample'
                samplefile = dictionary[step]['output']['pruned_samples']
                if os.path.isfile(samplefile):
                    pruned = pd.read_csv(samplefile, sep='\t')
                    if pruned.shape[0] > 0:
                        pruned.loc[:,'step'] = step
                        pruned_df = pd.concat([pruned_df, pruned[['#FID','IID','step']]], ignore_index=True)
            else:
                level = 'variant'

            for metric, value in dictionary[step]['metrics'].items():
                tmp_metrics_df = pd.DataFrame({'step':[qc_step], 'pruned_count':[value], 'metric':[metric], 'ancestry':[ancestry_label], 'level':[level], 'pass': [pf]})
                metrics_df = pd.concat([metrics_df, tmp_metrics_df], ignore_index=True)

    return metrics_df, pruned_df


if __name__=='__main__':
    # definte arg parse
    parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')

    # file i/o arguments
    parser.add_argument('--geno_path', type=str, nargs='?', default=None, const=None, help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].', required=True)
    parser.add_argument('--out_path', type=str, nargs='?', default=None, const=None, help='Prefix for output (including path)', required=True)
    parser.add_argument('--full_output', type=str, nargs='?', default='True', const='True', help='Output everything')
    parser.add_argument('--skip_fails', type=str, nargs='?', default='False', const='True', help='Skip up front check for fails')

    # ancerstry arguments
    parser.add_argument('--ancestry', type=str, nargs='?', default='False', const='True', help='Split by ancestry')
    parser.add_argument('--ref_panel', type=str, nargs='?', default=None, const=None, help='Genotype: (string file path). Path to PLINK format reference genotype file, everything before the *.bed/bim/fam.')
    parser.add_argument('--ref_labels', type=str, nargs='?', default=None, const=None, help='tab-separated plink-style IDs with ancestry label (FID  IID label) with no header')
    parser.add_argument('--model', type=str, nargs='?', default=None, const='path', help='Path to pickle file with trained ancestry model for passed reference panel')
    parser.add_argument('--container', type=str, nargs='?', default='False', const='True', help='Run predictions in container')
    parser.add_argument('--singularity', type=str, nargs='?', default='False', const='True', help='Run containerized precitions via singularity')
    parser.add_argument('--subset_ancestry', nargs='*', help='Subset to continue analysis for')

    # sample-level qc arguments
    parser.add_argument('--callrate', type=float, nargs='?', default=None, const=0.02, help='Minimum Callrate threshold for QC')
    parser.add_argument('--sex', nargs='*', help='Sex prune with cutoffs')
    parser.add_argument('--related', type=str, nargs='?', default='False', const='True', help='Relatedness prune')
    parser.add_argument('--related_cutoff', type=float, nargs='?', default=0.0884, const=0.0884, help='Relatedness cutoff')
    parser.add_argument('--duplicate_cutoff', type=float, nargs='?', default=0.354, const=0.354, help='Relatedness cutoff')
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
    parser.add_argument('--all_variant', type=str, nargs='?', default='False', const='True', help='Run all variant-level QC')

    # GWAS and PCA argument
    parser.add_argument('--pca', type=int, nargs='?', default=None, const=10, help='PCA and number of PCs')
    parser.add_argument('--build', type=str, nargs='?', default='hg38', const='hg38', help='Build for PCA')
    parser.add_argument('--gwas', type=str, nargs='?', default='False', const='True', help='Run GWAS')
    parser.add_argument('--covars', type=str, nargs='?', default=None, const=None, help='Path to external covars')
    parser.add_argument('--covar_names', type=str, nargs='?', default=None, const=None, help='Covar names to use from external file')


    # parse args and turn into dict
    args = parser.parse_args()
    args_dict = vars(args)
    print(args_dict)

    # initialize classes
    samp_qc = SampleQC()
    var_qc = VariantQC()
    ancestry = Ancestry()
    assoc = Assoc()

    # ordered steps with their methods to be called
    ordered_steps =  {'ancestry':ancestry.run_ancestry,'callrate':samp_qc.run_callrate_prune,'sex':samp_qc.run_sex_prune,
                    'related':samp_qc.run_related_prune,'het':samp_qc.run_het_prune,'case_control':var_qc.run_case_control_prune,
                    'haplotype':var_qc.run_haplotype_prune,'hwe':var_qc.run_hwe_prune,'geno':var_qc.run_geno_prune,
                    'assoc':assoc.run_association}

    # some up-front editing of pipeline arguments (booleans and lists)
    for step in args_dict:
        if (args_dict[step] == 'True') or (args_dict[step] == 'False'):
            args_dict[step] = bool(args_dict[step] == 'True')

    if args_dict['sex'] is not None:
        if len(args_dict['sex']) == 0:
            args_dict['sex'] = [0.25, 0.75]

        else:
            args_dict['sex'] = [float(i) for i in args_dict['sex']]

    if args_dict['het'] is not None:
        if len(args_dict['het']) == 0:
            args_dict['het'] = [-0.25, 0.25]

        else:
            args_dict['het'] = [float(i) for i in args_dict['sex']]

    # if all sample or all variant called, replace necessary items with defaults
    if args_dict['all_sample']:
        args_dict['callrate'] = 0.02
        args_dict['sex'] = [0.25, 0.75]
        args_dict['related'] = True
        args_dict['het'] = [-0.25, 0.25]
    
    if args_dict['all_variant']:
        args_dict['geno'] = 0.05
        args_dict['case_control'] = 1e-4
        args_dict['haplotype'] = 1e-4
        args_dict['hwe'] = 1e-4
        args_dict['filter_controls'] = True

    print(args_dict)

    # clear log files if repeated out path
    if os.path.exists(f"{args_dict['out_path']}_all_logs.log"):
        os.remove(f"{args_dict['out_path']}_all_logs.log")
    if os.path.exists(f"{args_dict['out_path']}_cleaned_logs.log"):
        os.remove(f"{args_dict['out_path']}_cleaned_logs.log")

    # create empty log files in output directory
    with open(f"{args_dict['out_path']}_all_logs.log", "w") as fp: 
        pass
    with open(f"{args_dict['out_path']}_cleaned_logs.log", "w") as fp: 
        pass
    
    # run data breakdows
    args_dict = upfront_check(args_dict['geno_path'], args_dict)

    # get correct order of steps to run
    run_steps_list = []
    for key in args_dict:
        if (args_dict[key] != None) and args_dict[key]:
            if key in ordered_steps.keys():
                run_steps_list.append(key)
            if ((key == 'pca') or (key == 'gwas')) and ('assoc' not in run_steps_list):
                run_steps_list.append('assoc')
    
    # check run steps and output step
    print(run_steps_list)
    print(f'Output steps: {run_steps_list[-1]}')

    # create tmp dir
    out_dir = os.path.dirname(args_dict['out_path'])
    tmp_dir = tempfile.TemporaryDirectory(suffix='_tmp', prefix='.', dir=out_dir)

    # run pipeline
    out_dict = execute_pipeline(run_steps_list, ordered_steps, args_dict['geno_path'], args_dict['out_path'], samp_qc=samp_qc, var_qc=var_qc, ancestry=ancestry, assoc=assoc, args=args_dict, tmp_dir=tmp_dir)
    
    # build output
    clean_out_dict = dict()
    metrics_df = pd.DataFrame()
    pruned_df = pd.DataFrame()

    if 'ancestry' in out_dict.keys():
        ancestry_counts_df = pd.DataFrame(out_dict['ancestry']['metrics']['predicted_counts']).reset_index()
        ancestry_counts_df.columns = ['label', 'count']
        clean_out_dict['ancestry_counts'] = ancestry_counts_df.to_dict()

        le = out_dict['ancestry']['data']['label_encoder']
        confusion_matrix = out_dict['ancestry']['data']['confusion_matrix']
        confusion_matrix = pd.DataFrame(confusion_matrix)
        confusion_matrix.columns = le.inverse_transform([i for i in range(10)])
        confusion_matrix.index = le.inverse_transform([i for i in range(10)])
        clean_out_dict['confusion_matrix'] = confusion_matrix.to_dict()

        clean_out_dict['ref_pcs'] = out_dict['ancestry']['data']['ref_pcs'].to_dict()
        clean_out_dict['projected_pcs'] = out_dict['ancestry']['data']['projected_pcs'].to_dict()
        clean_out_dict['total_umap'] = out_dict['ancestry']['data']['total_umap'].to_dict()
        clean_out_dict['ref_umap'] = out_dict['ancestry']['data']['ref_umap'].to_dict()
        clean_out_dict['new_samples_umap'] = out_dict['ancestry']['data']['new_samples_umap'].to_dict()
        clean_out_dict['pred_ancestry_labels'] = out_dict['ancestry']['data']['predict_data']['ids'].to_dict()
    
        for ancestry in ['AFR', 'SAS', 'EAS', 'EUR', 'AMR', 'AJ', 'CAS', 'MDE', 'FIN', 'AAC']:
            if ancestry in out_dict.keys():
                metrics_df, pruned_df = build_metrics_pruned_df(metrics_df=metrics_df, pruned_df=pruned_df, dictionary=out_dict[ancestry], ancestry=ancestry)

    else:
        metrics_df, pruned_df = build_metrics_pruned_df(metrics_df=metrics_df, pruned_df=pruned_df, dictionary=out_dict)

    # for weird error with the first sample in pruned file showing up twice when run in tmp file
    pruned_df = pruned_df.drop_duplicates(subset=['#FID','IID'], ignore_index=True)

    clean_out_dict['QC'] = metrics_df.to_dict()
    clean_out_dict['pruned_samples'] = pruned_df.to_dict()

    # dump output to json
    with open(f'{args_dict["out_path"]}.json', 'w') as f:
        json.dump(clean_out_dict, f)

    tmp_dir.cleanup() # to delete directory

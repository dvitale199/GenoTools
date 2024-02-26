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
import tempfile
import pandas as pd
import json
from genotools.pipeline import gt_argparse


def handle_main():
    
    args = gt_argparse()

    args_dict = vars(args)
    
    from numba.core.errors import NumbaDeprecationWarning
    import warnings
    warnings.simplefilter('ignore', category=NumbaDeprecationWarning)

    from umap import UMAP
    from genotools.utils import upfront_check, bfiles_to_pfiles, vcf_to_pfiles, gt_header
    from genotools.qc import SampleQC, VariantQC
    from genotools.ancestry import Ancestry
    from genotools.gwas import Assoc
    from genotools.pipeline import execute_pipeline, build_metrics_pruned_df

    # initialize classes
    samp_qc = SampleQC()
    var_qc = VariantQC()
    ancestry = Ancestry()
    assoc = Assoc()

    # ordered steps with their methods to be called
    ordered_steps =  {'ancestry':ancestry.run_ancestry,'callrate':samp_qc.run_callrate_prune,'sex':samp_qc.run_sex_prune,
                    'related':samp_qc.run_related_prune,'het':samp_qc.run_het_prune,'case_control':var_qc.run_case_control_prune,
                    'haplotype':var_qc.run_haplotype_prune,'hwe':var_qc.run_hwe_prune,'geno':var_qc.run_geno_prune,
                    'ld':var_qc.run_ld_prune,'assoc':assoc.run_association}

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
        
    if args_dict['ld'] is not None:
        if len(args_dict['ld']) == 0:
            args_dict['ld'] = [50, 5, 0.5]

        else:
            args_dict['ld'] = [int(args_dict['ld'][0]), int(args_dict['ld'][1]), float(args_dict['ld'][2])]

    # if all sample or all variant called, replace necessary items with defaults
    if args_dict['all_sample']:
        args_dict['callrate'] = 0.05
        args_dict['sex'] = [0.25, 0.75]
        args_dict['related'] = True
        args_dict['het'] = [-0.25, 0.25]
    
    if args_dict['all_variant']:
        args_dict['geno'] = 0.05
        args_dict['case_control'] = 1e-4
        args_dict['haplotype'] = 1e-4
        args_dict['hwe'] = 1e-4
        args_dict['filter_controls'] = True
        args_dict['ld'] = None
    
    # run data breakdows
    if (args_dict['bfile'] is None) and (args_dict['pfile'] is None) and (args_dict['vcf'] is None):
        raise KeyError('No bfiles, pfiles, or vcf genotypes were provided!')
    
    elif args_dict['bfile'] and (args_dict['pfile'] is None):
        bfiles_to_pfiles(bfile_path=args_dict['bfile'])
        args_dict['geno_path'] = args_dict['bfile']
    
    elif args_dict['vcf'] and (args_dict['pfile'] is None):
        vcf_to_pfiles(vcf_path=args_dict['vcf'])
        args_dict['geno_path'] = args_dict['vcf'].split('.vcf')[0]
    
    else:
        args_dict['geno_path'] = args_dict['pfile']
    
    # run model-side error catching
    if args_dict['model'] and args_dict['container']:
        raise KeyError('Cannot pass a pre-trained model and run predictions using the NeuroBooster array model in a container! Please choose one of the two.')
    elif args_dict['model'] and args_dict['cloud']:
        raise KeyError('Cannot pass a pre-trained model and run predictions in the cloud! Please choose one of the two.')
    elif args_dict['container'] and args_dict['cloud']:
        raise KeyError('Cannot run predictions using the NeuroBooster array model in a container and run predictions in the cloud! Please choose one of the two.')

    args_dict = upfront_check(args_dict['geno_path'], args_dict)

    # clear log files if repeated out path
    if os.path.exists(f"{args_dict['out']}_all_logs.log"):
        os.remove(f"{args_dict['out']}_all_logs.log")
    if os.path.exists(f"{args_dict['out']}_cleaned_logs.log"):
        os.remove(f"{args_dict['out']}_cleaned_logs.log")

    # create empty log files in output directory
    header = gt_header()
    with open(f"{args_dict['out']}_all_logs.log", "w") as fp: 
        fp.write(header)
        fp.write("\n")
    with open(f"{args_dict['out']}_cleaned_logs.log", "w") as fp: 
        pass

    # get correct order of steps to run
    run_steps_list = []
    for key in args_dict:
        if (args_dict[key] != None) and args_dict[key]:
            if key in ordered_steps.keys():
                run_steps_list.append(key)
            if ((key == 'pca') or (key == 'gwas')) and ('assoc' not in run_steps_list):
                run_steps_list.append('assoc')
    
    # check run steps and output step
    if len(run_steps_list) == 0:
        raise KeyError('No main Ancestry, QC, or GWAS flags were used.')
    else:
        print(f'Output steps: {run_steps_list[-1]}')

    # create tmp dir
    out_dir = os.path.dirname(args_dict['out'])
    tmp_dir = tempfile.TemporaryDirectory(suffix='_tmp', prefix='.', dir=out_dir)

    # run pipeline
    out_dict = execute_pipeline(run_steps_list, ordered_steps, args_dict['geno_path'], args_dict['out'], samp_qc=samp_qc, var_qc=var_qc, ancestry=ancestry, assoc=assoc, args=args_dict, tmp_dir=tmp_dir)
    
    # build output
    clean_out_dict = dict()
    metrics_df = pd.DataFrame()
    pruned_df = pd.DataFrame()
    gwas_df = pd.DataFrame()
    related_df = pd.DataFrame()

    # input psam in out_dict
    if os.path.isfile(f'{args_dict["geno_path"]}.psam'):
        psam = pd.read_csv(f'{args_dict["geno_path"]}.psam', sep='\s+')
        clean_out_dict['input_samples'] = psam.to_dict()

    if 'ancestry' in out_dict.keys():
        ancestry_counts_df = pd.DataFrame(out_dict['ancestry']['metrics']['predicted_counts']).reset_index()
        ancestry_counts_df.columns = ['label', 'count']
        clean_out_dict['ancestry_counts'] = ancestry_counts_df.to_dict()

        ancestry_labels_df = pd.DataFrame(out_dict['ancestry']['data']['predict_data']['ids'])
        clean_out_dict['ancestry_labels'] = ancestry_labels_df.to_dict()

        le = out_dict['ancestry']['data']['label_encoder']
        confusion_matrix = out_dict['ancestry']['data']['confusion_matrix']
        confusion_matrix = pd.DataFrame(confusion_matrix)
        confusion_matrix.columns = le.inverse_transform([i for i in range(10)])
        confusion_matrix.index = le.inverse_transform([i for i in range(10)])
        clean_out_dict['confusion_matrix'] = confusion_matrix.to_dict()

        clean_out_dict['test_accuracy'] = out_dict['ancestry']['metrics']['test_accuracy']

        clean_out_dict['ref_pcs'] = out_dict['ancestry']['data']['ref_pcs'].to_dict()
        clean_out_dict['projected_pcs'] = out_dict['ancestry']['data']['projected_pcs'].to_dict()
        clean_out_dict['total_umap'] = out_dict['ancestry']['data']['total_umap'].to_dict()
        clean_out_dict['ref_umap'] = out_dict['ancestry']['data']['ref_umap'].to_dict()
        clean_out_dict['new_samples_umap'] = out_dict['ancestry']['data']['new_samples_umap'].to_dict()
    
        for ancestry in ['AFR', 'SAS', 'EAS', 'EUR', 'AMR', 'AJ', 'CAS', 'MDE', 'FIN', 'AAC', 'CAH']:
            if ancestry in out_dict.keys():
                metrics_df, pruned_df, gwas_df, related_df = build_metrics_pruned_df(metrics_df=metrics_df, pruned_df=pruned_df, gwas_df=gwas_df, related_df=related_df, dictionary=out_dict[ancestry], out=args_dict['out'], ancestry=ancestry)

    else:
        metrics_df, pruned_df, gwas_df, related_df = build_metrics_pruned_df(metrics_df=metrics_df, pruned_df=pruned_df, gwas_df=gwas_df, related_df=related_df, dictionary=out_dict, out=args_dict['out'])

    # for weird error with the first sample in pruned file showing up twice when run in tmp file
    pruned_df = pruned_df.drop_duplicates(subset=['#FID','IID'], ignore_index=True)
    pruned_df = pruned_df.rename({'#FID':'FID'}, axis=1)

    # ensure no empty df is being output to JSON
    output_dfs = {'QC':metrics_df, 'GWAS':gwas_df, 'pruned_samples':pruned_df, 'related_samples':related_df}

    for df in output_dfs:
        if not output_dfs[df].empty:
            if df == 'pruned_samples':
                if 'ancestry_labels' in list(clean_out_dict.keys()):
                    labels = pd.DataFrame(clean_out_dict['ancestry_labels'])
                    labeled_pruned_df = output_dfs[df].merge(labels[['FID','IID','label']], how='left', on=['FID','IID'])

                    ancestry_pruned_df = out_dict['ancestry']['data']['pruned_samples']
                    full_labeled_pruned_df = pd.concat([ancestry_pruned_df, labeled_pruned_df], axis=0, ignore_index=True)
                    clean_out_dict[df] = full_labeled_pruned_df.to_dict()
                else:
                    clean_out_dict[df] = output_dfs[df].to_dict()
            
            else:
                clean_out_dict[df] = output_dfs[df].to_dict()

    # for df in [metrics_df, pruned_df, gwas_df, related_df]:
    #     if not df.empty:
    #         clean_out_dict['QC'] = metrics_df.to_dict()
    #         clean_out_dict['GWAS'] = gwas_df.to_dict()

    #         if 'ancestry_labels' in list(clean_out_dict.keys()):
    #             labels = pd.DataFrame(clean_out_dict['ancestry_labels'])
    #             labeled_pruned_df = pruned_df.merge(labels[['FID','IID','label']], how='left', on=['FID','IID'])

    #             ancestry_pruned_df = out_dict['ancestry']['data']['pruned_samples']
    #             full_labeled_pruned_df = pd.concat([ancestry_pruned_df, labeled_pruned_df], axis=0, ignore_index=True)
    #             clean_out_dict['pruned_samples'] = full_labeled_pruned_df.to_dict()
    #         else:
    #             clean_out_dict['pruned_samples'] = pruned_df.to_dict()

    # dump output to json
    with open(f'{args_dict["out"]}.json', 'w') as f:
        json.dump(clean_out_dict, f)

    tmp_dir.cleanup() # to delete directory


if __name__ == "__main__":
    handle_main()
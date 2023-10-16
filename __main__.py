import subprocess
import sys
import argparse
import os
import tempfile

from genotools.utils import shell_do
from genotools.qc import SampleQC, VariantQC
from genotools.ancestry import Ancestry

def execute_pipeline(steps, steps_dict, geno_path, out_path, samp_qc, var_qc, ancestry, args):
    # to know which class to call
    samp_steps = ['callrate','sex','related','het']
    var_steps = ['case_control','haplotype','hwe','geno']

    # if full output requested, go to out path
    if args['full_output']:
        step_paths = [out_path]

    # otherwise tmpdir
    else:
        out_dir = os.path.dirname(out_path)
        tmp_dir = tempfile.TemporaryDirectory(suffix='_tmp', prefix='.', dir=out_dir)
        step_paths = [f'{tmp_dir.name}/out']

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

        # keep track of paths
        step_paths.append(step_output)

        # ancestry setup and call
        if step == 'ancestry':
            ancestry.geno_path = step_input
            ancestry.out_path = step_output
            ancestry.ref_panel = args['ref_panel']
            ancestry.ref_labels = args['ref_labels']
            ancestry.model_path = args['model']
            ancestry.containerized = args['container']
            ancestry.singularity = args['singularity']
            ancestry.subset = args['subset']
            out_dict[step] = steps_dict[step]()

            # call ancestry specific steps within each group
            if len(steps_ancestry) > 0:
                for geno, label in zip(out_dict[step]['output']['split_paths'], out_dict[step]['data']['labels_list']):
                    out_dict[label] = execute_pipeline(steps_ancestry, steps_dict, geno, f'{out_path}_{label}', samp_qc, var_qc, ancestry, args)

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
    
    out_dict['paths'] = step_paths

    return out_dict


if __name__=='__main__':
    # definte arg parse
    parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')

    # file i/o arguments
    parser.add_argument('--geno_path', type=str, nargs='?', default=None, const=None, help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].', required=True)
    parser.add_argument('--out_path', type=str, nargs='?', default=None, const=None, help='Prefix for output (including path)', required=True)
    parser.add_argument('--full_output', type=str, nargs='?', default='True', const='True', help='Output everything')

    # ancerstry arguments
    parser.add_argument('--ancestry', type=str, nargs='?', default='False', const='True', help='Split by ancestry')
    parser.add_argument('--ref_panel', type=str, nargs='?', default=None, const=None, help='Genotype: (string file path). Path to PLINK format reference genotype file, everything before the *.bed/bim/fam.')
    parser.add_argument('--ref_labels', type=str, nargs='?', default=None, const=None, help='tab-separated plink-style IDs with ancestry label (FID  IID label) with no header')
    parser.add_argument('--model', type=str, nargs='?', default=None, const='path', help='Path to pickle file with trained ancestry model for passed reference panel')
    parser.add_argument('--container', type=str, nargs='?', default='False', const='True', help='Run predictions in container')
    parser.add_argument('--singularity', type=str, nargs='?', default='False', const='True', help='Run containerized precitions via singularity')
    parser.add_argument('--subset', nargs='*', help='Subset to continue analysis for')

    # sample-level qc arguments
    parser.add_argument('--callrate', type=float, nargs='?', default=None, const=0.05, help='Minimum Callrate threshold for QC')
    parser.add_argument('--sex', nargs='*', help='Sex prune with cutoffs')
    parser.add_argument('--related', type=str, nargs='?', default='False', const='True', help='Relatedness prune')
    parser.add_argument('--related_cutoff', type=float, nargs='?', default=0.0884, const=0.0884, help='Relatedness cutoff')
    parser.add_argument('--duplicate_cutoff', type=float, nargs='?', default=0.354, const=0.354, help='Relatedness cutoff')
    parser.add_argument('--prune_related', type=str, nargs='?', default='False', const='True', help='Relatedness prune')
    parser.add_argument('--prune_duplicated', type=str, nargs='?', default='True', const='True', help='Relatedness prune')
    parser.add_argument('--het', nargs='*', help='Het prune with cutoffs')

    # variant-level qc arguments
    parser.add_argument('--geno', type=float, nargs='?', default=None, const=0.05, help='Minimum Missingness threshold for QC')
    parser.add_argument('--case_control', type=float, nargs='?', default=None, const=1e-4, help='Case control prune')
    parser.add_argument('--haplotype', type=float, nargs='?', default=None, const=1e-4, help='Haplotype prune')
    parser.add_argument('--hwe', type=float, nargs='?', default=None, const=1e-4, help='HWE pruning')
    parser.add_argument('--filter_controls', type=str, nargs='?', default='False', const='True', help='Control filter for HWE prune')

    # parse args and turn into dict
    args = parser.parse_args()
    args_dict = vars(args)
    print(args_dict)

    # initialize classes
    samp_qc = SampleQC()
    var_qc = VariantQC()
    ancestry = Ancestry()

    # ordered steps with their methods to be called
    ordered_steps =  {'ancestry':ancestry.run_ancestry,'callrate':samp_qc.run_callrate_prune,'sex':samp_qc.run_sex_prune,
                    'related':samp_qc.run_related_prune,'het':samp_qc.run_het_prune,'case_control':var_qc.run_case_control_prune,
                    'haplotype':var_qc.run_haplotype_prune,'hwe':var_qc.run_hwe_prune,'geno':var_qc.run_geno_prune}
    
    # get correct order of steps to run
    run_steps_list = []
    for key in args_dict:
        if (args_dict[key] != None) and (args_dict[key] != 'False') and (key in ordered_steps.keys()):
            run_steps_list.append(key)

    # check run steps and output step
    print(run_steps_list)
    print(f'Output steps: {run_steps_list[-1]}')

    # some up-front editing of pipeline arguments (booleans and lists)
    args_dict['full_output'] = bool(args_dict['full_output'] == 'True')
    args_dict['ancestry'] = bool(args_dict['ancestry'] == 'True')
    args_dict['container'] = bool(args_dict['container'] == 'True')
    args_dict['singularity'] = bool(args_dict['singularity'] == 'True')
    args_dict['related'] = bool(args_dict['related'] == 'True')
    args_dict['prune_related'] = bool(args_dict['prune_related'] == 'True')
    args_dict['prune_duplicated'] = bool(args_dict['prune_duplicated'] == 'True')
    args_dict['filter_controls'] =  bool(args_dict['filter_controls'] == 'True')

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

    # run pipeline
    execute_pipeline(run_steps_list, ordered_steps, args_dict['geno_path'], args_dict['out_path'], samp_qc=samp_qc, var_qc=var_qc, ancestry=ancestry, args=args_dict)
            
    # tmp_dir.cleanup() # to delete directory

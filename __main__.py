import subprocess
import sys
import argparse
import os
import tempfile

from genotools.utils import shell_do
    

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
    parser.add_argument('--geno_path', type=str, nargs='?', default=None, const=None, help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].', required=True)
    parser.add_argument('--out_path', type=str, nargs='?', default=None, const=None, help='Prefix for output (including path)', required=True)

    parser.add_argument('--ancestry', type=str, nargs='?', default='False', const='True', help='Split by ancestry')
    parser.add_argument('--ref_panel', type=str, nargs='?', default=None, const=None, help='Genotype: (string file path). Path to PLINK format reference genotype file, everything before the *.bed/bim/fam.')
    parser.add_argument('--ref_labels', type=str, nargs='?', default=None, const=None, help='tab-separated plink-style IDs with ancestry label (FID  IID label) with no header')
    parser.add_argument('--model', type=str, nargs='?', default=None, const='path', help='Path to pickle file with trained ancestry model for passed reference panel')
    parser.add_argument('--container', type=str, nargs='?', default='False', const='True', help='Run predictions in container')
    parser.add_argument('--singularity', type=str, nargs='?', default='False', const='True', help='Run containerized precitions via singularity')

    parser.add_argument('--callrate', type=float, nargs='?', default=None, const=0.05, help='Minimum Callrate threshold for QC')
    parser.add_argument('--sex', nargs='*', help='Sex prune with cutoffs')
    parser.add_argument('--related', type=str, nargs='?', default='False', const='True', help='Relatedness prune')
    parser.add_argument('--related_cutoff', type=float, nargs='?', default=0.0884, const=0.0884, help='Relatedness cutoff')
    parser.add_argument('--duplicate_cutoff', type=float, nargs='?', default=0.354, const=0.354, help='Relatedness cutoff')
    parser.add_argument('--prune_related', type=str, nargs='?', default='False', const='True', help='Relatedness prune')
    parser.add_argument('--prune_duplicated', type=str, nargs='?', default='True', const='False', help='Relatedness prune')
    parser.add_argument('--het', nargs='*', help='Het prune with cutoffs')

    parser.add_argument('--geno', type=float, nargs='?', default=None, const=0.05, help='Minimum Missingness threshold for QC')
    parser.add_argument('--case_control', type=float, nargs='?', default=None, const=1e-4, help='Case control prune')
    parser.add_argument('--haplotype', type=float, nargs='?', default=None, const=1e-4, help='Haplotype prune')
    parser.add_argument('--hwe', type=float, nargs='?', default=None, const=1e-4, help='HWE pruning')


    # default_list

    args = parser.parse_args()
    print(args)
    passed_steps_dict = vars(args)

    # order this list based on passed steps
    ordered_steps_list = ['ancestry', 'callrate', 'sex', 'related', 'het', 'case_control', 'haplotype', 'hwe', 'geno']
    run_steps_list = []

    for key in passed_steps_dict:
        print(key)
        if (passed_steps_dict[key] != None) and (passed_steps_dict[key] != 'False') and (key in ordered_steps_list):
            run_steps_list.append(key)

    print(run_steps_list)
    print(f'Output steps: {run_steps_list[-1]}')


    geno_path = args.geno_path
    out_path = args.out_path

    ancestry = bool(args.ancestry == 'True')
    ref_panel = args.ref_panel
    ref_labels = args.ref_labels
    model = args.model
    container = bool(args.container == 'True')
    singularity = bool(args.singularity == 'True')

    callrate = args.callrate
    sex = args.sex
    related = bool(args.related == 'True')
    related_cutoff = args.related_cutoff
    duplicate_cutoff = args.duplicate_cutoff
    prune_related = bool(args.prune_related == 'True')
    prune_duplicated = bool(args.prune_duplicated == 'True')
    het = args.het

    geno = args.geno
    case_control = args.case_control
    haplotype = args.haplotype
    hwe = args.hwe

    # I/O workflow
    if not geno_path:
        raise Exception('A path to genotype files for processing must be provided!')

    if not out_path:
        raise Exception('A path for output must be provided!')

    out_dir = os.path.dirname(out_path)
    # tmp_dir = tempfile.TemporaryDirectory(suffix='_tmp', prefix='.', dir={out_dir})
    # tmp_dir.name # to feed to modules as out_path 


    # ancestry workflow
    if ancestry:
        if not (ref_panel and ref_labels):
            raise Exception('A reference panel and reference panel labels must be provided to predict ancestry!')
        
        if model and container:
            raise Warning('Model path provided and containerized predictions requested! Defaulting to containerized predictions!')
        
        ##### Can check if singularity should be run #####
        # if container and not singularity:
        #     shell_do('docker --version')

    else:
        if ref_panel or ref_labels or model or container or singularity:
            raise Warning('Ancestry-specific parameters passed, but ancestry was not called!')


    # sample-level qc workflow
    if callrate:
        print('CALL CALLRATE PRUNE')

    if sex is not None:
        if len(sex) == 0:
            print('CALL SEX PRUNE WITH DEFAULTS')
        
        elif (len(sex) > 0) and (len(sex) != 2):
            raise Exception('Need to pass two cutoff values when customizing sex prune!')

        else:
            sex = [float(i) for i in sex]
            print(sex)
            print('CALL SEX PRUNE')
    
    if related:
        print('CALL RELATED PRUNE')

    if het is not None:
        if len(het) == 0:
            print('CALL HET PRUNE WITH DEFAULTS')
        
        elif (len(het) > 0) and (len(het) != 2):
            raise Exception('Need to pass two cutoff values when customizing het prune!')

        else:
            het = [float(i) for i in het]
            print(het)
            print('CALL HET PRUNE')


    # variant-level qc workflow
    if geno:
        print('CALL GENO')

    if case_control:
        print('CALL CASE CONTROL')

    if haplotype:
        print('CALL HAPLOTYPE')

    if hwe:
        print('CALL HWE')
            
    # tmp_dir.cleanup() # to delete directory

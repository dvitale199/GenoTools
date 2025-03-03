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


# rough pipeline for sample-level QC for WGS data
def handle_wgs():

    args = gt_argparse()

    args_dict = vars(args)

    from numba.core.errors import NumbaDeprecationWarning
    import warnings
    warnings.simplefilter('ignore', category=NumbaDeprecationWarning)

    from umap import UMAP
    from genotools.utils import upfront_check, bfiles_to_pfiles, vcf_to_pfiles, gt_header
    from genotools.wgs_qc import WholeGenomeSeqQC
    from genotools.ancestry import Ancestry
    from genotools.gwas import Assoc
    from genotools.pipeline import execute_pipeline, build_metrics_pruned_df

    # some up-front editing of pipeline arguments (booleans and lists)
    for step in args_dict:
        if (args_dict[step] == 'True') or (args_dict[step] == 'False'):
            args_dict[step] = bool(args_dict[step] == 'True')

    if args_dict['wgs']:
        args_dict['freemix'] = True
        args_dict['coverage'] = True
        args_dict['titv'] = True
        args_dict['wgs_het'] = [-0.25, 0.25]
        args_dict['wgs_callrate'] = True
        args_dict['wgs_sex'] = [0.25, 0.75]
        args_dict['wgs_relatedness'] = True

    # currently just assuming that files in shard_dir is in plink2.0 format
    if args_dict['shards_dir'] is None:
        raise KeyError('No input directory was provided!')

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

    # initialize class
    wgs_qc = WholeGenomeSeqQC(shards_dir=args_dict['shards_dir'], out_path=args_dict['out'], \
                                keep_all=args_dict['keep_all'], slurm=args_dict['slurm'], slurm_user=args_dict['slurm_user'], \
                                shard_key_path=args_dict['shard_key'], \
                                preBqsr_path=args_dict['preBqsr'], \
                                wgs_metrics_path=args_dict['wgs_metrics'], \
                                variant_calling_summary_metrics_path=args_dict['variant_calling_summary_metrics'], \
                                ref_variants_path=args_dict['ref_variants'])

    # ordered steps with their methods to be called
    ordered_steps =  {'freemix':wgs_qc.run_freemix_check,'coverage':wgs_qc.run_coverage_check,
                    'titv':wgs_qc.run_titv_check,'wgs_het':wgs_qc.merge_sample_het,'wgs_callrate':wgs_qc.merge_sample_callrate,
                    'wgs_sex':wgs_qc.run_sex_check,'wgs_relatedness':wgs_qc.run_relatedness}

    # assuming all steps of wgs sample qc are being run
    steps = ['freemix', 'coverage', 'titv', 'wgs_het', 'wgs_callrate', 'wgs_sex', 'wgs_relatedness']

    # do freemix check
    ordered_steps['freemix']()

    # do coverage check
    ordered_steps['coverage']()

    # do titv ratio check
    ordered_steps['titv']()

    # do wgs het check
    ordered_steps['wgs_het']()

    # do wgs callrate check
    ordered_steps['wgs_callrate']()

    # do wgs sex check
    ordered_steps['wgs_sex']()

    # do wgs relatedness check
    ordered_steps['wgs_relatedness']()
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

    # initialize classes
    wgs_qc = WholeGenomeSeqQC()
    ancestry = Ancestry()
    assoc = Assoc()

    # ordered steps with their methods to be called
    ordered_steps =  {'ancestry':ancestry.run_ancestry,'freemix':wgs_qc.run_freemix_check,'coverage':wgs_qc.run_coverage_check,
                    'titv':wgs_qc.run_titv_check,'wgs_het':wgs_qc.merge_sample_het,'wgs_callrate':wgs_qc.merge_sample_callrate,
                    'wgs_sex':wgs_qc.run_sex_check,'wgs_relatedness':wgs_qc.run_relatedness,'assoc':assoc.run_association}

    # some up-front editing of pipeline arguments (booleans and lists)
    for step in args_dict:
        if (args_dict[step] == 'True') or (args_dict[step] == 'False'):
            args_dict[step] = bool(args_dict[step] == 'True')

    if args_dict['wgs_sex'] is not None:
        if len(args_dict['wgs_sex']) == 0:
            args_dict['wgs_sex'] = [0.25, 0.75]

        else:
            args_dict['wgs_sex'] = [float(i) for i in args_dict['wgs_sex']]

    if args_dict['wgs_het'] is not None:
        if len(args_dict['wgs_het']) == 0:
            args_dict['wgs_het'] = [-0.25, 0.25]

        else:
            args_dict['wgs_het'] = [float(i) for i in args_dict['wgs_het']]
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

def handle_wgs_main():

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
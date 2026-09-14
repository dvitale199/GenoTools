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


import argparse
import os
import sys
import hashlib
import requests
import zipfile
from tqdm import tqdm

# Downloadable artifacts, keyed by the name users pass to --ref / --model. The
# value is the md5 of `{url_base}/{kind}/{name}.zip`, checked after download.
#
# Ancestry models come in two mutually incompatible formats. A 1.x model is a
# bare sklearn Pipeline pickle; 2.x writes a directory (pipeline.pkl,
# common_snps.txt, metadata.json, requirements.txt) and rejects a 1.x model by
# design. The format is recorded here and reported at download time.
REF_PANELS = {
    "1kg_30x_hgdp_ashk_ref_panel": "6cf0764ae6e99f60127e42b12b4af5d7",
}
DEFAULT_REF = "1kg_30x_hgdp_ashk_ref_panel"

MODELS = {
    # name: (md5, format, description)
    "nba_gp2_r12": (
        "8972607acae64f3147959281cc097dc9",
        "2.x",
        "NeuroBooster array, trained on GP2 release 12 (43,173 SNPs, 10 labels)",
    ),
}
DEFAULT_MODEL = "nba_gp2_r12"

URL_BASE = "https://storage.googleapis.com/genotools_refs"
ARCHIVE_URL_BASE = f"{URL_BASE}/models/archive"

# Retired names, and why. `nba_v1`, `nba_v2` and `neurochip_v1` are 1.x
# pickles: `AncestryModel.load` rejects that format by design, so 2.x could
# never load them and serving them only offered users a download that could
# not work. They also predate `core.provenance`, so they record no library
# versions -- and sklearn can change an estimator's behaviour between versions
# without saying so, which is exactly the silent drift provenance exists to
# catch.
#
# The archives are not deleted. They remain under `ARCHIVE_URL_BASE` so an
# analysis pinned to a 1.x GenoTools can still be reproduced; they are simply
# no longer served by name here. Asking for one gets the message in
# `resolve_name` rather than a bare "Unknown model", matching how
# `_REMOVED_FLAGS` treats removed CLI flags.
RETIRED_MODELS = {
    "nba_v1": "1.x pickle (NeuroBooster array)",
    "nba_v2": "1.x pickle (NeuroBooster array)",
    "neurochip_v1": "1.x pickle (NeuroChip array)",
}


def resolve_name(requested, catalogue, default, kind):
    """Map a user-supplied artifact name to a key in `catalogue`.

    `default` (and the literal string "default", which the CLI help has always
    advertised) resolve to the default artifact. An unknown name raises with the
    available names, rather than the bare KeyError this used to produce.
    """
    name = default if requested in (None, "default") else requested
    if name not in catalogue:
        available = ", ".join(sorted(catalogue))
        if kind == "model" and name in RETIRED_MODELS:
            raise SystemExit(
                f"The {name!r} model has been retired -- it is a "
                f"{RETIRED_MODELS[name]}, and GenoTools 2.x cannot load that "
                f"format. Use {default!r} instead, or train your own with "
                f"--ref-panel/--ref-labels. The archive is still available at "
                f"{ARCHIVE_URL_BASE}/{name}.zip for reproducing an analysis "
                f"pinned to GenoTools 1.x."
            )
        raise SystemExit(
            f"Unknown {kind} {name!r}. Available: {available}."
        )
    return name


def compute_checksum(file_path):
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as file:
        for chunk in iter(lambda: file.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def validate_checksum(file_path, checksum):
    calculated_checksum = compute_checksum(file_path)
    return calculated_checksum == checksum

def download_data_from_gcs(url, destination_file_path, force=False):
    if os.path.exists(destination_file_path) and not force:
        print(f"File already exists at {destination_file_path}")
        return

    response = requests.get(url, stream=True)
    total_size_in_bytes = int(response.headers.get('content-length', 0))
    
    if response.status_code == 200:
        with open(destination_file_path, 'wb') as file, \
             tqdm(
                desc=destination_file_path,
                total=total_size_in_bytes,
                unit='iB',
                unit_scale=True,
                unit_divisor=1024
             ) as bar:
                for chunk in response.iter_content(chunk_size=1024):
                    file.write(chunk)
                    bar.update(len(chunk))
        print(f"File downloaded to {destination_file_path}")
    else:
        response.raise_for_status()

def unzip_file(zip_file_path, destination_dir):
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        zip_ref.extractall(destination_dir)
    print(f"Extracted {zip_file_path} to {destination_dir}")

def handle_download():
    parser = argparse.ArgumentParser(description="Download, validate, and unzip reference data")
    default_destination = os.path.expanduser("~/.genotools/ref")
    parser.add_argument('--destination', type=str, default=default_destination, help="Local destination directory for the download (default: ~/.genotools/refs)")
    parser.add_argument('--model', type=str, help=f"Ancestry model to download (default: {DEFAULT_MODEL}). Available: {', '.join(sorted(MODELS))}. Only 2.x-format models load in GenoTools 2.x")
    parser.add_argument('--ref', type=str, help=f"Reference panel to download (default: {DEFAULT_REF}). Available: {', '.join(sorted(REF_PANELS))}")

    args = parser.parse_args()

    url_base = URL_BASE
    download_ref = args.ref is not None or (args.model is None and args.ref is None)
    download_model = args.model is not None or (args.model is None and args.ref is None)

    if download_ref:
        ref = resolve_name(args.ref, REF_PANELS, DEFAULT_REF, "reference panel")
        url = f"{url_base}/ref_panel/{ref}.zip"
        checksum = REF_PANELS[ref]
        ref_panel_path = f'{args.destination}/ref_panel'
        print(f'Pulling reference panel {ref}')
        os.makedirs(ref_panel_path, exist_ok=True)
        destination_file_path = os.path.join(ref_panel_path, os.path.basename(url))

        if os.path.exists(destination_file_path) and validate_checksum(destination_file_path, checksum):
            print(f"Reference panel already downloaded and validated: {destination_file_path}")
        else:
            # force=True: reaching here means either there is no local copy or
            # the one there failed validation, and a stale copy that keeps
            # failing is not worth preserving. Without this, re-publishing an
            # artifact strands everyone holding the previous one.
            download_data_from_gcs(url, destination_file_path, force=True)
            if not validate_checksum(destination_file_path, checksum):
                print(
                    f"Error: checksum mismatch for {destination_file_path}. "
                    f"Expected {checksum}. Delete that file and retry; if it "
                    f"persists, the published artifact and this version of "
                    f"GenoTools disagree."
                )
                sys.exit(1)
            unzip_file(destination_file_path, ref_panel_path)

    if download_model:
        model = resolve_name(args.model, MODELS, DEFAULT_MODEL, "model")
        checksum, model_format, description = MODELS[model]
        url = f"{url_base}/models/{model}.zip"
        model_path = f'{args.destination}/models'
        # The 1.x warning that used to sit here is gone with the 1.x names:
        # every served model is 2.x now, so the branch could never fire, and a
        # guard that cannot fire reads as protection it is not providing
        # (REFACTOR.md round 22). The format is still reported, because a name
        # alone does not tell you whether a model will load.
        print(f'Pulling model: {model} ({model_format} format, {description})')
        os.makedirs(model_path, exist_ok=True)
        destination_file_path = os.path.join(model_path, os.path.basename(url))

        if os.path.exists(destination_file_path) and validate_checksum(destination_file_path, checksum):
            print(f"Model already downloaded and validated: {destination_file_path}")
        else:
            # force=True: reaching here means either there is no local copy or
            # the one there failed validation, and a stale copy that keeps
            # failing is not worth preserving. Without this, re-publishing an
            # artifact strands everyone holding the previous one.
            download_data_from_gcs(url, destination_file_path, force=True)
            if not validate_checksum(destination_file_path, checksum):
                print(
                    f"Error: checksum mismatch for {destination_file_path}. "
                    f"Expected {checksum}. Delete that file and retry; if it "
                    f"persists, the published artifact and this version of "
                    f"GenoTools disagree."
                )
                sys.exit(1)
            unzip_file(destination_file_path, model_path)

if __name__ == "__main__":
    handle_download()
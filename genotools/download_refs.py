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
# design. Serving both means a name alone cannot tell you whether the model
# will load, so the format is recorded here and reported at download time.
REF_PANELS = {
    "1kg_30x_hgdp_ashk_ref_panel": "6cf0764ae6e99f60127e42b12b4af5d7",
}
DEFAULT_REF = "1kg_30x_hgdp_ashk_ref_panel"

MODELS = {
    # name: (md5, format, description)
    "nba_v1": ("755042b6a1e600a06b10352d42c57d20", "1.x", "NeuroBooster array"),
    "nba_v2": ("7618cd9be74a6f8da96ae99016851cce", "1.x", "NeuroBooster array"),
    "neurochip_v1": ("8825d8b490bab62d91752ba64e960c2d", "1.x", "NeuroChip array"),
    "nba_gp2_r12": (
        "2ff0af5218cc7f93bad737b821924438",
        "2.x",
        "NeuroBooster array, trained on GP2 release 12 (43,173 SNPs, 10 labels)",
    ),
}
DEFAULT_MODEL = "nba_gp2_r12"


def resolve_name(requested, catalogue, default, kind):
    """Map a user-supplied artifact name to a key in `catalogue`.

    `default` (and the literal string "default", which the CLI help has always
    advertised) resolve to the default artifact. An unknown name raises with the
    available names, rather than the bare KeyError this used to produce.
    """
    name = default if requested in (None, "default") else requested
    if name not in catalogue:
        available = ", ".join(sorted(catalogue))
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

def download_data_from_gcs(url, destination_file_path):
    if os.path.exists(destination_file_path):
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

    url_base = "https://storage.googleapis.com/genotools_refs"
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
            download_data_from_gcs(url, destination_file_path)
            if not validate_checksum(destination_file_path, checksum):
                print("Error: Checksum validation failed for reference panel.")
                sys.exit(1)
            unzip_file(destination_file_path, ref_panel_path)

    if download_model:
        model = resolve_name(args.model, MODELS, DEFAULT_MODEL, "model")
        checksum, model_format, description = MODELS[model]
        url = f"{url_base}/models/{model}.zip"
        if model_format == "1.x":
            print(
                f"Warning: {model} is a GenoTools 1.x model ({description}) and "
                f"cannot be loaded by 2.x. Pass --model {DEFAULT_MODEL} for the "
                f"2.x model, or train your own with --ref-panel/--ref-labels."
            )
        model_path = f'{args.destination}/models'
        print(f'Pulling model: {model}')
        os.makedirs(model_path, exist_ok=True)
        destination_file_path = os.path.join(model_path, os.path.basename(url))

        if os.path.exists(destination_file_path) and validate_checksum(destination_file_path, checksum):
            print(f"Model already downloaded and validated: {destination_file_path}")
        else:
            download_data_from_gcs(url, destination_file_path)
            if not validate_checksum(destination_file_path, checksum):
                print("Error: Checksum validation failed for model.")
                sys.exit(1)
            unzip_file(destination_file_path, model_path)

if __name__ == "__main__":
    handle_download()
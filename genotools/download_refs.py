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
    parser.add_argument('--model', type=str, help="Version of the model to use (provide 'default' for default model)")
    parser.add_argument('--ref', type=str, help="Version of the reference panel to use (provide 'default' for default reference panel)")

    args = parser.parse_args()

    url_base = "https://storage.googleapis.com/genotools_refs"
    download_ref = args.ref is not None or (args.model is None and args.ref is None)
    download_model = args.model is not None or (args.model is None and args.ref is None)

    if download_ref:
        checksums_dict = {
            "1kg_30x_hgdp_ashk_ref_panel": "e97728c6cc0e672cf783cb4ca338b184"
        }

        ref = args.ref if args.ref else "1kg_30x_hgdp_ashk_ref_panel"
        url = f"{url_base}/ref_panel/{ref}.zip"
        checksum = checksums_dict[ref]
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
        checksums_dict = {
            'nba_v1': 'e167af2a192bd7460f23a0868be4d7bd',
            'neurochip_v1': '8f46bab3bb1226a36f8d3b9b7dab4b1a'
        }

        model = args.model if args.model else "nba_v1"
        url = f"{url_base}/models/{model}.zip"
        checksum = checksums_dict[model]
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
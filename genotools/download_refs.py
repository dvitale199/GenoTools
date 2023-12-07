import argparse
import os
import sys
import hashlib
import requests
import zipfile


def validate_checksum(file_path, checksum):
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    if hash_md5.hexdigest() != checksum:
        return False
    return True


def download_data_from_gcs(public_url, destination_file_name):
    """Downloads a file from a public GCS bucket."""
    if not os.path.exists(os.path.dirname(destination_file_name)):
        os.makedirs(os.path.dirname(destination_file_name), exist_ok=True)

    response = requests.get(public_url, stream=True)
    if response.status_code == 200:
        with open(destination_file_name, 'wb') as f:
            for chunk in response.iter_content(chunk_size=128):
                f.write(chunk)
        print(f"File downloaded to {destination_file_name}")

    else:
        response.raise_for_status()


def unzip_file(zip_file_path, destination_dir):
    """Unzips a file to a specified destination directory."""
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        zip_ref.extractall(destination_dir)
    print(f"Extracted {zip_file_path} to {destination_dir}")


# Main function to handle the download process
def handle_download():
    parser = argparse.ArgumentParser(description="Download, validate, and unzip reference data")
    parser.add_argument('--destination', type=str, required=True, help="Local destination directory for the download")
    parser.add_argument('--model', type=str, default="nba_v1", help="Version of the model to use")
    parser.add_argument('--ref', type=str, default="1kg_30x_hgdp_ashk_ref_panel", help="Version of the reference panel to use")

    args = parser.parse_args()

    url_base = "https://storage.googleapis.com/genotools_refs"

    if args.ref:
        url = f"{url_base}/ref_panel/{args.ref}.zip"
        checksum = "fd3d79b9e1c0054b10881fa130eb9b21"

    if args.model:
        url = f"{url_base}/models/{args.model}.zip"
        checksum = None

    # Download the file
    destination_file_path = os.path.join(args.destination, os.path.basename(args.url))
    download_data_from_gcs(args.url, destination_file_path)

    # Validate the checksum if provided
    if args.checksum and not validate_checksum(destination_file_path, checksum):
        print("Error: Checksum validation failed.")
        sys.exit(1)

    # Unzip the file
    unzip_file(destination_file_path, args.destination)

if __name__ == "__main__":
    handle_download()
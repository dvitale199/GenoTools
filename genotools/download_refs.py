import argparse
import os
import sys
import hashlib
import requests
import zipfile
from tqdm import tqdm


def validate_checksum(file_path, checksum):
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    if hash_md5.hexdigest() != checksum:
        return False
    return True


# def download_data_from_gcs(public_url, destination_file_name):
#     """Downloads a file from a public GCS bucket."""
#     if not os.path.exists(os.path.dirname(destination_file_name)):
#         os.makedirs(os.path.dirname(destination_file_name), exist_ok=True)

#     response = requests.get(public_url, stream=True)
#     if response.status_code == 200:
#         with open(destination_file_name, 'wb') as f:
#             for chunk in response.iter_content(chunk_size=128):
#                 f.write(chunk)
#         print(f"File downloaded to {destination_file_name}")

#     else:
#         response.raise_for_status()
def download_data_from_gcs(url, destination_file_path):
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
    """Unzips a file to a specified destination directory."""
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
        ref = args.ref if args.ref else "1kg_30x_hgdp_ashk_ref_panel"
        url = f"{url_base}/ref_panel/{ref}.zip"
        checksum = "fd3d79b9e1c0054b10881fa130eb9b21"
        print(f'Pulling reference panel {ref}')
        destination_file_path = os.path.join(args.destination, os.path.basename(url))
        download_data_from_gcs(url, destination_file_path)
        if not validate_checksum(destination_file_path, checksum):
            print("Error: Checksum validation failed for reference panel.")
            sys.exit(1)
        unzip_file(destination_file_path, args.destination)

    if download_model:
        model = args.model if args.model else "nba_v1"
        url = f"{url_base}/models/{model}.zip"
        # Add checksum for model if available
        print(f'Pulling model: {model}')
        destination_file_path = os.path.join(args.destination, os.path.basename(url))
        download_data_from_gcs(url, destination_file_path)
        # Add checksum validation for the model if you have a checksum
        unzip_file(destination_file_path, args.destination)

# # Main function to handle the download process
# def handle_download():
#     parser = argparse.ArgumentParser(description="Download, validate, and unzip reference data")
#     default_destination = os.path.expanduser("~/.genotools/ref")
#     parser.add_argument('--destination', type=str, default=default_destination, help="Local destination directory for the download (default: ~/.genotools/refs)")
#     parser.add_argument('--model', type=str, default="nba_v1", help="Version of the model to use (default: nba_v1)")
#     parser.add_argument('--ref', type=str, default="1kg_30x_hgdp_ashk_ref_panel", help="Version of the reference panel to use (default: 1kg_30x_hgdp_ashk_ref_panel)")

#     args = parser.parse_args()

#     url_base = "https://storage.googleapis.com/genotools_refs"

#     if args.ref:
#         url = f"{url_base}/ref_panel/{args.ref}.zip"
#         checksum = "fd3d79b9e1c0054b10881fa130eb9b21"
#         print(f'Pulling reference panel {args.ref}')
#         print(f"Downloading {url} to {args.destination}")
#         destination_file_path = os.path.join(args.destination, os.path.basename(args.url))
#         download_data_from_gcs(args.url, destination_file_path)
#         if not validate_checksum(destination_file_path, checksum):
#             print("Error: Checksum validation failed.")
#             sys.exit(1)
#         unzip_file(destination_file_path, args.destination)
#     if args.model:
#         url = f"{url_base}/models/{args.model}.zip"
#         checksum = None
#         print(f'Pulling model: {args.model}')
#         print(f"Downloading {url} to {args.destination}")
#         destination_file_path = os.path.join(args.destination, os.path.basename(args.url))
#         download_data_from_gcs(args.url, destination_file_path)
#             # Validate the checksum if provided
#         if not validate_checksum(destination_file_path, checksum):
#             print("Error: Checksum validation failed.")
#             sys.exit(1)
#         unzip_file(destination_file_path, args.destination)


if __name__ == "__main__":
    handle_download()
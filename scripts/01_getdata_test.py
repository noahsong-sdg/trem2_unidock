#!/usr/bin/env python3
"""
Test script for downloading PDBQT.gz files from ZINC database.
This version downloads only a few files for testing the workflow.
"""
import requests
import os
import gzip
import shutil
from pathlib import Path

# Get the absolute path of the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def download_zinc_subset(url, output_dir, filename=None):
    """
    Downloads a subset of the ZINC database.

    Args:
        url (str): URL to the ZINC subset (e.g., a direct download link to a PDBQT.gz file).
        output_dir (str): Directory to save the downloaded file.
        filename (str): Name of the file to save. If None, extracts from URL.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Extract filename from URL if not provided
    if filename is None:
        filename = url.split('/')[-1]

    filepath = os.path.join(output_dir, filename)

    print(f"Downloading ZINC file from {url} to {filepath}...")
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        # Get file size if available
        total_size = int(response.headers.get('content-length', 0))
        downloaded_size = 0
        
        with open(filepath, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                downloaded_size += len(chunk)
                if total_size > 0:
                    percent = (downloaded_size / total_size) * 100
                    print(f"\rDownloading {filename}: {percent:.1f}%", end='', flush=True)
        
        print(f"\nDownload complete: {filename}")
        file_size = os.path.getsize(filepath)
        
    except requests.exceptions.RequestException as e:
        print(f"Error downloading {filename} from {url}: {e}")
        return None
    return filepath

def extract_pdbqt_files(raw_dir, output_dir):
    """
    Extract .pdbqt.gz files to individual .pdbqt files ready for docking.
    
    Args:
        raw_dir (str): Directory containing downloaded .pdbqt.gz files
        output_dir (str): Directory to save extracted .pdbqt files
    
    Returns:
        tuple: (successful_extractions, failed_extractions)
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Find all .pdbqt.gz files
    gz_files = list(Path(raw_dir).glob("*.pdbqt.gz"))
    
    if not gz_files:
        print(f"No .pdbqt.gz files found in {raw_dir}")
        return 0, 0
    
    print(f"Extracting {len(gz_files)} PDBQT files...")
    
    successful = 0
    failed = 0
    
    for gz_file in gz_files:
        try:
            # Create output filename (remove .gz extension)
            output_file = Path(output_dir) / gz_file.stem
            
            # Extract the gzipped file
            with gzip.open(gz_file, 'rb') as f_in:
                with open(output_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            print(f"Extracted: {gz_file.name} → {output_file.name}")
            successful += 1
            
        except Exception as e:
            print(f"Error extracting {gz_file.name}: {e}")
            failed += 1
    
    print(f"Extraction complete! Successful: {successful}, Failed: {failed}")
    return successful, failed

def download_all_from_uri_file(uri_file_path, base_output_dir):
    """
    Reads a list of URLs from a .uri file and downloads data from each.

    Args:
        uri_file_path (str): Path to the .uri file containing URLs.
        base_output_dir (str): The base directory to save downloaded files.
    """
    if not os.path.exists(uri_file_path):
        print(f"Error: URI file not found at {uri_file_path}")
        return

    # Count URLs first for progress tracking
    urls = []
    with open(uri_file_path, 'r') as f:
        for line in f:
            url = line.strip()
            if url and not url.startswith("#"):
                urls.append(url)
    
    downloaded_files_count = 0
    failed_downloads_count = 0
    
    for i, url in enumerate(urls, 1):
        try:
            # Extract filename from URL (e.g., http://files.docking.org/3D/AC/AAML/ACAAML.xaa.pdbqt.gz)
            filename = url.split('/')[-1]
            if not filename: # Handle cases like trailing slash
                filename = f"downloaded_ligand_{i}.pdbqt.gz"
            
            print(f"Processing file {i}/{len(urls)}: {filename}")
            downloaded_file = download_zinc_subset(url, base_output_dir, filename=filename)
            if downloaded_file:
                downloaded_files_count += 1
            else:
                failed_downloads_count += 1
                            
        except Exception as e:
            print(f"Failed to process URL {url}: {e}")
            failed_downloads_count += 1
        
    print(f"Successfully downloaded files: {downloaded_files_count}")
    print(f"Failed downloads: {failed_downloads_count}")
    
    return downloaded_files_count, failed_downloads_count

if __name__ == "__main__":
    
    try:
        RAW_LIGANDS_DIR = os.path.join(SCRIPT_DIR, "../data/ligands_raw_test")
        URI_FILE = os.path.join(SCRIPT_DIR, "../data/test_pdbqt.uri") # Using test PDBQT.gz file URLs
        
        # Check if URI file exists
        if not os.path.exists(URI_FILE):
            print(f"Error: URI file not found at {URI_FILE}")
            print("Please ensure the URI file exists with URLs to download.")
            exit(1)
        
        # Step 1: Download PDBQT.gz files
        print("=== Step 1: Downloading PDBQT.gz files ===")
        successful_downloads, failed_downloads = download_all_from_uri_file(URI_FILE, RAW_LIGANDS_DIR)
        
        # Step 2: Extract PDBQT files to the standard location
        print("\n=== Step 2: Extracting PDBQT files ===")
        pdbqt_dir = os.path.join(SCRIPT_DIR, "../data/ligands_pdbqt_from_zinc")
        successful_extractions, failed_extractions = extract_pdbqt_files(RAW_LIGANDS_DIR, pdbqt_dir)
        
        print(f"\n=== DOWNLOAD & EXTRACTION SUMMARY ===")
        print(f"✓ Successful downloads: {successful_downloads}")
        print(f"✗ Failed downloads: {failed_downloads}")
        print(f"✓ Successful extractions: {successful_extractions}")
        print(f"✗ Failed extractions: {failed_extractions}")
        print(f"📁 Raw .gz files: {RAW_LIGANDS_DIR}")
        print(f"📁 PDBQT files ready for docking: {pdbqt_dir}")
        
        if successful_extractions > 0:
            print(f"\n🎉 SUCCESS: {successful_extractions} PDBQT files are ready for docking!")
            print("You can now skip the entire 02_prep.py script and go directly to 03_dock.py")
            print("Just point your docking script to the PDBQT directory above.")
        else:
            print("\n⚠️  WARNING: No PDBQT files were successfully extracted. Check the download and extraction logs.")
                        
    except Exception as e:
        print(f"Error during data download and extraction: {e}")
        exit(1)

#!/usr/bin/env python3
"""
Data download script with comprehensive timing functionality and parallel downloads.
Downloads ligand data and tracks performance metrics.
"""
import requests
import os
import gzip
import shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# Get the absolute path of the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Thread-safe progress tracking
download_lock = threading.Lock()
progress_counter = {'completed': 0, 'failed': 0}

def download_zinc_subset(url, output_dir, filename=None):
    """
    Downloads a subset of the ZINC database.

    Args:
        url (str): URL to the ZINC subset (e.g., a direct download link to a PDBQT.gz file).
        output_dir (str): Directory to save the downloaded file.
        filename (str): Name of the file to save. If None, extracts from URL.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Extract filename from URL if not provided
    if filename is None:
        filename = url.split('/')[-1]

    filepath = os.path.join(output_dir, filename)

    try:
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        # Get file size if available
        total_size = int(response.headers.get('content-length', 0))
        downloaded_size = 0
        
        with open(filepath, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                downloaded_size += len(chunk)
        
        file_size = os.path.getsize(filepath)
        
        # Thread-safe progress update
        with download_lock:
            progress_counter['completed'] += 1
            print(f"✓ Downloaded ({progress_counter['completed']}) {filename} ({file_size:,} bytes)")
        
        return filepath
        
    except requests.exceptions.RequestException as e:
        with download_lock:
            progress_counter['failed'] += 1
            print(f"✗ Failed ({progress_counter['failed']}) {filename}: {e}")
        return None
    except Exception as e:
        with download_lock:
            progress_counter['failed'] += 1
            print(f"✗ Error ({progress_counter['failed']}) {filename}: {e}")
        return None

def download_single_file(args):
    """Helper function for parallel downloads"""
    url, output_dir, filename = args
    return download_zinc_subset(url, output_dir, filename)

def download_all_from_uri_file(uri_file_path, base_output_dir, max_workers=4):
    """
    Reads a list of URLs from a .uri file and downloads data from each using parallel processing.

    Args:
        uri_file_path (str): Path to the .uri file containing URLs.
        base_output_dir (str): The base directory to save downloaded files.
        max_workers (int): Number of parallel download threads.
    """
    if not os.path.exists(uri_file_path):
        print(f"Error: URI file not found at {uri_file_path}")
        return 0, 0

    # Read all URLs
    urls = []
    with open(uri_file_path, 'r') as f:
        for line in f:
            url = line.strip()
            if url and not url.startswith("#"):
                urls.append(url)
    
    if not urls:
        print("No valid URLs found in URI file")
        return 0, 0
    
    print(f"Starting parallel download of {len(urls)} files using {max_workers} workers...")
    
    # Reset progress counters
    with download_lock:
        progress_counter['completed'] = 0
        progress_counter['failed'] = 0
    
    # Prepare download arguments
    download_args = []
    for i, url in enumerate(urls, 1):
        filename = url.split('/')[-1]
        if not filename:
            filename = f"downloaded_ligand_{i}.pdbqt.gz"
        download_args.append((url, base_output_dir, filename))
    
    # Execute parallel downloads
    downloaded_files_count = 0
    failed_downloads_count = 0
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all download tasks
        future_to_url = {executor.submit(download_single_file, args): args[0] 
                        for args in download_args}
        
        # Process completed downloads
        for future in as_completed(future_to_url):
            url = future_to_url[future]
            try:
                result = future.result()
                if result:
                    downloaded_files_count += 1
                else:
                    failed_downloads_count += 1
            except Exception as e:
                failed_downloads_count += 1
                print(f"✗ Exception downloading {url}: {e}")
    
    print(f"\n=== PARALLEL DOWNLOAD COMPLETE ===")
    print(f"✓ Successful downloads: {downloaded_files_count}")
    print(f"✗ Failed downloads: {failed_downloads_count}")
    
    return downloaded_files_count, failed_downloads_count

def extract_pdbqt_files(raw_dir, output_dir, max_workers=4):
    """
    Extract .pdbqt.gz files to individual .pdbqt files ready for docking with parallel processing.
    
    Args:
        raw_dir (str): Directory containing downloaded .pdbqt.gz files
        output_dir (str): Directory to save extracted .pdbqt files
        max_workers (int): Number of parallel extraction threads
    
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
    
    print(f"Extracting {len(gz_files)} PDBQT files using {max_workers} workers...")
    
    def extract_single_file(gz_file):
        """Helper function for parallel extraction"""
        try:
            # Create output filename (remove .gz extension)
            output_file = Path(output_dir) / gz_file.stem
            
            # Extract the gzipped file
            with gzip.open(gz_file, 'rb') as f_in:
                with open(output_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            return gz_file.name, True, None
            
        except Exception as e:
            return gz_file.name, False, str(e)
    
    successful = 0
    failed = 0
    
    # Execute parallel extractions
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all extraction tasks
        future_to_file = {executor.submit(extract_single_file, gz_file): gz_file 
                         for gz_file in gz_files}
        
        # Process completed extractions
        for future in as_completed(future_to_file):
            gz_file = future_to_file[future]
            try:
                filename, success, error = future.result()
                if success:
                    successful += 1
                    print(f"✓ Extracted ({successful}): {filename}")
                else:
                    failed += 1
                    print(f"✗ Failed ({failed}): {filename} - {error}")
            except Exception as e:
                failed += 1
                print(f"✗ Exception extracting {gz_file.name}: {e}")
    
    print(f"Extraction complete! Successful: {successful}, Failed: {failed}")
    return successful, failed

def get_tranche_name_from_filename(filename):
    """Extract tranche name from ZINC filename for organization."""
    # ZINC files are typically named like: BFEDMM.xaa.pdbqt.gz
    # Extract the tranche part (e.g., BFEDMM.xaa)
    if '.pdbqt' in filename:
        base_name = filename.replace('.pdbqt.gz', '').replace('.pdbqt', '')
        if '.' in base_name:
            return base_name
    return "unknown_tranche"

def _save_molecule(molecule_lines, molecule_name, output_dir, tranche_name, molecule_index):
    """Save a single molecule to its own PDBQT file in appropriate tranche directory."""
    if not molecule_lines:
        return
    
    # Create tranche-specific directory
    tranche_dir = os.path.join(output_dir, tranche_name)
    os.makedirs(tranche_dir, exist_ok=True)
    
    # Clean molecule lines - remove MODEL and ENDMDL lines that cause UniDock errors
    clean_lines = []
    for line in molecule_lines:
        stripped_line = line.strip()
        if not (stripped_line.startswith('MODEL') or stripped_line.startswith('ENDMDL')):
            clean_lines.append(line)
    
    # Generate filename
    if molecule_name:
        # Use molecule name if available
        filename = f"{molecule_name}.pdbqt"
        # Clean filename - remove invalid characters
        filename = "".join(c for c in filename if c.isalnum() or c in "._-")
    else:
        # Use ZINC ID from molecule content if available
        zinc_id = None
        for line in clean_lines:
            if line.startswith('REMARK  Name ='):
                zinc_id = line.split('=', 1)[1].strip()
                break
        
        if zinc_id:
            filename = f"{zinc_id}.pdbqt"
        else:
            filename = f"molecule_{molecule_index:06d}.pdbqt"
    
    output_path = os.path.join(tranche_dir, filename)
    
    with open(output_path, 'w') as f:
        for line in clean_lines:
            f.write(line + '\n')

def split_pdbqt_files(input_dir, output_dir, max_workers=4):
    """
    Split multi-molecule PDBQT files into individual single-molecule files organized by tranche.
    This is required for UniDock and mcdock which expect one molecule per file.
    
    Args:
        input_dir (str): Directory containing multi-molecule PDBQT files
        output_dir (str): Directory to save individual PDBQT files (organized by tranche)
        max_workers (int): Number of parallel splitting threads
    
    Returns:
        tuple: (total_molecules_split, failed_files, tranche_count)
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Find all .pdbqt files
    pdbqt_files = list(Path(input_dir).glob("*.pdbqt"))
    
    if not pdbqt_files:
        print(f"No .pdbqt files found in {input_dir}")
        return 0, 0, 0
    
    print(f"Splitting {len(pdbqt_files)} PDBQT files into tranche-organized individual molecules using {max_workers} workers...")
    
    def split_single_pdbqt(pdbqt_file):
        """Helper function for parallel PDBQT splitting"""
        try:
            molecule_count = 0
            current_molecule = []
            current_name = None
            tranche_name = get_tranche_name_from_filename(pdbqt_file.name)
            
            with open(pdbqt_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    
                    if line.startswith('MODEL'):
                        # Start of new molecule
                        if current_molecule:
                            # Save previous molecule
                            _save_molecule(current_molecule, current_name, output_dir, 
                                         tranche_name, molecule_count)
                            molecule_count += 1
                        
                        # Reset for new molecule
                        current_molecule = [line]
                        current_name = None
                        
                    elif line.startswith('REMARK  Name ='):
                        # Extract molecule name
                        current_name = line.split('=', 1)[1].strip()
                        current_molecule.append(line)
                        
                    elif line.startswith('ENDMDL'):
                        # End of molecule
                        current_molecule.append(line)
                        _save_molecule(current_molecule, current_name, output_dir, 
                                     tranche_name, molecule_count)
                        molecule_count += 1
                        current_molecule = []
                        current_name = None
                        
                    else:
                        # Regular line - add to current molecule
                        if current_molecule:  # Only if we're inside a MODEL
                            current_molecule.append(line)
            
            # Handle case where file doesn't end with ENDMDL
            if current_molecule:
                _save_molecule(current_molecule, current_name, output_dir, 
                             tranche_name, molecule_count)
                molecule_count += 1
            
            return pdbqt_file.name, molecule_count, tranche_name, None
            
        except Exception as e:
            return pdbqt_file.name, 0, "unknown", str(e)
    
    total_molecules = 0
    failed_files = 0
    tranches_created = set()
    
    # Execute parallel splitting
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all splitting tasks
        future_to_file = {executor.submit(split_single_pdbqt, pdbqt_file): pdbqt_file 
                         for pdbqt_file in pdbqt_files}
        
        # Process completed splits
        for future in as_completed(future_to_file):
            pdbqt_file = future_to_file[future]
            try:
                filename, molecule_count, tranche_name, error = future.result()
                if error:
                    failed_files += 1
                    print(f"✗ Failed to split {filename}: {error}")
                else:
                    total_molecules += molecule_count
                    tranches_created.add(tranche_name)
                    print(f"✓ Split {filename}: {molecule_count} molecules → tranche {tranche_name}")
            except Exception as e:
                failed_files += 1
                print(f"✗ Exception splitting {pdbqt_file.name}: {e}")
    
    print(f"Splitting complete! Total molecules: {total_molecules}, Failed files: {failed_files}, Tranches: {len(tranches_created)}")
    return total_molecules, failed_files, len(tranches_created)

if __name__ == "__main__":
    
    try:
        RAW_LIGANDS_DIR = os.path.join(SCRIPT_DIR, "../data/450_3/ligands_raw")
        URI_FILE = os.path.join(SCRIPT_DIR, "../data/450_3.uri") # Using the PDBQT.gz file URLs
        
        # Configuration for parallel processing
        DOWNLOAD_WORKERS = 8  # Number of parallel download threads
        EXTRACTION_WORKERS = 4  # Number of parallel extraction threads
        
        print(f"=== PARALLEL DATA DOWNLOAD SCRIPT ===")
        print(f"Download workers: {DOWNLOAD_WORKERS}")
        print(f"Extraction workers: {EXTRACTION_WORKERS}")
        
        # Check if URI file exists
        if not os.path.exists(URI_FILE):
            print(f"Error: URI file not found at {URI_FILE}")
            print("Please ensure the URI file exists with URLs to download.")
            exit(1)
        
        # Download all files with parallel processing
        successful_downloads, failed_downloads = download_all_from_uri_file(
            URI_FILE, RAW_LIGANDS_DIR, max_workers=DOWNLOAD_WORKERS)
                
        print(f"\n=== DOWNLOAD SUMMARY ===")
        print(f"✓ Successful downloads: {successful_downloads}")
        print(f"✗ Failed downloads: {failed_downloads}")
        print(f"📁 Files saved to: {RAW_LIGANDS_DIR}")
        
        # Extract PDBQT files with parallel processing
        pdbqt_dir = os.path.join(SCRIPT_DIR, "../data/450_3/ligands_pdbqt")
        successful_extractions, failed_extractions = extract_pdbqt_files(
            RAW_LIGANDS_DIR, pdbqt_dir, max_workers=EXTRACTION_WORKERS)
        print(f"\n=== EXTRACTION SUMMARY ===")
        print(f"✓ Successful extractions: {successful_extractions}")
        print(f"✗ Failed extractions: {failed_extractions}")
        print(f"📁 PDBQT files ready for splitting: {pdbqt_dir}")
        
        # Check if data has already been processed
        split_dir = os.path.join(SCRIPT_DIR, "../data/450_3/ligands_pdbqt_split")
        
        if os.path.exists(split_dir) and os.listdir(split_dir):
            # Count existing processed molecules
            existing_tranches = 0
            existing_molecules = 0
            for item in os.listdir(split_dir):
                item_path = os.path.join(split_dir, item)
                if os.path.isdir(item_path):
                    pdbqt_files = list(Path(item_path).glob("*.pdbqt"))
                    if pdbqt_files:
                        existing_tranches += 1
                        existing_molecules += len(pdbqt_files)
            
            print(f"\n=== EXISTING PROCESSED DATA DETECTED ===")
            print(f"✓ Found {existing_tranches} tranches with {existing_molecules:,} individual molecules")
            print(f"📁 Location: {split_dir}")
            print(f"\n🚀 Data is already processed and ready for docking!")
            print(f"\n📋 Next Steps:")
            print(f"   • For tranche-aware UniDock: LIGAND_DIR is already set to {split_dir}")
            print(f"   • Run: python scripts/mcdock.py")
            
        elif successful_extractions > 0:
            # Split multi-molecule PDBQT files into individual molecules organized by tranche
            print(f"\n--- Splitting extracted PDBQT files into individual molecules ---")
            total_molecules, failed_splits, tranche_count = split_pdbqt_files(
                pdbqt_dir, split_dir, max_workers=EXTRACTION_WORKERS)
            print(f"\n=== SPLITTING SUMMARY ===")
            print(f"✓ Total individual molecules created: {total_molecules}")
            print(f"✗ Failed file splits: {failed_splits}")
            print(f"📁 Tranches created: {tranche_count}")
            print(f"📁 Individual PDBQT files organized by tranche: {split_dir}")
            
            if total_molecules > 0:
                print(f"\n🚀 SUCCESS: {total_molecules} individual PDBQT molecules are ready for docking!")
                print(f"   Organized across {tranche_count} tranches for efficient HTS processing")
                print(f"\n📋 Next Steps:")
                print(f"   • For tranche-aware UniDock: Update LIGAND_DIR to {split_dir}")
                print(f"   • For mcdock: Update LIGAND_DIR to {split_dir}")
                print(f"   • Run: python scripts/mcdock.py")
            else:
                print("\n⚠️  WARNING: No molecules were successfully split.")
        else:
            print("\n⚠️  WARNING: No PDBQT files were successfully extracted. Check the download and extraction logs.")
                        
    except Exception as e:
        print(f"Error during data download: {e}")
        exit(1)

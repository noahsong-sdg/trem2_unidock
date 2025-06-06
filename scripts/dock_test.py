#!/usr/bin/env python3
"""
Test version of UniDock docking script for subset testing.
This script tests molecular docking on a small subset of prepared ligands.
"""
import os
import subprocess
import shutil
from pathlib import Path

# Import timing utilities
from timing_utils import TimingTracker

# Get the absolute path of the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# --- Test Configuration ---
RECEPTOR_FILE = os.path.join(SCRIPT_DIR, "../data/receptor/cluster1_receptor.pdbqt")
LIGAND_DIR = os.path.join(SCRIPT_DIR, "../data/ligands_pdbqt/")
TEST_LIGAND_DIR = os.path.join(SCRIPT_DIR, "../data/ligands_pdbqt_test/")
DOCKING_OUTPUT_DIR = os.path.join(SCRIPT_DIR, "../results/docking_outputs_test/")
CONFIG_FILE = os.path.join(SCRIPT_DIR, "../configs/docking_config_test.txt")

# Test with only 3 ligands for quick validation
MAX_TEST_LIGANDS = 3

def create_test_ligand_subset():
    """Create a test subset of ligands for testing."""
    # Check if test ligands already exist
    if os.path.exists(TEST_LIGAND_DIR):
        existing_ligands = [f for f in os.listdir(TEST_LIGAND_DIR) if f.endswith('.pdbqt')]
        if existing_ligands:
            print(f"Using existing test ligands: {existing_ligands}")
            return len(existing_ligands)
    
    os.makedirs(TEST_LIGAND_DIR, exist_ok=True)
    
    # Check if main ligand directory exists
    if not os.path.exists(LIGAND_DIR):
        print(f"Warning: Main ligand directory not found: {LIGAND_DIR}")
        print("Will proceed with any existing test ligands")
        existing_ligands = [f for f in os.listdir(TEST_LIGAND_DIR) if f.endswith('.pdbqt')]
        return len(existing_ligands)
    
    # Get list of available ligands
    ligand_files = [f for f in os.listdir(LIGAND_DIR) if f.endswith('.pdbqt')]
    
    # Copy first few ligands to test directory
    test_ligands = ligand_files[:MAX_TEST_LIGANDS]
    
    for ligand_file in test_ligands:
        src = os.path.join(LIGAND_DIR, ligand_file)
        dst = os.path.join(TEST_LIGAND_DIR, ligand_file)
        shutil.copy2(src, dst)
        print(f"Copied {ligand_file} to test directory")
    
    print(f"Created test subset with {len(test_ligands)} ligands")
    return len(test_ligands)

def create_unidock_config(config_filepath, receptor_filepath, ligand_dir, 
                          center_x, center_y, center_z, 
                          size_x, size_y, size_z, 
                          num_modes=9, search_mode="balance", scoring_function="vinardo"):
    """
    Creates a configuration file for Uni-Dock.
    """
    # Ensure the directory for the config file exists
    os.makedirs(os.path.dirname(config_filepath), exist_ok=True)

    config_content = f"""
receptor = {os.path.abspath(receptor_filepath)}

# Search space configuration (TREM2 binding site coordinates)
center_x = {center_x}
center_y = {center_y}
center_z = {center_z}

size_x = {size_x}
size_y = {size_y}
size_z = {size_z}

# Docking parameters
num_modes = {num_modes}
search_mode = {search_mode}
scoring = {scoring_function}
"""
    with open(config_filepath, 'w') as f:
        f.write(config_content)
    print(f"Uni-Dock configuration file created at: {config_filepath}")

def run_unidock_test(unidock_executable, receptor_file, ligand_input, output_dir, 
                     center_x, center_y, center_z, size_x, size_y, size_z, 
                     scoring_function="vinardo", num_modes=9, timer=None):
    """
    Runs Uni-Dock for a test set of ligands against a receptor.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Get ligand files
    ligand_files = []
    if os.path.isdir(ligand_input):
        for file in os.listdir(ligand_input):
            if file.endswith('.pdbqt'):
                ligand_files.append(os.path.join(ligand_input, file))
    
    if not ligand_files:
        print(f"No valid ligand files found in {ligand_input}")
        return 0, 0
    
    if timer:
        timer.start_step("Molecular docking with UniDock", len(ligand_files))
    
    print(f"Found {len(ligand_files)} ligand file(s) to dock")
    
    successful_dockings = 0
    failed_dockings = 0
    
    # For small test set, do individual docking for better error handling
    for i, ligand_file in enumerate(ligand_files):
        ligand_name = os.path.splitext(os.path.basename(ligand_file))[0]
        output_name = f"{ligand_name}_docked.pdbqt"
        output_path = os.path.join(output_dir, output_name)
        
        print(f"\nDocking ligand {i+1}/{len(ligand_files)}: {ligand_name}")
        
        command = [
            unidock_executable,
            "--receptor", os.path.abspath(receptor_file),
            "--ligand", os.path.abspath(ligand_file),
            "--center_x", str(center_x),
            "--center_y", str(center_y),
            "--center_z", str(center_z),
            "--size_x", str(size_x),
            "--size_y", str(size_y),
            "--size_z", str(size_z),
            "--scoring", scoring_function,
            "--num_modes", str(num_modes),
            "--out", output_path
        ]
        
        print(f"Command: {' '.join(command)}")
        
        try:
            result = subprocess.run(command, check=True, text=True, capture_output=True)
            
            if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
                successful_dockings += 1
                print(f"✓ Docking successful: {output_path}")
                print(f"  Output size: {os.path.getsize(output_path)} bytes")
            else:
                failed_dockings += 1
                print(f"✗ Docking failed - no valid output generated")
            
            if timer:
                timer.update_progress(1)
                
        except subprocess.CalledProcessError as e:
            print(f"✗ Error during UniDock execution:")
            print(f"  Command: {' '.join(e.cmd)}")
            print(f"  Return code: {e.returncode}")
            print(f"  Error output: {e.stderr}")
            failed_dockings += 1
        except FileNotFoundError:
            print(f"✗ Error: {unidock_executable} not found. Please ensure it's in your PATH.")
            failed_dockings += 1
        except Exception as e:
            print(f"✗ Unexpected error during UniDock execution: {e}")
            failed_dockings += 1
    
    if timer:
        timer.end_step()
    
    return successful_dockings, failed_dockings

if __name__ == "__main__":
    # Initialize timing tracker
    timer = TimingTracker("03_dock_test")
    
    try:
        # --- TREM2 Binding Site Coordinates ---
        # These coordinates are for the TREM2 binding site
        CENTER_X, CENTER_Y, CENTER_Z = 42.328, 28.604, 21.648 
        SIZE_X, SIZE_Y, SIZE_Z = 30.0, 30.0, 30.0  # In Angstroms

        # UniDock executable
        UNIDOCK_EXECUTABLE = "unidock"

        print("=== UniDock Test Workflow ===")
        
        # 1. Create test ligand subset
        print("\n--- Step 1: Creating Test Ligand Subset ---")
        timer.start_step("Create test ligand subset")
        ligand_count = create_test_ligand_subset()
        timer.end_step()

        # 2. Create Uni-Dock configuration file
        print("\n--- Step 2: Creating Uni-Dock Configuration ---")
        timer.start_step("Create UniDock configuration")
        create_unidock_config(CONFIG_FILE, RECEPTOR_FILE, TEST_LIGAND_DIR, 
                              CENTER_X, CENTER_Y, CENTER_Z, 
                              SIZE_X, SIZE_Y, SIZE_Z,
                              scoring_function="vinardo")
        timer.end_step()

        # 3. Validate input files
        print("\n--- Step 3: Validating Input Files ---")
        timer.start_step("Validate input files")
        
        if not os.path.exists(RECEPTOR_FILE):
            print(f"\nError: Receptor file not found at {RECEPTOR_FILE}")
            timer.finish()
            exit(1)
        else:
            print(f"✓ Receptor file found: {RECEPTOR_FILE}")

        if not os.path.exists(TEST_LIGAND_DIR) or not os.listdir(TEST_LIGAND_DIR):
            print(f"\nError: No test ligand files found in {TEST_LIGAND_DIR}")
            timer.finish()
            exit(1)
        else:
            test_ligand_count = len([f for f in os.listdir(TEST_LIGAND_DIR) if f.endswith('.pdbqt')])
            print(f"✓ Found {test_ligand_count} test ligands in: {TEST_LIGAND_DIR}")
        
        timer.end_step()
        
        # 4. Run Uni-Dock
        print("\n--- Step 4: Running UniDock Molecular Docking ---")
        print("Testing docking workflow on small subset of ligands...")
        
        successful_dockings, failed_dockings = run_unidock_test(
            UNIDOCK_EXECUTABLE, RECEPTOR_FILE, TEST_LIGAND_DIR, DOCKING_OUTPUT_DIR, 
            CENTER_X, CENTER_Y, CENTER_Z, SIZE_X, SIZE_Y, SIZE_Z, timer=timer
        )
        
        # Set final ligand count for performance metrics
        timer.set_final_ligand_count(successful_dockings)
        
        # Generate final timing report
        report = timer.finish()
        
        print(f"\n=== TEST DOCKING WORKFLOW SUMMARY ===")
        print(f"✓ Successful dockings: {successful_dockings}")
        print(f"✗ Failed dockings: {failed_dockings}")
        print(f" Test docking outputs saved to: {DOCKING_OUTPUT_DIR}")
        
        # Performance summary
        if "performance_metrics" in report:
            metrics = report["performance_metrics"]
            print(f"\n Test Docking Performance Metrics:")
            print(f"   Docking rate: {metrics['ligands_per_minute']:.1f} ligands/minute")
            print(f"   Average time per ligand: {metrics['average_seconds_per_ligand']:.3f} seconds")
            print(f"\n Scale Estimates:")
            print(f"   1M ligands would take: {metrics['estimated_time_for_1M_ligands']}")
            print(f"   10M ligands would take: {metrics['estimated_time_for_10M_ligands']}")
            
        # Show sample output files
        if successful_dockings > 0:
            print(f"\n Sample docking output files:")
            output_files = [f for f in os.listdir(DOCKING_OUTPUT_DIR) if f.endswith('.pdbqt')]
            for output_file in output_files[:3]:
                file_path = os.path.join(DOCKING_OUTPUT_DIR, output_file)
                file_size = os.path.getsize(file_path)
                print(f"   {output_file} ({file_size} bytes)")
    
    except Exception as e:
        print(f"Error during test docking workflow: {e}")
        timer.finish()
        exit(1)

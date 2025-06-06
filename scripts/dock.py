#!/usr/bin/env python3
"""
UniDock docking script with comprehensive timing functionality.
This script performs molecular docking using UniDock and tracks performance metrics.
"""
import os
import subprocess
from pathlib import Path

# Import timing utilities
from timing_utils import TimingTracker

# Get the absolute path of the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# --- Configuration ---
RECEPTOR_FILE = os.path.join(SCRIPT_DIR, "../data/receptor/cluster1_receptor.pdbqt")  # Prepared receptor file
LIGAND_DIR = os.path.join(SCRIPT_DIR, "../data/ligands_pdbqt/")       # Directory containing individual ligand PDBQT files
DOCKING_OUTPUT_DIR = os.path.join(SCRIPT_DIR, "../results/docking_outputs/")
CONFIG_FILE = os.path.join(SCRIPT_DIR, "../configs/docking_config.txt") # Uni-Dock configuration file

# --- Create Uni-Dock Configuration File ---
def create_unidock_config(config_filepath, receptor_filepath, ligand_dir, 
                          center_x, center_y, center_z, 
                          size_x, size_y, size_z, 
                          num_modes=5, search_mode="balance", scoring_function="vinardo"):
    """
    Creates a configuration file for Uni-Dock.

    Args:
        config_filepath (str): Path to save the configuration file.
        receptor_filepath (str): Path to the prepared receptor PDBQT file.
        ligand_filepath_template (str): Path to the ligand directory for batch processing.
        center_x, center_y, center_z (float): Coordinates of the search space center.
        size_x, size_y, size_z (float): Dimensions of the search space (Angstroms).
        num_modes (int): Number of binding modes to generate.
        exhaustiveness (int): Exhaustiveness of the search.
        scoring_function (str): Scoring function to use (e.g., "vina", "vinardo", "ad4").
    """
    # Ensure the directory for the config file exists
    os.makedirs(os.path.dirname(config_filepath), exist_ok=True)

    config_content = f"""
receptor = {os.path.abspath(receptor_filepath)}

# For docking a single ligand file, use 'ligand =' path/to/ligand.pdbqt
# For batch docking of all ligands in a directory, UniDock typically handles this via command line args
# or by specifying a directory. We will use command line for batch.

# Search space configuration (replace with actual TREM2 binding site coordinates)
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

# Output format (optional, UniDock defaults usually work)
# out_flex = # Path for flexible part of receptor (if any)
# log = docking.log # Log file name (UniDock usually names this based on output)
"""
    with open(config_filepath, 'w') as f:
        f.write(config_content) # Added missing write
    print(f"Uni-Dock configuration file created at: {config_filepath}")

# --- Run Uni-Dock ---
def run_unidock(unidock_executable, receptor_file, ligand_input, output_dir, 
                center_x, center_y, center_z, size_x, size_y, size_z, 
                scoring_function="vinardo", num_modes=5, timer=None):
    """
    Runs Uni-Dock for a set of ligands against a receptor.

    Args:
        unidock_executable (str): Path to the Uni-Dock executable.
        receptor_file (str): Path to the prepared receptor PDBQT file.
        ligand_input (str): Path to the directory containing prepared ligand files (e.g., SDF, PDBQT) 
                            OR path to a single ligand file.
        output_dir (str): Directory to save docking results.
        center_x, center_y, center_z (float): Coordinates of the search space center.
        size_x, size_y, size_z (float): Dimensions of the search space (Angstroms).
        scoring_function (str): Scoring function to use.
        num_modes (int): Number of binding modes to generate.
        timer (TimingTracker): Timing tracker instance.
    
    Returns:
        tuple: (successful_dockings, failed_dockings)
    """
    os.makedirs(output_dir, exist_ok=True)

    # For single ligand file, use --ligand
    # For batch processing, create a ligand index file and use --gpu_batch
    ligand_files = []
    
    if os.path.isfile(ligand_input):
        # Single ligand file
        ligand_files = [ligand_input]
    elif os.path.isdir(ligand_input):
        # Directory of ligand files
        for file in os.listdir(ligand_input):
            if file.endswith(('.pdbqt', '.sdf')):
                ligand_files.append(os.path.join(ligand_input, file))
    
    if not ligand_files:
        print(f"No valid ligand files found in {ligand_input}")
        return 0, 0
    
    if timer:
        timer.start_step("Molecular docking with UniDock", len(ligand_files))
    
    print(f"Found {len(ligand_files)} ligand file(s) to dock")
    
    successful_dockings = 0
    failed_dockings = 0
    
    # For batch processing, create ligand index file
    if len(ligand_files) > 1:
        ligand_index_file = os.path.join(output_dir, "ligand_index.txt")
        with open(ligand_index_file, 'w') as f:
            for ligand_file in ligand_files:
                f.write(f"{os.path.abspath(ligand_file)}\n")
        
        command = [
            unidock_executable,
            "--receptor", os.path.abspath(receptor_file),
            "--ligand_index", ligand_index_file,
            "--center_x", str(center_x),
            "--center_y", str(center_y), 
            "--center_z", str(center_z),
            "--size_x", str(size_x),
            "--size_y", str(size_y),
            "--size_z", str(size_z),
            "--scoring", scoring_function,
            "--num_modes", str(num_modes),
            "--dir", os.path.abspath(output_dir)
        ]
        
        print(f"Running batch docking command: {' '.join(command)}")
        try:
            result = subprocess.run(command, check=True, text=True, capture_output=True)
            
            # Count successful outputs
            for ligand_file in ligand_files:
                ligand_name = os.path.splitext(os.path.basename(ligand_file))[0]
                expected_output = os.path.join(output_dir, f"{ligand_name}_out.pdbqt")
                if os.path.exists(expected_output) and os.path.getsize(expected_output) > 0:
                    successful_dockings += 1
                else:
                    failed_dockings += 1
            
            if timer:
                timer.update_progress(len(ligand_files))
            
            print(f"Batch docking completed: {successful_dockings} successful, {failed_dockings} failed")
            
        except subprocess.CalledProcessError as e:
            print(f"Error during batch UniDock execution:")
            print(f"Command: {' '.join(e.cmd)}")
            print(f"Return code: {e.returncode}")
            print(f"Error output: {e.stderr}")
            failed_dockings = len(ligand_files)
            
    else:
        # Single ligand docking
        ligand_file = ligand_files[0]
        output_name = os.path.splitext(os.path.basename(ligand_file))[0] + "_docked.pdbqt"
        output_path = os.path.join(output_dir, output_name)
        
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

        print(f"Running single ligand docking: {' '.join(command)}")
        try:
            result = subprocess.run(command, check=True, text=True, capture_output=True)
            
            if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
                successful_dockings = 1
                print(f"Single ligand docking successful: {output_path}")
            else:
                failed_dockings = 1
                print(f"Docking failed - no valid output generated")
            
            if timer:
                timer.update_progress(1)
                
        except subprocess.CalledProcessError as e:
            print(f"Error during single ligand UniDock execution:")
            print(f"Command: {' '.join(e.cmd)}")
            print(f"Return code: {e.returncode}")
            print(f"Error output: {e.stderr}")
            failed_dockings = 1
        except FileNotFoundError:
            print(f"Error: {unidock_executable} not found. Please ensure it's in your PATH or provide the full path.")
            failed_dockings = 1
        except Exception as e:
            print(f"An unexpected error occurred during UniDock execution: {e}")
            failed_dockings = 1
    
    if timer:
        timer.end_step()
    
    return successful_dockings, failed_dockings

if __name__ == "__main__":
    # Initialize timing tracker
    timer = TimingTracker("03_dock")
    
    try:
        # --- !!! IMPORTANT: Define Receptor Binding Site !!! ---
        # These are placeholder values. Replace with actual coordinates for TREM2.
        # You can get these from literature, by analyzing the receptor structure in a
        # molecular viewer (like PyMOL, ChimeraX) with a known ligand, or using pocket detection tools.
        CENTER_X, CENTER_Y, CENTER_Z = 42.328, 28.604, 21.648 
        SIZE_X, SIZE_Y, SIZE_Z = 30.0, 30.0, 30.0 # In Angstroms

        # --- Path to Uni-Dock executable ---
        # If unidock is in your PATH, you can just use "unidock".
        # Otherwise, provide the full path, e.g., "/path/to/unidock_1.1.1_Linux_x86_64/unidock"
        UNIDOCK_EXECUTABLE = "unidock" # Or "/path/to/your/unidock/executable"

        print("=== UniDock HTS Workflow with Performance Tracking ===")

        # 1. Create Uni-Dock configuration file
        print("\n--- Step 1: Creating Uni-Dock Configuration ---")
        timer.start_step("Create UniDock configuration")
        create_unidock_config(CONFIG_FILE, RECEPTOR_FILE, LIGAND_DIR, 
                              CENTER_X, CENTER_Y, CENTER_Z, 
                              SIZE_X, SIZE_Y, SIZE_Z,
                              scoring_function="vinardo") # Using vinardo as an example
        timer.end_step()

        # 2. Check for receptor and ligands
        print("\n--- Step 2: Validating Input Files ---")
        timer.start_step("Validate input files")
        
        if not os.path.exists(RECEPTOR_FILE):
            print(f"\nError: Receptor file not found at {RECEPTOR_FILE}")
            print("Please prepare your TREM2 receptor, convert it to PDBQT, and place it there.")
            print("Example using Open Babel: obabel trem2.pdb -O ../data/receptor/cluster1_receptor.pdbqt -xr")
            timer.finish()
            exit(1)

        if not os.path.exists(LIGAND_DIR) or not os.listdir(LIGAND_DIR):
            print(f"\nError: No prepared ligand files found in {LIGAND_DIR}")
            print("Please run the `02_prep.py` script first or place prepared ligands there.")
            timer.finish()
            exit(1)
        else:
            ligand_count = len([f for f in os.listdir(LIGAND_DIR) if f.endswith(('.pdbqt', '.sdf'))])
            print(f"\nFound {ligand_count} prepared ligands in: {LIGAND_DIR}")
        
        timer.end_step()
        
        # 3. Run Uni-Dock
        print("\n--- Step 3: Running UniDock Molecular Docking ---")
        print("This step will process all ligands in the specified directory against the receptor.")
        print("Ensure UniDock executable is correctly specified and receptor/ligand paths are correct.")
        
        successful_dockings, failed_dockings = run_unidock(
            UNIDOCK_EXECUTABLE, RECEPTOR_FILE, LIGAND_DIR, DOCKING_OUTPUT_DIR, 
            CENTER_X, CENTER_Y, CENTER_Z, SIZE_X, SIZE_Y, SIZE_Z, 
            scoring_function="vinardo", num_modes=5, timer=timer
        )
        
        # Set final ligand count for performance metrics
        timer.set_final_ligand_count(successful_dockings)
        
        # Generate final timing report
        report = timer.finish()
        
        print(f"\n=== DOCKING WORKFLOW SUMMARY ===")
        print(f"✓ Successful dockings: {successful_dockings}")
        print(f"✗ Failed dockings: {failed_dockings}")
        print(f"📁 Docking outputs saved to: {DOCKING_OUTPUT_DIR}")
        
        # Performance summary
        if "performance_metrics" in report:
            metrics = report["performance_metrics"]
            print(f"\n📊 Docking Performance Metrics:")
            print(f"   Docking rate: {metrics['ligands_per_minute']:.1f} ligands/minute")
            print(f"   Average time per ligand: {metrics['average_seconds_per_ligand']:.3f} seconds")
            print(f"\n🚀 Scale Estimates:")
            print(f"   1M ligands would take: {metrics['estimated_time_for_1M_ligands']}")
            print(f"   10M ligands would take: {metrics['estimated_time_for_10M_ligands']}")
    
    except Exception as e:
        print(f"Error during docking workflow: {e}")
        timer.finish()
        exit(1)

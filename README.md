# A High-Throughput Screening (HTS) Workflow Using UniDock

A set of Python scripts to perform a high-throughput virtual screening workflow using Uni-Dock, targeting TREM2. I've tried to make this workflow operable for possibly other datasets and/or receptors, but there are some pecularities. There should be markers in scripts that need to be modified otherwise

Naming convention for the results:
[Loop type]_[config type]\_[column numbers].dat
Ie, `CL_1_3.dat`; the closed loop configuration, using the first docking box config, and the third column 

## Workflow Steps

### 0. Prepare the Receptor (Manual Step)

1.  Obtain the PDB structure of your target protein (e.g., TREM2).
2.  Prepare the receptor for docking. This usually involves:
    *   Removing water molecules and other non-essential heteroatoms.
    *   Adding hydrogen atoms (if missing).
    *   Repairing missing residues or side chains (if necessary).
    *   Converting the receptor to PDBQT format.
        *   **In this example, we use Open Babel**:
            ```bash
            obabel TREM2.pdb -O TREM2.pdbqt -xr
            ```
            (The `-xr` flag tells Open Babel to prepare the receptor by removing non-standard atoms, adding hydrogens, and assigning atom types suitable for AutoDock Vina/Uni-Dock.)
3.  Place the prepared receptor file into the `data/receptor/` directory.

### 1. Download Ligand Dataset (01_getdata.py)

This script downloads ligand datasets. Its currently only set up for ZINC22 `.uri` file types, which can be obtained from the [official Zinc website](https://zinc.docking.org/tranches/home/).

Select the 3D representation and the cells that have the desired properties. For TREM2, we scanned all molecules between 1 and 3 LogP with up to 500 Daltons in molecular weight.

Press the download button near the top right, and you should see an option an option to change the file type. Make sure it downloads as `.pdbqt.gz`. When you press download, the result will be a series of URLs. The next step will unpack this URL.

### 3. Configure and Run Docking (dock.py)

This script sets up the Uni-Dock configuration and runs the docking calculations.

*   **Modify `scripts/03_run_docking.py`**:
    *   **Receptor File**: Ensure `RECEPTOR_FILE` points to your prepared receptor in `data/receptor/`.
    *   **Ligand Directory**: `LIGAND_DIR` should point to `data/ligands_prepared/`.
    *   **Uni-Dock Executable**: Set `UNIDOCK_EXECUTABLE` to `unidock` or full path of your Uni-Dock executable.
    *   **Binding Site Definition**: Update `CENTER_X, CENTER_Y, CENTER_Z` and `SIZE_X, SIZE_Y, SIZE_Z` with the coordinates and dimensions of the binding pocket on your receptor.
    *   You can also adjust docking parameters like `num_modes`, `exhaustiveness`, and `scoring_function`.
*   **Run the script**:
    ```bash
    python 03_run_docking.py
    ```
    This script will:
    1.  Create a `docking_config.txt` in the `configs/` directory.
    2.  Run Uni-Dock using the specified receptor, ligands, and configuration.
    3.  Save docking outputs (e.g., docked poses in PDBQT format) into `results/docking_outputs/. Each ligand will have its own file.

### 4. Analyze Results

This script parses the docking output files, extracts binding scores, and generates a summary.

*   **Modify `scripts/04_analyze_results.py` (if needed)**:
    *   The script currently expects Uni-Dock to produce PDBQT files (Vina-like format) where scores are in `REMARK VINA RESULT:` lines.
    *   The `search_pattern` for finding output files might need adjustment based on how Uni-Dock structures its output directory (e.g., if it creates subfolders for each ligand or names files differently).
*   **Run the script**:
    ```bash
    python 04_analyze_results.py
    ```
    This will:
    1.  Scan the `results/docking_outputs/` directory for PDBQT files.
    2.  Parse each file to extract binding affinities (scores).
    3.  Save a `docking_summary.csv` file in the `results/` directory, sorted by binding affinity.
    4.  Print the top N results to the console.

## Extending the Workflow
Note: Running MCDock also requires AmberTools.
*   **Other Datasets**: To use a different ligand dataset, place the raw ligand files (preferably in SMILES format) into `data/ligands_raw/` and re-run the pipeline from step 3 (or step 2 if you modify the download script).
*   **Other Targets**: 
    1.  Prepare your new target receptor and place its PDBQT file in `data/receptor/`.
    2.  Update the `RECEPTOR_FILE` variable in `scripts/03_run_docking.py`.
    3.  **Crucially**, update the binding site coordinates (`CENTER_X,Y,Z` and `SIZE_X,Y,Z`) in `scripts/03_run_docking.py` for the new target.
*   **Different Ligand Preparation Tools**: Modify `scripts/02_prepare_ligands.py` to use different tools or steps for ligand preparation if needed.
*   **Uni-Dock Parameters**: Adjust Uni-Dock parameters in `scripts/03_run_docking.py` (e.g., `scoring_function`, `exhaustiveness`, `num_modes`) as required for your specific study.
*   **Output Analysis**: Enhance `scripts/04_analyze_results.py` for more sophisticated analysis, visualization, or filtering of results.

## Important Notes



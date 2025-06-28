#!/bin/bash
#SBATCH --job-name=zinc_download
#SBATCH --partition=amilan
#SBATCH --account=ucb-general
#SBATCH --qos=normal
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --output=zinc_download_%j.out
#SBATCH --error=zinc_download_%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=noso3320@colorado.edu

# Load necessary modules
module purge
module load gcc
module load cmake

# Set working directory paths to match project structure
SCRIPT_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock/scripts"
DATA_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock/data"
WORK_DIR="/scratch/alpine/noso3320/trem2/trem2_unidock"
cd $WORK_DIR

# Log start time
echo "Starting ZINC database download at $(date)"
echo "Working directory: $WORK_DIR"
echo "Script directory: $SCRIPT_DIR"
echo "Data directory: $DATA_DIR"

# Function to download with retry logic
download_with_retry() {
    local url=$1
    local filename=$(basename $url)
    local max_attempts=3
    local attempt=1
    
    while [ $attempt -le $max_attempts ]; do
        echo "Attempt $attempt of $max_attempts for $filename"
        
        # Try wget first
        if wget --continue --timeout=600 --connect-timeout=60 --read-timeout=300 --tries=3 --retry-connrefused --waitretry=10 --progress=bar:force "$url"; then
            echo "Successfully downloaded $filename with wget"
            return 0
        else
            echo "wget failed for $filename, trying curl..."
            # Try curl as fallback
            if curl -L -C - --connect-timeout 60 --max-time 600 --retry 3 --retry-delay 10 -o "$filename" "$url"; then
                echo "Successfully downloaded $filename with curl"
                return 0
            else
                echo "Both wget and curl failed for $filename (attempt $attempt)"
                attempt=$((attempt + 1))
                sleep 60  # Wait 1 minute before retry
            fi
        fi
    done
    
    echo "ERROR: Failed to download $filename after $max_attempts attempts"
    return 1
}

# Create log file
LOG_FILE="$WORK_DIR/download_log_$(date +%Y%m%d_%H%M%S).txt"
exec 1> >(tee -a "$LOG_FILE")
exec 2> >(tee -a "$LOG_FILE" >&2)

# Download specific ZINC subsets (modify URLs as needed)
echo "Downloading ZINC database subsets..."

# Read URLs from column_one.uri file
URI_FILE="$DATA_DIR/column_one.uri"

# Check if URI file exists
if [ ! -f "$URI_FILE" ]; then
    echo "ERROR: $URI_FILE not found! Please ensure the file is in the data directory."
    exit 1
fi

# Test connectivity to files.docking.org
echo "Testing connectivity to files.docking.org..."
if ping -c 3 files.docking.org > /dev/null 2>&1; then
    echo "✓ Can ping files.docking.org"
else
    echo "⚠️  Cannot ping files.docking.org - may still work via HTTP"
fi

# Test basic HTTP connectivity
echo "Testing HTTP connectivity..."
if wget --spider --timeout=30 http://files.docking.org > /dev/null 2>&1; then
    echo "✓ Can reach files.docking.org via HTTP"
else
    echo "⚠️  Cannot reach files.docking.org via HTTP - downloads may fail"
fi

# Count total URLs
TOTAL_URLS=$(wc -l < "$URI_FILE")
echo "Found $TOTAL_URLS URLs to download from $URI_FILE"

# TESTING: Limit to first 3 URLs for quick test
TEST_LIMIT=3
echo "🧪 TESTING MODE: Processing only first $TEST_LIMIT URLs (out of $TOTAL_URLS total)"

# Download each file from URI list
FAILED_DOWNLOADS=()
CURRENT_COUNT=0

while IFS= read -r url; do
    # Skip empty lines
    if [ -z "$url" ]; then
        continue
    fi
    
    CURRENT_COUNT=$((CURRENT_COUNT + 1))
    echo "[$CURRENT_COUNT/$TEST_LIMIT] Processing: $url"
    
    if ! download_with_retry "$url"; then
        FAILED_DOWNLOADS+=("$url")
    fi
    
    # TESTING: Break after processing the test limit
    if [ $CURRENT_COUNT -ge $TEST_LIMIT ]; then
        echo "🧪 TESTING: Reached limit of $TEST_LIMIT files, stopping early"
        break
    fi
    
    # Progress update every 100 files
    if [ $((CURRENT_COUNT % 100)) -eq 0 ]; then
        echo "Progress: $CURRENT_COUNT/$TOTAL_URLS files processed"
        echo "Failed so far: ${#FAILED_DOWNLOADS[@]}"
        echo "Current time: $(date)"
    fi
    
done < "$URI_FILE"

# Don't extract .pdbqt.gz files - keep them compressed
echo "Download completed - PDBQT files kept in compressed format"

# Generate summary
echo "========================================="
echo "Download Summary:"
echo "========================================="
echo "Completed at: $(date)"
echo "Total files processed: $TOTAL_URLS"
echo "Successfully downloaded: $((TOTAL_URLS - ${#FAILED_DOWNLOADS[@]}))"
echo "Failed downloads: ${#FAILED_DOWNLOADS[@]}"

if [ ${#FAILED_DOWNLOADS[@]} -gt 0 ]; then
    echo "Failed URLs:"
    for failed_url in "${FAILED_DOWNLOADS[@]}"; do
        echo "  - $failed_url"
    done
fi

# Show disk usage
echo "Disk usage:"
du -sh $WORK_DIR

# Check available space
echo "Available space in scratch:"
df -h /scratch/summit/$USER

echo "ZINC database download job completed!"

# Optional: Copy to data directory for use by other scripts
echo "Copying downloaded files to data directory..."
LIGANDS_RAW_DIR="$DATA_DIR/column_one/ligands_raw"
mkdir -p $LIGANDS_RAW_DIR
rsync -av --progress $WORK_DIR/*.pdbqt.gz $LIGANDS_RAW_DIR/
echo "Copy to data directory completed!"

exit 0

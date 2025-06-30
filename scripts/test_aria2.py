#!/usr/bin/env python3
"""
Test aria2c as alternative download method for ZINC data
"""
import subprocess
import os
import sys

def test_aria2c():
    """Test if aria2c is available and works better than requests"""
    print("=== TESTING ARIA2C ALTERNATIVE ===")
    
    # Check if aria2c is installed
    try:
        result = subprocess.run(['aria2c', '--version'], capture_output=True, text=True)
        if result.returncode == 0:
            print("✓ aria2c is installed")
            print(f"   Version: {result.stdout.split()[2]}")
        else:
            print("✗ aria2c not found - installing...")
            subprocess.run(['conda', 'install', '-c', 'bioconda', 'aria2', '-y'])
    except FileNotFoundError:
        print("✗ aria2c not found - installing...")
        subprocess.run(['conda', 'install', '-c', 'bioconda', 'aria2', '-y'])
    
    # Create test URI file with just a few URLs
    test_urls = [
        "http://files.docking.org/3D/AC/AAMP/ACAAMP.xaa.pdbqt.gz",
        "http://files.docking.org/3D/AC/AAMM/ACAAMM.xaa.pdbqt.gz",
        "http://files.docking.org/3D/AC/AAMN/ACAAMN.xaa.pdbqt.gz"
    ]
    
    with open('test_aria2.uri', 'w') as f:
        for url in test_urls:
            f.write(url + '\n')
    
    print(f"\n📁 Created test URI file with {len(test_urls)} URLs")
    
    # Test aria2c download
    print("\n🚀 Testing aria2c download...")
    try:
        cmd = [
            'aria2c',
            '-i', 'test_aria2.uri',  # Input file
            '-x', '4',               # 4 connections per file
            '-j', '3',               # 3 parallel downloads
            '--continue=true',       # Resume interrupted downloads
            '--max-tries=3',         # Retry 3 times
            '--retry-wait=5',        # Wait 5s between retries
            '--timeout=300',         # 5 minute timeout
            '--connect-timeout=30',  # 30s connection timeout
            '--dir=test_download'    # Output directory
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        
        if result.returncode == 0:
            print("✓ aria2c download successful!")
            
            # Check downloaded files
            if os.path.exists('test_download'):
                files = os.listdir('test_download')
                total_size = 0
                for file in files:
                    size = os.path.getsize(os.path.join('test_download', file))
                    total_size += size
                    print(f"   Downloaded: {file} ({size:,} bytes)")
                print(f"   Total: {len(files)} files, {total_size:,} bytes")
                return True
        else:
            print("✗ aria2c download failed")
            print(f"   stdout: {result.stdout}")
            print(f"   stderr: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ aria2c download timed out")
        return False
    except Exception as e:
        print(f"✗ aria2c error: {e}")
        return False

def suggest_alternatives():
    """Suggest other alternatives if aria2c fails"""
    print("\n=== OTHER ALTERNATIVES ===")
    print("1. 🌐 ZINC API + curl:")
    print("   curl -o molecules.mol2 --data-urlencode zinc.ids='ZINC123456' -d page.format=mol2 http://zinc.docking.org/results")
    
    print("\n2. 📦 Pre-built libraries:")
    print("   - ChEMBL: https://www.ebi.ac.uk/chembl/")
    print("   - PubChem: https://pubchem.ncbi.nlm.nih.gov/")
    print("   - GitHub mirrors: https://github.com/quantaosun/Zinc-Million")
    
    print("\n3. 🏠 Local download + transfer:")
    print("   - Download on local machine with good internet")
    print("   - Transfer via rsync/scp to HPC")
    
    print("\n4. 🎯 Smaller subsets:")
    print("   - Download specific tranches instead of everything")
    print("   - Use drug-like/lead-like subsets only")

if __name__ == "__main__":
    success = test_aria2c()
    
    if success:
        print("\n🎉 SUCCESS: aria2c works! Use this instead of your Python script:")
        print("   aria2c -i data/column_one.uri -x 8 -j 8 --continue=true --max-tries=5")
    else:
        suggest_alternatives()
        
    # Cleanup
    if os.path.exists('test_aria2.uri'):
        os.remove('test_aria2.uri') 

#!/usr/bin/env python3
"""
Debug script to test connectivity to files.docking.org
"""
import requests
import time
import sys

def test_connectivity():
    """Test basic connectivity to files.docking.org"""
    print("=== DEBUGGING DOWNLOAD CONNECTIVITY ===")
    
    # Test 1: Basic connectivity
    print("\n1. Testing basic connectivity...")
    test_url = "http://files.docking.org"
    try:
        response = requests.get(test_url, timeout=30)
        print(f"✓ Can reach {test_url} - Status: {response.status_code}")
    except Exception as e:
        print(f"✗ Cannot reach {test_url}: {e}")
        return False
    
    # Test 2: Try downloading a small file
    print("\n2. Testing small file download...")
    # Use first URL from your URI file
    test_urls = [
        "http://files.docking.org/3D/AC/AAMP/ACAAMP.xaa.pdbqt.gz",
        "http://files.docking.org/3D/AC/AAMM/ACAAMM.xaa.pdbqt.gz"
    ]
    
    for i, url in enumerate(test_urls, 1):
        print(f"\nTest {i+1}: Downloading {url.split('/')[-1]}...")
        try:
            start_time = time.time()
            response = requests.get(url, stream=True, timeout=300)
            
            # Get file size
            total_size = int(response.headers.get('content-length', 0))
            print(f"   File size: {total_size:,} bytes")
            
            # Download with progress
            downloaded = 0
            chunk_count = 0
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    downloaded += len(chunk)
                    chunk_count += 1
                    if chunk_count % 100 == 0:  # Progress every 100 chunks
                        elapsed = time.time() - start_time
                        speed = downloaded / elapsed / 1024  # KB/s
                        print(f"   Progress: {downloaded:,}/{total_size:,} bytes ({speed:.1f} KB/s)")
            
            elapsed = time.time() - start_time
            speed = downloaded / elapsed / 1024  # KB/s
            print(f"✓ Download completed: {downloaded:,} bytes in {elapsed:.1f}s ({speed:.1f} KB/s)")
            return True
            
        except requests.exceptions.Timeout as e:
            elapsed = time.time() - start_time
            print(f"✗ TIMEOUT after {elapsed:.1f}s: {e}")
        except Exception as e:
            elapsed = time.time() - start_time
            print(f"✗ ERROR after {elapsed:.1f}s: {e}")
    
    return False

def test_network_diagnostics():
    """Run network diagnostics"""
    print("\n=== NETWORK DIAGNOSTICS ===")
    
    import subprocess
    import os
    
    # Test ping
    print("\n1. Testing ping to files.docking.org...")
    try:
        result = subprocess.run(['ping', '-c', '3', 'files.docking.org'], 
                              capture_output=True, text=True, timeout=30)
        if result.returncode == 0:
            print("✓ Ping successful")
            print(result.stdout)
        else:
            print("✗ Ping failed")
            print(result.stderr)
    except Exception as e:
        print(f"✗ Ping test failed: {e}")
    
    # Test DNS resolution
    print("\n2. Testing DNS resolution...")
    try:
        import socket
        ip = socket.gethostbyname('files.docking.org')
        print(f"✓ DNS resolved: files.docking.org -> {ip}")
    except Exception as e:
        print(f"✗ DNS resolution failed: {e}")
    
    # Check environment
    print("\n3. Environment info...")
    print(f"   Python: {sys.version}")
    print(f"   Working directory: {os.getcwd()}")
    
    # Check for proxy settings
    http_proxy = os.environ.get('HTTP_PROXY') or os.environ.get('http_proxy')
    https_proxy = os.environ.get('HTTPS_PROXY') or os.environ.get('https_proxy')
    if http_proxy or https_proxy:
        print(f"   HTTP proxy: {http_proxy}")
        print(f"   HTTPS proxy: {https_proxy}")
    else:
        print("   No proxy settings detected")

if __name__ == "__main__":
    success = test_connectivity()
    test_network_diagnostics()
    
    if not success:
        print("\n🚨 DIAGNOSIS: Download connectivity issues detected!")
        print("\nPossible solutions:")
        print("1. Check if your HPC has outbound internet restrictions")
        print("2. Ask HPC support about accessing external websites")
        print("3. Try downloading from a login node instead of compute node")
        print("4. Check if files.docking.org is accessible from your location")
        print("5. Consider downloading files locally and transferring to HPC")
    else:
        print("\n✓ DIAGNOSIS: Downloads should work - check getdata.py logic") 

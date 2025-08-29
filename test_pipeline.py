#!/usr/bin/env python3
"""
Test script to verify TRACE pipeline components are properly installed.

Authors: Michal Izydorczyk, Nicola Wong, John Adedeji, Julian Chiu, and Thomas X. Garcia
"""

import subprocess
import sys
import os
from pathlib import Path

def check_python_module(module_name):
    """Check if a Python module can be imported."""
    try:
        __import__(module_name)
        return True, "OK"
    except ImportError as e:
        return False, str(e)

def check_script_exists(script_name):
    """Check if a pipeline script exists."""
    script_path = Path(__file__).parent / script_name
    if script_path.exists():
        return True, "OK"
    else:
        return False, f"Not found at {script_path}"

def check_command(command):
    """Check if a system command is available."""
    try:
        result = subprocess.run(
            command.split() + ['--version'],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            return True, "OK"
        else:
            return True, "Found (version check failed)"
    except FileNotFoundError:
        return False, "Not found in PATH"
    except subprocess.TimeoutExpired:
        return True, "Found (timeout on version check)"
    except Exception as e:
        return False, str(e)

def main():
    """Run all checks."""
    print("=" * 60)
    print("TRACE Pipeline Component Check")
    print("=" * 60)
    
    # Check Python version
    print(f"\nPython Version: {sys.version}")
    if sys.version_info < (3, 10):
        print("WARNING: Python 3.10+ is recommended")
    
    # Check pipeline scripts
    print("\n--- Pipeline Scripts ---")
    scripts = [
        "TRACE.py",
        "vcf2fasta.py",
        "out2bed.py",
        "squash_repeat_masker.py",
        "vcf2bed.py",
        "bed_processor.py",
        "squash_intersect.py",
        "annotate_vcf.py"
    ]
    
    all_scripts_ok = True
    for script in scripts:
        exists, msg = check_script_exists(script)
        status = "✓" if exists else "✗"
        print(f"{status} {script}: {msg}")
        if not exists:
            all_scripts_ok = False
    
    # Check Python modules
    print("\n--- Python Modules ---")
    modules = ["argparse", "pathlib", "logging", "subprocess", "os", "sys"]
    
    for module in modules:
        exists, msg = check_python_module(module)
        status = "✓" if exists else "✗"
        print(f"{status} {module}: {msg}")
    
    # Check external tools
    print("\n--- External Tools ---")
    print("NOTE: These tools may be in separate conda environments")
    
    tools = {
        "bedtools": "Required for genomic intersections",
        "RepeatMasker": "Required for transposon annotation (usually in separate env)",
        "sniffles": "Optional - for initial SV calling (usually in separate env)"
    }
    
    for tool, description in tools.items():
        exists, msg = check_command(tool)
        status = "✓" if exists else "✗"
        print(f"{status} {tool}: {msg}")
        if not exists:
            print(f"  → {description}")
    
    # Check for reference files
    print("\n--- Reference Files ---")
    data_dir = Path(__file__).parent / "data"
    
    if data_dir.exists():
        print(f"✓ Data directory exists: {data_dir}")
    else:
        print(f"✗ Data directory not found: {data_dir}")
        print("  → Create with: mkdir data")
    
    # Check for specific reference files
    ref_files = {
        "hg38.fa.bed": "RepeatMasker annotations for human genome",
        "Dfam-RepeatMasker.lib": "Repeat library for RepeatMasker"
    }
    
    for filename, description in ref_files.items():
        file_path = data_dir / filename
        if file_path.exists():
            size_mb = file_path.stat().st_size / (1024 * 1024)
            print(f"✓ {filename}: {size_mb:.1f} MB")
        else:
            # Check in parent directory
            parent_path = Path("../") / filename
            if parent_path.exists():
                size_mb = parent_path.stat().st_size / (1024 * 1024)
                print(f"⚠ {filename}: Found in parent directory ({size_mb:.1f} MB)")
                print(f"  → Consider copying/linking to data/ directory")
            else:
                print(f"✗ {filename}: Not found")
                print(f"  → {description}")
    
    # Summary
    print("\n" + "=" * 60)
    if all_scripts_ok:
        print("✓ All pipeline scripts are present")
    else:
        print("✗ Some pipeline scripts are missing")
    
    print("\nTo run the pipeline:")
    print("  python TRACE.py input.vcf \\")
    print("    --intersect data/hg38.fa.bed \\")
    print("    --lib data/Dfam-RepeatMasker.lib \\")
    print("    --threads 24")
    
    print("\nFor tools in separate environments:")
    print("  conda activate repeatmasker  # For RepeatMasker")
    print("  conda activate sniffles      # For Sniffles")
    
    print("=" * 60)

if __name__ == "__main__":
    main()
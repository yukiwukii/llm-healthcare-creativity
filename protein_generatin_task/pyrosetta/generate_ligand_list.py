#!/usr/bin/env python3
"""
generate_ligand_list.py

This script scans a directory for PDB files and creates a ligands.txt file
that can be used with the workaround_docking.py script.

Usage:
  python generate_ligand_list.py --dir pdb_files/llama --output ligands.txt [--pattern "*.pdb"]
"""

import os
import glob
import argparse
from pathlib import Path

def find_pdb_files(directory, pattern="*.pdb", recursive=False):
    """
    Find all PDB files in the given directory
    
    Parameters:
    directory (str): Directory to search
    pattern (str): File pattern to match (default: "*.pdb")
    recursive (bool): Whether to search recursively in subdirectories
    
    Returns:
    list: List of full paths to PDB files
    """
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory not found: {directory}")
    
    # Normalize the directory path
    directory = os.path.abspath(directory)
    
    # Find files using glob
    if recursive:
        search_pattern = os.path.join(directory, "**", pattern)
        pdb_files = glob.glob(search_pattern, recursive=True)
    else:
        search_pattern = os.path.join(directory, pattern)
        pdb_files = glob.glob(search_pattern)
    
    # Sort files for consistent ordering
    pdb_files.sort()
    
    return pdb_files

def generate_ligand_list(directory, output_file, pattern="*.pdb", recursive=False, absolute_paths=True):
    """
    Generate a text file with paths to PDB files
    
    Parameters:
    directory (str): Directory to search
    output_file (str): Path to output file
    pattern (str): File pattern to match (default: "*.pdb")
    recursive (bool): Whether to search recursively in subdirectories
    absolute_paths (bool): Whether to use absolute paths in the output file
    
    Returns:
    int: Number of PDB files found
    """
    try:
        # Find PDB files
        pdb_files = find_pdb_files(directory, pattern, recursive)
        
        if not pdb_files:
            print(f"No PDB files found in {directory} matching pattern '{pattern}'")
            return 0
        
        # Convert to relative paths if requested
        if not absolute_paths:
            base_dir = os.path.abspath(os.path.dirname(output_file))
            pdb_files = [os.path.relpath(pdb_path, base_dir) for pdb_path in pdb_files]
        
        # Write to output file
        with open(output_file, 'w') as f:
            f.write("# Automatically generated ligand list\n")
            f.write(f"# Source directory: {os.path.abspath(directory)}\n")
            f.write(f"# Pattern: {pattern}\n")
            f.write(f"# Generated on: {os.popen('date').read().strip()}\n")
            f.write(f"# Total files: {len(pdb_files)}\n\n")
            
            for pdb_file in pdb_files:
                f.write(f"{pdb_file}\n")
        
        print(f"Successfully wrote {len(pdb_files)} PDB file paths to {output_file}")
        return len(pdb_files)
        
    except Exception as e:
        print(f"Error generating ligand list: {e}")
        return 0

def main():
    parser = argparse.ArgumentParser(description="Generate a list of PDB files for docking")
    parser.add_argument("--dir", required=True, help="Directory containing PDB files")
    parser.add_argument("--output", required=True, help="Output file path")
    parser.add_argument("--pattern", default="*.pdb", help="File pattern (default: *.pdb)")
    parser.add_argument("--recursive", action="store_true", help="Search recursively in subdirectories")
    parser.add_argument("--relative", action="store_true", help="Use relative paths instead of absolute paths")
    
    args = parser.parse_args()
    
    # Ensure output directory exists
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Generate the ligand list
    num_files = generate_ligand_list(
        args.dir,
        args.output,
        args.pattern,
        args.recursive,
        not args.relative  # Invert the flag - we want absolute paths by default
    )
    
    if num_files > 0:
        # Print example command
        print("\nExample command for docking:")
        print(f"python workaround_docking.py --receptor your_receptor.pdb --ligand_list {args.output} --output_dir docking_results")

if __name__ == "__main__":
    main()

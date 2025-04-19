#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
HSSP-Based Protein Novelty Analysis Script

This script takes existing BLAST results and a directory of FASTA files
to apply HSSP curve analysis for determining protein novelty with
length-dependent thresholds.

Configuration is hardcoded in the script for simplicity.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import re
from pathlib import Path
import logging
import sys
from Bio import SeqIO
import time

# ====== CONFIGURATION - MODIFY THESE VALUES ======
# Path to your BLAST results CSV file
INPUT_CSV = './novelty/protein_novelty/gemma/novelty_report_20250417_181729.csv'

# Directory where results will be saved
OUTPUT_DIR = './novelty/novelty_results/gemma/'

# Directory containing your FASTA files
FASTA_DIR = './novelty/fasta_files/gemma/'
# ===============================================

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("hssp_analysis.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def calculate_hssp_threshold(length):
    """
    Calculate the identity threshold for structural homology based on the HSSP curve.
    
    The HSSP curve is a length-dependent threshold where shorter sequences need
    higher identity to be considered similar.
    
    Using the Rost (1999) formulation which produces values in the valid 0-100% range.
    
    References:
    - Rost, B. (1999). Twilight zone of protein sequence alignments. Protein Engineering, 12(2), 85-94
    """
    if length <= 11:
        # For very short sequences, require very high identity
        return 90.0
    
    # Rost (1999) formulation: threshold = 290.15 * L^(-0.562) + 5
    threshold = 290.15 * (length ** -0.562) + 5
    
    # Cap at 100% since identity can't be higher than that
    return min(100.0, threshold)

def extract_length_from_id(seq_id):
    """Try to extract sequence length from ID based on common naming patterns."""
    # Pattern 1: seq_L20 (L followed by number)
    match = re.search(r'[_-]L(\d+)', seq_id)
    if match:
        return int(match.group(1))
    
    # Pattern 2: seq_20aa (number followed by 'aa')
    match = re.search(r'[_-](\d+)aa', seq_id)
    if match:
        return int(match.group(1))
    
    # Pattern 3: just look for a number after the last underscore or dash
    match = re.search(r'[_-](\d+)$', seq_id)
    if match:
        return int(match.group(1))
    
    # Default if no length can be extracted
    return None

def find_fasta_files(directory):
    """Find all FASTA files in directory and subdirectories."""
    fasta_extensions = ('*.fasta', '*.fa', '*.fna', '*.faa')
    found_files = []
    
    # Search for FASTA files in the directory and subdirectories
    for ext in fasta_extensions:
        found_files.extend(list(Path(directory).glob(ext)))
        found_files.extend(list(Path(directory).glob(f'**/{ext}')))
    
    return found_files

def get_sequence_lengths(blast_results, fasta_dir=None):
    """Get sequence lengths either from FASTA files or by estimating from sequence IDs."""
    seq_lengths = {}
    sequence_ids = set(blast_results['Sequence ID'])
    
    # If FASTA directory is provided, read lengths from files
    if fasta_dir and os.path.exists(fasta_dir):
        logger.info(f"Reading sequence lengths from FASTA files in {fasta_dir}")
        
        # Find all FASTA files
        fasta_files = find_fasta_files(fasta_dir)
        logger.info(f"Found {len(fasta_files)} FASTA files to scan")
        
        # Track matched sequence IDs to avoid duplicate processing
        matched_ids = set()
        
        # Process all found FASTA files
        for i, file_path in enumerate(fasta_files):
            try:
                if (i+1) % 10 == 0 or i+1 == len(fasta_files):
                    logger.info(f"Processing file {i+1}/{len(fasta_files)}: {file_path}")
                
                # Read all sequences from this file
                for record in SeqIO.parse(file_path, 'fasta'):
                    # Try exact match first
                    if record.id in sequence_ids and record.id not in matched_ids:
                        seq_lengths[record.id] = len(record.seq)
                        matched_ids.add(record.id)
                    
                    # Try partial matches for records not yet matched
                    # (BLAST sometimes truncates IDs or formats them differently)
                    for seq_id in sequence_ids:
                        if seq_id not in matched_ids:
                            # Check if one is prefix of the other
                            if seq_id.startswith(record.id) or record.id.startswith(seq_id):
                                seq_lengths[seq_id] = len(record.seq)
                                matched_ids.add(seq_id)
                                break
                            
                            # Try removing common prefix formats (gi|, lcl|, etc.)
                            clean_record_id = re.sub(r'^[a-z]+\|', '', record.id)
                            clean_seq_id = re.sub(r'^[a-z]+\|', '', seq_id)
                            
                            if clean_record_id == clean_seq_id:
                                seq_lengths[seq_id] = len(record.seq)
                                matched_ids.add(seq_id)
                                break
                
                # If we've matched all sequences, we can stop searching
                if len(matched_ids) == len(sequence_ids):
                    logger.info("Found all sequence lengths, stopping file search")
                    break
                    
            except Exception as e:
                logger.warning(f"Error parsing FASTA file {file_path}: {e}")
        
        logger.info(f"Successfully found lengths for {len(matched_ids)}/{len(sequence_ids)} sequences from FASTA files")
    
    # For any missing lengths, try to extract from sequence IDs
    missing_lengths = []
    for seq_id in sequence_ids:
        if seq_id not in seq_lengths:
            length = extract_length_from_id(seq_id)
            if length:
                seq_lengths[seq_id] = length
            else:
                missing_lengths.append(seq_id)
    
    if len(missing_lengths) < len(sequence_ids):
        logger.info(f"Extracted {len(sequence_ids) - len(missing_lengths)} sequence lengths from IDs")
    
    # For remaining sequences with unknown length, use average length or a reasonable default
    if missing_lengths:
        logger.warning(f"Could not determine length for {len(missing_lengths)} sequences")
        
        # Use average length from known sequences or default to 30
        if seq_lengths:
            avg_length = sum(seq_lengths.values()) / len(seq_lengths)
            default_length = round(avg_length)
        else:
            default_length = 30  # Reasonable default for short proteins
        
        logger.warning(f"Using default length of {default_length} for sequences with unknown length")
        
        for seq_id in missing_lengths:
            seq_lengths[seq_id] = default_length
    
    return seq_lengths

def plot_hssp_curve(output_dir):
    """Generate a plot of the HSSP curve for reference."""
    lengths = list(range(10, 51))  # Ensure we plot from 10 to 50 inclusive
    thresholds = [calculate_hssp_threshold(length) for length in lengths]
    
    plt.figure(figsize=(10, 6))
    plt.plot(lengths, thresholds, 'b-', linewidth=2, label='HSSP threshold')
    
    plt.xlabel('Sequence Length')
    plt.ylabel('Identity Threshold (%)')
    plt.title('HSSP Curve: Identity Threshold vs Sequence Length')
    plt.xlim(10, 50)  # Set x-axis limits explicitly
    plt.ylim(0, 100)  # Set y-axis limits to standard percentage range
    plt.grid(True)
    plt.legend()
    
    # Save the plot
    os.makedirs(output_dir, exist_ok=True)
    plot_path = os.path.join(output_dir, 'hssp_curve.png')
    plt.savefig(plot_path)
    logger.info(f"HSSP curve plot saved to {plot_path}")
    
    # Close the plot to free memory
    plt.close()

def analyze_with_hssp(blast_results, seq_lengths, output_dir):
    """Analyze novelty using the HSSP curve."""
    # Create a copy of the results and add columns for HSSP analysis
    results_df = blast_results.copy()
    
    # Add sequence length column
    results_df['Sequence Length'] = results_df['Sequence ID'].map(seq_lengths)
    
    # Calculate HSSP threshold for each sequence
    results_df['HSSP Threshold'] = results_df['Sequence Length'].apply(calculate_hssp_threshold)
    
    # Determine novelty status using HSSP threshold
    def determine_hssp_novelty(row):
        if row['Novelty Status'] == 'Unknown (BLAST failed)':
            return 'Unknown (BLAST failed)'
        if pd.isna(row['Top Identity (%)']):
            return 'Novel'  # No hits
        return 'Novel' if row['Top Identity (%)'] < row['HSSP Threshold'] else 'Known'
    
    results_df['HSSP Novelty Status'] = results_df.apply(determine_hssp_novelty, axis=1)
    
    # Calculate percentages above and below the HSSP curve
    valid_results = results_df[results_df['Top Identity (%)'].notna()]
    
    # Split data into categories for plotting
    novel_points = valid_results[valid_results['HSSP Novelty Status'] == 'Novel']
    known_points = valid_results[valid_results['HSSP Novelty Status'] == 'Known']
    
    # Calculate percentages based on total number of points that will be plotted
    total_plotted = len(novel_points) + len(known_points)
    novel_percent = len(novel_points) / total_plotted * 100 if total_plotted > 0 else 0
    known_percent = len(known_points) / total_plotted * 100 if total_plotted > 0 else 0
    
    logger.info(f"HSSP threshold analysis:")
    logger.info(f"  Total sequences with valid identity plotted: {total_plotted}")
    logger.info(f"  Sequences above HSSP curve (Known): {len(known_points)} ({known_percent:.1f}%)")
    logger.info(f"  Sequences below HSSP curve (Novel): {len(novel_points)} ({novel_percent:.1f}%)")
    
    # Create a visualization of the sequences
    plt.figure(figsize=(12, 8))
    
    # Plot HSSP curve
    lengths = list(range(10, 51))
    thresholds = [calculate_hssp_threshold(length) for length in lengths]
    plt.plot(lengths, thresholds, 'g-', linewidth=2, label='HSSP threshold')
    
    # Scatter plot of novel and known sequences
    plt.scatter(
        novel_points['Sequence Length'], 
        novel_points['Top Identity (%)'],
        c='blue', alpha=0.6, label=f'Novel ({len(novel_points)} sequences, {novel_percent:.1f}%)'
    )
    
    plt.scatter(
        known_points['Sequence Length'], 
        known_points['Top Identity (%)'],
        c='red', alpha=0.6, label=f'Known ({len(known_points)} sequences, {known_percent:.1f}%)'
    )
    
    plt.xlabel('Sequence Length')
    plt.ylabel('Identity (%)')
    plt.title('Protein Sequence Novelty Analysis using HSSP Threshold')
    plt.xlim(9, 51)  # Set x-axis limits with slight padding
    plt.ylim(0, 100)  # Set y-axis to standard percentage range
    plt.grid(True)
    plt.legend()
    
    # Save the plot
    os.makedirs(output_dir, exist_ok=True)
    comparison_plot_path = os.path.join(output_dir, 'hssp_analysis.png')
    plt.savefig(comparison_plot_path)
    logger.info(f"Analysis plot saved to {comparison_plot_path}")
    
    # Save the updated results
    results_path = os.path.join(output_dir, 'novelty_hssp_results.csv')
    results_df.to_csv(results_path, index=False)
    logger.info(f"Updated results with HSSP analysis saved to {results_path}")
    
    # Generate summary text file
    summary_path = os.path.join(output_dir, 'hssp_summary.txt')
    with open(summary_path, 'w') as f:
        f.write(f"HSSP Curve Analysis Summary\n")
        f.write(f"=========================\n\n")
        f.write(f"Total sequences analyzed: {len(results_df)}\n")
        f.write(f"Sequences plotted with valid identity values: {total_plotted}\n\n")
        
        f.write(f"HSSP threshold results:\n")
        f.write(f"  Sequences above curve (Known): {len(known_points)} ({known_percent:.1f}%)\n")
        f.write(f"  Sequences below curve (Novel): {len(novel_points)} ({novel_percent:.1f}%)\n")
        f.write(f"  Unknown (BLAST failed): {len(results_df) - total_plotted} ({(len(results_df) - total_plotted)/len(results_df)*100:.1f}%)\n\n")
        
        # Add information about sequence lengths
        lengths = list(seq_lengths.values())
        f.write(f"Sequence Length Statistics:\n")
        f.write(f"  Minimum length: {min(lengths)}\n")
        f.write(f"  Maximum length: {max(lengths)}\n")
        f.write(f"  Average length: {sum(lengths)/len(lengths):.1f}\n\n")
        
        # Add details of the HSSP curve
        f.write(f"HSSP Curve Information:\n")
        f.write(f"  Formula: Threshold = 290.15 * (length^-0.562) + 5 (Rost, 1999)\n")
        f.write(f"  For length=10: Threshold = {calculate_hssp_threshold(10):.1f}%\n")
        f.write(f"  For length=20: Threshold = {calculate_hssp_threshold(20):.1f}%\n")
        f.write(f"  For length=30: Threshold = {calculate_hssp_threshold(30):.1f}%\n")
        f.write(f"  For length=40: Threshold = {calculate_hssp_threshold(40):.1f}%\n")
        f.write(f"  For length=50: Threshold = {calculate_hssp_threshold(50):.1f}%\n\n")
        
    logger.info(f"Summary saved to {summary_path}")
    
    # Close the plot to free memory
    plt.close()
    
    return results_df

def main():
    start_time = time.time()
    
    # Log configuration
    logger.info(f"HSSP Analysis Configuration:")
    logger.info(f"  Input CSV: {INPUT_CSV}")
    logger.info(f"  Output Directory: {OUTPUT_DIR}")
    logger.info(f"  FASTA Directory: {FASTA_DIR}")
    
    # Read BLAST results
    try:
        blast_results = pd.read_csv(INPUT_CSV)
        logger.info(f"Loaded {len(blast_results)} BLAST results from {INPUT_CSV}")
    except Exception as e:
        logger.error(f"Error loading BLAST results: {e}")
        return
    
    # Get sequence lengths
    seq_lengths = get_sequence_lengths(blast_results, FASTA_DIR)
    
    # Generate HSSP curve plot
    plot_hssp_curve(OUTPUT_DIR)
    
    # Analyze with HSSP curve
    analyze_with_hssp(blast_results, seq_lengths, OUTPUT_DIR)
    
    elapsed_time = time.time() - start_time
    logger.info(f"Analysis completed in {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.exception(f"An unexpected error occurred: {e}")
        sys.exit(1)
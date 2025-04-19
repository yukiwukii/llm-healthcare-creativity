#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Short Protein Novelty Analysis using BLAST

This script analyzes FASTA files to determine the novelty of short protein sequences (10-50 amino acids)
"""

import os
import time
import argparse
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import logging
import random
import sys

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("protein_analysis.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Analyze short protein sequence novelty using BLAST')
    parser.add_argument('--input_dir', type=str, required=True, 
                        help='Directory containing FASTA files')
    parser.add_argument('--output_dir', type=str, required=True, 
                        help='Directory to save results')
    parser.add_argument('--min_length', type=int, default=10,
                        help='Minimum sequence length to analyze (default: 10)')
    parser.add_argument('--max_length', type=int, default=50,
                        help='Maximum sequence length to analyze (default: 50)')
    parser.add_argument('--timeout', type=int, default=450,
                        help='Timeout in seconds for each BLAST request (default: 300)')
    parser.add_argument('--database', type=str, default='swissprot',
                        help='BLAST database to search against (default: swissprot)')
    parser.add_argument('--resume', type=str, default=None,
                        help='Path to previous results CSV file to resume from')
    return parser.parse_args()

def read_fasta_sequences(input_dir, min_length, max_length):
    """Read protein sequences from FASTA files in the input directory."""
    sequences = []
    skipped = 0
    
    for file_path in Path(input_dir).glob('*.fa*'):
        for record in SeqIO.parse(file_path, 'fasta'):
            sequence = str(record.seq)
            seq_length = len(sequence)
            
            # Check if length is within desired range
            if seq_length < min_length or seq_length > max_length:
                skipped += 1
                continue
                
            sequences.append((record.id, sequence, file_path.name))
            
    logger.info(f"Found {len(sequences)} valid protein sequences in {input_dir}")
    logger.info(f"Skipped {skipped} sequences outside length range {min_length}-{max_length}")
    
    return sequences

def run_blast(sequence, sequence_id, timeout=300, database="swissprot"):
    """Submit a sequence to NCBI BLAST with simplified parameters and timeout."""
    logger.info(f"Running BLAST for sequence {sequence_id} (length: {len(sequence)})")
    
    # Import required for timeout functionality
    import signal
    import platform
    
    # Timeout handling is platform-specific
    if platform.system() != "Windows":  # signal.SIGALRM not available on Windows
        class TimeoutError(Exception):
            pass
        
        def timeout_handler(signum, frame):
            raise TimeoutError("BLAST request timed out")
        
        # Set the timeout (not on Windows)
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(timeout)  # Timeout in seconds
    else:
        logger.warning("Timeout functionality not available on Windows. BLAST requests may hang indefinitely.")
    
    try:
        # Use the most basic BLAST parameters with specified database
        result_handle = NCBIWWW.qblast(
            program="blastp",
            database=database,  # Use the database specified by user
            sequence=sequence,
            hitlist_size=10  # Limit results for faster response
        )
        
        # Parse the BLAST results
        blast_records = list(NCBIXML.parse(result_handle))
        result_handle.close()  # Close the handle
        
        # Cancel the timeout (not on Windows)
        if platform.system() != "Windows":
            signal.alarm(0)
        
        if len(blast_records) == 0:
            logger.warning(f"No BLAST records returned for {sequence_id}")
            return None
            
        return blast_records[0]  # Return the first BLAST record
        
    except Exception as e:
        if platform.system() != "Windows":
            # Cancel the timeout
            signal.alarm(0)
            
        if str(e) == "BLAST request timed out":
            logger.error(f"BLAST request for {sequence_id} timed out after {timeout} seconds")
        else:
            logger.error(f"BLAST failed for {sequence_id}: {e}")
        return None

def analyze_novelty(blast_record, sequence_length):
    """Analyze BLAST results to determine sequence novelty."""
    if not blast_record or not blast_record.alignments:
        return {'novel': True, 'top_hit': None}
    
    # Extract information from the top hit
    top_alignment = blast_record.alignments[0]
    top_hsp = top_alignment.hsps[0]
    
    # Calculate percent identity for aligned region
    identical_residues = sum(1 for q, s in zip(top_hsp.query, top_hsp.sbjct) if q == s)
    aligned_length = len(top_hsp.query)
    aligned_identity = (identical_residues / aligned_length) * 100
    
    # Use 85% identity as cutoff - more sensitive to small differences
    is_novel = aligned_identity < 85
    
    return {
        'novel': is_novel,
        'top_hit': top_alignment.title,
        'top_identity': aligned_identity
    }

def generate_report(results, output_dir):
    """Generate a simple report of novelty analysis results."""
    report_data = []
    
    for seq_id, novelty_data in results.items():
        # Define the novelty status text
        if novelty_data['novel'] is None:
            novelty_status = "Unknown (BLAST failed)"
        elif novelty_data['novel']:
            novelty_status = "Novel"
        else:
            novelty_status = "Known"
            
        report_data.append({
            'Sequence ID': seq_id,
            'Novelty Status': novelty_status,
            'Top Hit': novelty_data.get('top_hit', None),
            'Top Identity (%)': round(novelty_data.get('top_identity', 0), 2)
        })
    
    # Create a DataFrame and save to CSV
    df = pd.DataFrame(report_data)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Add timestamp to filename
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Save the report
    report_path = os.path.join(output_dir, f'novelty_report_{timestamp}.csv')
    df.to_csv(report_path, index=False)
    
    # Create a summary file
    novel_count = sum(1 for data in results.values() if data.get('novel') is True)
    known_count = sum(1 for data in results.values() if data.get('novel') is False)
    unknown_count = sum(1 for data in results.values() if data.get('novel') is None)
    total_count = len(results)
    
    summary_path = os.path.join(output_dir, f'summary_{timestamp}.txt')
    with open(summary_path, 'w') as f:
        f.write(f"Short Protein Novelty Analysis Summary\n")
        f.write(f"===================================\n\n")
        f.write(f"Total sequences analyzed: {total_count}\n")
        f.write(f"Novel sequences: {novel_count} ({round(novel_count/total_count*100, 2)}%)\n")
        f.write(f"Known sequences: {known_count} ({round(known_count/total_count*100, 2)}%)\n")
        f.write(f"Unknown (BLAST failed): {unknown_count} ({round(unknown_count/total_count*100, 2)}%)\n\n")
        f.write(f"Identity threshold: 85%\n")
        f.write(f"Detailed results saved to: {report_path}\n")
    
    logger.info(f"Report generated and saved to {report_path}")
    logger.info(f"Summary saved to {summary_path}")
    return df

def main():
    args = parse_arguments()
    
    # Read sequences from FASTA files
    sequences = read_fasta_sequences(args.input_dir, args.min_length, args.max_length)
    
    if not sequences:
        logger.error(f"No valid sequences found in {args.input_dir}")
        return
    
    # Initialize results dictionary
    results = {}
    
    # Check for resumption
    already_processed = set()
    unknown_sequences = {}  # Will store sequences that were unknown and need reprocessing
    
    if args.resume and os.path.exists(args.resume):
        try:
            # Load previous results
            prev_results = pd.read_csv(args.resume)
            for _, row in prev_results.iterrows():
                seq_id = row['Sequence ID']
                
                # Convert back to our internal format
                if row['Novelty Status'] == "Novel":
                    novel = True
                    results[seq_id] = {
                        'novel': novel,
                        'top_hit': row['Top Hit'],
                        'top_identity': row['Top Identity (%)']
                    }
                    already_processed.add(seq_id)
                elif row['Novelty Status'] == "Known":
                    novel = False
                    results[seq_id] = {
                        'novel': novel,
                        'top_hit': row['Top Hit'],
                        'top_identity': row['Top Identity (%)']
                    }
                    already_processed.add(seq_id)
                else:
                    # This sequence was unknown - we'll reprocess it
                    # We need to find its sequence data first
                    for s_id, seq, file_name in sequences:
                        if s_id == seq_id:
                            unknown_sequences[seq_id] = (seq, file_name)
                            break
            
            logger.info(f"Resuming from previous run. Loaded {len(already_processed)} successfully processed sequences.")
            logger.info(f"Found {len(unknown_sequences)} previous unknown results that will be retried.")
        except Exception as e:
            logger.error(f"Error loading previous results for resumption: {e}")
            logger.info("Starting from scratch.")
            unknown_sequences = {}
    
    # Run BLAST and analyze novelty for each sequence
    timeout_count = 0
    
    # Filter sequences to only process new ones (not in already_processed)
    sequences_to_process = [(seq_id, seq, file_name) for seq_id, seq, file_name in sequences 
                           if seq_id not in already_processed]
    
    logger.info(f"Total sequences: {len(sequences)}")
    logger.info(f"Already successfully processed: {len(already_processed)}")
    logger.info(f"Unknown to retry: {len(unknown_sequences)}")
    logger.info(f"New to process: {len(sequences_to_process)}")
    
    # First retry the unknown sequences
    if unknown_sequences:
        logger.info("First retrying sequences that were unknown in previous run...")
        
        for i, (seq_id, (seq, file_name)) in enumerate(unknown_sequences.items()):
            logger.info(f"Retrying unknown sequence {i+1}/{len(unknown_sequences)}: {seq_id}")
            
            # Run BLAST with user-specified timeout
            blast_record = run_blast(seq, seq_id, timeout=args.timeout, database=args.database)
            
            if blast_record is None:
                # Still unknown
                logger.warning(f"Still could not determine novelty for {seq_id} - marking as unknown")
                results[seq_id] = {'novel': None, 'top_hit': 'BLAST analysis failed', 'top_identity': 0}
                timeout_count += 1
            else:
                # Analyze novelty
                novelty_data = analyze_novelty(blast_record, len(seq))
                results[seq_id] = novelty_data
                
                # Log the result
                status = 'Novel' if novelty_data['novel'] else 'Known'
                logger.info(f"Retry successful! Sequence {seq_id} from {file_name}: {status}")
            
            # Take a pause between requests
            time.sleep(1)
            
            # Save intermediate results periodically
            if (i + 1) % 5 == 0 or i == len(unknown_sequences) - 1:
                generate_report(results, args.output_dir)
                logger.info(f"Saved intermediate results after retrying {i+1}/{len(unknown_sequences)} previously unknown sequences")
    
    # Now process the new sequences
    for i, (seq_id, seq, file_name) in enumerate(sequences_to_process):
        logger.info(f"Processing new sequence {i+1}/{len(sequences_to_process)}: {seq_id}")
        
        # Run BLAST with user-specified timeout
        blast_record = run_blast(seq, seq_id, timeout=args.timeout, database=args.database)
        
        if blast_record is None:
            # If BLAST timed out or failed, mark as "unknown" but continue processing
            logger.warning(f"Could not determine novelty for {seq_id} - marking as unknown")
            results[seq_id] = {'novel': None, 'top_hit': 'BLAST analysis failed', 'top_identity': 0}
            timeout_count += 1
        else:
            # Analyze novelty
            novelty_data = analyze_novelty(blast_record, len(seq))
            results[seq_id] = novelty_data
            
            # Log the result
            status = 'Novel' if novelty_data['novel'] else 'Known'
            logger.info(f"Sequence {seq_id} from {file_name}: {status}")
        
        # Take a pause between requests (NCBI recommends no more than 3 requests per second)
        time.sleep(1)
        
        # Save intermediate results every 5 sequences or when finished
        if (i + 1) % 5 == 0 or i == len(sequences_to_process) - 1:
            # Generate report with all results (previously loaded + new ones)
            generate_report(results, args.output_dir)
            logger.info(f"Saved intermediate results after processing {i+1}/{len(sequences_to_process)} new sequences")
        
    if timeout_count > 0:
        logger.warning(f"{timeout_count} sequences timed out or failed during BLAST analysis")
        
    # Final report generation
    generate_report(results, args.output_dir)
    logger.info("Analysis complete. Final results saved.")
    
    # Return final results count
    novel_count = sum(1 for data in results.values() if data.get('novel') is True)
    known_count = sum(1 for data in results.values() if data.get('novel') is False)
    unknown_count = sum(1 for data in results.values() if data.get('novel') is None)
    logger.info(f"Final counts - Novel: {novel_count}, Known: {known_count}, Unknown: {unknown_count}")
    
    # Generate and save the report
    report_df = generate_report(results, args.output_dir)
    
    # Print a brief summary
    novel_count = sum(1 for data in results.values() if data['novel'])
    logger.info(f"\nAnalysis complete: {novel_count} novel sequences out of {len(sequences)}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.exception(f"An unexpected error occurred: {e}")
        sys.exit(1)
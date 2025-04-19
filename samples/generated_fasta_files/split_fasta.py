import os
import argparse
from pathlib import Path

def parse_fasta(fasta_file):
    """Parse a FASTA file and return a list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_sequence = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Save the previous sequence if it exists
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_sequence)))
                
                # Start a new sequence
                current_header = line[1:]  # Remove the '>' character
                current_sequence = []
            else:
                # Add to the current sequence
                current_sequence.append(line)
    
    # Add the last sequence if it exists
    if current_header is not None and current_sequence:
        sequences.append((current_header, ''.join(current_sequence)))
    
    return sequences

def split_fasta(input_file, output_dir):
    """Split a FASTA file into individual files with one sequence per file."""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse the FASTA file
    sequences = parse_fasta(input_file)
    
    # Write each sequence to its own file
    for header, sequence in sequences:
        # Create a valid filename from the header
        # Replace any invalid filename characters with underscores
        safe_filename = ''.join(c if c.isalnum() or c in '_-.' else '_' for c in header)
        output_file = os.path.join(output_dir, f"{safe_filename}.fasta")
        
        with open(output_file, 'w') as f:
            f.write(f">{header}\n{sequence}\n")
    
    return len(sequences)

def main():
    parser = argparse.ArgumentParser(description='Split a FASTA file into individual files with one sequence per file.')
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('--output_dir', '-o', default='split_output',
                        help='Output directory for individual FASTA files (default: split_output)')
    
    args = parser.parse_args()
    
    # Get absolute paths
    input_path = Path(args.input_file).absolute()
    output_path = Path(args.output_dir).absolute()
    
    if not input_path.exists():
        print(f"Error: Input file '{input_path}' does not exist.")
        return 1
    
    # Split the FASTA file
    num_sequences = split_fasta(input_path, output_path)
    
    print(f"Successfully split {input_path.name} into {num_sequences} individual FASTA files.")
    print(f"Output files saved to: {output_path}")
    
    return 0

if __name__ == "__main__":
    exit(main())

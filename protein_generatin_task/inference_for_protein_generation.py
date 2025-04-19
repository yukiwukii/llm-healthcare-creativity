import ollama
import os
from tqdm import tqdm

# List of models to use
models = ['deepseek-r1:8b']

# Number of proteins to generate
num_proteins = 100

# Protein generation prompt with OmpT binding specificity
prompt = '''Design a novel protein sequence specifically designed to bind to OmpT (Outer Membrane Protease T), a bacterial outer membrane protease from E. coli.

Important considerations for OmpT binding:
1. OmpT is a 35.5 kDa aspartyl protease that cleaves between basic amino acids (Arg-Arg, Lys-Lys, Arg-Lys, Lys-Arg)
2. OmpT has a negative electrostatic potential around its active site (acidic residues)
3. Target positively charged residues (Arg, Lys) at the binding interface to interact with OmpT's negatively charged surface
4. Include hydrophobic residues for interacting with OmpT's membrane-embedded regions
5. Consider beta-sheet structures as OmpT is a beta-barrel protein
6. Aim for a sequence of 10-50 amino acids, which is standard for antimicrobial proteins

The sequence should only use the 20 standard amino acids (A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V).
Format the output as a valid protein sequence using the one-letter amino acid codes with no spaces or other characters. Make sure that the sequence generated is novel, and has not been discovered before.
Do not include any explanations, just return the sequence.

Here are some examples of protein sequences that strongly binds to OmpT [RQIKIWFQNRRMKWKK, RQIKIWFQWRRWKWKK, RWIKIQFQIRRWKNKK. These are all mutations of the protein penetratin.
'''

# Function to create FASTA format
def create_fasta_entry(id, sequence):
    return f">{id}\n{sequence}\n"

# Process each model
for model in models:
    print(f"Generating OmpT-binding proteins with {model}...")
    model_name = model.split(':')[0]
    fasta_file = f"{model_name}_ompT_binding_proteins.fasta"
    
    # Create or clear the output file
    with open(fasta_file, 'w') as f:
        f.write("")
    
    # Generate proteins
    for i in tqdm(range(1, num_proteins + 1)):
        try:
            # Call the Ollama chat model
            response = ollama.chat(
                model=model,
                messages=[
                    {"role": "system", "content": prompt},
                    {"role": "user", "content": f"Generate OmpT-binding protein number {i}"}
                ]
            )
            
            # Extract protein sequence, cleaning any non-amino acid characters
            protein_seq = response['message']['content'].strip()
            
            # Basic validation (keep only valid amino acid letters)
            valid_aa = "ARNDCEQGHILKMFPSTWYV"
            protein_seq = ''.join([aa for aa in protein_seq if aa in valid_aa])
            
            # Skip if sequence is too short
            if len(protein_seq) < 10:
                print(f"Generated sequence too short ({len(protein_seq)}), retrying...")
                i -= 1
                continue
            elif len(protein_seq) > 50:
                print(f"Generated sequence too long ({len(protein_seq)}), retrying...")
                i -= 1
                continue
                
            # Create FASTA entry and append to file
            protein_id = f"{model_name}_ompT_binder_{i}"
            
            with open(fasta_file, 'a') as f:
                f.write(create_fasta_entry(protein_id, protein_seq))
                
        except Exception as e:
            print(f"Error with {model} on protein {i}: {str(e)}")
            # Try again
            i -= 1
    
    print(f"Completed {model}. Output saved to {fasta_file}")

print("All OmpT-binding protein generation complete.")

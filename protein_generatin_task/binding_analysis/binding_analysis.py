#!/usr/bin/env python3
"""
Modified version of the binding energy and pLDDT analysis script
with focus on negative binding energy, fixed pLDDT threshold,
and terminology change from 'binder' to 'sequence'.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import os
from scipy.stats import spearmanr, pearsonr  # Added import for correlation coefficients

# Define key configuration parameters at the top of the script
# CHANGE: Make total_sequences a variable that can be easily changed
total_sequences = 100  # Change this value to adjust the expected total number of sequences
plddt_threshold = 70   # Fixed pLDDT threshold
MODEL = 'DEEPSEEK'
binding_energy_threshold = -13

# Create output directory
output_dir = "binding_analysis"
os.makedirs(output_dir, exist_ok=True)

# 1. Load and prepare docking data
print("Loading docking data...")
docking_df = pd.read_csv(f"{MODEL}/docking_summary_fixed.tsv", sep='\t')

# Check if we have a duplicate header row and remove it
if "Receptor" in str(docking_df.iloc[0].values):
    docking_df = docking_df.iloc[1:].reset_index(drop=True)

# Convert Score to numeric
docking_df['Score'] = pd.to_numeric(docking_df['Score'])

# 2. Extract sequence ID and rank from Ligand names
print("\nExtracting sequence ID and rank from model names...")
sequence_ids = []
ranks = []
for model_name in docking_df['Ligand']:
    # Extract sequence ID
    seq_match = re.search(r'binder_(\d+)', model_name)
    seq_id = int(seq_match.group(1)) if seq_match else None
    sequence_ids.append(seq_id)
    
    # Extract rank
    rank_match = re.search(r'rank_00(\d)', model_name)
    rank = int(rank_match.group(1)) if rank_match else None
    ranks.append(rank)

docking_df['sequence_id'] = sequence_ids
docking_df['docking_rank'] = ranks

print(f"Sample docking data with extracted info:")
print(docking_df[['Ligand', 'sequence_id', 'docking_rank', 'Score']].head())

# 3. Find models with negative binding energy
negative_models = docking_df[docking_df['Score'] <= 0].copy()
negative_models = negative_models.sort_values('Score')

print(f"\nFound {len(negative_models)} models with negative binding energy:")
print(negative_models[['Ligand', 'sequence_id', 'docking_rank', 'Score']].head(10))

# Save to CSV
negative_models.to_csv(f"{output_dir}/negative_energy_models.csv", index=False)

# 4. Load pLDDT data
print("\nLoading pLDDT data...")
plddt_df = pd.read_csv(f"{MODEL}/log.csv")
print(f"pLDDT file columns: {plddt_df.columns.tolist()}")
print(f"Sample pLDDT data:\n{plddt_df.head()}")

# 5. Extract model information from pLDDT data
print("\nExtracting model information from pLDDT data...")
plddt_ranks = []
for model_name in plddt_df['model']:
    # Extract rank 
    if 'model_' in model_name:
        rank_match = re.search(r'model_(\d+)', model_name)
        rank = int(rank_match.group(1)) if rank_match else None
        plddt_ranks.append(rank)
    else:
        plddt_ranks.append(None)

plddt_df['plddt_rank'] = plddt_ranks

# Extract sequence ID from sequence_id
sequence_ids = []
for seq_id in plddt_df['sequence_id']:
    seq_match = re.search(r'binder_(\d+)', seq_id)
    seq_id_num = int(seq_match.group(1)) if seq_match else None
    sequence_ids.append(seq_id_num)

plddt_df['sequence_id_num'] = sequence_ids

print(f"Sample pLDDT data with extracted info:")
print(plddt_df[['sequence_id', 'sequence_id_num', 'model', 'plddt_rank', 'pLDDT']].head())

# 6. Find high pLDDT models
# Set fixed pLDDT threshold as defined at the top of the script
high_plddt = plddt_df[plddt_df['pLDDT'] >= plddt_threshold].copy()
high_plddt = high_plddt.sort_values('pLDDT', ascending=False)

unique_high_plddt_sequences = high_plddt['sequence_id_num'].nunique()
print(f"\nFound {len(high_plddt)} models with high pLDDT (≥ {plddt_threshold:.2f}):")
print(f"These models represent {unique_high_plddt_sequences} unique sequences out of {total_sequences} total sequences")
print(high_plddt[['sequence_id', 'model', 'pLDDT']].head(10))

# Save to CSV
high_plddt.to_csv(f"{output_dir}/high_plddt_models.csv", index=False)

# 7. Match negative energy models with their pLDDT scores
print("\nMatching negative energy models with pLDDT scores...")

# Create a map from sequence_id and model rank to pLDDT score
plddt_map = {}

# Group pLDDT data by sequence_id and rank
for _, row in plddt_df.iterrows():
    if row['sequence_id_num'] is not None and row['plddt_rank'] is not None:
        key = (row['sequence_id_num'], row['plddt_rank'])
        plddt_map[key] = row['pLDDT']

# Match negative energy models with pLDDT scores
negative_plddt = []
for _, row in negative_models.iterrows():
    sequence_id = row['sequence_id']
    
    # For matching the rank, we need to check if docking_rank matches plddt_rank
    # or if we should try other matching approaches
    docking_rank = row['docking_rank']
    
    # First, try direct match (sequence_id, docking_rank)
    if (sequence_id, docking_rank) in plddt_map:
        plddt = plddt_map[(sequence_id, docking_rank)]
    else:
        # Try alternative rankings
        plddt = None
        for rank in range(1, 6):  # Check all 5 ranks
            if (sequence_id, rank) in plddt_map:
                plddt = plddt_map[(sequence_id, rank)]
                break
    
    negative_plddt.append({
        'Ligand': row['Ligand'],
        'Sequence_ID': sequence_id,
        'Docking_Rank': docking_rank,
        'Binding_Energy': row['Score'],
        'pLDDT': plddt
    })

# Create DataFrame
negative_plddt_df = pd.DataFrame(negative_plddt)
    
print("\nNegative energy models with pLDDT values:")
print(negative_plddt_df[['Sequence_ID', 'Docking_Rank', 'Binding_Energy', 'pLDDT']].head(10))

# Save to CSV
negative_plddt_df.to_csv(f"{output_dir}/negative_energy_with_plddt.csv", index=False)

# 8. Find high pLDDT models with their binding energy
print("\nFinding binding energy for high pLDDT models...")

# Create a map from sequence_id and model rank to binding energy
energy_map = {}

# Group docking data by sequence_id and rank
for _, row in docking_df.iterrows():
    if row['sequence_id'] is not None and row['docking_rank'] is not None:
        key = (row['sequence_id'], row['docking_rank'])
        energy_map[key] = row['Score']

# Match high pLDDT models with binding energy
high_plddt_energy = []
for _, row in high_plddt.iterrows():
    sequence_id = row['sequence_id_num']
    
    # For matching the rank, try a few approaches
    plddt_rank = row['plddt_rank']
    
    # First, try direct match (sequence_id, plddt_rank)
    if (sequence_id, plddt_rank) in energy_map:
        energy = energy_map[(sequence_id, plddt_rank)]
    else:
        # Try alternative rankings
        energy = None
        for rank in range(1, 6):  # Check all 5 ranks
            if (sequence_id, rank) in energy_map:
                energy = energy_map[(sequence_id, rank)]
                break
    
    high_plddt_energy.append({
        'Sequence_ID_Full': row['sequence_id'],
        'Sequence_ID': sequence_id,
        'Model': row['model'],
        'PLDDT_Rank': plddt_rank,
        'pLDDT': row['pLDDT'],
        'Binding_Energy': energy
    })

# Create DataFrame
high_plddt_energy_df = pd.DataFrame(high_plddt_energy)

print("\nHigh pLDDT models with binding energy:")
print(high_plddt_energy_df[['Sequence_ID', 'PLDDT_Rank', 'pLDDT', 'Binding_Energy']].head(10))

# Save to CSV
high_plddt_energy_df.to_csv(f"{output_dir}/high_plddt_with_energy.csv", index=False)

# 9. Find optimal models (negative energy AND high pLDDT)
print("\nFinding optimal models (negative energy AND high pLDDT)...")

optimal_models = negative_plddt_df[negative_plddt_df['pLDDT'] >= plddt_threshold].copy()
optimal_models = optimal_models.sort_values('Binding_Energy')

# Count unique sequences that have both negative energy and high pLDDT
unique_optimal_sequences = optimal_models['Sequence_ID'].nunique()
optimal_accuracy = (unique_optimal_sequences / total_sequences) * 100

if len(optimal_models) > 0:
    print(f"Found {len(optimal_models)} optimal models from {unique_optimal_sequences} unique sequences out of {total_sequences} total:")
    print(optimal_models[['Sequence_ID', 'Docking_Rank', 'Binding_Energy', 'pLDDT']].to_string(index=False))
    
    # Save to CSV
    optimal_models.to_csv(f"{output_dir}/optimal_models.csv", index=False)
else:
    print("No optimal models found (no overlap between negative energy and high pLDDT).")

# 10. Create visualization - Distribution of binding energies
# CHANGE: Focus on negative binding energy and use binding_energy_threshold
print("\nCreating visualizations...")

plt.figure(figsize=(10, 6))
# Only use negative energies for the histogram
neg_energy_df = docking_df[docking_df['Score'] < 0]

# Count unique sequences with negative binding energy
unique_seqs_with_neg_energy = neg_energy_df['sequence_id'].nunique()
seq_percentage = (unique_seqs_with_neg_energy / total_sequences) * 100

# Count unique sequences with binding energy more negative or equal to binding_energy_threshold
threshold_energy_df = docking_df[docking_df['Score'] <= binding_energy_threshold]
unique_seqs_with_threshold_energy = threshold_energy_df['sequence_id'].nunique()
threshold_seq_percentage = (unique_seqs_with_threshold_energy / total_sequences) * 100

sns.histplot(data=neg_energy_df, x='Score', bins=30, kde=True)
plt.axvline(x=binding_energy_threshold, color='red', linestyle='--', 
           label=f'Binding Energy Threshold ({binding_energy_threshold})')

# Add count and percentage to legend
plt.plot([], [], ' ', label=f'Sequences with negative energy: {unique_seqs_with_neg_energy} of {total_sequences}')
plt.plot([], [], ' ', label=f'Sequences with energy ≤ {binding_energy_threshold}: {unique_seqs_with_threshold_energy} of {total_sequences}')

plt.title("Distribution of Negative Binding Energies", fontsize=14)
plt.xlabel("Binding Energy", fontsize=12)
plt.ylabel("Count of Models", fontsize=12)
plt.grid(alpha=0.3)
plt.legend()
plt.savefig(f"{output_dir}/negative_binding_energy_distribution.png", dpi=300, bbox_inches='tight')

# 11. Create visualization - Distribution of pLDDT scores
# CHANGE: Keep mean but remove std dev, add percentage of unique sequences with high pLDDT
plt.figure(figsize=(10, 6))
sns.histplot(data=plddt_df, x='pLDDT', bins=30, kde=True)

# Calculate and add mean to the legend
plddt_mean = plddt_df['pLDDT'].mean()
plt.axvline(x=plddt_mean, color='green', linestyle='-', 
           label=f'Mean pLDDT = {plddt_mean:.2f}')

# Add threshold line
plt.axvline(x=plddt_threshold, color='blue', linestyle='--', 
           label=f'High pLDDT threshold ({plddt_threshold:.2f})')

# Verify our counts with extensive diagnostics
print("\n--- DIAGNOSTIC INFORMATION ---")
# Count total models and high pLDDT models
total_models = len(plddt_df)
high_plddt_models = len(plddt_df[plddt_df['pLDDT'] >= plddt_threshold])
high_plddt_model_percentage = (high_plddt_models / total_models) * 100
print(f"Total models: {total_models}")
print(f"Models with pLDDT >= {plddt_threshold}: {high_plddt_models}")
print(f"Percentage of models with high pLDDT: {high_plddt_model_percentage:.1f}%")

# Check for NaN values in sequence_id_num
null_count = plddt_df['sequence_id_num'].isnull().sum()
print(f"Number of rows with null sequence_id_num: {null_count}")

# Count unique sequences in entire dataset
total_sequences = plddt_df['sequence_id_num'].nunique()
print(f"Total unique sequences detected: {total_sequences}")

# Count unique sequences with high pLDDT
high_plddt_rows = plddt_df[plddt_df['pLDDT'] >= plddt_threshold]
high_plddt_sequences = high_plddt_rows['sequence_id_num'].nunique()
high_plddt_percentage = (high_plddt_sequences / total_sequences) * 100
print(f"Unique sequences with high pLDDT: {high_plddt_sequences} out of {total_sequences}")
print(f"Percentage of sequences with high pLDDT: {high_plddt_percentage:.1f}%")

# Distribution of models per sequence
models_per_seq = plddt_df.groupby('sequence_id_num').size()
print(f"Average models per sequence: {models_per_seq.mean():.1f}")
print(f"Min models per sequence: {models_per_seq.min()}")
print(f"Max models per sequence: {models_per_seq.max()}")

# Add info to legend
plt.plot([], [], ' ', label=f'Sequences with high pLDDT: {high_plddt_sequences} of {total_sequences}')
plt.plot([], [], ' ', label=f'Percentage: {high_plddt_percentage:.1f}%')

plt.title("Distribution of pLDDT Scores", fontsize=14)
plt.xlabel("pLDDT Score", fontsize=12)
plt.ylabel("Count of Models", fontsize=12)
plt.grid(alpha=0.3)
plt.legend()
plt.savefig(f"{output_dir}/plddt_distribution.png", dpi=300, bbox_inches='tight')

# 12. REMOVED: Individual model scatter plot (as requested, focusing only on sequence-level data)

# 13. Create sequence-level analysis
print("\nCreating sequence-level analysis...")

# Group by sequence_id and find best binding energy
sequence_stats = docking_df.groupby('sequence_id')['Score'].agg(['min', 'mean']).reset_index()
sequence_stats.columns = ['Sequence_ID', 'Best_Energy', 'Avg_Energy']

# Count negative energy models per sequence
neg_counts = docking_df[docking_df['Score'] < 0].groupby('sequence_id').size().reset_index()
neg_counts.columns = ['Sequence_ID', 'Negative_Energy_Count']

# Merge with sequence_stats
sequence_stats = pd.merge(sequence_stats, neg_counts, on='Sequence_ID', how='left')
sequence_stats['Negative_Energy_Count'] = sequence_stats['Negative_Energy_Count'].fillna(0).astype(int)

# Find best pLDDT score per sequence
plddt_stats = plddt_df.groupby('sequence_id_num')['pLDDT'].max().reset_index()
plddt_stats.columns = ['Sequence_ID', 'Best_pLDDT']

# Merge with sequence_stats
sequence_stats = pd.merge(sequence_stats, plddt_stats, on='Sequence_ID', how='left')

# Sort by best energy
sequence_stats = sequence_stats.sort_values('Best_Energy')

print("\nSequence performance (top 15):")
print(sequence_stats.head(15).to_string(index=False))

# Save to CSV
sequence_stats.to_csv(f"{output_dir}/sequence_performance.csv", index=False)

# Create visualization - Best energy by sequence ID
plt.figure(figsize=(12, 6))
top_sequences = sequence_stats.head(15)  # Show top 15 sequences
plt.bar(
    range(len(top_sequences)), 
    top_sequences['Best_Energy'], 
    color=['green' if x < 0 else 'red' for x in top_sequences['Best_Energy']]
)
plt.axhline(y=0, color='black', linestyle='--')
plt.xticks(range(len(top_sequences)), top_sequences['Sequence_ID'], rotation=45)
plt.title("Best Binding Energy by Sequence ID (Top 15)", fontsize=14)
plt.xlabel("Sequence ID", fontsize=12)
plt.ylabel("Best Binding Energy", fontsize=12)
plt.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig(f"{output_dir}/best_energy_by_sequence.png", dpi=300, bbox_inches='tight')

# 14. Create visualization - Top sequences with their best pLDDT scores and best binding energy
plt.figure(figsize=(12, 8))

# Get sequences with negative energy
neg_energy_sequences = sequence_stats[sequence_stats['Best_Energy'] < 0].copy()

# Count the unique sequences with negative binding energy (regardless of pLDDT)
unique_sequences_neg_energy = len(neg_energy_sequences)

# Calculate correlation coefficients between Best_pLDDT and Best_Energy
# Using only sequences with negative energy for the correlation
spearman_corr, spearman_p = spearmanr(neg_energy_sequences['Best_pLDDT'], neg_energy_sequences['Best_Energy'])
pearson_corr, pearson_p = pearsonr(neg_energy_sequences['Best_pLDDT'], neg_energy_sequences['Best_Energy'])

# Create scatter plot
plt.scatter(
    neg_energy_sequences['Best_pLDDT'],
    neg_energy_sequences['Best_Energy'],
    s=neg_energy_sequences['Negative_Energy_Count'] * 50 + 30,
    alpha=0.7,
    label=f'Sequences with negative energy: {unique_sequences_neg_energy}%'
)

# Add labels for all sequences
for _, row in neg_energy_sequences.iterrows():
    plt.annotate(
        f"Sequence {row['Sequence_ID']}",
        (row['Best_pLDDT'], row['Best_Energy']),
        xytext=(5, 0),
        textcoords='offset points',
        fontsize=10,
        bbox=dict(boxstyle='round,pad=0.3', fc='yellow', alpha=0.5)
    )

# Add reference lines
plt.axhline(y=0, color='red', linestyle='--')
plt.axvline(x=plddt_threshold, color='blue', linestyle='--', 
           label=f'High pLDDT threshold ({plddt_threshold:.2f})')

# Add correlation coefficients to the legend
plt.plot([], [], ' ', label=f'Spearman correlation: {spearman_corr:.3f} (p={spearman_p:.3f})')
plt.plot([], [], ' ', label=f'Pearson correlation: {pearson_corr:.3f} (p={pearson_p:.3f})')
# Keep only the accuracy percentage, remove the count of optimal sequences
plt.plot([], [], ' ', label=f'Optimal Model Accuracy: {optimal_accuracy:.1f}%')

plt.xlabel('Best pLDDT Score', fontsize=14)
plt.ylabel('Best Binding Energy', fontsize=14)
plt.title('Best Binding Energy vs. Best pLDDT by Sequence ID', fontsize=16)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=11)

plt.savefig(f"{output_dir}/negative_energy_sequences_plddt.png", dpi=300, bbox_inches='tight')

print(f"\nAnalysis complete. Results saved to {output_dir}/")
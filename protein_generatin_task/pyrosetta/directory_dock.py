# directory_docking.py
import os
import sys
import traceback
import argparse
import time
import json
import glob
from datetime import datetime

def checkpoint(msg, logfile=None):
    """Print checkpoint message and flush to ensure it's written"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    message = f"[{timestamp}] {msg}"
    print(message)
    sys.stdout.flush()
    if logfile:
        with open(logfile, 'a') as f:
            f.write(message + "\n")
            f.flush()

def dock_proteins(receptor_path, ligand_path, output_dir="docking_results", log_path=None, skip_if_exists=True):
    """
    Perform Monte Carlo docking between receptor and ligand
    
    Parameters:
    receptor_path (str): Path to receptor PDB file
    ligand_path (str): Path to ligand PDB file
    output_dir (str): Directory to save results
    log_path (str): Path to log file
    skip_if_exists (bool): Skip docking if output files already exist
    
    Returns:
    dict: Dictionary with docking results
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # Get base names for files
    receptor_name = os.path.basename(receptor_path).split('.')[0]
    ligand_name = os.path.basename(ligand_path).split('.')[0]
    complex_name = f"{receptor_name}_{ligand_name}"
    
    # Define output paths
    docked_path = os.path.join(output_dir, f"{complex_name}_docked.pdb")
    report_path = os.path.join(output_dir, f"{complex_name}_report.txt")
    
    # Check if docking was already completed
    if skip_if_exists and os.path.exists(docked_path) and os.path.exists(report_path):
        print(f"Docking results for {complex_name} already exist. Skipping.")
        
        # Try to read binding energy from existing report
        try:
            with open(report_path, 'r') as f:
                for line in f:
                    if "Binding energy:" in line:
                        binding_energy = float(line.split(':')[1].split()[0])
                        break
            
            return {
                "receptor": receptor_path,
                "ligand": ligand_path,
                "complex_name": complex_name,
                "success": True,
                "binding_energy": binding_energy,
                "docked_pdb": docked_path,
                "report_path": report_path,
                "skipped": True
            }
        except Exception:
            pass  # If we can't read the report, just continue with docking
    
    # If no log path provided, create one in the output directory
    if not log_path:
        log_path = os.path.join(output_dir, f"{complex_name}_docking.log")
    
    # Clear log file if it exists
    open(log_path, 'w').close()
    
    # Dictionary to store results
    results = {
        "receptor": receptor_path,
        "ligand": ligand_path,
        "complex_name": complex_name,
        "success": False
    }
    
    try:
        checkpoint(f"Starting docking of {receptor_name} with {ligand_name}", log_path)
        checkpoint("Importing PyRosetta modules", log_path)
        import pyrosetta
        from pyrosetta import Pose
        from pyrosetta.rosetta.protocols.docking import DockMCMProtocol
        from pyrosetta.rosetta.core.scoring import ScoreFunction, ScoreType
        
        checkpoint("Initializing PyRosetta", log_path)
        pyrosetta.init(options="-mute basic -mute core -ex1 -ex2aro")
        
        # Check if files exist
        if not os.path.exists(receptor_path):
            raise FileNotFoundError(f"Receptor file not found: {receptor_path}")
        if not os.path.exists(ligand_path):
            raise FileNotFoundError(f"Ligand file not found: {ligand_path}")
        
        checkpoint("Loading receptor and ligand", log_path)
        receptor = pyrosetta.pose_from_pdb(receptor_path)
        ligand = pyrosetta.pose_from_pdb(ligand_path)
        checkpoint(f"Loaded receptor ({receptor.total_residue()} residues) and ligand ({ligand.total_residue()} residues)", log_path)
        
        # Create a score function
        checkpoint("Creating score function", log_path)
        scorefxn = ScoreFunction()
        # Add standard scoring terms
        scorefxn.set_weight(ScoreType.fa_atr, 1.0)
        scorefxn.set_weight(ScoreType.fa_rep, 0.55)
        scorefxn.set_weight(ScoreType.fa_sol, 0.9)
        scorefxn.set_weight(ScoreType.fa_elec, 0.5)
        scorefxn.set_weight(ScoreType.hbond_sc, 1.0)
        scorefxn.set_weight(ScoreType.hbond_bb_sc, 1.0)
        checkpoint("Score function created", log_path)
        
        # Create merged PDB file for manual complex creation
        checkpoint("Creating merged PDB file", log_path)
        # Place ligand at an appropriate distance from receptor
        # Calculate centroids
        rec_atoms = 0
        rec_com_x, rec_com_y, rec_com_z = 0, 0, 0
        
        for i in range(1, receptor.total_residue() + 1):
            for j in range(1, receptor.residue(i).natoms() + 1):
                atom_xyz = receptor.residue(i).atom(j).xyz()
                rec_com_x += atom_xyz[0]
                rec_com_y += atom_xyz[1]
                rec_com_z += atom_xyz[2]
                rec_atoms += 1
        
        rec_com_x /= rec_atoms
        rec_com_y /= rec_atoms
        rec_com_z /= rec_atoms
        
        lig_atoms = 0
        lig_com_x, lig_com_y, lig_com_z = 0, 0, 0
        
        for i in range(1, ligand.total_residue() + 1):
            for j in range(1, ligand.residue(i).natoms() + 1):
                atom_xyz = ligand.residue(i).atom(j).xyz()
                lig_com_x += atom_xyz[0]
                lig_com_y += atom_xyz[1]
                lig_com_z += atom_xyz[2]
                lig_atoms += 1
        
        lig_com_x /= lig_atoms
        lig_com_y /= lig_atoms
        lig_com_z /= lig_atoms
        
        # Calculate translation to place ligand 25Å away from receptor
        trans_x = rec_com_x + 25.0 - lig_com_x
        trans_y = rec_com_y - lig_com_y
        trans_z = rec_com_z - lig_com_z
        
        # Create a transformed ligand
        ligand_transformed = ligand.clone()
        for i in range(1, ligand_transformed.total_residue() + 1):
            for j in range(1, ligand_transformed.residue(i).natoms() + 1):
                old_xyz = ligand_transformed.residue(i).atom(j).xyz()
                # Update coordinates
                from pyrosetta.rosetta.numeric import xyzVector_double_t
                new_xyz = xyzVector_double_t(old_xyz[0] + trans_x, 
                                             old_xyz[1] + trans_y,
                                             old_xyz[2] + trans_z)
                ligand_transformed.residue(i).atom(j).xyz(new_xyz)
        
        # Use unique file paths in the output directory
        receptor_prep_path = os.path.join(output_dir, f"{receptor_name}_prep.pdb")
        ligand_prep_path = os.path.join(output_dir, f"{ligand_name}_prep.pdb")
        merged_path = os.path.join(output_dir, f"{complex_name}_merged.pdb")
        
        # Save receptor and transformed ligand to separate PDBs
        receptor.dump_pdb(receptor_prep_path)
        ligand_transformed.dump_pdb(ligand_prep_path)
        
        # Create merged PDB with labeled chains
        with open(merged_path, "w") as outfile:
            with open(receptor_prep_path, "r") as infile:
                for line in infile:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        # Ensure chain ID is set to A
                        outfile.write(line[:21] + "A" + line[22:])
                        
            with open(ligand_prep_path, "r") as infile:
                for line in infile:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        # Set chain ID to B for ligand
                        outfile.write(line[:21] + "B" + line[22:])
                        
            outfile.write("END\n")
        
        # Load the merged structure
        checkpoint("Loading merged structure", log_path)
        docking_pose = pyrosetta.pose_from_pdb(merged_path)
        checkpoint(f"Merged structure loaded with {docking_pose.total_residue()} residues", log_path)
        
        # Score initial complex
        initial_score = scorefxn(docking_pose)
        checkpoint(f"Initial complex score: {initial_score:.2f} REU", log_path)
        
        # Perform Monte Carlo docking
        checkpoint("Setting up Monte Carlo docking protocol", log_path)
        dock_protocol = DockMCMProtocol()
        dock_protocol.set_scorefxn(scorefxn)
        
        checkpoint("Running Monte Carlo docking - this may take some time", log_path)
        dock_protocol.apply(docking_pose)
        checkpoint("Monte Carlo docking completed", log_path)
        
        # Score final complex
        final_score = scorefxn(docking_pose)
        checkpoint(f"Final complex score: {final_score:.2f} REU", log_path)
        checkpoint(f"Score improvement: {initial_score - final_score:.2f} REU", log_path)
        
        # Save the docked structure
        docking_pose.dump_pdb(docked_path)
        checkpoint(f"Saved docked structure as {docked_path}", log_path)
        
        # Calculate binding energy
        receptor = pyrosetta.pose_from_pdb(receptor_prep_path)
        ligand = pyrosetta.pose_from_pdb(ligand_prep_path)
        
        rec_score = scorefxn(receptor)
        lig_score = scorefxn(ligand)
        binding_energy = final_score - (rec_score + lig_score)
        
        checkpoint(f"Receptor score: {rec_score:.2f} REU", log_path)
        checkpoint(f"Ligand score: {lig_score:.2f} REU", log_path)
        checkpoint(f"Binding energy: {binding_energy:.2f} REU", log_path)
        
        # Analyze interface
        checkpoint("Analyzing interface", log_path)
        
        # Count interface contacts
        interface_contacts = []
        for i in range(1, docking_pose.total_residue() + 1):
            chain_i = docking_pose.pdb_info().chain(i)
            if chain_i != "A":
                continue
            
            for j in range(1, docking_pose.total_residue() + 1):
                chain_j = docking_pose.pdb_info().chain(j)
                if chain_j != "B":
                    continue
                    
                # Check distance between residues
                min_dist = float('inf')
                for iatom in range(1, docking_pose.residue(i).natoms() + 1):
                    for jatom in range(1, docking_pose.residue(j).natoms() + 1):
                        dist = docking_pose.residue(i).atom(iatom).xyz().distance(
                               docking_pose.residue(j).atom(jatom).xyz())
                        min_dist = min(min_dist, dist)
                
                if min_dist < 8.0:  # 8Å cutoff
                    interface_contacts.append((i, j, min_dist))
        
        unique_rec_residues = set(c[0] for c in interface_contacts)
        unique_lig_residues = set(c[1] for c in interface_contacts)
        
        checkpoint(f"Interface contacts: {len(interface_contacts)}", log_path)
        checkpoint(f"Receptor residues in interface: {len(unique_rec_residues)}", log_path)
        checkpoint(f"Ligand residues in interface: {len(unique_lig_residues)}", log_path)
        
        # Create a detailed report
        with open(report_path, "w") as report:
            report.write("Monte Carlo Docking Results\n")
            report.write("==========================\n\n")
            report.write(f"Receptor: {receptor_path}\n")
            report.write(f"Ligand: {ligand_path}\n\n")
            report.write(f"Initial score: {initial_score:.2f} REU\n")
            report.write(f"Final score: {final_score:.2f} REU\n")
            report.write(f"Score improvement: {initial_score - final_score:.2f} REU\n\n")
            report.write(f"Receptor score: {rec_score:.2f} REU\n")
            report.write(f"Ligand score: {lig_score:.2f} REU\n")
            report.write(f"Binding energy: {binding_energy:.2f} REU\n\n")
            
            # Binding energy interpretation
            report.write("Binding Energy Interpretation:\n")
            if binding_energy < -10.0:
                report.write("Very strong binding (< -10 REU)\n")
            elif binding_energy < -5.0:
                report.write("Strong binding (-5 to -10 REU)\n")
            elif binding_energy < 0:
                report.write("Moderate binding (0 to -5 REU)\n")
            elif binding_energy < 10:
                report.write("Weak binding (0 to 10 REU)\n")
            else:
                report.write("Unfavorable binding (> 10 REU)\n")
            
            report.write("\nInterface Analysis:\n")
            report.write(f"Number of interface contacts: {len(interface_contacts)}\n")
            report.write(f"Receptor residues in interface: {len(unique_rec_residues)}\n")
            report.write(f"Ligand residues in interface: {len(unique_lig_residues)}\n\n")
            
            # List interface residues
            report.write("Interface Residues:\n")
            report.write("Receptor: ")
            report.write(", ".join([f"{docking_pose.pdb_info().number(r)}" for r in sorted(unique_rec_residues)]))
            report.write("\n")
            report.write("Ligand: ")
            report.write(", ".join([f"{docking_pose.pdb_info().number(r)}" for r in sorted(unique_lig_residues)]))
            report.write("\n\n")
            
            # List closest contacts
            report.write("Closest Contacts:\n")
            sorted_contacts = sorted(interface_contacts, key=lambda x: x[2])
            for i, (rec_res, lig_res, dist) in enumerate(sorted_contacts[:10]):  # Top 10 contacts
                rec_num = docking_pose.pdb_info().number(rec_res)
                lig_num = docking_pose.pdb_info().number(lig_res)
                rec_name = docking_pose.residue(rec_res).name()
                lig_name = docking_pose.residue(lig_res).name()
                report.write(f"{i+1}. {rec_name} {rec_num} - {lig_name} {lig_num}: {dist:.2f} Å\n")
        
        checkpoint(f"Docking report written to {report_path}", log_path)
        checkpoint("Docking analysis complete", log_path)
        
        # Store results in dictionary
        results["success"] = True
        results["initial_score"] = initial_score
        results["final_score"] = final_score
        results["score_improvement"] = initial_score - final_score
        results["receptor_score"] = rec_score
        results["ligand_score"] = lig_score
        results["binding_energy"] = binding_energy
        results["interface_contacts"] = len(interface_contacts)
        results["receptor_interface_residues"] = len(unique_rec_residues)
        results["ligand_interface_residues"] = len(unique_lig_residues)
        results["docked_pdb"] = docked_path
        results["report_path"] = report_path
        
    except Exception as e:
        error_msg = f"ERROR: {e}"
        checkpoint(error_msg, log_path)
        traceback.print_exc()
        
        # Store error information
        results["error"] = str(e)
        results["traceback"] = traceback.format_exc()
    
    return results

def batch_dock_ligands(receptor_path, ligand_paths, output_dir="docking_results", start_idx=0, end_idx=None, 
                       skip_completed=True, progress_file=None):
    """
    Dock multiple ligands against the same receptor
    
    Parameters:
    receptor_path (str): Path to receptor PDB file
    ligand_paths (list): List of paths to ligand PDB files
    output_dir (str): Directory to save results
    start_idx (int): Index to start docking from (for resuming)
    end_idx (int): Index to end docking at (for batching)
    skip_completed (bool): Skip ligands that have already been docked
    progress_file (str): Path to file tracking docking progress
    
    Returns:
    dict: Dictionary with results for each ligand
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Limit the range if specified
    if end_idx is None:
        end_idx = len(ligand_paths)
    
    ligand_paths = ligand_paths[start_idx:end_idx]
    total_ligands = len(ligand_paths)
    
    # Dictionary to store results
    batch_results = {}
    
    # Create or update progress file
    if not progress_file:
        progress_file = os.path.join(output_dir, "docking_progress.json")
    
    # Initialize or load progress
    if os.path.exists(progress_file):
        with open(progress_file, "r") as f:
            progress_data = json.load(f)
    else:
        progress_data = {
            "receptor": receptor_path,
            "total_ligands": total_ligands,
            "completed": 0,
            "successful": 0,
            "failed": 0,
            "skipped": 0,
            "started_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "last_updated": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "completed_ligands": [],
            "failed_ligands": []
        }
    
    # Summary file
    summary_path = os.path.join(output_dir, "docking_summary.tsv")
    if not os.path.exists(summary_path):
        with open(summary_path, "w") as summary:
            summary.write("Receptor\tLigand\tBinding Energy (REU)\tInterface Contacts\tReceptor Interface Residues\tLigand Interface Residues\tStatus\tTimestamp\n")
    
    # Process each ligand
    for i, ligand_path in enumerate(ligand_paths):
        ligand_name = os.path.basename(ligand_path).split('.')[0]
        current_idx = start_idx + i
        
        print(f"\n=== Docking ligand {current_idx+1}/{start_idx+total_ligands} ({i+1}/{total_ligands}): {ligand_name} ===\n")
        
        # Check if already completed
        if skip_completed and ligand_name in progress_data.get("completed_ligands", []):
            print(f"Ligand {ligand_name} already completed. Skipping.")
            progress_data["skipped"] += 1
            continue
        
        start_time = time.time()
        
        # Perform docking
        result = dock_proteins(receptor_path, ligand_path, output_dir, skip_if_exists=skip_completed)
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        
        # Store result
        batch_results[ligand_name] = result
        
        # Update progress data
        progress_data["last_updated"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        if result.get("skipped", False):
            progress_data["skipped"] += 1
        elif result["success"]:
            progress_data["completed"] += 1
            progress_data["successful"] += 1
            if ligand_name not in progress_data["completed_ligands"]:
                progress_data["completed_ligands"].append(ligand_name)
        else:
            progress_data["completed"] += 1
            progress_data["failed"] += 1
            if ligand_name not in progress_data["failed_ligands"]:
                progress_data["failed_ligands"].append(ligand_name)
        
        # Save progress
        with open(progress_file, "w") as f:
            json.dump(progress_data, f, indent=2)
        
        # Update summary file
        with open(summary_path, "a") as summary:
            receptor_name = os.path.basename(receptor_path).split('.')[0]
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            
            if result.get("skipped", False):
                summary.write(f"{receptor_name}\t{ligand_name}\t{result.get('binding_energy', '-')}\t-\t-\t-\tSkipped\t{timestamp}\n")
            elif result["success"]:
                summary.write(f"{receptor_name}\t{ligand_name}\t{result['binding_energy']:.2f}\t{result['interface_contacts']}\t{result['receptor_interface_residues']}\t{result['ligand_interface_residues']}\tSuccess\t{timestamp}\n")
            else:
                summary.write(f"{receptor_name}\t{ligand_name}\t-\t-\t-\t-\tFailed: {result.get('error', 'Unknown error')}\t{timestamp}\n")
        
        print(f"Docking completed in {elapsed_time:.1f} seconds")
        print(f"Progress: {progress_data['completed']}/{total_ligands} completed ({progress_data['successful']} successful, {progress_data['failed']} failed, {progress_data['skipped']} skipped)")
    
    print(f"\nBatch docking completed. Summary written to {summary_path}")
    return batch_results

def get_ligand_files(ligand_dir, file_pattern="*.pdb"):
    """
    Get all ligand PDB files from a directory
    
    Parameters:
    ligand_dir (str): Directory containing ligand PDB files
    file_pattern (str): Pattern to match ligand files (default: "*.pdb")
    
    Returns:
    list: List of paths to ligand PDB files
    """
    if not os.path.exists(ligand_dir):
        raise FileNotFoundError(f"Ligand directory not found: {ligand_dir}")
    
    # Get all PDB files in the directory
    ligand_paths = glob.glob(os.path.join(ligand_dir, file_pattern))
    
    # Sort for consistent ordering
    ligand_paths.sort()
    
    return ligand_paths

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Dock ligands to a receptor using PyRosetta")
    parser.add_argument("--receptor", required=True, help="Path to receptor PDB file")
    parser.add_argument("--ligand", help="Path to ligand PDB file (for single docking)")
    parser.add_argument("--ligand_dir", help="Directory containing ligand PDB files")
    parser.add_argument("--ligand_list", help="Path to text file with list of ligand PDBs (one per line)")
    parser.add_argument("--output_dir", default="docking_results", help="Directory to save results")
    parser.add_argument("--file_pattern", default="*.pdb", help="Pattern to match ligand files in directory (default: *.pdb)")
    
    # Additional parameters for large-scale docking
    parser.add_argument("--start", type=int, default=0, help="Starting index for batch processing")
    parser.add_argument("--end", type=int, help="Ending index for batch processing")
    parser.add_argument("--no_skip", action="store_true", help="Don't skip already completed ligands")
    parser.add_argument("--job_id", help="SLURM job ID for distributed processing")
    
    args = parser.parse_args()
    
    # Count the options provided - only one source of ligands should be specified
    ligand_sources = sum([bool(args.ligand), bool(args.ligand_list), bool(args.ligand_dir)])
    
    # Validate arguments
    if ligand_sources == 0:
        parser.error("Either --ligand, --ligand_list, or --ligand_dir must be provided")
    
    if ligand_sources > 1:
        parser.error("Only one of --ligand, --ligand_list, or --ligand_dir should be provided")
    
    # Create job-specific output directory if job_id is provided
    output_dir = args.output_dir
    if args.job_id:
        output_dir = os.path.join(args.output_dir, f"job_{args.job_id}")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    if args.ligand:
        # Single docking
        print(f"Starting docking of {args.receptor} with {args.ligand}")
        result = dock_proteins(args.receptor, args.ligand, output_dir)
        
        if result["success"]:
            print(f"\nDocking completed successfully!")
            print(f"Binding energy: {result['binding_energy']:.2f} REU")
            print(f"Interface contacts: {result['interface_contacts']}")
            print(f"Docked structure saved to: {result['docked_pdb']}")
            print(f"Full report: {result['report_path']}")
        else:
            print("\nDocking failed!")
            print(f"Error: {result.get('error', 'Unknown error')}")
    
    elif args.ligand_dir:
        # Directory-based batch docking
        print(f"Starting batch docking with receptor {args.receptor}")
        print(f"Looking for ligand PDBs in directory: {args.ligand_dir}")
        
        # Get all ligand PDB files from directory
        try:
            ligand_paths = get_ligand_files(args.ligand_dir, args.file_pattern)
            print(f"Found {len(ligand_paths)} ligand files matching pattern '{args.file_pattern}'")
            
            if len(ligand_paths) == 0:
                print(f"No ligand files found in {args.ligand_dir} matching pattern '{args.file_pattern}'")
                sys.exit(1)
                
            # Perform batch docking
            batch_results = batch_dock_ligands(
                args.receptor, 
                ligand_paths, 
                output_dir, 
                start_idx=args.start,
                end_idx=args.end,
                skip_completed=not args.no_skip
            )
            
            # Count successes
            success_count = sum(1 for result in batch_results.values() if result["success"])
            skip_count = sum(1 for result in batch_results.values() if result.get("skipped", False))
            
            print(f"\nBatch docking completed: {success_count}/{len(batch_results)} successful, {skip_count} skipped")
            print(f"Results saved to {output_dir}")
            
        except FileNotFoundError as e:
            print(f"Error: {e}")
            sys.exit(1)
    
    else:
        # Text file-based batch docking
        print(f"Starting batch docking with receptor {args.receptor}")
        
        # Read ligand list
        with open(args.ligand_list, "r") as f:
            ligand_paths = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        
        print(f"Found {len(ligand_paths)} ligands in list file")
        print(f"Will process ligands {args.start+1} to {args.end if args.end else len(ligand_paths)}")
        
        # Perform batch docking
        batch_results = batch_dock_ligands(
            args.receptor, 
            ligand_paths, 
            output_dir, 
            start_idx=args.start,
            end_idx=args.end,
            skip_completed=not args.no_skip
        )
        
        # Count successes
        success_count = sum(1 for result in batch_results.values() if result["success"])
        skip_count = sum(1 for result in batch_results.values() if result.get("skipped", False))
        
        print(f"\nBatch docking completed: {success_count}/{len(batch_results)} successful, {skip_count} skipped")
        print(f"Results saved to {output_dir}")

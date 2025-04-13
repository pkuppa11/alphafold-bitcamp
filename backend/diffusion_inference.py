import os
import subprocess

def run_rfdiffusion_binder_design(
    target_pdb,
    binder_length=80,
    binder_chain='B',
    interface_residues='A:100-110',
    output_dir='rfdiff_output',
    rf_script_path='RFdiffusion/scripts/run_inference.py',
    num_designs=5
):
    """
    Run RFdiffusion to design a binder protein to a target site.
    Args:
        target_pdb (str): Path to target PDB structure.
        binder_length (int): Length of the binder to generate.
        binder_chain (str): Chain ID to assign to binder.
        interface_residues (str): Target residues to bind to (e.g. "A:100-110").
        output_dir (str): Where to store output.
        rf_script_path (str): Path to RFdiffusion script.
        num_designs (int): Number of designs to generate.
    """

    os.makedirs(output_dir, exist_ok=True)

    # Format contig string: binder (new chain) binds to a target region (interface)
    contig = f"{binder_chain}:0-{binder_length-1} {interface_residues}"

    # Command to run RFdiffusion
    cmd = [
        "python", rf_script_path,
        "--pdb", target_pdb,
        "--contig", contig,
        "--binder_design",
        "--out_dir", output_dir,
        "--num_designs", str(num_designs),
    ]

    print(f"Running RFdiffusion with:\n{' '.join(cmd)}\n")

    subprocess.run(cmd, check=True)

    print(f"\n✅ Binder design complete! Results in: {output_dir}/")

# Example usage:
run_rfdiffusion_binder_design(
    target_pdb="receptor.pdb",
    binder_length=60,
    interface_residues="A:45-55",
    output_dir="binders_out",
    num_designs=3
)

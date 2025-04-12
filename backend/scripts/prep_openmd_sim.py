import os
import subprocess
import sys
from pathlib import Path
from datetime import datetime

# -----------------------------
# User Configuration Section
# -----------------------------
FORCE_FIELD = "oplsaa"
ENSEMBLE = "NVT"
NSTEPS = 100000
DT = 1.0
TEMPERATURE = 300.0
NMOL = 1
COMPONENT_TYPE = "Polymer"

# -----------------------------
# Helper Functions
# -----------------------------
def run_command(command, check=True):
    print(f"Running: {' '.join(command)}")
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check and result.returncode != 0:
        print(f"Error running command: {result.stderr}")
        sys.exit(1)
    return result.stdout.strip()

def convert_pdb_to_xyz(pdb_file, xyz_file):
    run_command(["obabel", "-ipdb", pdb_file, "-oxyz", "-O", xyz_file])

def generate_omd_file(xyz_file, omd_file, system_name):
    content = f"""################################################
# Auto-generated OpenMD input file
# Created: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
################################################

component {{
    type = {COMPONENT_TYPE}
    nMol = {NMOL}
    structure = "{xyz_file}"
}}

forceField = "{FORCE_FIELD}"
ensemble = {ENSEMBLE}
nSteps = {NSTEPS}
dt = {DT}

temperature = {TEMPERATURE}
"""
    with open(omd_file, 'w') as f:
        f.write(content)
    print(f"Generated .omd file: {omd_file}")

def validate_with_build_atom_types(omd_file):
    print("\n--- Validating Atom Types ---")
    run_command(["buildAtomTypes", "-f", omd_file], check=False)

# -----------------------------
# Main Workflow
# -----------------------------
def main():
    if len(sys.argv) < 2:
        print("Usage: python prepare_openmd_sim.py <input.pdb>")
        sys.exit(1)

    input_path = Path(sys.argv[1])
    if not input_path.exists() or input_path.suffix.lower() != ".pdb":
        print("Error: Input file must be a .pdb file and must exist.")
        sys.exit(1)

    base_name = input_path.stem
    xyz_file = f"{base_name}.xyz"
    omd_file = f"{base_name}.omd"

    print(f"\n--- Starting OpenMD Prep for: {input_path.name} ---")

    # Step 1: Convert to XYZ
    convert_pdb_to_xyz(str(input_path), xyz_file)

    # Step 2: Generate OMD
    generate_omd_file(xyz_file, omd_file, base_name)

    # Step 3: Validate
    validate_with_build_atom_types(omd_file)

    print("\n✅ Done. You can now run your simulation with:")
    print(f"   openmd {omd_file}")

if __name__ == "__main__":
    main()

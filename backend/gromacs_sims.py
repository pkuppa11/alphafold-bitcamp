from Bio.PDB import PDBParser
import subprocess
import MDAnalysis as mda
from MDAnalysis.analysis import rms, polymer


parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", "pdb_outputs/P24941.pdb")

# Optionally, remove water, heteroatoms, etc.
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":  # Not a standard residue
                chain.detach_child(residue.id)

def run_gmx(command, input=None):
    full_cmd = f"gmx {command}"
    print(f"Running: {full_cmd}")
    result = subprocess.run(full_cmd.split(), input=input, text=True, capture_output=True)
    print(result.stdout)
    if result.stderr:
        print("Errors:", result.stderr)
    return result

# Convert PDB to GRO and generate topology
run_gmx("pdb2gmx -f pdb_outputs/P24941.pdb -o processed.gro -water tip3p")

# Create box
run_gmx("editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic")

# Solvate
run_gmx("solvate -cp newbox.gro -cs spc216.gro -o solvated.gro -p topol.top")

# Prepare for ion addition
run_gmx("grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr")

# Add ions
run_gmx("genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral")



# Energy minimization
run_gmx("grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr")
run_gmx("mdrun -deffnm em")

# NVT
run_gmx("grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr")
run_gmx("mdrun -deffnm nvt")

# NPT
run_gmx("grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr")
run_gmx("mdrun -deffnm npt")

# Production MD
run_gmx("grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr")
run_gmx("mdrun -deffnm md")




u = mda.Universe("topol.top", "md.xtc")

rmsd = rms.RMSD(u, select="protein")
rmsd.run()

rg = polymer.RadiusOfGyration(u, select="protein")
rg.run()

from Bio.PDB import PDBParser, PDBIO, is_aa
from Bio.Data import IUPACData

def one_to_three(letter):
    letter = letter.upper()
    if letter not in IUPACData.protein_letters_1to3:
        raise ValueError(f"Invalid amino acid code: '{letter}'")
    return IUPACData.protein_letters_1to3[letter].upper()

def replace_residues(pdb_input_path, pdb_output_path, replacements):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_input_path)

    for model in structure:
        for chain in model:
            for residue in chain:
                res_id = residue.get_id()
                if not is_aa(residue):
                    continue
                res_pos = res_id[1]
                if res_pos in replacements:
                    one_letter = replacements[res_pos].upper()
                    try:
                        new_res_name = one_to_three(one_letter)
                        print(f"Replacing residue {residue.resname} at position {res_pos} with {new_res_name}")
                        residue.resname = new_res_name
                    except ValueError as e:
                        print(e)

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_output_path)
    print(f"Modified PDB saved to {pdb_output_path}")

# Example usage
if __name__ == "__main__":
    pdb_input = "pdb_outputs/P24941.pdb"
    pdb_output = "pdb_outputs/P24941_mut.pdb"
    
    # Format: {position: 'NEW_RES_FASTA'}, e.g., 'A' = Alanine, 'G' = Glycine
    replacements = {
        10: 'A',  # ALA
        25: 'G',  # GLY
        42: 'S',   # SER
        59: 'T',   # SER
        63: 'S',   # SER
        70: 'S',   # SER
        72: 'S',   # SER
        200: 'T',   # SER
        42: 'S',   # SER
        123: 'G',   # SER
    }

    replace_residues(pdb_input, pdb_output, replacements)

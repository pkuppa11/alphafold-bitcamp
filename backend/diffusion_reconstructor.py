from Bio import SeqIO
from Bio.Seq import Seq
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.PDBParser import PDBParser

def extract_fasta_from_seqres(pdb_path):
    three_to_one_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z', 'XLE': 'J',
        'UNK': 'X'
    }

    fasta_seqs = {}
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
    
    seqres = {}
    for line in lines:
        if line.startswith("SEQRES"):
            chain_id = line[11]
            residues = line[19:].split()
            if chain_id not in seqres:
                seqres[chain_id] = []
            seqres[chain_id].extend(residues)

    for chain, reslist in seqres.items():
        seq = ""
        for res in reslist:
            aa = three_to_one_map.get(res.upper(), 'X')  # fallback to X for unknowns
            seq += aa
        fasta_seqs[chain] = seq

    return fasta_seqs

sequences = extract_fasta_from_seqres("pdb_outputs/P24941.pdb")
for chain, seq in sequences.items():
    print(f">Chain_{chain}\n{seq}")

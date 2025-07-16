import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pymsaviz import MsaViz
from Bio import AlignIO, PDB
from Bio.PDB import PDBParser, PPBuilder
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def extract_sequence_from_pdb(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", pdb_path)
    ppb = PPBuilder()

    chains = list(structure[0])  # First model only

    chain = chains[0]
    peptides = ppb.build_peptides(chain)

    sequence = ''.join([str(peptide.get_sequence()) for peptide in peptides])

    pdb_id = os.path.splitext(os.path.basename(pdb_path))[0]
    return SeqRecord(Seq(sequence), id=pdb_id, description="")


def frequency(align):
    amino_acids = ['-']
    df = pd.DataFrame(columns=amino_acids, data=np.zeros((len(align[0]), len(amino_acids))))

    for a in range(len(align[0])):
        for x in align[:, a]:
            if x in amino_acids:
                df.at[a, x] += 1
            else:
                amino_acids.append(x)
                df[x] = 0
                df.at[a, x] = 1
    df.index = df.index + 1
    return df


def main(msa_path, pdb_file, output_dir):
    if not os.path.exists(msa_path):
        print(f"MSA file not found: {msa_path}")
        sys.exit(1)

    if not os.path.exists(pdb_file):
        print(f"PDB file not found: {pdb_file}")
        sys.exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read MSA
    align = AlignIO.read(msa_path, "fasta")

    # Frequency table
    freq_table = frequency(align)

    most_frequent_amino_acids = []
    identity = []

    for _, row in freq_table.iterrows():
        most_freq = row.idxmax()
        max_count = row.max()
        total_count = row.sum()

        identity_percentage = 0 if most_freq == "-" else (max_count / total_count) * 100
        most_frequent_amino_acids.append(most_freq)
        identity.append(identity_percentage)

    freq_table['MostFrequentAminoAcid'] = most_frequent_amino_acids
    freq_table['FrequencyPercentage'] = identity

    print(freq_table)

    # Load and modify PDB
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("my_protein", pdb_file)

    query_record = extract_sequence_from_pdb(pdb_file)

    if len(query_record) != len(identity):
        print(f"WARNING: Number of residues in PDB ({len(query_record)}) "
              f"does not match number of entries in identity list ({len(identity)}).")

    chain = list(structure[0])[0]  # First model, first chain
    residues = [res for res in chain if res.id[0] == ' ']

    if len(residues) != len(identity):
        print(f"⛔ Chain has {len(residues)} residues, expected {len(identity)}")
        sys.exit(1)

    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
    output_pdb = os.path.join(output_dir, f"{pdb_name}_with_frequencies.pdb")

    print(f"✅ Processing file {pdb_file} chain {chain.id} with {len(residues)} residues")

    for i, residue in enumerate(residues):
        for atom in residue:
            atom.set_bfactor(round(identity[i], 2))


    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

    mv = MsaViz(open(msa_path), wrap_length=60, show_grid=True, show_consensus=True)
    msa_image_path = os.path.join(output_dir, "msa.png")
    mv.savefig(msa_image_path)

    print(f"Saved modified PDB to {output_pdb}")
    print(f"Saved MSA image to {msa_image_path}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python conservationPDB.py <msa_fasta_file> <pdb_file> <output_directory>")
        sys.exit(1)

    msa_input = sys.argv[1]
    pdb_file = sys.argv[2]
    output_dir = sys.argv[3]
    main(msa_input, pdb_file, output_dir)

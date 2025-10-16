import sys 
if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python prueba2.py <msa_fasta_file> <pdb_file> <output_directory>")
        sys.exit(1)

    msa_input = sys.argv[1]
    pdb_file = sys.argv[2]
    output_dir = sys.argv[3]
    main(msa_input, pdb_file, output_dir)

def main(msa_input, pdb_file, output_dir):
    if not os.path.exists(msa_input):
        print(f"MSA file not found: {msa_input}")
        sys.exit(1)

    if not os.path.exists(pdb_file):
        print(f"PDB file not found: {pdb_file}")
        sys.exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

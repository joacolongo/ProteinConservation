import sys 
if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python prueba2.py <msa_fasta_file> <pdb_file> <output_directory>")
        sys.exit(1)

    msa_input = sys.argv[1]
    pdb_file = sys.argv[2]
    output_dir = sys.argv[3]
    main(msa_input, pdb_file, output_dir)


import os
import requests
from Bio.PDB import PDBParser, PPBuilder
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
import subprocess

def download_pdb(pdb_id, save_dir):
    pdb_id = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_path = os.path.join(save_dir, f"{pdb_id}.pdb")
    if os.path.exists(output_path):
        print(f"PDB file {output_path} already exists, skipping download.")
        return output_path
    print(f"Downloading PDB {pdb_id} to {output_path} ...")
    response = requests.get(url)
    if response.status_code == 200:
        with open(output_path, 'w') as f:
            f.write(response.text)
        print(f"Downloaded {pdb_id}.pdb to {output_path}")
        return output_path
    else:
        raise Exception(f"Could not download PDB file: {pdb_id}")

def extract_sequence_from_pdb(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", pdb_path)
    ppb = PPBuilder()
    records = []
    seen_chains = set()

    for model in structure:
        for chain in model:
            chain_id = chain.id
            if chain_id in seen_chains:
                continue
            seen_chains.add(chain_id)
            peptides = ppb.build_peptides(chain)
            if not peptides:
                continue
            sequence = ''.join([str(peptide.get_sequence()) for peptide in peptides])
            record = SeqRecord(Seq(sequence), id=f"Chain_{chain_id}", description="")
            records.append(record)
    if not records:
        raise ValueError("No peptide sequences found in PDB.")
    return records[0]  # Use first chain as query

def run_blastp(query_seq_record, blast_output_xml, max_hits=100):
    print("Running BLAST search against nr...")
    result_handle = NCBIWWW.qblast(
        program="blastp",
        database="nr",
        sequence=query_seq_record.seq,
        hitlist_size=max_hits,
        format_type="XML"
    )
    with open(blast_output_xml, "w") as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()
    print(f"BLAST search complete. Results saved to {blast_output_xml}")

def parse_blast_results(blast_xml_file, min_identity=30.0):
    records = []
    with open(blast_xml_file) as result_handle:
        blast_record = NCBIXML.read(result_handle)
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                identity = (hsp.identities / hsp.align_length) * 100
                if identity >= min_identity:
                    seq = hsp.sbjct.replace('-', '')  # Remove gaps
                    record = SeqRecord(Seq(seq), id=alignment.accession, description=alignment.title)
                    records.append(record)
                    break  # only one HSP per alignment
    return records

def write_fasta(records, filename):
    SeqIO.write(records, filename, "fasta")
    print(f"Wrote {len(records)} sequences to {filename}")

def run_clustal_omega(fasta_in, fasta_out):
    subprocess.run([
        "clustalo", "-i", fasta_in, "-o", fasta_out, "--force", "--outfmt=fasta"
    ], check=True)
    print(f"Multiple Sequence Alignment saved to {fasta_out}")



def main(pdb_source, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Check if pdb_source looks like a PDB code (4 letters, no extension)
    if len(pdb_source) == 4 and pdb_source.isalnum() and not os.path.isfile(pdb_source):
        pdb_file = download_pdb(pdb_source, output_dir)
    else:
        pdb_file = pdb_source
        if not os.path.isfile(pdb_file):
            raise FileNotFoundError(f"File not found: {pdb_file}")

    query_record = extract_sequence_from_pdb(pdb_file)
    print(f"Extracted query sequence {query_record.id} length {len(query_record.seq)}")

    blast_xml = os.path.join(output_dir, "blast_results.xml")
    run_blastp(query_record, blast_xml, max_hits=100)

    hits = parse_blast_results(blast_xml, min_identity=30.0)
    print(f"Parsed {len(hits)} BLAST hits with >=30% identity")

    all_seqs = [query_record] + hits

    fasta_file = os.path.join(output_dir, "all_sequences.fasta")
    write_fasta(all_seqs, fasta_file)

    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]

    aligned_file = os.path.join(output_dir, f"{pdb_name}_aligned.fasta")
    run_clustal_omega(fasta_file, aligned_file)


# === Example usage ===
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage: python create_MSA.py <pdb_code_or_file> <output_directory>")
        sys.exit(1)
    pdb_source = sys.argv[1]
    output_dir = sys.argv[2]
    main(pdb_source, output_dir)

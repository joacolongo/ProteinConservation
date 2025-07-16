# Protein Conservation

From a multiple sequence alignment (MSA) in FASTA format and a corresponding PDB protein structure, sequence identity at each alignment position is computed and map it onto the PDB file in the B-factor position.

If a MSA of the protein is not available, create_MSA.py script can be used.

WARNING: Only works when the protein is a monomer and the multifasta sequence have the same length as the PDB file

## Run example  

Create the MSA
```bash
python create_MSA.py example/AF-P9WKK7-F1-model_v4.pdb example
```
Create the PDB with identity values

```bash
python conservationPDB.py  example/AF-P9WKK7-F1-model_v4_aligned.fasta example/AF-P9WKK7-F1-model_v4.pdb  example
```

import sys
import subprocess
import pandas as pd
from Bio import Entrez, SeqIO
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

Entrez.email = "your_email@example.com"  # Replace with your actual email

# --- Configurable inputs ---
TARGET_GENE = sys.argv[1] if len(sys.argv) > 1 else "Laccase"
TAXA_FILTER = sys.argv[2] if len(sys.argv) > 2 else "Bacteria[Organism] OR Fungi[Organism]"
MAX_SEQS = 20

query = f"{TARGET_GENE}[Gene/Protein Name] AND ({TAXA_FILTER})"

print(f"Searching NCBI for: {query}")
handle = Entrez.esearch(db="protein", term=query, retmax=MAX_SEQS, idtype="acc")
record = Entrez.read(handle)
ids = record["IdList"]

if not ids:
    print("No sequences found — try a broader gene name or taxa filter.")
    sys.exit(1)

print(f"Fetching {len(ids)} sequences...")
fetch_handle = Entrez.efetch(db="protein", id=ids, rettype="gb", retmode="text")
records = list(SeqIO.parse(fetch_handle, "genbank"))

fasta_path = "data/raw_sequences.fasta"
metadata = []
with open(fasta_path, "w") as f:
    for rec in records:
        f.write(f">{rec.id}\n{rec.seq}\n")
        metadata.append({
            "Accession": rec.id,
            "Organism": rec.annotations.get("organism", "Unknown"),
            "Length": len(rec.seq)
        })

pd.DataFrame(metadata).to_csv("output/metadata.csv", index=False)
print("Metadata saved to output/metadata.csv")

# --- Alignment via subprocess ---
print("Aligning sequences with MAFFT...")
aligned_path = "data/aligned.fasta"

with open(aligned_path, "w") as out_f:
    result = subprocess.run(["mafft", fasta_path], stdout=out_f, stderr=subprocess.PIPE, text=True)

if result.returncode != 0:
    print(f"MAFFT failed with error:\n{result.stderr}")
    sys.exit(1)

# Verify alignment
alignment = AlignIO.read(aligned_path, "fasta")
if len(alignment) < 2:
    print("Alignment failed or produced too few sequences. Stopping.")
    sys.exit(1)

print(f"Alignment complete: {len(alignment)} sequences, {alignment.get_alignment_length()} columns")

# --- Build the tree ---
print("Calculating distances and building tree...")
calculator = DistanceCalculator("identity")
dm = calculator.get_distance(alignment)

constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

tree_path = "output/tree.nwk"
Phylo.write(tree, tree_path, "newick")
print(f"Tree saved to {tree_path}")

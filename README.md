# Phylogenomics-Pipeline

Automated tool that takes a gene/protein name, pulls homologous
sequences from NCBI across specified taxa, aligns them with MAFFT,
builds a neighbor-joining phylogenetic tree, and outputs an
interactive HTML report.

## Usage
```bash
conda env create -f environment.yml
conda activate bio_env
python run_pipeline.py "YourGeneName" "Bacteria[Organism] OR Fungi[Organism]"
python generate_report.py
```

## Output
- `output/metadata.csv` — sequence accessions, organisms, lengths
- `output/tree.nwk` — phylogenetic tree (Newick format)
- `output/phylogenetic_report.html` — interactive report

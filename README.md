# GeomDrug: A Bioinformatics Utility Toolbox

A modular collection of Python tools for cheminformatics, structural biology, and sequence analysis. Designed to accelerate daily bioinformatics research workflows.

---

## Overview

This repository provides a set of reusable utilities spanning small-molecule cheminformatics (SMILES handling, conformer generation, CCD database queries), macromolecular structure parsing (PDB/CIF I/O, RNA structure extraction, ligand analysis), sequence processing (FASTA I/O, pairwise identity computation), and computational resource management (GPU allocation). All functions are self-contained and can be imported individually into analysis scripts or computational pipelines.

---

## Repository Organization

| File | Description |
|------|-------------|
| `bio_utils.py` | Core module: 35+ functions covering SMILES/cheminformatics, PDB/CIF structure parsing, RNA structure analysis, FASTA I/O, distance matrix computation, and RCSB CCD queries |
| `data_process.py` | Data processing utilities: true positive rate under thresholds, Gaussian noise augmentation |
| `device.py` | GPU management: automatic GPU selection via `nvidia-smi`, SLURM-aware device allocation |
| `file_process.py` | File utilities: PNG-to-PDF conversion, white-to-transparent image processing |
| `draw.ipynb` | Visualization recipe notebook with publication-quality plotting examples (matplotlib, seaborn, py3Dmol) |
| `bio.sh` | Shell utility: extract sequence IDs from FASTA files |
| `Some Scripts/sdf_more_conformers.py` | Standalone workflow: conformer generation, RMSD-based hierarchical clustering, filtering |

---

## Documentation

Complete function-level documentation with signatures, parameters, return values, and usage notes is available in **[DOCUMENTATION.md](DOCUMENTATION.md)**.

---

## Quick Start

```python
# Cheminformatics
from bio_utils import smiles_to_ecfp1024, get_smiles_from_ccd
fps = smiles_to_ecfp1024(["c1ccccc1", "CC(=O)O"])
smiles = get_smiles_from_ccd("ATP")

# Structural biology
from bio_utils import load_rna_bases, extract_ligand_info_from_cif
base_atoms, rna_seqs = load_rna_bases("structure.pdb", return_seq=True)
ligand_mol = extract_ligand_info_from_cif("structure.cif", "LIG", chain_id="A")

# Sequence analysis
from bio_utils import get_fasta_ids_and_sequences, calc_identity
ids, seqs = get_fasta_ids_and_sequences("sequences.fasta")
identity = calc_identity(seqs[0], seqs[1])
```

---

## Requirements

Core dependencies: `rdkit`, `biotite`, `biopython`, `prody`, `gemmi`, `numpy`, `pandas`, `matplotlib`, `seaborn`, `pillow`, `requests`, `GPUtil`

Optional: `py3Dmol` (for Jupyter 3D visualization), `scipy` (for conformer clustering), `reportlab` (for PDF generation)

External CLI tools: [Open Babel](https://openbabel.org/) (`obabel`), NVIDIA `nvidia-smi`

---

## Contact

**Xia Yuping** — [xiayp99@gmail.com](mailto:xiayp99@gmail.com)

For inquiries regarding custom pipeline development, multi-omics data integration, or structural biology analysis, please reach out via email.

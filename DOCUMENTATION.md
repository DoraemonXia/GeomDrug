# Bioinfo-Toolbox Documentation

## Overview

Bioinfo-Toolbox is a modular collection of Python utilities for cheminformatics, structural biology, and bioinformatics data processing. The tools cover small-molecule handling (SMILES, conformer generation, CCD queries), macromolecular structure parsing (PDB/CIF I/O, RNA structure extraction), sequence analysis (FASTA I/O, pairwise identity), GPU resource management, and image/file processing.

All functions are designed to be imported individually and used interactively or within pipeline scripts.

---

## Dependencies

| Package | Purpose |
|---------|---------|
| `rdkit` | Cheminformatics: SMILES, fingerprints, conformers, MOL2/SDF I/O |
| `biotite` | PDB file I/O, structure filtering |
| `biopython` | PDB/CIF parsing, FASTA I/O, pairwise sequence alignment |
| `prody` | PDB parsing, backbone extraction |
| `gemmi` | CIF structure parsing |
| `numpy` | Numerical arrays and linear algebra |
| `pandas` | Tabular data (imported but used sparingly) |
| `matplotlib` | Publication-quality plotting |
| `seaborn` | Statistical visualizations |
| `py3Dmol` | Interactive 3D molecular visualization in Jupyter |
| `pillow` (PIL) | Image I/O and pixel manipulation |
| `reportlab` | PDF generation from images |
| `scipy` | Hierarchical clustering (`sdf_more_conformers.py`) |
| `requests` | HTTP downloads, GraphQL API queries |
| `GPUtil` | GPU memory querying (`device.py`) |

**External CLI tools**: `obabel` (Open Babel) for PDB-to-SDF conversion, `nvidia-smi` for GPU status.

---

## Module Reference

### 1. `bio_utils.py` — Core Utilities

The central module. Organised into six sub-domains.

#### 1.1 SMILES and Cheminformatics

##### `contains_functional_group(molecule_smiles, functional_group_smarts) -> bool`

Determine whether a molecule contains a specified functional group.

- **molecule_smiles** (`str`): SMILES representation of the molecule.
- **functional_group_smarts** (`str`): SMARTS representation of the functional group.
- **Returns**: `True` if the substructure is found, `False` otherwise.

##### `smiles_to_ecfp1024(smiles_list) -> np.ndarray`

Convert a list of SMILES strings to a 1024-bit ECFP (Morgan fingerprint, radius 2) matrix.

- **smiles_list** (`list[str]`): Input SMILES strings.
- **Returns**: NumPy array of shape `(N, 1024)`. Invalid SMILES are replaced with an all-zero row.

##### `smiles_to_png(smiles, output_file) -> None`

Render a SMILES string as a 300×300 PNG image with a transparent background.

- **smiles** (`str`): SMILES string to render.
- **output_file** (`str`): Path for the output PNG file.
- **Raises**: `ValueError` if the SMILES string is invalid.

##### `save_multiple_conformers_to_sdf(smiles, num_confs=100, output_file='output.sdf') -> None`

Generate multiple 3D conformers for a SMILES molecule using ETKDG, UFF-optimize each, and write to an SDF file with energy annotations.

- **smiles** (`str`): Input SMILES.
- **num_confs** (`int`): Number of conformers to generate. Default 100.
- **output_file** (`str`): Output SDF path.

##### `translate_molecule(sdf_file, output_file, dx, dy, dz) -> None`

Translate all atom coordinates in all molecules of an SDF file by a given displacement vector.

- **sdf_file** (`str`): Path to the input SDF file.
- **output_file** (`str`): Path for the output SDF file.
- **dx, dy, dz** (`float`): Translation along each axis.

##### `get_smiles_from_ccd(ccd_id: str) -> str | None`

Query the RCSB GraphQL API for the canonical SMILES of a CCD (Chemical Component Dictionary) identifier. Uses a three-tier priority fallback: (1) OpenEye canonical, (2) any canonical, (3) any SMILES. The returned SMILES is re-standardized through RDKit.

- **ccd_id** (`str`): CCD code (e.g., `"ATP"`).
- **Returns**: Canonical SMILES string, or `None` if no SMILES is found.

---

#### 1.2 PDB / CIF Structural Biology

##### `get_pdb_file_information(file_name, if_multi_struc=False, atoms=["N","CA","C","O"]) -> tuple`

Parse a PDB file and return backbone atom objects and the amino acid sequence for each chain.

- **file_name** (`str`): Path to the PDB file.
- **if_multi_struc** (`bool`): If `True`, return lists for all chains; otherwise return the first chain only.
- **atoms** (`list[str]`): Atom names to select. Default `["N", "CA", "C", "O"]`.
- **Returns**: `(backbone_list, pdb_seq_list)` — lists of backbone atom arrays and sequences.

##### `write_coords_to_pdb(coords: np.ndarray, out_fname: str, ...) -> str`

Write a NumPy array of backbone coordinates (N, CA, C, O per residue) back to a PDB file using biotite.

- **coords** (`np.ndarray`): Coordinate array of shape `(4 * N_residues, 3)`.
- **out_fname** (`str`): Output PDB filename.
- **if_bonds** (`bool`): If `True`, add covalent bonds between consecutive atoms.
- **select_atoms** (`list[str]`): Atom names per residue. Default `["N", "CA", "C", "O"]`.
- **Returns**: The output filename.

##### `extract_backbone_coords(fname, atoms=["N","CA","C","O"]) -> np.ndarray | None`

Extract backbone atom coordinates from a PDB file (supports `.gz` compression).

- **fname** (`str`): Path to the PDB or `.pdb.gz` file.
- **atoms** (`Collection`): Atom names to extract.
- **Returns**: NumPy array of shape `(N_atoms, 3)`.

##### `extract_backbone_atoms(pdb_file, output_dir) -> None`

Use ProDy to parse a PDB, select the longest chain (max 1500 residues), extract N/CA/C/O backbone atoms, and save to a new PDB file.

- **pdb_file** (`str`): Path to the input PDB file.
- **output_dir** (`str`): Directory for the output PDB file.

##### `extract_ligand_info_from_cif(cif_file, ligand_name, chain_id="A") -> rdkit.Chem.Mol`

Extract a specific non-RNA residue (ligand) from a CIF file as an RDKit Mol with 3D coordinates. Residues named A, C, G, U, T are skipped (treated as RNA).

- **cif_file** (`str`): Path to the CIF file.
- **ligand_name** (`str`): Three-letter residue name of the ligand.
- **chain_id** (`str`): Chain to search. Default `"A"`.
- **Returns**: RDKit Mol object with a 3D conformer.

##### `extract_rna_atoms_from_cif(cif_file, chain_id=-1) -> list`

Extract all RNA base (A/U/C/G) atoms from a CIF file.

- **cif_file** (`str`): Path to the CIF file.
- **chain_id** (`int/str`): Chain identifier. If `-1`, extract from all chains.
- **Returns**: List of BioPython Atom objects.

##### `save_atoms_to_pdb_simple(atoms, output_file) -> None`

Write extracted atoms to a minimal PDB file in ATOM record format.

- **atoms** (`list`): List of Atom objects.
- **output_file** (`str`): Output PDB path.

##### `download_pdb_files(pdb_ids, output_folder, timeout=600) -> None`

Batch download PDB files from RCSB. Falls back to CIF download if a PDB file is unavailable or the request times out.

- **pdb_ids** (`list[str]`): List of PDB IDs.
- **output_folder** (`str`): Output directory.
- **timeout** (`int`): Download timeout in seconds. Default 600.

##### `download_cif_file(url, pdb_id, output_folder, timeout) -> None`

Internal helper for downloading a single CIF file.

##### `detect_chain_break(cif_path: str) -> bool`

Check whether a CIF file contains discontinuous chain breaks by catching `PDBConstructionWarning` during parsing.

- **cif_path** (`str`): Path to the CIF file.
- **Returns**: `True` if a break is detected, `False` otherwise.

##### `get_hetatm_residues(cif_file) -> set`

Return a set of HETATM residue names from a CIF file. Useful for identifying ligands, cofactors, and solvent.

- **cif_file** (`str`): Path to the CIF file.
- **Returns**: Set of residue name strings.

##### `convert_pdb_to_sdf(input_folder) -> list`

Batch convert all `.pdb` files in a directory to `.sdf` using Open Babel (`obabel`). Two versions exist in the module; the effective version (line ~1084) uses the `-d` flag to remove hydrogen atoms.

- **input_folder** (`str`): Path to the directory containing `.pdb` files.
- **Returns**: List of filenames that failed to convert.

---

#### 1.3 RNA Structure

##### `parse_RNAfold(file_path) -> dict`

Parse a Vienna RNAfold output file and return a dictionary mapping sequence names to dot-bracket secondary structure strings.

- **file_path** (`str`): Path to the RNAfold output file.
- **Returns**: `dict` of `{seq_name: dot_bracket_structure}`.

##### `parse_rna_structure(structure) -> list[tuple]`

Convert a dot-bracket RNA secondary structure string into a list of directed edge pairs. Base-pairing edges come from parentheses; covalent backbone edges come from sequential adjacency.

- **structure** (`str`): Dot-bracket notation (e.g., `"(((...)))"`).
- **Returns**: List of `(from_idx, to_idx)` tuples.

##### `load_rna_bases_from_mol2(file_path) -> list[list[int]]`

Parse a MOL2 file containing multiple `@<TRIPOS>MOLECULE` blocks and extract atom index groups per RNA base using the residue ID column.

- **file_path** (`str`): Path to the MOL2 file.
- **Returns**: List of lists, each sublist containing atom indices for one base.

##### `load_rna_bases_from_pdb(file_path, return_seq=False, chain_id=-1) -> dict | tuple`

Parse a PDB file and extract RNA base atom coordinates grouped by residue.

- **file_path** (`str`): Path to the PDB file.
- **return_seq** (`bool`): If `True`, also return sequence strings.
- **chain_id** (`int/str`): Chain identifier. `-1` for all chains.
- **Returns**: `base_atoms` dict (`{base_id: [coord, ...]}`) and optionally `rna_sequences` list.

##### `load_rna_bases(file_path, return_seq=False, chain_id=-1) -> dict | tuple`

Unified function for PDB and CIF files. Uses BioPython's `PDBParser` for PDB files and `gemmi` for CIF files. Auto-detects format by file extension.

- **file_path** (`str`): Path to the PDB or CIF file.
- **return_seq** (`bool`): If `True`, also return sequences.
- **chain_id** (`int/str`): Chain identifier. `-1` for all chains.
- **Returns**: `base_atoms` dict of NumPy arrays, optionally with `rna_sequences`.
- **Raises**: `ValueError` on unsupported file formats.

##### `convert_base_atoms_to_matrix(base_atoms) -> dict`

Convert a dict of base atom lists into a dict of NumPy matrices (`n_atoms × 3`).

- **base_atoms** (`dict`): Mapping from base IDs to atom coordinate lists.
- **Returns**: Dict of `{base_key: np.ndarray}`.

##### `calculate_min_distance_matrix(base_atoms, mol_coords) -> np.ndarray`

Compute a distance matrix: for each RNA base and each ligand atom, the minimum Euclidean distance between any atom in that base and that ligand atom.

- **base_atoms** (`dict`): Base atom coordinate matrices.
- **mol_coords** (`np.ndarray`): Ligand atom coordinates of shape `(M, 3)`.
- **Returns**: NumPy array of shape `(N_bases, M_atoms)`.

##### `extract_atom_coordinates(pdb_file, atoms_to_extract=['C4\'']) -> tuple`

Extract specific named atoms (e.g., C4', N1, P) from all standard residues in a PDB file.

- **pdb_file** (`str`): Path to the PDB file.
- **atoms_to_extract** (`list[str]`): Atom names to select.
- **Returns**: `(atom_ids_array, coordinates_array)`.

> **Note**: This function references an undefined variable `RNA_atom_id`. A user-provided mapping from atom names to numeric IDs must be defined before calling this function.

---

#### 1.4 MOL2 File Handling

##### `load_multiple_mol2(file_path) -> list`

Load a MOL2 file containing multiple `@<TRIPOS>MOLECULE` blocks into a list of RDKit Mol objects (without sanitization).

- **file_path** (`str`): Path to the MOL2 file.
- **Returns**: List of RDKit Mol objects.
- **Raises**: `ValueError` if no molecules could be loaded.

---

#### 1.5 FASTA and Sequence Utilities

##### `get_fasta_ids_and_sequences(fasta_file) -> tuple`

Parse a FASTA file and return two lists: sequence IDs and sequences.

- **fasta_file** (`str`): Path to the FASTA file.
- **Returns**: `(ids: list[str], sequences: list[str])`.

##### `generate_fasta(rna_sequences, names=None, fasta_file_path='output.fasta', reverse=False) -> None`

Write a list of RNA sequences to a FASTA file. If `reverse=True`, replaces T with U in all sequences.

- **rna_sequences** (`list[str]`): RNA sequences.
- **names** (`list[str]`, optional): Sequence names. Autogenerated if `None`.
- **fasta_file_path** (`str`): Output file path.
- **reverse** (`bool`): T→U substitution flag.
- **Raises**: `ValueError` if lengths of sequences and names do not match.

##### `extract_fasta_ids(input_folder, output_file, sequence_index=None) -> None`

Merge FASTA files from a directory into a single file. Optionally extract only a specific sequence index from each file (useful when generating sequences with ProteinMPNN/LigandMPNN).

- **input_folder** (`str`): Directory containing `.fasta`/`.fa` files.
- **output_file** (`str`): Output FASTA path.
- **sequence_index** (`int`, optional): If specified, extract only the `sequence_index`-th sequence from each file.

##### `split_fasta(input_fasta, output_dir, chunk_size=100) -> None`

Split a large FASTA file into chunks of `chunk_size` records each.

- **input_fasta** (`str`): Path to the input FASTA file.
- **output_dir** (`str`): Output directory.
- **chunk_size** (`int`): Records per chunk. Default 100.

##### `merge_fasta_files(fasta_files, output_filename) -> None`

Concatenate multiple FASTA files into one output file, preserving blank line separation between source files.

- **fasta_files** (`list[str]`): List of FASTA file paths.
- **output_filename** (`str`): Output file path.

---

#### 1.6 Sequence Analysis

##### `calc_identity(seq1, seq2) -> float`

Calculate pairwise sequence identity (percentage) using BioPython's global alignment (`pairwise2.align.globalxx`).

- **seq1** (`str`): First sequence.
- **seq2** (`str`): Second sequence.
- **Returns**: Identity percentage as a float (0–100).

---

#### 1.7 Miscellaneous

##### `save_mol_to_pdb(mol, filename) -> None`

Save an RDKit Mol object to a PDB format file.

- **mol** (`rdkit.Chem.Mol`): RDKit molecule object.
- **filename** (`str`): Output filename including path.

---

### 2. `data_process.py` — Data Processing

##### `get_true_positive_rate_under_thresholds(probs, label) -> list`

Compute the true positive rate at 100 evenly spaced thresholds from 0 to 1 (step 0.01).

- **probs** (`array-like`): Predicted probabilities.
- **label** (`array-like`): Binary ground-truth labels.
- **Returns**: List of TPR values at each threshold.

##### `add_gaussian_noise(matrix, scale_factor=0.35) -> np.ndarray`

Add scaled Gaussian noise to a matrix. The noise standard deviation at each element is `scale_factor * |matrix[i,j]|`.

- **matrix** (`np.ndarray`): Input matrix.
- **scale_factor** (`float`): Noise scaling factor. Default 0.35.
- **Returns**: Matrix with added noise.

---

### 3. `device.py` — GPU Management

##### `find_max_gpu(file_path) -> str`

Query `nvidia-smi` to find the GPU with the most free memory and set the `CUDA_VISIBLE_DEVICES` environment variable to that GPU's index.

- **file_path** (unused): Reserved parameter.
- **Returns**: GPU index string.

##### `get_available_device_slurm(num_device, only_empty=True) -> list[int]`

Retrieve available GPU local indices based on `CUDA_VISIBLE_DEVICES`. Optionally filter to GPUs with low memory usage (<1000 MB used).

- **num_device** (`int`): Maximum number of GPU indices to return.
- **only_empty** (`bool`): If `True`, exclude GPUs using ≥1000 MB.
- **Returns**: List of local GPU indices.

---

### 4. `file_process.py` — File Utilities

##### `images_to_pdf(folder_path, output_pdf, include_title=False) -> None`

Convert all PNG images in a folder into a multi-page PDF.

- **folder_path** (`str`): Path to the folder containing PNG images.
- **output_pdf** (`str`): Output PDF file path.
- **include_title** (`bool`): If `True`, render each filename as a page title.

##### `convert_white_to_transparent(input_image_path, output_image_path) -> None`

Replace white and near-white pixels (RGB all >200) with full transparency (RGBA alpha=0) in a PNG image. Useful for cleaning figure backgrounds.

- **input_image_path** (`str`): Path to the input PNG.
- **output_image_path** (`str`): Path for the output PNG.

---

### 5. `bio.sh` — FASTA ID Extraction

A shell utility that reads a FASTA file and prints all sequence IDs (header lines stripped of the `>` character).

**Usage**:
```bash
./bio.sh <path_to_fasta_file>
```

---

### 6. `Some Scripts/sdf_more_conformers.py` — Conformer Clustering

A standalone workflow script that:

1. Reads a molecule from `PC4.sdf`
2. Generates 500 conformers using ETKDG
3. Computes the pairwise RMSD matrix across all conformers
4. Performs hierarchical clustering (`scipy.cluster.hierarchy.linkage`, average linkage) into 5 clusters
5. Filters conformers with RMSD ≤ 1.5 Å from the reference (conf_id=0)
6. Saves filtered conformers as SDF files into per-cluster subdirectories (`Cluster_1/` through `Cluster_5/`)

**Dependencies**: RDKit, NumPy, SciPy.

---

### 7. `draw.ipynb` — Visualization Recipes

A Jupyter notebook containing 7 visualization cells suitable for publication-quality figures:

| Cell | Type | Description |
|------|------|-------------|
| 0 | Grouped bar chart | Residue-type-matched vs. regardless comparisons at 1/2/4 Å cutoffs |
| 1 | Custom heatmap | Green palette, no color bar, no tick labels |
| 2 | 2D histogram | Helix vs. strand composition density |
| 3 | Score distribution | Histogram with customizable bins and range clipping |
| 4 | Box + strip plot | Decoy rank scores by method with significance annotations |
| 5 | Joint plot + KDE | Scatter density with KL divergence annotation |
| 6 | 3D molecule animation | py3Dmol protein trajectory with cartoon style and ligand sticks |

All cells use inline English comments. Copy individual cells to your own notebooks for reuse.

---

## Appendix: External Dependencies

**Python packages** (installable via `pip`):

| Package | Typical Import Name |
|---------|---------------------|
| numpy | `import numpy as np` |
| pandas | `import pandas as pd` |
| rdkit | `from rdkit import Chem` |
| biotite | `import biotite.structure as struc` |
| biopython | `from Bio import SeqIO` |
| prody | `import prody as pr` |
| gemmi | `import gemmi` |
| scipy | `from scipy.cluster.hierarchy import linkage` |
| matplotlib | `import matplotlib.pyplot as plt` |
| seaborn | `import seaborn as sns` |
| py3Dmol | `import py3Dmol` |
| pillow | `from PIL import Image` |
| reportlab | `from reportlab.pdfgen import canvas` |
| requests | `import requests` |
| GPUtil | `import GPUtil` |

**Command-line tools**:

| Tool | Required By | Purpose |
|------|-------------|---------|
| `obabel` (Open Babel) | `bio_utils.py::convert_pdb_to_sdf` | PDB-to-SDF format conversion |
| `nvidia-smi` | `device.py::find_max_gpu` | GPU memory status query |

---

*Documentation maintained for the Bioinfo-Toolbox project. Contributions and extensions are welcome.*

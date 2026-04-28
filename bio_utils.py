import gzip
import numpy as np
import pandas as pd
from typing import Collection, Optional, Literal

from rdkit import Chem
from rdkit.Chem import AllChem

from biotite.structure.io.pdb import PDBFile
from biotite.sequence import ProteinSequence
import biotite.structure as struc
from Bio import SeqIO



from biotite.structure import AtomArray, filter_backbone,filter_amino_acids


def contains_functional_group(molecule_smiles, functional_group_smarts):
    """
    Determine whether a given molecule contains a specified functional group.

    :param molecule_smiles: SMILES representation of the molecule.
    :param functional_group_smarts: SMARTS representation of the functional group.
    :return: True if the molecule contains the functional group, False otherwise.
    """
    # Create molecule and functional group objects
    molecule = Chem.MolFromSmiles(molecule_smiles)
    functional_group = Chem.MolFromSmarts(functional_group_smarts)

    # Use substructure matching to test for the functional group
    if molecule.HasSubstructMatch(functional_group):
        return True
    else:
        return False

"""
how to use
"""
# count  = 0
# for i in range(len(Smiles)):
#     if contains_functional_group(Smiles[i], "[*]CC([*])(CC(=O)OCCCC)C(=O)OCCCC"):
#         count+=1
# print(count)



def get_fasta_ids_and_sequences(fasta_file):
    ids, sequences = [], []
    for record in SeqIO.parse(fasta_file, "fasta"):
        ids.append(record.id)
        sequences.append(str(record.seq))
    return ids, sequences

"""
Example Usage
fasta_file = 'path_to_your_fasta_file.fasta'
sequence_ids, sequences = get_fasta_ids_and_sequences(fasta_file)
print(sequence_ids)
print(sequences)
"""



def generate_fasta(rna_sequences, names=None, fasta_file_path='output.fasta', reverse=False):
    """
    Generate a FASTA file from RNA sequences.

    Parameters:
    - rna_sequences (list): List of RNA sequences.
    - names (list, optional): List of names corresponding to RNA sequences. If None, default names will be used.
    - fasta_file_path (str, optional): Path to save the generated FASTA file.
    - reverse (bool, optional): If True, replace 'T' with 'U' in the sequences.

    Returns:
    - None
    """
    if not names:
        # If names are not provided, use default names (seq_0, seq_1, ...)
        names = [f"seq_{i}" for i in range(len(rna_sequences))]
    
    if len(rna_sequences) != len(names):
        raise ValueError("Number of RNA sequences must match the number of names.")

    with open(fasta_file_path, 'w') as fasta_file:
        for i in range(len(rna_sequences)):
            sequence = rna_sequences[i]
            name = names[i]
            
            if reverse:
                # If reverse is True, replace 'T' with 'U'
                sequence = sequence.replace('T', 'U')

            fasta_file.write(f">{name}\n{sequence}\n")

'''
# Example usage:
rna_sequences = ["AUGCUAGUAC", "CCGUAGCGUA", "UAGCUAGCUA"]
names = ["seq_1", "seq_2", "seq_3"]
generate_fasta(rna_sequences, names, fasta_file_path='output.fasta', reverse=True)
'''


# Merge FASTA files from a directory. sequence_index can be None or the index of a specific
# sequence to extract from each FASTA file. (Useful when generating sequences with ProteinMPNN/LigandMPNN.)
import os
from Bio import SeqIO

def extract_fasta_ids(input_folder, output_file, sequence_index=None):
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(input_folder):
            if filename.endswith('.fasta') or filename.endswith('.fa'):
                file_path = os.path.join(input_folder, filename)
                sequences = list(SeqIO.parse(file_path, "fasta"))

                if sequence_index is not None:
                    if len(sequences) > sequence_index:
                        SeqIO.write(sequences[sequence_index], outfile, "fasta")
                else:
                    SeqIO.write(sequences, outfile, "fasta")

'''
# Example usage:
input_folder = "/xcfhome/ypxia/Workspace/LigandMPNN/outputs/ligandmpnn_temperature_01_gn_01/seqs/"
output_file = "/xcfhome/ypxia/Workspace/esm/esm/results/fasta_ids.fasta"
extract_fasta_ids(input_folder, output_file, sequence_index=1)
'''

def smiles_to_ecfp1024(smiles_list):
    fingerprints = []
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                ecfp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                fingerprints.append(np.array(ecfp))
            else:
                print(f"Invalid SMILES string: {smiles}")
        except Exception as e:
            print(f"Error processing SMILES string: {smiles}\nException: {e}")
            fingerprints.append(np.zeros((1, 1024), dtype=np.bool))

    fingerprint_matrix = np.vstack(fingerprints)
    return fingerprint_matrix

'''
# Example usage:
smiles_list = ["","",""]
smiles_to_ecfp1024(smiles_list)
'''


def get_pdb_file_information( file_name, if_multi_struc = False, atoms = ["N", "CA", "C", "O"] ):
    """
    Get PDB seq and backbone atom from .pdb file
    """
    structure = PDBFile.read(file_name)
    backbone_list = []
    pdb_seq_list = []
    
    for i in range(len(structure.get_structure())):
        chain = structure.get_structure()[i]  #There has a problem about multi model or multi chain

        # backbone = chain[struc.filter_backbone(chain)]   #this is backbone including N, CA, C atom
        
        amino_acids_atoms = chain[struc.filter_amino_acids(chain)]  #filter atom in amino acids
        
        selected_atoms = [c for c in amino_acids_atoms if c.atom_name in atoms] #this is backbone including N, CA, C, and O atom
        
        amino_acids = [ProteinSequence.convert_letter_3to1( selected_atoms[len(atoms)*i].res_name ) for i in range(int(len(selected_atoms)/len(atoms)))]
    
        pdb_seq = ""
        for j in amino_acids:
            pdb_seq+=j
            
        backbone_list.append(selected_atoms)
        pdb_seq_list.append(pdb_seq)
            
    if if_multi_struc:
        return backbone_list, pdb_seq_list
    else:
        return backbone_list[0], pdb_seq_list[0]
 
'''
# Example usage:
file_name = "../../../Workspace/Combs/data/LigandMPNN/Protein_Metal/1dwh.pdb"
backbone, pdb_seq = get_pdb_file_information(file_name)
backbone_coord = np.array( [ backbone[i].coord for i in range(len(backbone)) ] )
'''

def write_coords_to_pdb(coords: np.ndarray, out_fname: str, if_bonds = False, select_atoms = ["N", "CA", "C", "O"] ) -> str:
    """
    Write the coordinates to the given pdb fname
    """
    # Create a new PDB file using biotite
    # https://www.biotite-python.org/tutorial/target/index.html#creating-structures
    assert len(coords) % len(select_atoms) == 0, f"Expected "+str(len(select_atoms))+"N coords, got {len(coords)}"
    atoms = []
    for i, (n_coord, ca_coord, c_coord, o_coord) in enumerate(
        (coords[j : j + len(select_atoms)] for j in range(0, len(coords), len(select_atoms) ))
    ):
        atom1 = struc.Atom(
            n_coord,
            chain_id="A",
            res_id=i + 1,
            atom_id=i * len(select_atoms) + 1,
            res_name="GLY",
            atom_name="N",
            element="N",
            occupancy=1.0,
            hetero=False,
            b_factor=5.0,
        )
        atom2 = struc.Atom(
            ca_coord,
            chain_id="A",
            res_id=i + 1,
            atom_id=i * len(select_atoms) + 2,
            res_name="GLY",
            atom_name="CA",
            element="C",
            occupancy=1.0,
            hetero=False,
            b_factor=5.0,
        )
        atom3 = struc.Atom(
            c_coord,
            chain_id="A",
            res_id=i + 1,
            atom_id=i * len(select_atoms) + 3,
            res_name="GLY",
            atom_name="C",
            element="C",
            occupancy=1.0,
            hetero=False,
            b_factor=5.0,
        )
        atom4 = struc.Atom(
            o_coord,
            chain_id="A",
            res_id=i + 1,
            atom_id=i * len(select_atoms) + 4,
            res_name="GLY",
            atom_name="O",
            element="O",
            occupancy=1.0,
            hetero=False,
            b_factor=5.0,
        )
        atoms.extend([atom1, atom2, atom3, atom4])
    full_structure = struc.array(atoms)

    # Add bonds
    if if_bonds:
        full_structure.bonds = struc.BondList(full_structure.array_length())
        indices = list(range(full_structure.array_length()))
        for a, b in zip(indices[:-1], indices[1:]):
            full_structure.bonds.add_bond(a, b, bond_type=struc.BondType.SINGLE)

    # Annotate secondary structure using CA coordinates
    # https://www.biotite-python.org/apidoc/biotite.structure.annotate_sse.html
    # https://academic.oup.com/bioinformatics/article/13/3/291/423201
    # a = alpha helix, b = beta sheet, c = coil
    # ss = struc.annotate_sse(full_structure, "A")
    # full_structure.set_annotation("secondary_structure_psea", ss)

    sink = PDBFile()
    sink.set_structure(full_structure)
    sink.write(out_fname)
    return out_fname

'''
# Example usage:
origin_seq = []
file_list = []
for i in range( 1500 ):
    try:
        backbone, pdb_seq = get_pdb_file_information( files[i] )
        #backbone_coord = extract_backbone_coords(files[i])
        assert len(backbone)==4*len(pdb_seq)
        backbone_coord = np.array( [ backbone[i].coord for i in range(len(backbone)) ] )
        write_coords_to_pdb( backbone_coord, "backbone/"+all_files_and_folders[i] )
        origin_seq.append( pdb_seq )
        file_list.append(all_files_and_folders[i])
    except:
        print(all_files_and_folders[i])
'''

import gzip
from typing import Collection, Optional, Literal
import numpy as np
from biotite.structure.io.pdb import PDBFile
from biotite.structure import AtomArray, filter_backbone,filter_amino_acids

def extract_backbone_coords(
    fname: str, atoms: Collection[Literal["N", "CA", "C", "O"]] = ["N", "CA", "C", "O"]
) -> Optional[np.ndarray]:
    """Extract the coordinates of specified backbone atoms"""
    opener = gzip.open if fname.endswith(".gz") else open
    with opener(str(fname), "rt") as f:
        structure = PDBFile.read(f)
    # if structure.get_model_count() > 1:
    #     return None
    chain = structure.get_structure()[0]
    #backbone = chain
    amino_acids = chain[filter_amino_acids(chain)]  #filter atom in amino acids
    selected_atoms = [c for c in amino_acids if c.atom_name in atoms]
    coords = np.vstack([c.coord for c in selected_atoms])
    return coords

import prody as pr
import os

def extract_backbone_atoms(pdb_file, output_dir):
    structure = pr.parsePDB(pdb_file)

    chains = list(structure.getHierView().iterChains())
    
    for chain in chains:
        chain_length = len(chain)
        chain_sequence = chain.getSequence() if chain.getSequence() else "No sequence available"

        if chain_length > 1500:
            print(f"Chain {chain.getChid()} is too long ({chain_length} residues). Skipping.")
            return

    first_chain = max(chains, key=lambda chain: len(chain))
    chain_id = first_chain.getChid()

    backbone_atoms = first_chain.select('name N CA C O')

    if backbone_atoms is None:
        print(f"No backbone atoms found in chain {chain_id}.")
        return

    pdb_base_name = os.path.basename(pdb_file).split('.')[0]  # Extract filename without extension
    output_pdb = os.path.join(output_dir, f"{pdb_base_name}_{chain_id}.pdb")

    pr.writePDB(output_pdb, backbone_atoms)
    print(f"Backbone atoms saved to {output_pdb}")

'''
Example Usage
pdb_file = "/xcfhome/ypxia/Workspace/Combs/data/PDB_Download/1a2p.pdb"  # Path to the input PDB file
output_dir = '/xcfhome/ypxia/Workspace/esm/esm/Origin_backbone/'  # Path to the output directory
extract_backbone_atoms(pdb_file, output_dir)
'''



from Bio import SeqIO
import os

def split_fasta(input_fasta, output_dir, chunk_size=100):
    """
    Split fasta files
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)  # Create the output directory

    base_name = os.path.basename(input_fasta)
    name, ext = os.path.splitext(base_name)

    records = list(SeqIO.parse(input_fasta, "fasta"))

    for i in range(0, len(records), chunk_size):
        chunk_records = records[i:i+chunk_size]
        output_file = os.path.join(output_dir, f"{name}_{i//chunk_size}{ext}")
        SeqIO.write(chunk_records, output_file, "fasta")
        print(f"Saved {len(chunk_records)} records to {output_file}")

'''
# Example usage:
input_fasta = "output.fasta"  # Path to the input FASTA file
output_dir = "output/"  # Output directory
split_fasta(input_fasta, output_dir)
'''


def parse_RNAfold(file_path):
    '''
    Transfer a RNA secondary structure file into secondary structure dict. 
    '''
    sequences_dict = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith('>'):
                seq_name = line[1:]
                i += 2
                sequence = ''
                while i < len(lines) and not lines[i].startswith('>'):
                    sequence += lines[i].split()[0]
                    i += 1
                sequences_dict[seq_name] = sequence
            else:
                i += 1
    return sequences_dict

'''
# Example usage:
Robin_path = 'Robin/RNAfold.out'  # replace with your filepath
Robin_RNA_ss = parse_RNAfold(Robin_path)
'''


def parse_rna_structure(structure):
    edges = []
    stack = []
    for i, symbol in enumerate(structure):
        if symbol == '(':
            stack.append(i)
        elif symbol == ')':
            if stack:
                opening_bracket_index = stack.pop()
                edges.append((opening_bracket_index, i))
                edges.append((i, opening_bracket_index))
    
    for i in range(len(structure)-1):
        edges.append( (i, i+1) )
        edges.append( (i+1, i) )
    return edges

'''
# Example usage:
Robin_RNA_Graph = {}
for i in Robin_sequences.keys():
    x = torch.tensor( Robin_RNA_FM_repre[i] , dtype=torch.float)  #modify torch.long=>torch.float  #node features
    edge_index = torch.tensor(np.array(parse_rna_structure(Robin_RNA_ss[i]) ).T, dtype=torch.long)  #transfer rna_structure into edge_index
    data = Data(x=x, edge_index=edge_index)
    Robin_RNA_Graph[i]=data
'''

from rdkit import Chem

def load_multiple_mol2(file_path):
    """
    Load a MOL2 file containing multiple molecules and return a list of molecule objects.
    
    Parameters:
    - file_path (str): Path to the MOL2 file.
    
    Returns:
    - mol_list (list of rdkit.Chem.rdchem.Mol): List of molecule objects.
    """
    mol_list = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    mol_block = []
    for line in lines:
        if line.startswith('@<TRIPOS>MOLECULE'):
            if mol_block:
                mol = Chem.MolFromMol2Block(''.join(mol_block), sanitize=False)
                if mol is not None:
                    mol_list.append(mol)
            mol_block = [line]
        else:
            mol_block.append(line)
    
    # Add the last molecule
    if mol_block:
        mol = Chem.MolFromMol2Block(''.join(mol_block), sanitize=False)
        if mol is not None:
            mol_list.append(mol)
    
    if not mol_list:
        raise ValueError(f"Failed to load any molecules from {file_path}")
    
    return mol_list

'''
# Example usage:
small_molecules = load_multiple_mol2("RLDock/RLDock/mol_0/_cluster.mol2")
'''


def load_rna_bases_from_mol2(file_path):
    """
    Load RNA bases from a MOL2 file and return a list of atom indices for each base.
    
    Parameters:
    - file_path (str): Path to the MOL2 file.
    
    Returns:
    - base_atoms (list of list of int): List of atom indices for each base in the RNA.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    base_atoms = {}
    
    in_atom_section = False

    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom_section = True
        elif line.startswith("@<TRIPOS>"):
            in_atom_section = False
        elif in_atom_section:
            tokens = line.split()
            atom_index = int(tokens[0]) - 1  # Atom index starts from 1 in MOL2
            base_index = int(tokens[6])      # Base index from the 8th column

            if base_index not in base_atoms:
                base_atoms[base_index] = []
            
            base_atoms[base_index].append(atom_index)
    
    return [base_atoms[key] for key in sorted(base_atoms.keys())]

'''
# Example usage:
rna_base_atoms = load_rna_bases_from_mol2(rna_file_path)
'''

from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

def smiles_to_png(smiles, output_file):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    img = Draw.MolToImage(mol, size=(300, 300), kekulize=True, wedgeBonds=True, includeAtomNumbers=False)
    img = img.convert("RGBA")

    datas = img.getdata()
    new_data = []
    for item in datas:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            new_data.append((255, 255, 255, 0))  # Set white background to transparent
        else:
            new_data.append(item)
    img.putdata(new_data)

    img.save(output_file, "PNG")

'''
# Example usage:
smiles = "CN1C(=O)N(C)c2nc[nH]c2C1=O"
output_file = "1eht.png"
smiles_to_png(smiles, output_file)
'''


from rdkit import Chem
from rdkit.Chem import AllChem

def save_multiple_conformers_to_sdf(smiles, num_confs=100, output_file='output.sdf'):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDG()
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
    
    # Optimize each conformer
    results = AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=4)

    writer = Chem.SDWriter(output_file)
    for conf_id, (energy, _) in zip(conf_ids, results):
        mol.SetProp("_Name", f'Conformer_{conf_id}_Energy_{energy}')
        writer.write(mol, confId=conf_id)
        print(f"Conformer ID: {conf_id}, Energy: {energy}")
    writer.close()

'''
# Example usage:
smiles = "Cc1ccc(CNC(c2nccn2C)c2ccc(F)cc2)cc1C"  # Example SMILES
save_multiple_conformers_to_sdf(smiles, num_confs=100, output_file='mol_0.sdf')
'''

from rdkit import Chem
from rdkit.Chem import AllChem

def translate_molecule(sdf_file, output_file, dx, dy, dz):
    # read sdf files
    suppl = Chem.SDMolSupplier(sdf_file)
    writer = Chem.SDWriter(output_file)
    
    for mol in suppl:
        if mol is None:
            continue
        # get all coordinates from mol.
        conf = mol.GetConformer()
        for atom_idx in range(mol.GetNumAtoms()):
            pos = list(conf.GetAtomPosition(atom_idx))
            # move coords
            pos[0] += dx
            pos[1] += dy
            pos[2] += dz
            conf.SetAtomPosition(atom_idx, pos)
        # write in new files.
        writer.write(mol)
    writer.close()

'''
# Example usage:
translate_molecule('mol/mol_0.sdf', 'mol/output.sdf', -4, 15, 15)
'''



import os
def convert_pdb_to_sdf(input_folder):
    '''
    Leverage obabel to transform mol.pdb to sdf.
    '''
    failed_conversions = []  # Store names of files that failed to convert

    # Iterate over all files in the directory
    for filename in os.listdir(input_folder):
        # Check if the file has a .pdb extension
        if filename.endswith(".pdb"):
            pdb_file = os.path.join(input_folder, filename)
            sdf_file = os.path.join(input_folder, filename.replace(".pdb", ".sdf"))

            # Build the Open Babel command
            command = f"obabel -i pdb {pdb_file} -o sdf -O {sdf_file}"

            # Execute the command and check the return value
            result = os.system(command)
            if result != 0:
                # If the command failed, record the filename
                failed_conversions.append(filename)
                print(f"Failed to convert {pdb_file}")
            else:
                print(f"Converted {pdb_file} to {sdf_file}")

    return failed_conversions

'''
# Example usage:
input_folder = "source_data/Hariboss/mol/"  # Replace with your PDB folder path
failed_files = convert_pdb_to_sdf(input_folder)
'''


'''
Then we will display how to get ligand from .cif file and transform it into pdb.
I think its' mean just use obabel to add bonds.
'''

from rdkit.Chem import rdmolfiles
def save_mol_to_pdb(mol, filename):
    """
    Save an RDKit Mol object to a PDB format file.

    :param mol: RDKit Mol object.
    :param filename: Output filename (including path).
    """
    if mol is None:
        raise ValueError("Invalid molecule object.")

    # Write to PDB file
    with open(filename, 'w') as pdb_file:
        pdb_block = Chem.MolToPDBBlock(mol)
        pdb_file.write(pdb_block)


from Bio.PDB import MMCIFParser

def extract_ligand_info_from_cif(cif_file, ligand_name, chain_id="A" ):
    """
    Extract atom information for a specified ligand from a CIF file.

    Parameters:
    - cif_file: Path to the CIF file.
    - ligand_name: Name of the ligand.
    - chain_id: Chain identifier to search within. Default is "A".

    Returns:
    - mol: RDKit Mol object with 3D coordinates of the ligand atoms.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_file)[0]

    ligand_info = {}
    atoms = []
    for chain in structure:
        if chain.id!=chain_id:
            continue
        not_found_residue = True
        for residue in chain:
            #print(residue.resname)
            if residue.resname not in ['A', 'C', 'G', 'U', 'T']:
                if residue.resname == ligand_name and not_found_residue:
                    not_found_residue = False
                    for atom in residue.get_atoms():
                        atom_name = atom.get_id()  # Extract atom name
                        residue_name = residue.get_resname()  # Extract residue name
                        x, y, z = atom.get_coord()  # Extract coordinates
                        atoms.append((atom_name, residue_name, x, y, z))


    mol = Chem.RWMol()

    # Add atoms
    atom_indices = {}
    for idx, (atom_name, residue_name, x, y, z) in enumerate(atoms):
        atom_idx = mol.AddAtom(Chem.Atom(atom_name[0]))
        atom_indices[idx] = atom_idx

    # Add coordinates
    conf = Chem.Conformer(mol.GetNumAtoms())
    for idx, (atom_name, residue_name, x, y, z) in enumerate(atoms):
        conf.SetAtomPosition(atom_indices[idx], (float(x), float(y), float(z)))

    mol.AddConformer(conf)
    return mol

'''
Example Usage
for i in range(len(rna_list)):
    cif_file = "/home/raojh/datasets/HARIBOSS/cifs/"+rna_list[i]+".cif"
    mol = extract_ligand_info_from_cif(cif_file, smile_name[i], smile_chain[i] )
    save_mol_to_pdb(mol, 'source_data/Hariboss/mol/'+rna_list[i]+'_mol_'+str(i)+'.pdb')
'''


#read atom in a chain from cif and save it into pdb
from Bio.PDB import MMCIFParser, PDBIO, Select

def extract_rna_atoms_from_cif(cif_file, chain_id=-1):
    """
    Extract RNA atoms from a specified chain (or all chains) in a CIF file.

    Parameters:
    - cif_file: Path to the CIF file.
    - chain_id: Chain identifier to extract. If -1, extract RNA atoms from all chains.

    Returns:
    - rna_atoms: A list of RNA atom objects from the specified chain(s).
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_file)[0]

    rna_atoms = []
    rna_bases = {"A", "U", "C", "G"}  # The four RNA bases

    for chain in structure:
        if chain_id != -1 and chain.id != chain_id:
            continue

        for residue in chain:
            if residue.resname in rna_bases:
                for atom in residue.get_atoms():
                    rna_atoms.append(atom)

    return rna_atoms

def save_atoms_to_pdb_simple(atoms, output_file):
    """
    Save extracted RNA atoms to a PDB file (simplified format).

    Parameters:
    - atoms: List of RNA atom objects.
    - output_file: Path for the output PDB file.
    """
    with open(output_file, 'w') as pdb_file:
        for i, atom in enumerate(atoms, start=1):
            atom_line = (
                f"ATOM  {i:5d} {atom.get_name():>4} {atom.get_parent().get_resname():>3} "
                f"{atom.get_parent().get_parent().id:>1}{atom.get_parent().id[1]:>4}    "
                f"{atom.coord[0]:8.3f}{atom.coord[1]:8.3f}{atom.coord[2]:8.3f}"
                f"  1.00  0.00           {atom.element:>2}\n"
            )
            pdb_file.write(atom_line)
        pdb_file.write("END\n")

'''
Example Usage
cif_file = "/home/raojh/datasets/HARIBOSS/cifs/8hba.cif"
rna_atoms = extract_rna_atoms_from_cif(cif_file, chain_id=-1)
output_pdb_file = "rna_chain_A.pdb"
save_atoms_to_pdb_simple(rna_atoms, output_pdb_file)

# Note: extract_rna_atoms_from_cif may occasionally fail to read atoms. If this occurs,
# re-download the PDB file and re-parse as PDB instead.

'''


import os
import requests
def download_pdb_files(pdb_ids, output_folder, timeout=600):
    """
    Batch download PDB files for the given PDB IDs to the specified folder.
    If a PDB download fails, fall back to downloading the CIF file.

    Parameters:
    - pdb_ids (list of str): List of PDB IDs.
    - output_folder (str): Path to the directory for saving downloaded files.
    - timeout (int): Download timeout in seconds. Default is 600.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    pdb_base_url = "https://files.rcsb.org/download/"
    cif_base_url = "https://files.rcsb.org/download/"

    for pdb_id in pdb_ids:
        pdb_id = pdb_id.upper()
        pdb_url = f"{pdb_base_url}{pdb_id}.pdb"
        cif_url = f"{cif_base_url}{pdb_id}.cif"

        try:
            response = requests.get(pdb_url, timeout=timeout)

            if response.status_code == 200:
                file_path = os.path.join(output_folder, f"{pdb_id}.pdb")
                with open(file_path, "wb") as f:
                    f.write(response.content)
                print(f"{pdb_id}.pdb downloaded and saved to {file_path}")
            else:
                print(f"Download of {pdb_id}.pdb failed, status code: {response.status_code}")
                # Try downloading the CIF file instead
                download_cif_file(cif_url, pdb_id, output_folder, timeout)

        except requests.exceptions.Timeout:
            print(f"Download of {pdb_id}.pdb timed out. Attempting CIF download.")
            download_cif_file(cif_url, pdb_id, output_folder, timeout)
        except requests.exceptions.RequestException as e:
            print(f"Error downloading {pdb_id}.pdb: {e}. Attempting CIF download.")
            download_cif_file(cif_url, pdb_id, output_folder, timeout)

def download_cif_file(url, pdb_id, output_folder, timeout):
    """
    Attempt to download a CIF file.

    Parameters:
    - url (str): Download URL for the CIF file.
    - pdb_id (str): PDB ID used for naming the file.
    - output_folder (str): Path to the directory for saving downloaded files.
    - timeout (int): Download timeout in seconds.
    """
    try:
        response = requests.get(url, timeout=timeout)
        if response.status_code == 200:
            file_path = os.path.join(output_folder, f"{pdb_id}.cif")
            with open(file_path, "wb") as f:
                f.write(response.content)
            print(f"{pdb_id}.cif downloaded and saved to {file_path}")
        else:
            print(f"Download of {pdb_id}.cif failed, status code: {response.status_code}")
    except requests.exceptions.Timeout:
        print(f"Download of {pdb_id}.cif timed out. Skipping this file.")
    except requests.exceptions.RequestException as e:
        print(f"Error downloading {pdb_id}.cif: {e}")

'''
Example Usage
download_pdb_files(["5afi"], "source_data/Hariboss/pdb/")  # Downloads PDB by default; falls back to CIF if PDB is unavailable
'''

# After saving the simplified RNA molecule, the next step is to convert it into coordinate lists.
from Bio.PDB import PDBParser

def load_rna_bases_from_pdb(file_path, return_seq=False, chain_id=-1):
    """
    Load RNA bases from a PDB file and return the atom list for each base.

    Parameters:
    - file_path (str): Path to the PDB file.
    - return_seq (bool): If True, also return sequence information.
    - chain_id (str/int): Chain identifier to extract. If -1, extract RNA atoms from all chains.

    Returns:
    - base_atoms (dict): Atom coordinates per base, keyed by base ID (e.g., "17_G").
    - rna_sequences (list): Sequence information per chain (only if return_seq is True).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", file_path)

    base_atoms = {}
    rna_sequences = []  # Store sequence information per chain
    rna_bases = {"A", "U", "C", "G"}  # The four RNA bases
    index = 1  # Sequential index counter

    for model in structure:
        for chain in model:
            if chain_id != -1 and chain.id != chain_id:
                continue

            chain_seq = []
            for residue in chain:
                resname = residue.get_resname().strip()
                if resname in rna_bases:
                    base_id = f"{index}_{resname}"  # Construct base ID as index_resname
                    atoms = [atom.get_coord() for atom in residue]
                    base_atoms[base_id] = atoms
                    if return_seq:
                        chain_seq.append(resname)
                    index += 1

            if return_seq and chain_seq:
                rna_sequences.append("".join(chain_seq))

    if return_seq:
        return base_atoms, rna_sequences
    else:
        return base_atoms

import numpy as np
from Bio.PDB import PDBParser
import gemmi

def load_rna_bases(file_path, return_seq=False, chain_id=-1):
    """
    Load RNA base information from a PDB or CIF file, returning an atom coordinate matrix per base.

    Parameters:
    - file_path (str): Path to the PDB or CIF file.
    - return_seq (bool): If True, also return sequence information.
    - chain_id (str/int): Chain identifier to extract. If -1, extract RNA bases from all chains.

    Returns:
    - base_atoms (dict): Atom coordinate matrix (NumPy array) per base, keyed by base ID (e.g., "17_G").
    - rna_sequences (dict): RNA base sequence per chain (optional; keyed by chain ID).
    """
    def parse_pdb(file_path):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", file_path)

        base_atoms = {}
        rna_sequences = {}
        rna_bases = {"A", "U", "C", "G"}  # The four RNA bases
        index = 1  # Sequential index counter

        for model in structure:
            for chain in model:
                if chain_id != -1 and chain.id != chain_id:
                    continue

                chain_seq = []
                for residue in chain:
                    resname = residue.get_resname().strip()
                    if resname in rna_bases:
                        base_id = f"{index}_{resname}"  # Construct base ID as index_resname
                        coordinates = [
                            atom.get_coord() for atom in residue
                        ]
                        base_atoms[base_id] = np.array(coordinates)  # Convert to NumPy array
                        chain_seq.append(resname)
                        index += 1

                if return_seq and chain_seq:
                    rna_sequences[chain.id] = "".join(chain_seq)

        if return_seq:
            return base_atoms, rna_sequences
        else:
            return base_atoms

    def parse_cif(file_path):
        doc = gemmi.cif.read_file(file_path)
        block = doc[0]  # Use the first block by default
        structure = gemmi.make_structure_from_block(block)

        base_atoms = {}
        rna_sequences = {}
        rna_bases = {"A", "U", "C", "G"}  # The four RNA bases
        index = 1  # Sequential index counter

        for model in structure:
            for chain in model:
                if chain_id != -1 and chain.name != chain_id:
                    continue

                chain_seq = []
                for residue in chain:
                    if residue.name in rna_bases:
                        base_id = f"{index}_{residue.name}"  # Construct base ID as index_resname
                        coordinates = [
                            [atom.pos.x, atom.pos.y, atom.pos.z] for atom in residue
                        ]
                        base_atoms[base_id] = np.array(coordinates)  # Convert to NumPy array
                        chain_seq.append(residue.name)
                        index += 1

                if return_seq and chain_seq:
                    rna_sequences[chain.name] = "".join(chain_seq)
        if return_seq:
            return base_atoms, rna_sequences
        else:
            return base_atoms

    if file_path.endswith(".pdb"):
        return parse_pdb(file_path)
    elif file_path.endswith(".cif"):
        return parse_cif(file_path)
    else:
        raise ValueError("Unsupported file format. Only PDB and CIF formats are supported.")

# Corresponding utility: convert the base_atoms dict to a matrix representation.
def convert_base_atoms_to_matrix(base_atoms):
    # For each base, extract atom coordinates and convert to a matrix
    base_matrices = {}
    for base_key, atom_coords in base_atoms.items():
        # Convert atom_coords list to an n*3 matrix
        base_matrix = np.array(atom_coords)
        base_matrices[base_key] = base_matrix
    return base_matrices

# Corresponding function: compute the contact (distance) matrix.
import numpy as np
def calculate_min_distance_matrix(base_atoms, mol_coords):
    # Initialize the distance matrix
    num_bases = len(base_atoms)
    num_mol_atoms = len(mol_coords)
    distance_matrix = np.zeros((num_bases, num_mol_atoms))

    # Iterate over each RNA base
    for i, (base_key, base_coords) in enumerate(base_atoms.items()):
        # Iterate over each atom in the molecule
        for j, mol_coord in enumerate(mol_coords):
            # Compute the distance between the molecule atom and all atoms in the current base
            distances = np.linalg.norm(base_coords - mol_coord, axis=1)
            # Take the minimum distance
            min_distance = np.min(distances)
            distance_matrix[i, j] = min_distance

    return distance_matrix

'''
Example Usage

rna_seq_list = []
pairwise_list = []

for i in range(len(rna_list)):
    pdb_file = "source_data/Hariboss/rna/"+rna_list[i]+"_"+rna_chain[i]+".pdb"
    base_atoms, rna_seq = load_rna_bases_from_pdb(pdb_file, return_seq=True, chain_id=rna_chain[i])

    rna_seq_list.append(rna_seq[0])
    supplier = Chem.SDMolSupplier("source_data/Hariboss/mol/"+rna_list[i]+"_mol_"+str(i)+".sdf",sanitize=False)
    for mol in supplier:
        coords = []
        for atom in mol.GetAtoms():
            pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            coords.append({"atom_id": atom.GetIdx(), "position": (pos.x, pos.y, pos.z)})

    base_matrices = convert_base_atoms_to_matrix(base_atoms)
    coord_matrix = convert_atom_dict_to_matrix(coords)

    distance_matrix = calculate_min_distance_matrix(base_matrices, coord_matrix)
    
    pairwise_list.append(distance_matrix)
'''


def extract_atom_coordinates(pdb_file, atoms_to_extract = ['C4\''] ):
    '''
    Extract the atoms which you want from pdb files.
    '''
    #atoms_to_extract = ['C4\'', 'N1', 'P']
    #atoms_to_extract = ['C4\'']
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

    coordinates_list = []
    atom_list = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0] == ' ':
                    atoms = {atom.name: atom for atom in residue.get_atoms()}
                    coordinates = [atoms[atom_name].get_coord() for atom_name in atoms if atom_name in atoms_to_extract]
                    atom_name = [atoms[atom_name].name for atom_name in atoms if atom_name in atoms_to_extract]
                    atom_list.extend( [RNA_atom_id[i] for i in atom_name ] )
                    coordinates_list.extend(coordinates)
    return np.array(atom_list),np.array(coordinates_list)

'''
Example Usage
C4_coords = extract_atom_coordinates(pdb_file_path)
'''


import os

def convert_pdb_to_sdf(input_folder):
    '''
    Transform .pdb files into .sdf files in a specific folder.
    '''
    failed_conversions = []  # Store names of files that failed to convert

    # Iterate over all files in the directory
    for filename in os.listdir(input_folder):
        # Check if the file has a .pdb extension
        if filename.endswith(".pdb"):
            pdb_file = os.path.join(input_folder, filename)
            sdf_file = os.path.join(input_folder, filename.replace(".pdb", ".sdf"))

            # Build the Open Babel command
            command = f"obabel -i pdb {pdb_file} -o sdf -O {sdf_file} -d"   # -d is used to remove H atoms.

            # Execute the command and check the return value
            result = os.system(command)
            if result != 0:
                # If the command failed, record the filename
                failed_conversions.append(filename)
                print(f"Failed to convert {pdb_file}")
            else:
                print(f"Converted {pdb_file} to {sdf_file}")

    return failed_conversions

'''
Example Usage
input_folder = "source_data/Hariboss/mol/"  # Replace with your PDB folder path
failed_files = convert_pdb_to_sdf(input_folder)
'''

from Bio import pairwise2
def calc_identity(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]
    seq1_aligned, seq2_aligned = best_alignment[:2]
    matches = sum(res1 == res2 for res1, res2 in zip(seq1_aligned, seq2_aligned))
    identity = matches / max(len(seq1), len(seq2))
    return identity * 100

'''
Example Usage
seq1 = "ACGTACGTACGT"
seq2 = "ACGTTCGTACGT"
print(f"Sequence identity: {calc_identity(seq1, seq2):.2f}%")
'''


from Bio.PDB.PDBExceptions import PDBConstructionWarning  
# Ref: biopython 1.79, due to the limitation of prody==2.4.1, there usually has a limitation for biopython version.
from Bio.PDB import MMCIFParser
import warnings

def detect_chain_break(cif_path: str) -> bool:
    """check break in cif file.
    Args:
        cif_path: CIF file path
    Returns:
        bool: return True if check a break，else return False
    """
    with warnings.catch_warnings(record=True) as warn_log:
        warnings.simplefilter("always", PDBConstructionWarning)
        
        parser = MMCIFParser()
        structure = parser.get_structure('', cif_path)
        list(structure.get_chains())
        
    return any(
        "discontinuous" in str(w.message).lower()
        and issubclass(w.category, PDBConstructionWarning)
        for w in warn_log
    )


from Bio.PDB import MMCIFParser

def get_hetatm_residues(cif_file):
    '''
    give a cif file, and return a set of relevant hetatm residues.
    '''
    parser = MMCIFParser()
    structure = parser.get_structure('structure_id', cif_file)
    het_residues = set()
    
    # Take only the first model (0-indexed)
    model = structure[0]

    for chain in model:
        for residue in chain:
            # Determine whether the residue is a hetero-atom (HETATM)
            if residue.id[0].strip() != '':  # Non-empty hetero flag indicates HETATM
                het_residues.add(residue.resname.strip())
    
    return het_residues

def merge_fasta_files(fasta_files, output_filename):
    '''
    merge fasta files list to output_filename
    '''
    with open(output_filename, 'w') as output_file:
        for i, fasta_file in enumerate(fasta_files):
            with open(fasta_file, 'r') as f:
                lines = f.readlines()
                cleaned_lines = []
                for line in lines:
                    line = line.strip()
                    if line:
                        cleaned_lines.append(line)
                if i < len(fasta_files) - 1:
                    output_file.write("\n".join(cleaned_lines) + "\n")
                else:
                    output_file.write("\n".join(cleaned_lines))


import requests
from rdkit import Chem

def get_smiles_from_ccd(ccd_id: str) -> str | None:
    query = """
    query ($id: String!) {
      chem_comp(comp_id: $id) {
        pdbx_chem_comp_descriptor {
          type
          program
          descriptor
        }
      }
    }
    """

    resp = requests.post(
        "https://data.rcsb.org/graphql",
        json={"query": query, "variables": {"id": ccd_id.upper()}},
        timeout=30
    )
    resp.raise_for_status()
    data = resp.json()

    descs = data["data"]["chem_comp"]["pdbx_chem_comp_descriptor"]

    smiles = None

    # (1) OpenEye canonical SMILES (highest priority)
    for x in descs:
        if x.get("type") == "SMILES_CANONICAL" and x.get("program") == "OpenEye OEToolkits":
            smiles = x["descriptor"]
            break

    # (2) Fallback: any canonical SMILES
    if smiles is None:
        for x in descs:
            if x.get("type") == "SMILES_CANONICAL":
                smiles = x["descriptor"]
                break

    # (3) Fallback: generic SMILES
    if smiles is None:
        for x in descs:
            if x.get("type") == "SMILES":
                smiles = x["descriptor"]
                break

    if smiles is None:
        return None

    # Standardize to canonical SMILES (critical step)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return Chem.MolToSmiles(mol, canonical=True)
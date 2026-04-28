"""
Read a small molecule, generate multiple conformers, perform clustering and RMSD
alignment, filter conformers whose RMSD exceeds 1.5 A from the reference, and
save results into separate directories.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina
import numpy as np

# Read the SDF file
sdf_file = 'PC4.sdf'
suppl = Chem.SDMolSupplier(sdf_file)
mol = next(suppl)

# Generate 500 conformers
conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=500)

# Compute the pairwise RMSD matrix across all conformers
rmsd_matrix = []
for conf_id1 in conf_ids:
    row = []
    for conf_id2 in conf_ids:
        rmsd = AllChem.GetBestRMS(mol, mol, prbId=conf_id1, refId=conf_id2)
        row.append(rmsd)
    rmsd_matrix.append(row)


"""
Perform clustering.
"""
from scipy.cluster.hierarchy import linkage, fcluster
import numpy as np

# rmsd_matrix is a 2D array of pairwise RMSD values across all conformers
rmsd_matrix = np.array(rmsd_matrix)

# Perform hierarchical clustering
Z = linkage(rmsd_matrix, method='average')  # Additional linkage methods: 'single', 'complete', 'weighted', 'centroid', 'median', 'ward'

# Set the number of clusters
num_clusters = 5  # Adjust the number of clusters as needed

# Partition the dendrogram Z into the specified number of clusters
clusters = fcluster(Z, num_clusters, criterion='maxclust')

# Print the cluster assignment for each conformer
print("Cluster assignments:", clusters)


"""
Filter out conformers with excessive deviation.
"""

from rdkit import Chem
from rdkit.Chem import AllChem


# Retain conformers with RMSD <= 1.5 A
filtered_conf_ids = []
for conf_id in conf_ids:
    rmsd = AllChem.GetBestRMS(mol, mol, prbId=conf_id, refId=0)
    if rmsd <= 1.5:
        filtered_conf_ids.append(conf_id)

print("Filtered conformation IDs:", filtered_conf_ids)


"""
Save to separate directories.
"""

# Create five directories for storing molecules from each cluster
import os
num_clusters = 5
cluster_folders = ['Cluster_{}'.format(i+1) for i in range(num_clusters)]

# Create the directories
for folder in cluster_folders:
    os.makedirs(folder, exist_ok=True)

# Save molecules to their respective cluster directories
for molecule_id, cluster_id in enumerate(clusters):
    if molecule_id not in filtered_conf_ids:
        sdf_filename = 'conf_{}.sdf'.format(molecule_id)
        cluster_folder = cluster_folders[cluster_id - 1]  # cluster_id is 1-indexed, list indices are 0-indexed
        sdf_filepath = os.path.join(cluster_folder, sdf_filename)

        sdf_block = Chem.MolToMolBlock(mol, confId=molecule_id)

        with open(sdf_filepath, 'w') as f:
            f.write(sdf_block)

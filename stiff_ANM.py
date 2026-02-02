import mdtraj as md
import numpy as np
from prody import *
import matplotlib.pyplot as plt

# Load trajectory with MDAnalysis
structure = parsePDB("last_frame_pull2_high.pdb").select('index 1 to 30187 and name CA')

# Perform Normal Mode Analysis using Elastic Network Model (ENM)
anm = ANM('Protein ANM')
anm.buildHessian(structure)  # Create Hessian matrix (stiffness representation)
anm.calcModes()  # Calculate vibrational modes

# Compute stiffness as the mean of eigenvalues (higher = stiffer protein)
stiffness = np.mean(anm.getEigvals())
print(f'Protein Stiffness: {stiffness:.3f}')
np.save('stiff_last_pull2_high_protein_index.npy', stiffness)

fluctuations = calcSqFlucts(anm[:10])
np.save("last_pull2_high_fluctuations_index.npy", fluctuations)


# Load trajectory with MDAnalysis
structure = parsePDB("last_frame_pull2_high_replica2.pdb").select('resnum 1 to 30187 and name CA')

# Perform Normal Mode Analysis using Elastic Network Model (ENM)
anm = ANM('Protein ANM')
anm.buildHessian(structure)  # Create Hessian matrix (stiffness representation)
anm.calcModes()  # Calculate vibrational modes

# Compute stiffness as the mean of eigenvalues (higher = stiffer protein)
stiffness = np.mean(anm.getEigvals())
print(f'Protein Stiffness: {stiffness:.3f}')
np.save('stiff_last_pull2_high_protein_replica2_index.npy', stiffness)
fluctuations = calcSqFlucts(anm[:10])
np.save("last_pull2_high_fluctuations_replica2_index.npy", fluctuations)




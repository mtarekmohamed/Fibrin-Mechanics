import MDAnalysis as mda
import numpy as np
from MDAnalysisTests.datafiles import GRO, XTC
from MDAnalysis.analysis.dihedrals import Dihedral
import matplotlib.pyplot as plt
import numpy as np



u = mda.Universe('last_frame_pull1.gro')

# Assuming `u` is your Universe object


#Protein

selected_residues = u.residues[0:1913]

#Gamma1
#selected_residues = u.residues[696:956]

#Gamma2
#selected_residues = u.residues[1652:1913]

#beta1
#selected_residues = u.residues[314:575]

#beta2
#selected_residues = u.residues[1270:1531]


#ABC
#selected_residues_1 = u.residues[20:174]
#selected_residues_2 = u.residues[193:314]
#selected_residues_3 = u.residues[583:696]

#DEF
#selected_residues_1 = u.residues[976:1130]
#selected_residues_2 = u.residues[1149:1270]
#selected_residues_3 = u.residues[1539:1652]

#Central
#selected_residues_1 = u.residues[0:20]
#selected_residues_2 = u.residues[174:193]
#selected_residues_3 = u.residues[575:583]
#selected_residues_4 = u.residues[956:976]
#selected_residues_5 = u.residues[1130:1149]
#selected_residues_6 = u.residues[1531:1539]

# Concatenate the two residue selections
#selected_residues = selected_residues_1 + selected_residues_2 + selected_residues_3 +  selected_residues_4 + selected_residues_5 + selected_residues_6


# Concatenate the two residue selections
#selected_residues = selected_residues_1 + selected_residues_2 + selected_residues_3



rama_central = mda.analysis.dihedrals.Ramachandran(selected_residues)
rama_central.run()
angles_def = rama_central.results.angles  # This will have shape (n_frames, n_residues, 2)

# Alpha-helix regions (lower-left and upper-right quadrants)
alpha_helix_range_phi_left = (-120, -30)  # Lower-left quadrant: Phi for left-handed alpha-helix
alpha_helix_range_psi_left = (-60, -30)   # Lower-left quadrant: Psi for left-handed alpha-helix

alpha_helix_range_phi_right = (60, 90)  # Upper-right quadrant: Phi for right-handed alpha-helix
alpha_helix_range_psi_right = (0, 60)  # Upper-right quadrant: Psi for right-handed alpha-helix

# Beta-sheet region (upper-left quadrant)
beta_sheet_range_phi_upper_right = (-180, -45)  # Upper-left quadrant: Phi for beta-sheet
beta_sheet_range_psi_upper_right = (60, 180)    # Upper-left quadrant: Psi for beta-sheet


frame_ratios = []

# Loop over all frames and calculate the ratio for each frame
for frame in angles_def:
    alpha_helix_count = 0
    beta_sheet_count = 0

    for residue_angles in frame:
        phi, psi = residue_angles  # Unpack the phi and psi angles for this residue

        # Check if it belongs to alpha-helix (left-handed or right-handed)
        # Left-handed alpha-helix (lower-left quadrant)
        if alpha_helix_range_phi_left[0] <= phi <= alpha_helix_range_phi_left[1] and alpha_helix_range_psi_left[0] <= psi <= alpha_helix_range_psi_left[1]:
            alpha_helix_count += 1
        # Right-handed alpha-helix (upper-right quadrant)
        elif alpha_helix_range_phi_right[0] <= phi <= alpha_helix_range_phi_right[1] and alpha_helix_range_psi_right[0] <= psi <= alpha_helix_range_psi_right[1]:
            alpha_helix_count += 1

        # Check if it belongs to beta-sheet (upper-left quadrant)
        elif beta_sheet_range_phi_upper_right[0] <= phi <= beta_sheet_range_phi_upper_right[1] and beta_sheet_range_psi_upper_right[0] <= psi <= beta_sheet_range_psi_upper_right[1]:
            beta_sheet_count += 1

    # Calculate the ratio for the frame, avoid division by zero
    if beta_sheet_count == 0:  # To avoid division by zero
        frame_ratios.append(float('inf'))  # If there are no beta-sheet residues, we set the ratio to infinity
    else:
        frame_ratios.append(alpha_helix_count / beta_sheet_count)

# Plot the box plot of the alpha-helix to beta-sheet ratio
plt.figure(figsize=(6, 4))
plt.boxplot(frame_ratios, vert=True, patch_artist=True, boxprops=dict(facecolor="lightblue"))
plt.xlabel('Frames')
plt.ylabel('Alpha-helix to Beta-sheet Ratio')
plt.title('Alpha-helix to Beta-sheet Ratio Across Frames')

# Save the plot
plt.savefig('alpha_beta_ratio_boxplot.png')

# Save the frame ratios as a numpy array
np.save('pull1_box_protein_new', frame_ratios)

# Optionally, print the frame ratios
print(f"Alpha-helix to Beta-sheet ratios per frame: {frame_ratios}")


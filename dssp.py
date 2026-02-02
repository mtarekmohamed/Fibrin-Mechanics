import mdtraj as md

# Load trajectory and topology
traj = md.load("pull_strain_atom_indicies_1000.dcd", top="fibrin_solv_ions.gro")

# Save the last frame to a .gro file
dssp = md.compute_dssp(traj)
np.save(dssp, 'dssp')


# Get the secondary structure types: e.g., 'H' for alpha-helix, 'E' for beta-sheet
secondary_structures = dssp[1]  # The second element of the returned tuple contains the structure types

# Initialize counters for alpha-helix and beta-sheet
alpha_helix_count = 0
beta_sheet_count = 0

# Count occurrences of alpha-helix ('H') and beta-sheet ('E')
for frame in secondary_structures:
    for res in frame:
        if res == 'H':
            alpha_helix_count += 1
        elif res == 'E':
            beta_sheet_count += 1

# Plotting the bar chart for alpha-helix and beta-sheet counts
labels = ['Alpha-Helix (H)', 'Beta-Sheet (E)']
values = [alpha_helix_count, beta_sheet_count]

plt.figure(figsize=(8, 6))
plt.bar(labels, values, color=['blue', 'green'])
plt.xlabel('Secondary Structure Type')
plt.ylabel('Count')
plt.title('Alpha-Helix and Beta-Sheet Components')
plt.savefig('alphavsbeta')


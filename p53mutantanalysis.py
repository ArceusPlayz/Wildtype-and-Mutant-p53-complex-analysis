import math
import numpy as np
import pandas as pd
import matplotlib as plt 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
def read_file(filename):
    #reading the pdb file and getting the information regarding the atoms in the p53 
    atoms = []
    with open (filename, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atoms.append(line)
    return atoms
def extract_residue(atoms, chain, residue_num):
    #extracting the information about the atoms on a specific residue
    residue_atoms = []
    for atom_line in atoms:
        chain_id = atom_line[21]
        res_num = int(atom_line[22:26].strip())
        if chain_id == chain and res_num == residue_num:
            residue_atoms.append(atom_line)
    return(residue_atoms)
wildtype_file = r"C:\Users\Jason Vo\Desktop\Coding\3TS8.pdb"
atoms_wildtype = read_file(wildtype_file)
print("Wild-type residue 121:")
res_121 = extract_residue(atoms_wildtype, 'A', 121)
for line in res_121:
    print(line.strip())
print("\nWild-type residue 122:")
res_122 = extract_residue(atoms_wildtype, 'A', 122)
for line in res_122:
    print(line.strip())
#Retreiving 121 and 122 residues from mutant p53    
mutant_file = r"C:\Users\Jason Vo\Desktop\Coding\4MZR.pdb"
atoms_mutant = read_file(mutant_file)
print("Mutant-type residue 121:")
res_121mutant = extract_residue(atoms_mutant, 'A', 121)
for line in res_121mutant:
    print(line.strip())
res_122mutant = extract_residue(atoms_mutant, 'A', 122)
for line in res_122mutant:
    print(line.strip())
def get_coords(atom_line):
    #Extract x, y, z coordinates from a PDB ATOM line
    x = float(atom_line[30:38].strip())
    y = float(atom_line[38:46].strip())
    z = float(atom_line[46:54].strip())
    return (x, y, z)
def distance(coord1, coord2):
    #Calculate  distance between two 3D points
    x1, y1, z1 = coord1
    x2, y2, z2 = coord2
    return math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
def extract_dna_atoms(atoms):
    #Extracting DNA atoms from the pdb file
    dna_atoms = []
    for atom_line in atoms:
        chain_id = atom_line[21]
        if chain_id not in ['A', 'B', 'C', 'D']:
            dna_atoms.append(atom_line)
    return dna_atoms
dna_wildtype = extract_dna_atoms(atoms_wildtype)
dna_mutant = extract_dna_atoms(atoms_mutant)
print(f"Wild-type DNA atoms: {len(dna_wildtype)}")
print(f"Mutant DNA atoms: {len(dna_mutant)}")
def find_nearby_dna(protein_atoms, dna_atoms, distance_cutoff=5.0):
    """Find DNA atoms within distance_cutoff of any protein atom"""
    nearby = []
    for prot_line in protein_atoms:
        prot_coords = get_coords(prot_line)
        for dna_line in dna_atoms:
            dna_coords = get_coords(dna_line)
            dist = distance(prot_coords, dna_coords)
            if dist <= distance_cutoff:
                nearby.append({
                    'protein_atom': prot_line.strip(),
                    'dna_atom': dna_line.strip(),
                    'distance': dist
                })
    return nearby
# Wild-type: find DNA atoms near residues 121 and 122
nearby_dna_wildtype = find_nearby_dna(res_121 + res_122, dna_wildtype, distance_cutoff=5.0)
# Mutant: find DNA atoms near residues 121 and 122
nearby_dna_mutant = find_nearby_dna(res_121mutant + res_122mutant, dna_mutant, distance_cutoff=5.0)
print(f"Wild-type: {len(nearby_dna_wildtype)} protein-DNA contacts within 5 Angstroms")
print(f"Mutant: {len(nearby_dna_mutant)} protein-DNA contacts within 5 Angstroms")
print("\nFirst 5 wild-type contacts:")
for contact in nearby_dna_wildtype[:5]:
    print(f"  {contact['distance']:.2f} A")
def find_closest_dna(protein_atoms, dna_atoms):
    """Find the closest DNA atom to any protein atom"""
    min_dist = float('inf')
    closest_pair = None
    for prot_line in protein_atoms:
        prot_coords = get_coords(prot_line)
        for dna_line in dna_atoms:
            dna_coords = get_coords(dna_line)
            dist = distance(prot_coords, dna_coords)
            if dist < min_dist:
                min_dist = dist
                closest_pair = (prot_line.strip(), dna_line.strip())
    return min_dist, closest_pair
# Wild-type
closest_wt, pair_wt = find_closest_dna(res_121 + res_122, dna_wildtype)
print(f"Wild-type closest contact: {closest_wt:.2f} angstroms")
print(f"  Protein atom line:\n    {pair_wt[0]}")
print(f"  DNA atom line:\n    {pair_wt[1]}")
# Extract and print coords explicitly
prot_coords = get_coords(pair_wt[0])
dna_coords = get_coords(pair_wt[1])
print(f"  Protein coords: {prot_coords}")
print(f"  DNA coords: {dna_coords}")    
# Mutant
closest_mut, pair_mut = find_closest_dna(res_121mutant + res_122mutant, dna_mutant)
print(f"\nMutant closest contact: {closest_mut:.2f} angstroms")
print(f"  Protein atom line:\n    {pair_mut[0]}")
print(f"  DNA atom line:\n    {pair_mut[1]}")
# Extract and print coords explicitly
prot_coords = get_coords(pair_mut[0])
dna_coords = get_coords(pair_mut[1])
print(f"  Protein coords: {prot_coords}")
print(f"  DNA coords: {dna_coords}")
print(f"\nChange in closest contact: {closest_mut - closest_wt:.2f} angstroms")
# Find closest contact specifically for residue 121
closest_wt_121, pair_wt_121 = find_closest_dna(res_121, dna_wildtype)
closest_mut_121, pair_mut_121 = find_closest_dna(res_121mutant, dna_mutant)
print(f"Residue 121 wild-type closest: {closest_wt_121:.2f} angstroms")
print(f"Residue 121 mutant closest: {closest_mut_121:.2f} angstroms")
print(f"Change: {closest_mut_121 - closest_wt_121:.2f} angstroms")
# And residue 122
closest_wt_122, pair_wt_122 = find_closest_dna(res_122, dna_wildtype)
closest_mut_122, pair_mut_122 = find_closest_dna(res_122mutant, dna_mutant)
print(f"\nResidue 122 wild-type closest: {closest_wt_122:.2f} angstroms")
print(f"Residue 122 mutant closest: {closest_mut_122:.2f} angstroms")
print(f"Change: {closest_mut_122 - closest_wt_122:.2f} angstroms")
# Check residues around the mutation
for res_num in [119, 120, 121, 122, 123, 124]:
    res_wt = extract_residue(atoms_wildtype, 'A', res_num)
    res_mut = extract_residue(atoms_mutant, 'A', res_num)

    if res_wt and res_mut:
        closest_wt, _ = find_closest_dna(res_wt, dna_wildtype)
        closest_mut, _ = find_closest_dna(res_mut, dna_mutant)
        change = closest_mut - closest_wt

        print(f"Residue {res_num}: {closest_wt:.2f} → {closest_mut:.2f} angstroms (Δ = {change:+.2f})")
# Find actual DNA contacts (hydrogen bonding range)
def find_hbond_contacts(protein_atoms, dna_atoms, distance_cutoff=3.5):
    """Find DNA atoms within hydrogen bonding distance"""
    contacts = []
    for prot_line in protein_atoms:
        prot_coords = get_coords(prot_line)
        atom_name = prot_line[12:16].strip()
        for dna_line in dna_atoms:
            dna_coords = get_coords(dna_line)
            dist = distance(prot_coords, dna_coords)
            if dist <= distance_cutoff:
                contacts.append({
                    'protein_atom': atom_name,
                    'distance': dist,
                    'dna_atom': dna_line[12:16].strip()
                })
    return contacts
# Check for hydrogen bonding contacts
print("Wild-type residues 119-124: hydrogen bonding contacts with DNA (≤3.5 Å)")
for res_num in [119, 120, 121, 122, 123, 124]:
    res_wt = extract_residue(atoms_wildtype, 'A', res_num)
    contacts_wt = find_hbond_contacts(res_wt, dna_wildtype)
    if contacts_wt:
        print(f"  Residue {res_num}: {len(contacts_wt)} contacts")
        for c in contacts_wt[:2]:  # Show first 2
            print(f"    {c['protein_atom']} - {c['dna_atom']}: {c['distance']:.2f} Å")
    else:
        print(f"  Residue {res_num}: NO contacts")
print("\nMutant residues 119-124: hydrogen bonding contacts with DNA (≤3.5 Å)")
for res_num in [119, 120, 121, 122, 123, 124]:
    res_mut = extract_residue(atoms_mutant, 'A', res_num)
    contacts_mut = find_hbond_contacts(res_mut, dna_mutant)
    if contacts_mut:
        print(f"  Residue {res_num}: {len(contacts_mut)} contacts")
        for c in contacts_mut[:2]:  # Show first 2
            print(f"    {c['protein_atom']} - {c['dna_atom']}: {c['distance']:.2f} Å")
    else:
        print(f"  Residue {res_num}: NO contacts")
def find_contacting_residues(atoms, dna_atoms, distance_cutoff=5.0):
    """Find which residues contact DNA and at what distance"""
    residue_contacts = {}

    for atom_line in atoms:
        chain_id = atom_line[21]
        if chain_id != 'A':  # Only chain A (protein)
            continue

        res_num = int(atom_line[22:26].strip())
        prot_coords = get_coords(atom_line)

        # Skip if not in DNA-binding domain (roughly 102-292)
        if res_num < 102 or res_num > 292:
            continue

        for dna_line in dna_atoms:
            dna_coords = get_coords(dna_line)
            dist = distance(prot_coords, dna_coords)

            if dist <= distance_cutoff:
                if res_num not in residue_contacts:
                    residue_contacts[res_num] = {'min_dist': dist, 'count': 0}
                residue_contacts[res_num]['count'] += 1
                residue_contacts[res_num]['min_dist'] = min(residue_contacts[res_num]['min_dist'], dist)

    return residue_contacts
# Run on both structures
print("Scanning wild-type DNA-binding domain (residues 102-292)...")
contacts_wt = find_contacting_residues(atoms_wildtype, dna_wildtype)
print("Scanning mutant DNA-binding domain (residues 102-292)...")
contacts_mut = find_contacting_residues(atoms_mutant, dna_mutant)
print(f"\nWild-type: {len(contacts_wt)} residues contact DNA (≤5 Å)")
print(f"Mutant: {len(contacts_mut)} residues contact DNA (≤5 Å)")
# Show top 10 closest residues in wild-type
print("\nTop 10 closest residues in wild-type:")
sorted_wt = sorted(contacts_wt.items(), key=lambda x: x[1]['min_dist'])
for res_num, data in sorted_wt[:10]:
    print(f"  Residue {res_num}: {data['min_dist']:.2f} Å ({data['count']} atoms)")
print("Residue-by-residue comparison (wild-type → mutant):\n")
print("Residue | Wild-type | Mutant | Change")
print("--------|-----------|--------|--------")
all_residues = set(contacts_wt.keys()) | set(contacts_mut.keys())
changes = []
for res_num in sorted(all_residues):
    wt_dist = contacts_wt.get(res_num, {}).get('min_dist', None)
    mut_dist = contacts_mut.get(res_num, {}).get('min_dist', None)

    if wt_dist and mut_dist:
        change = mut_dist - wt_dist
        changes.append((res_num, wt_dist, mut_dist, change))
        print(f"{res_num:7d} | {wt_dist:9.2f} | {mut_dist:6.2f} | {change:+6.2f}")
# Find the biggest changes
print("\n\nBiggest IMPROVEMENTS (negative = closer to DNA):")
sorted_changes = sorted(changes, key=lambda x: x[3])
for res_num, wt, mut, change in sorted_changes[:5]:
    print(f"  Residue {res_num}: {change:+.2f} Å (now {mut:.2f} Å)")
print("\nBiggest DISRUPTIONS (positive = farther from DNA):")
for res_num, wt, mut, change in sorted(changes, key=lambda x: x[3], reverse=True)[:5]:
    print(f"  Residue {res_num}: {change:+.2f} Å (now {mut:.2f} Å)")
# Get full residue information for all 14 binding residues
def get_residue_centroid(atoms, res_num):
    """Calculate the center of mass of a residue"""
    coords = []
    for atom_line in atoms:
        chain_id = atom_line[21]
        res_n = int(atom_line[22:26].strip())
        if chain_id == 'A' and res_n == res_num:
            coords.append(get_coords(atom_line))

    if not coords:
        return None

    # Average position
    x = sum(c[0] for c in coords) / len(coords)
    y = sum(c[1] for c in coords) / len(coords)
    z = sum(c[2] for c in coords) / len(coords)
    return (x, y, z)
# Get centroids for all binding residues
binding_residues = sorted(contacts_wt.keys())
wt_positions = {}
mut_positions = {}
for res_num in binding_residues:
    wt_pos = get_residue_centroid(atoms_wildtype, res_num)
    mut_pos = get_residue_centroid(atoms_mutant, res_num)

    if wt_pos and mut_pos:
        wt_positions[res_num] = wt_pos
        mut_positions[res_num] = mut_pos
print("Residue centroids collected:")
print(f"  Wild-type: {len(wt_positions)} residues")
print(f"  Mutant: {len(mut_positions)} residues")
# Prepare data for plotting
res_nums = list(wt_positions.keys())
wt_x = [wt_positions[r][0] for r in res_nums]
wt_y = [wt_positions[r][1] for r in res_nums]
wt_z = [wt_positions[r][2] for r in res_nums]
mut_x = [mut_positions[r][0] for r in res_nums]
mut_y = [mut_positions[r][1] for r in res_nums]
mut_z = [mut_positions[r][2] for r in res_nums]
# Get distance changes for coloring
# Only keep residues that exist in BOTH structures
valid_residues = []
valid_wt_x, valid_wt_y, valid_wt_z = [], [], []
valid_mut_x, valid_mut_y, valid_mut_z = [], [], []
distance_changes = []
for res_num in res_nums:
    if res_num in contacts_wt and res_num in contacts_mut:
        wt_dna_dist = contacts_wt[res_num]['min_dist']
        mut_dna_dist = contacts_mut[res_num]['min_dist']
        change = mut_dna_dist - wt_dna_dist

        valid_residues.append(res_num)
        valid_wt_x.append(wt_positions[res_num][0])
        valid_wt_y.append(wt_positions[res_num][1])
        valid_wt_z.append(wt_positions[res_num][2])
        valid_mut_x.append(mut_positions[res_num][0])
        valid_mut_y.append(mut_positions[res_num][1])
        valid_mut_z.append(mut_positions[res_num][2])
        distance_changes.append(change)
# Use valid data instead
res_nums = valid_residues
wt_x, wt_y, wt_z = valid_wt_x, valid_wt_y, valid_wt_z
mut_x, mut_y, mut_z = valid_mut_x, valid_mut_y, valid_mut_z
# Extract CA atoms (alpha carbons) from both proteins
def extract_ca_atoms(atoms, chain='A'):
    """Extract only CA (alpha carbon) atoms from a specific chain"""
    ca_atoms = []
    for atom_line in atoms:
        atom_name = atom_line[12:16].strip()
        chain_id = atom_line[21]
        if atom_name == 'CA' and chain_id == chain:
            ca_atoms.append(atom_line)
    return ca_atoms
ca_wt = extract_ca_atoms(atoms_wildtype, 'A')
ca_mut = extract_ca_atoms(atoms_mutant, 'A')
print(f"Wild-type CA atoms: {len(ca_wt)}")
print(f"Mutant CA atoms: {len(ca_mut)}")
# Extract coordinates
wt_ca_coords = [get_coords(line) for line in ca_wt]
mut_ca_coords = [get_coords(line) for line in ca_mut]
def extract_phosphate_backbone(atoms):
    """Extract P (phosphate) atoms from DNA chains (not A/B/C/D)"""
    p_atoms = []
    for atom_line in atoms:
        atom_name = atom_line[12:16].strip()
        chain_id = atom_line[21]
        # DNA chains are everything NOT A, B, C, D
        if atom_name == 'P' and chain_id not in ['A', 'B', 'C', 'D']:
            p_atoms.append(atom_line)
    return p_atoms
p_wt = extract_phosphate_backbone(atoms_wildtype)
p_mut = extract_phosphate_backbone(atoms_mutant)
print(f"Wild-type DNA P atoms: {len(p_wt)}")
print(f"Mutant DNA P atoms: {len(p_mut)}")
# Extract coordinates
wt_p_coords = [get_coords(line) for line in p_wt]
mut_p_coords = [get_coords(line) for line in p_mut]
# Extract coordinates we already got
wt_ca_coords = [get_coords(line) for line in ca_wt]
mut_ca_coords = [get_coords(line) for line in ca_mut]
wt_p_coords = [get_coords(line) for line in p_wt]
mut_p_coords = [get_coords(line) for line in p_mut]
# Separate into x, y, z for plotting
wt_ca_x = [c[0] for c in wt_ca_coords]
wt_ca_y = [c[1] for c in wt_ca_coords]
wt_ca_z = [c[2] for c in wt_ca_coords]
mut_ca_x = [c[0] for c in mut_ca_coords]
mut_ca_y = [c[1] for c in mut_ca_coords]
mut_ca_z = [c[2] for c in mut_ca_coords]
wt_p_x = [c[0] for c in wt_p_coords]
wt_p_y = [c[1] for c in wt_p_coords]
wt_p_z = [c[2] for c in wt_p_coords]
mut_p_x = [c[0] for c in mut_p_coords]
mut_p_y = [c[1] for c in mut_p_coords]
mut_p_z = [c[2] for c in mut_p_coords]
# Get positions of residues 121 and 122
res_121_wt = get_residue_centroid(atoms_wildtype, 121)
res_121_mut = get_residue_centroid(atoms_mutant, 121)
res_122_wt = get_residue_centroid(atoms_wildtype, 122)
res_122_mut = get_residue_centroid(atoms_mutant, 122)
print(f"Residue 121: wt {res_121_wt} -> mut {res_121_mut}")
print(f"Residue 122: wt {res_122_wt} -> mut {res_122_mut}")
# Create figure
fig = plt.figure(figsize=(14, 10))
ax = fig.add_subplot(111, projection='3d')
# Plot protein and DNA backbones
ax.plot(wt_ca_x, wt_ca_y, wt_ca_z, color='lightblue', linewidth=2, label='Wild-type p53')
ax.plot(mut_ca_x, mut_ca_y, mut_ca_z, color='orange', linewidth=2, label='Mutant p53')
ax.plot(wt_p_x, wt_p_y, wt_p_z, color='darkblue', linewidth=2, label='Wild-type DNA')
ax.plot(mut_p_x, mut_p_y, mut_p_z, color='red', linewidth=2, label='Mutant DNA')
# Add arrows for residue 121 and 122 shifts
ax.quiver(res_121_wt[0], res_121_wt[1], res_121_wt[2],
          res_121_mut[0] - res_121_wt[0], 
          res_121_mut[1] - res_121_wt[1], 
          res_121_mut[2] - res_121_wt[2],
          color='purple', arrow_length_ratio=0.2, linewidth=2.5)
ax.quiver(res_122_wt[0], res_122_wt[1], res_122_wt[2],
          res_122_mut[0] - res_122_wt[0], 
          res_122_mut[1] - res_122_wt[1], 
          res_122_mut[2] - res_122_wt[2],
          color='darkgreen', arrow_length_ratio=0.2, linewidth=2.5)
# Plot binding residues color-coded by distance change
# Get mutant positions for binding residues
binding_mut_x = [mut_positions[r][0] for r in valid_residues]
binding_mut_y = [mut_positions[r][1] for r in valid_residues]
binding_mut_z = [mut_positions[r][2] for r in valid_residues]
# Create colormap: green (negative/improved) to red (positive/worsened)
norm = Normalize(vmin=min(distance_changes), vmax=max(distance_changes))
cmap = plt.cm.RdYlGn_r  # Red-Yellow-Green (reversed so green=good)
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
# Scatter plot binding residues
colors = cmap(norm(np.array(distance_changes)))
scatter = ax.scatter(binding_mut_x, binding_mut_y, binding_mut_z, 
                     c=distance_changes, cmap='RdYlGn_r', s=150, 
                     alpha=0.8, edgecolors='black', linewidth=1.5,
                     vmin=min(distance_changes), vmax=max(distance_changes))
# Add colorbar
cbar = plt.colorbar(scatter, ax=ax, pad=0.1, shrink=0.8)
cbar.set_label('Distance Change (Å)\n← Closer to DNA | Farther from DNA →', fontsize=10)
# Set view
ax.view_init(elev=20, azim=290)
# Labels and legend
ax.set_xlabel('X (Å)', fontsize=11)
ax.set_ylabel('Y (Å)', fontsize=11)
ax.set_zlabel('Z (Å)', fontsize=11)
ax.set_title('p53 DNA-Binding Domain: Mutation Effects on Binding Residues\nWild-type vs S121F/V122G Mutant', 
             fontsize=13, fontweight='bold')
# Create custom legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='lightblue', linewidth=2, label='Wild-type p53'),
    Line2D([0], [0], color='orange', linewidth=2, label='Mutant p53'),
    Line2D([0], [0], color='darkblue', linewidth=2, label='Wild-type DNA'),
    Line2D([0], [0], color='red', linewidth=2, label='Mutant DNA'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='purple', markersize=8, label='Residue 121 shift'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='darkgreen', markersize=8, label='Residue 122 shift'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=10, 
           markeredgecolor='black', markeredgewidth=1.5, label='DNA-binding residues')
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=10)
plt.tight_layout()
plt.savefig('p53_mutation_binding_analysis.png', dpi=300, bbox_inches='tight')
plt.show()
print("Plot saved as 'p53_mutation_binding_analysis.png'")
print(f"\nBinding residues colored by distance change:")
print(f"  Green: residues that got closer to DNA (improved contact)")
print(f"  Red: residues that got farther from DNA (lost contact)")
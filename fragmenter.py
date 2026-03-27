import argparse
from pymatgen.core import Structure, Molecule
from collections import deque
import numpy as np

def build_distance_graph(struct, max_bond=2.0):
    all_neighbors = struct.get_all_neighbors(r=max_bond)
    adjacency = {}
    for i, neighs in enumerate(all_neighbors):
        adjacency[i] = [(n.index, n.image) for n in neighs]
    return adjacency

def extract_fragment(mof_path, center_idx=-1, radius=6.0, output_path="fragment.xyz"):
    if radius > 4.5:
        print("\n===========================================================")
        print(f"WARNING: You are using a radius of {radius} A.")
        print("Using a radius > 4.5 A may extract disjoint or asymmetric")
        print("linkers from adjacent pores. For a clean symmetrical cluster")
        print("model (e.g. 3-Metal core), radius=4.0 A is recommended.")
        print("===========================================================\n")
        
    print(f"Loading '{mof_path}'...")
    struct = Structure.from_file(mof_path)
    
    metals = {"Mg", "Zn", "Cu", "Fe", "Co", "Ni", "Mn", "Zr", "Ti", "V", "Cr"}
    if center_idx == -1:
        for i, site in enumerate(struct):
            if site.species_string in metals:
                center_idx = i
                break
                
    center_site = struct[center_idx]
    neighbors = struct.get_neighbors(center_site, radius)
    print("Creating 3x3x3 supercell to handle periodic boundaries naturally...")
    supercell = struct * [3, 3, 3]
    
    center_frac = [0.5, 0.5, 0.5]
    best_dist = float('inf')
    sc_center_idx = -1
    for i, site in enumerate(supercell):
        if site.species_string in metals:
            dist, _ = supercell.lattice.get_distance_and_image(center_frac, site.frac_coords)
            if dist < best_dist:
                best_dist = dist
                sc_center_idx = i
                
    if sc_center_idx == -1:
        raise ValueError("No metal found in the structure to use as center!")
        
    sc_center_site = supercell[sc_center_idx]
    print(f"Center atom is {sc_center_site.species_string} at {sc_center_site.coords}")
    
    sc_neighbors = supercell.get_neighbors(sc_center_site, radius)
    initial_indices = {n.index for n in sc_neighbors}
    initial_indices.add(sc_center_idx)
    
    print("Building bond graph for supercell...")
    sc_all_neighbors = supercell.get_all_neighbors(r=3.0) 
    
    final_indices = set(initial_indices)
    queue = deque(initial_indices)
    visited = set(initial_indices)
    
    broken_bonds = []
    
    while queue:
        curr_idx = queue.popleft()
        curr_site = supercell[curr_idx]
        
        if curr_site.species_string in metals and curr_idx != sc_center_idx:
            continue
            
        for n in sc_all_neighbors[curr_idx]:
            n_idx = n.index
            n_species = n.species_string
            n_dist = getattr(n, 'nn_distance', supercell.lattice.get_distance_and_image(curr_site.frac_coords, n.frac_coords)[0])
            
            is_valid_bond = False
            if n_species in metals and n_dist < 2.6:
                is_valid_bond = True
            elif n_species not in metals and n_dist < 1.8:
                is_valid_bond = True
                
            if not is_valid_bond:
                continue
            
            if n_species in metals and n_idx not in initial_indices:
                broken_bonds.append((curr_idx, n))
            elif n_idx not in visited:
                final_indices.add(n_idx)
                visited.add(n_idx)
                queue.append(n_idx)

    print(f"Extracted {len(final_indices)} raw atoms from MOF.")
    
    sites = [supercell[idx] for idx in final_indices]
    species = [s.species_string for s in sites]
    coords = [s.coords for s in sites]
    
    bonds_by_ligand = {}
    for ligand_idx, missing_metal in broken_bonds:
        if ligand_idx not in bonds_by_ligand:
            bonds_by_ligand[ligand_idx] = []
        bonds_by_ligand[ligand_idx].append(missing_metal.coords)
        
    capped_count = 0
    for ligand_idx, missing_coords in bonds_by_ligand.items():
        ligand_site = supercell[ligand_idx]
        
        kept_bonds = 0
        for n in sc_all_neighbors[ligand_idx]:
            if n.index in final_indices:
                n_species2 = n.species_string
                n_dist2 = getattr(n, 'nn_distance', supercell.lattice.get_distance_and_image(ligand_site.frac_coords, n.frac_coords)[0])
                if n_species2 in metals and n_dist2 < 2.6:
                    kept_bonds += 1
                elif n_species2 not in metals and n_dist2 < 1.8:
                    kept_bonds += 1
                    
        if ligand_site.species_string == "O" and kept_bonds >= 2:
            continue
            
        capped_count += 1
        avg_vec = np.zeros(3)
        for coords_ in missing_coords:
            avg_vec += (coords_ - ligand_site.coords)
            
        dist = np.linalg.norm(avg_vec)
        if dist > 0:
            unit_vec = avg_vec / dist
            bond_length = 0.96
            if ligand_site.species_string == "C": bond_length = 1.09
            elif ligand_site.species_string == "N": bond_length = 1.01
            
            h_coords = ligand_site.coords + unit_vec * bond_length
            species.append("H")
            coords.append(h_coords)
            
    print(f"Capped {capped_count} exposed atoms with 1 Hydrogen each...")
    print(f"Final fragment size: {len(species)} atoms.")
    mol = Molecule(species, coords)
    
    print(f"Saving fragment to {output_path}...")
    mol.to(filename=output_path, fmt="xyz")
    print("Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract MOF fragment.")
    parser.add_argument("cif_path", help="Path to the MOF CIF file.")
    parser.add_argument("--center", type=int, default=-1, help="Index of the center metal atom (-1 to auto-detect).")
    parser.add_argument("--radius", type=float, default=6.0, help="Initial extraction radius in Angstroms.")
    parser.add_argument("--output", default="extracted_fragment.xyz", help="Output file path.")
    
    args = parser.parse_args()
    extract_fragment(args.cif_path, args.center, args.radius, args.output)

import argparse
from pymatgen.core import Structure, Molecule
from collections import deque
import numpy as np

METALS = {"Mg", "Zn", "Cu", "Fe", "Co", "Ni", "Mn", "Zr", "Ti", "V", "Cr"}
LARGE_NON_METALS = {"Br", "I", "S", "P", "Cl"}

def is_valid_bond(s1_str, s2_str, dist):
    if s1_str in METALS or s2_str in METALS:
        return dist < 2.6
    if "H" in (s1_str, s2_str):
        return dist < 1.2
    if s1_str in LARGE_NON_METALS or s2_str in LARGE_NON_METALS:
        return dist < 2.2
    return dist < 1.8

def extract_fragment(mof_path, center_idx=-1, radius=6.0, nmetals=3, output_path="fragment.xyz", minimize=False):
    print(f"Loading '{mof_path}'...")
    struct = Structure.from_file(mof_path)
    if nmetals < 1:
        raise ValueError(f"--nmetals must be >= 1 (got {nmetals}).")
    
    metals = METALS
    user_center = center_idx != -1
    if user_center:
        if center_idx < 0 or center_idx >= len(struct):
            raise IndexError(f"--center index {center_idx} is out of range [0, {len(struct)-1}]")
        if struct[center_idx].species_string not in metals:
            raise ValueError(
                f"--center index {center_idx} is {struct[center_idx].species_string}, not a metal in {sorted(metals)}"
            )
                
    print("Creating supercell...")
    # Dynamically size the supercell to avoid massive memory usage on giant MOFs (MIL-100)
    dims = [max(1, int(np.ceil(30.0 / a))) for a in struct.lattice.abc]
    # Ensure a 3x3x3 minimum ONLY for very small cells (<15 A) to prevent self-intersection
    dims = [max(3, d) if a < 15.0 else d for d, a in zip(dims, struct.lattice.abc)]
    supercell = struct * dims
    
    sc_center_cart = supercell.lattice.get_cartesian_coords([0.5, 0.5, 0.5])
    best_dist = float("inf")
    sc_center_idx = -1

    if not user_center:
        # Default behavior: auto-detect the most central metal in the supercell.
        for i, site in enumerate(supercell):
            if site.species_string not in metals:
                continue
            d = np.linalg.norm(site.coords - sc_center_cart)
            if d < best_dist:
                best_dist = d
                sc_center_idx = i
    else:
        # Explicit behavior: use the user-selected unit-cell site and pick its
        # nearest periodic image to the supercell center.
        target_site = struct[center_idx]
        target_frac = np.array(target_site.frac_coords)
        dims_arr = np.array(dims, dtype=float)
        for i, site in enumerate(supercell):
            if site.species_string != target_site.species_string:
                continue
            parent_frac = np.mod(np.array(site.frac_coords) * dims_arr, 1.0)
            frac_delta = np.abs(np.mod(parent_frac - target_frac + 0.5, 1.0) - 0.5)
            if np.max(frac_delta) > 1e-3:
                continue
            d = np.linalg.norm(site.coords - sc_center_cart)
            if d < best_dist:
                best_dist = d
                sc_center_idx = i
    
    if sc_center_idx == -1:
        if user_center:
            raise ValueError(f"Could not map --center index {center_idx} into the generated supercell.")
        raise ValueError("No metal found in the input structure.")
    sc_center_site = supercell[sc_center_idx]
    
    # 1. Topology Detection
    print("Detecting topology...")
    sbu_all_neighs = supercell.get_all_neighbors(r=3.6)
    sbu_metals = {sc_center_idx}
    sbu_queue = deque([sc_center_idx])
    is_infinite_sbu = False
    while sbu_queue:
        curr = sbu_queue.popleft()
        for n in sbu_all_neighs[curr]:
            if n.species_string in metals and n.index not in sbu_metals:
                sbu_metals.add(n.index)
                sbu_queue.append(n.index)
        if len(sbu_metals) > 20:
            is_infinite_sbu = True
            break
            
    if not is_infinite_sbu and len(sbu_metals) > 1:
        coords = np.array([supercell[i].coords for i in sbu_metals])
        dists = np.sqrt(((coords[:, None, :] - coords[None, :, :])**2).sum(axis=-1))
        if dists.max() > min(struct.lattice.abc) * 0.5:
            is_infinite_sbu = True

    # 2. Seeding & Core Selection
    # The 'Old Version' behavior for Path B: take all atoms in radius.
    sc_neighbors = supercell.get_neighbors(sc_center_site, radius)
    initial_indices = {sc_center_idx}
    for n in sc_neighbors:
        if not is_infinite_sbu:
            initial_indices.add(n.index)
        else:
            # For Path B, only add non-metals from the spatial cut
            if n.species_string not in metals:
                initial_indices.add(n.index)
    
    if is_infinite_sbu:
        # Path B: Select exactly nmetals closest to center
        print(f"  -> Path B (Infinite). Metals: {nmetals}")
        all_sc_m = []
        for i, site in enumerate(supercell):
            if site.species_string in metals:
                all_sc_m.append((np.linalg.norm(site.coords - sc_center_site.coords), i))
        all_sc_m.sort()
        core_metals = {idx for _, idx in all_sc_m[:nmetals]}
        if len(core_metals) < nmetals:
            raise ValueError(
                f"Requested --nmetals={nmetals}, but only found {len(core_metals)} metals in generated supercell."
            )
        # Ensure our selected nmetals are in the initial set
        initial_indices.update(core_metals)
        is_path_c = False
    elif len(sbu_metals) == 2 and not minimize:
        # Path C: Paddlewheel / Small SBU -> extract 2 SBUs to get the full linker between them
        print(f"  -> Path C (Discrete, 2 SBUs). Auto-detected small SBU size: {len(sbu_metals)}")
        
        c0 = np.mean([supercell[m].coords for m in sbu_metals], axis=0)
        q = deque([(list(sbu_metals)[0], 0)])
        visited = set(sbu_metals)
        # Use 3.6 A here so metal-metal SBU reconstruction is consistent with
        # topology detection; 3.0 A can miss adjacent metals in some PCNs.
        struct_all_neighs = supercell.get_all_neighbors(r=3.6)
        
        found_sbus = []  # entries: (sbu_metals_set, path_length, center_coords)
        
        while q:  # completely explore local organic connected components
            curr_idx, dist = q.popleft()
            curr_site = supercell[curr_idx]
            
            for n in struct_all_neighs[curr_idx]:
                n_idx = n.index
                n_dist = np.linalg.norm(n.coords - curr_site.coords)
                if not is_valid_bond(n.species_string, curr_site.species_string, n_dist):
                    continue
                
                if n.species_string in metals and n_idx not in sbu_metals:
                    m_comp = {n_idx}
                    mq = deque([n_idx])
                    mv = {n_idx}
                    while mq:
                        mc = mq.popleft()
                        for mn in struct_all_neighs[mc]:
                            if mn.species_string in metals and mn.index not in mv:
                                md = np.linalg.norm(mn.coords - supercell[mc].coords)
                                if md < 3.6:
                                    mv.add(mn.index)
                                    m_comp.add(mn.index)
                                    mq.append(mn.index)
                    if len(m_comp) == len(sbu_metals):
                        c1 = np.mean([supercell[m].coords for m in m_comp], axis=0)
                        # ensure not already found
                        if not any(np.linalg.norm(x[2] - c1) < 1.0 for x in found_sbus):
                            found_sbus.append((m_comp, dist + 1, c1))
                            
                elif n_idx not in visited:
                    visited.add(n_idx)
                    q.append((n_idx, dist + 1))

        # Fallback for Path C: if BFS discovery misses adjacent SBU candidates,
        # build metal connected-components directly and pick the nearest SBU
        # with matching metal count.
        if not found_sbus:
            all_metal_indices = [i for i, s in enumerate(supercell) if s.species_string in metals]
            seen_m = set()
            metal_components = []
            for m in all_metal_indices:
                if m in seen_m:
                    continue
                comp = {m}
                mq = deque([m])
                seen_m.add(m)
                while mq:
                    mc = mq.popleft()
                    for mn in struct_all_neighs[mc]:
                        if mn.species_string in metals and mn.index not in seen_m:
                            md = np.linalg.norm(mn.coords - supercell[mc].coords)
                            if md < 3.6:
                                seen_m.add(mn.index)
                                comp.add(mn.index)
                                mq.append(mn.index)
                metal_components.append(comp)

            c_ref = np.mean([supercell[m].coords for m in sbu_metals], axis=0)
            fallback_candidates = []
            for comp in metal_components:
                if len(comp) != len(sbu_metals):
                    continue
                c_comp = np.mean([supercell[m].coords for m in comp], axis=0)
                if np.linalg.norm(c_comp - c_ref) < 1.0:
                    continue
                fallback_candidates.append((np.linalg.norm(c_comp - c_ref), comp, c_comp))

            if fallback_candidates:
                fallback_candidates.sort(key=lambda x: x[0])
                best_dist, best_comp, best_center = fallback_candidates[0]
                print(f"     Fallback Path C: selected nearest matching SBU at {best_dist:.2f} A")
                found_sbus.append((best_comp, 1, best_center))
        
        # Pick the second SBU by testing multiple candidates and selecting
        # the one that yields the smallest final fragment.
        if found_sbus:
            c1 = np.mean([supercell[m].coords for m in sbu_metals], axis=0)
            found_sbus.sort(key=lambda x: np.linalg.norm(x[2] - c1))
            
            max_candidates = min(20, len(found_sbus))
            candidates_to_test = found_sbus[:max_candidates]
            print(f"     Testing {len(candidates_to_test)} SBU candidates to minimize atom count...")
            
            best_size = float('inf')
            best_result = None
            best_dist = float('inf')
            
            for idx, candidate in enumerate(candidates_to_test):
                test_core = set(sbu_metals) | candidate[0]
                test_init = set(initial_indices) | test_core
                sp, co = get_fragment(test_core, test_init, supercell, sc_center_idx, metals, is_infinite_sbu, nmetals, minimize)
                
                sz = len(sp)
                dist_val = np.linalg.norm(candidate[2] - c1)
                print(f"       Candidate {idx+1} (dist {dist_val:.2f} A) -> {sz} atoms")
                if sz < best_size or (sz == best_size and dist_val < best_dist):
                    best_size = sz
                    best_result = (sp, co)
                    best_dist = dist_val
                    
            print(f"     Selected SBU candidate yielding {best_size} atoms.")
            final_species, final_coords = best_result
        else:
            print("     Could not find adjacent SBU. Reverting to 1 SBU.")
            core_metals = set(sbu_metals)
            initial_indices.update(core_metals)
            final_species, final_coords = get_fragment(core_metals, initial_indices, supercell, sc_center_idx, metals, is_infinite_sbu, nmetals, minimize)
            
        is_path_c = True
        
    else:
        # Path A: Discrete SBU
        print(f"  -> Path A (Discrete). SBU size: {len(sbu_metals)}")
        core_metals = sbu_metals
        initial_indices.update(core_metals)
        is_path_c = False

    if not is_path_c:
        final_species, final_coords = get_fragment(core_metals, initial_indices, supercell, sc_center_idx, metals, is_infinite_sbu, nmetals, minimize)

    print(f"Final size: {len(final_species)} atoms.")
    mol = Molecule(final_species, final_coords)
    mol.to(filename=output_path, fmt="xyz")


def get_fragment(core_metals, initial_indices, supercell, sc_center_idx, metals, is_infinite_sbu, nmetals, minimize):
    # 3. Traversal
    sc_all_neighbors = supercell.get_all_neighbors(r=3.0)
    final_indices = set(initial_indices)
    queue = deque(initial_indices)
    visited = set(initial_indices)
    broken_bonds = []
    
    unwrapped_coords = {idx: supercell[idx].coords for idx in initial_indices}
    
    while queue:
        curr_idx = queue.popleft()
        curr_site = supercell[curr_idx]
        
        # Path B Traversal Limit: only traverse outward from absolute center
        if is_infinite_sbu and curr_site.species_string in metals:
            if curr_idx != sc_center_idx:
                continue
        # Path A / Path C Traversal Limit: only traverse outward from SBU metals
        if not is_infinite_sbu and curr_site.species_string in metals:
            if curr_idx not in core_metals:
                continue
            
        for n in sc_all_neighbors[curr_idx]:
            n_idx = n.index
            n_dist = np.linalg.norm(n.coords - curr_site.coords)
            
            is_bond = is_valid_bond(n.species_string, curr_site.species_string, n_dist)
            if not is_bond: continue
            
            if n.species_string in metals:
                if n_idx not in core_metals:
                    unwrapped_nb_pos = unwrapped_coords[curr_idx] + (n.coords - curr_site.coords)
                    broken_bonds.append((curr_idx, unwrapped_nb_pos))
                elif n_idx not in visited:
                    final_indices.add(n_idx); visited.add(n_idx); queue.append(n_idx)
                    unwrapped_coords[n_idx] = unwrapped_coords[curr_idx] + (n.coords - curr_site.coords)
            elif n_idx not in visited:
                final_indices.add(n_idx); visited.add(n_idx); queue.append(n_idx)
                unwrapped_coords[n_idx] = unwrapped_coords[curr_idx] + (n.coords - curr_site.coords)

    # 3.5 Coordination Completion for Edge Metals (Path B only)
    # The main BFS only traverses from center metal, so edge metals
    # miss linkers on their far side. Do a FULL BFS from each edge
    # metal, stopping only at external metal boundaries.
    if is_infinite_sbu and nmetals > 1:
        print("Completing coordination for edge metals...")
        edge_metals = core_metals - {sc_center_idx}
        comp_queue = deque(edge_metals)
        
        while comp_queue:
            idx = comp_queue.popleft()
            site = supercell[idx]
            
            # Don't traverse outward from metals that aren't in our core
            if site.species_string in metals and idx not in core_metals:
                continue
                
            for n in sc_all_neighbors[idx]:
                n_idx = n.index
                n_dist = np.linalg.norm(n.coords - site.coords)
                
                is_bond = is_valid_bond(n.species_string, site.species_string, n_dist)
                if not is_bond: continue
                
                if n.species_string in metals:
                    if n_idx not in core_metals:
                        unwrapped_nb_pos = unwrapped_coords[idx] + (n.coords - site.coords)
                        broken_bonds.append((idx, unwrapped_nb_pos))
                elif n_idx not in visited:
                    final_indices.add(n_idx)
                    visited.add(n_idx)
                    comp_queue.append(n_idx)
                    unwrapped_coords[n_idx] = unwrapped_coords[idx] + (n.coords - site.coords)

    # 3.6 Prune Partial Linkers
    # A "complete" linker bridges 2+ core metals. A "partial" linker only
    # touches 1 core metal (dangling chain). Remove partial linkers and
    # cap the metal with H.
    print("Pruning partial linkers...")
    
    # Build local bond adjacency among final_indices
    organic_indices = {idx for idx in final_indices if supercell[idx].species_string not in metals}
    local_adj = {idx: [] for idx in final_indices}
    for idx in final_indices:
        s1 = supercell[idx]
        for n in sc_all_neighbors[idx]:
            if n.index in final_indices:
                nd = np.linalg.norm(n.coords - s1.coords)
                is_b = is_valid_bond(s1.species_string, n.species_string, nd)
                if is_b:
                    local_adj[idx].append(n.index)
    
    # Find organic connected components (BFS ignoring metal nodes)
    org_visited = set()
    components = []
    for seed in organic_indices:
        if seed in org_visited: continue
        component = set([seed])
        org_visited.add(seed)
        q = deque([seed])
        while q:
            cur = q.popleft()
            for nb in local_adj[cur]:
                if nb not in org_visited and supercell[nb].species_string not in metals:
                    org_visited.add(nb)
                    component.add(nb)
                    q.append(nb)
        if component:
            components.append(component)
    
    # For each component, count how many distinct core metals it touches
    if minimize:
        # ==========================================
        # ADVANCED MINIMIZED PRUNING
        # ==========================================
        partial_to_remove = set()
        bridge_atoms_to_cap = []
        linker_evaluations = []

        for comp in components:
            touching_metals = set()
            bridge_atoms = set()
            for atom_idx in comp:
                for nb in local_adj[atom_idx]:
                    if nb in core_metals:
                        touching_metals.add(nb)
                        bridge_atoms.add(atom_idx)

            # Strict dangling-linker removal:
            # if a component touches <2 distinct core metals, keep only the
            # first-shell bridge atoms (typically O) and cap them.
            if len(touching_metals) < 2:
                atoms_to_cut = comp - bridge_atoms
                partial_to_remove.update(atoms_to_cut)
                for ba in bridge_atoms:
                    removed_neighbor_coords = []
                    for nb in local_adj[ba]:
                        if nb in atoms_to_cut:
                            removed_neighbor_coords.append(unwrapped_coords[nb])
                    if removed_neighbor_coords:
                        bridge_atoms_to_cap.append((ba, removed_neighbor_coords))
                continue

            keep_linker = False
            touched_coords = [supercell[m].coords for m in touching_metals]
            max_dist = 0
            for i in range(len(touched_coords)):
                for j in range(i+1, len(touched_coords)):
                    d = np.linalg.norm(touched_coords[i] - touched_coords[j])
                    if d > max_dist: max_dist = d
            if max_dist > 4.5:
                keep_linker = True

            linker_evaluations.append((comp, bridge_atoms, keep_linker))
            
        kept_any_full_linker = any(keep for _, _, keep in linker_evaluations)
        if not kept_any_full_linker and linker_evaluations:
            linker_evaluations.sort(key=lambda x: len(x[0]), reverse=True)
            best_comp, best_ba, _ = linker_evaluations[0]
            linker_evaluations[0] = (best_comp, best_ba, True)
        
        for comp, bridge_atoms, keep_linker in linker_evaluations:
            if not keep_linker:
                keep_atoms = set(bridge_atoms)
                q = deque([(ba, 0) for ba in bridge_atoms])
                while q:
                    curr, depth = q.popleft()
                    if depth < 5:
                        for nb in local_adj[curr]:
                            if nb in comp and nb not in keep_atoms:
                                keep_atoms.add(nb)
                                q.append((nb, depth + 1))
                
                # Trim dangling boundary atoms: if a kept C atom has only 1 
                # neighbor inside keep_atoms, it's a stub from a fused ring 
                # (e.g. naphthalene junction C). Remove it iteratively.
                changed = True
                while changed:
                    changed = False
                    to_remove = set()
                    for ka in keep_atoms - bridge_atoms:
                        if supercell[ka].species_string == "C":
                            internal_bonds = sum(1 for nb in local_adj[ka] if nb in keep_atoms)
                            if internal_bonds <= 1:
                                to_remove.add(ka)
                    if to_remove:
                        keep_atoms -= to_remove
                        changed = True
                                
                atoms_to_cut = comp - keep_atoms
                partial_to_remove.update(atoms_to_cut)
                
                for ka in keep_atoms:
                    removed_neighbor_coords = []
                    for nb in local_adj[ka]:
                        if nb in atoms_to_cut:
                            removed_neighbor_coords.append(unwrapped_coords[nb])
                    if removed_neighbor_coords:
                        bridge_atoms_to_cap.append((ka, removed_neighbor_coords))
    else:
        # ==========================================
        # NORMAL/ORIGINAL PRUNING
        # ==========================================
        partial_to_remove = set()
        bridge_atoms_to_cap = []
        for comp in components:
            touching_metals = set()
            bridge_atoms = set()
            for atom_idx in comp:
                for nb in local_adj[atom_idx]:
                    if nb in core_metals:
                        touching_metals.add(nb)
                        bridge_atoms.add(atom_idx)
            # Eliminate components connected to only one (or zero) core metals.
            # Keep first-shell bridge atoms and cap them.
            if len(touching_metals) < 2:
                atoms_to_cut = comp - bridge_atoms
                partial_to_remove.update(atoms_to_cut)
                for ba in bridge_atoms:
                    removed_neighbor_coords = []
                    for nb in local_adj[ba]:
                        if nb in atoms_to_cut:
                            removed_neighbor_coords.append(unwrapped_coords[nb])
                    if removed_neighbor_coords:
                        bridge_atoms_to_cap.append((ba, removed_neighbor_coords))
    
    if partial_to_remove:
        final_indices -= partial_to_remove
        broken_bonds = [(l, n) for l, n in broken_bonds if l not in partial_to_remove]

    # 4. Capping & Save
    ordered_indices = sorted(final_indices)
    sites = [supercell[idx] for idx in ordered_indices]
    species = [s.species_string for s in sites]
    coords = [unwrapped_coords[idx] for idx in ordered_indices]
    local_index = {sc_idx: i for i, sc_idx in enumerate(ordered_indices)}
    
    bonds_by_ligand = {}
    for lid, unwrapped_m_pos in broken_bonds:
        if lid not in bonds_by_ligand: bonds_by_ligand[lid] = []
        bonds_by_ligand[lid].append(unwrapped_m_pos)
    
    def cap_bond_length(site_species):
        if site_species == "C":
            return 1.09
        if site_species == "N":
            return 1.01
        return 0.96

    def _orthonormal_basis(u):
        trial = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(u, trial)) > 0.9:
            trial = np.array([0.0, 1.0, 0.0])
        e1 = np.cross(u, trial)
        n1 = np.linalg.norm(e1)
        if n1 < 1e-12:
            e1 = np.array([0.0, 0.0, 1.0])
            n1 = np.linalg.norm(e1)
        e1 = e1 / n1
        e2 = np.cross(u, e1)
        e2 = e2 / np.linalg.norm(e2)
        return e1, e2

    def _score_candidate_h(candidate_pos, parent_idx):
        min_h = float("inf")
        min_heavy = float("inf")
        for i, s in enumerate(species):
            d = np.linalg.norm(candidate_pos - np.array(coords[i]))
            if s == "H":
                if d < min_h:
                    min_h = d
            else:
                if i != parent_idx and d < min_heavy:
                    min_heavy = d
        return min_h, min_heavy

    def place_capping_h(parent_idx, base_vec, bl, min_hh=1.5, min_heavy=0.9):
        parent_pos = np.array(coords[parent_idx])
        vnorm = np.linalg.norm(base_vec)
        if vnorm < 1e-12:
            return
        u = base_vec / vnorm
        e1, e2 = _orthonormal_basis(u)

        theta_list = [0.0, 20.0, 35.0, 50.0, 65.0, 80.0, 110.0, 140.0, 170.0]
        phi_list = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
        directions = []
        for th in theta_list:
            th_r = np.deg2rad(th)
            ct, st = np.cos(th_r), np.sin(th_r)
            if th == 0.0:
                directions.append(u)
            else:
                for ph in phi_list:
                    ph_r = np.deg2rad(ph)
                    dir_vec = ct * u + st * (np.cos(ph_r) * e1 + np.sin(ph_r) * e2)
                    directions.append(dir_vec / np.linalg.norm(dir_vec))

        best_pos = None
        best_tuple = (-1.0, -1.0)
        for dvec in directions:
            cand = parent_pos + dvec * bl
            mh, mheavy = _score_candidate_h(cand, parent_idx)
            if mh >= min_hh and mheavy >= min_heavy:
                species.append("H")
                coords.append(cand)
                return
            rank = (mh, mheavy)
            if rank > best_tuple:
                best_tuple = rank
                best_pos = cand

        # Fallback: place best available direction rather than dropping H.
        if best_pos is not None:
            species.append("H")
            coords.append(best_pos)

    for ba_idx, removed_coords in bridge_atoms_to_cap:
        ba_site = supercell[ba_idx]
        ba_pos = np.array(unwrapped_coords[ba_idx])
        bl = cap_bond_length(ba_site.species_string)
        parent_local_idx = local_index.get(ba_idx)
        if parent_local_idx is None:
            # Parent atom may not be in ordered_indices if pruned unexpectedly.
            continue
        for rc in removed_coords:
            vec = rc - ba_pos
            place_capping_h(parent_local_idx, vec, bl, min_hh=1.5)
    
    already_capped = {ba_idx for ba_idx, _ in bridge_atoms_to_cap}
        
    for lid, missing_coords in bonds_by_ligand.items():
        if lid in already_capped: continue
        lsite = supercell[lid]
        lpos = unwrapped_coords[lid]
        kept = 0
        for n in sc_all_neighbors[lid]:
            if n.index in final_indices:
                d = np.linalg.norm(n.coords - lsite.coords)
                if is_valid_bond(n.species_string, lsite.species_string, d):
                    kept = kept + 1
        if lsite.species_string == "O" and kept >= 2: continue
        avg_vec = np.zeros(3)
        for c in missing_coords:
            avg_vec = avg_vec + (c - lpos)
        dist = np.linalg.norm(avg_vec)
        if dist > 0:
            bl = cap_bond_length(lsite.species_string)
            parent_local_idx = local_index.get(lid)
            if parent_local_idx is not None:
                place_capping_h(parent_local_idx, avg_vec, bl, min_hh=1.5)
            
    return species, coords

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cif_path")
    parser.add_argument("--center", type=int, default=-1)
    parser.add_argument("--radius", type=float, default=6.0)
    parser.add_argument("--nmetals", type=int, default=3)
    parser.add_argument("--output", default="extracted_fragment.xyz")
    parser.add_argument("--minimize", action="store_true", help="Enable advanced structure minimization by restricting partials cleanly")
    args = parser.parse_args()
    extract_fragment(args.cif_path, args.center, args.radius, args.nmetals, args.output, args.minimize)

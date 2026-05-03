import argparse
from collections import deque

import numpy as np
from pymatgen.core import Molecule, Structure

# Conservative covalent radii (Angstrom) for common COF elements.
COV_RAD = {
    "H": 0.31,
    "B": 0.84,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "Si": 1.11,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Br": 1.20,
    "I": 1.39,
}

CAP_BOND = {
    "C": 1.09,
    "N": 1.01,
    "O": 0.96,
    "B": 1.19,
    "Si": 1.48,
    "P": 1.42,
    "S": 1.34,
}


def _rad(sym):
    return COV_RAD.get(sym, 0.77)


def is_valid_bond(s1, s2, dist):
    if s1 == "H" and s2 == "H":
        return dist < 0.9
    cutoff = 1.25 * (_rad(s1) + _rad(s2))
    cutoff = min(2.2, max(1.1, cutoff))
    return dist <= cutoff


def make_supercell(struct, target_span=28.0):
    dims = [max(1, int(np.ceil(target_span / a))) for a in struct.lattice.abc]
    dims = [max(3, d) if a < 15.0 else d for d, a in zip(dims, struct.lattice.abc)]
    return struct * dims, dims


def pick_center_atom(supercell, center_idx=-1):
    if center_idx >= 0:
        if center_idx >= len(supercell):
            raise IndexError(f"--center {center_idx} out of range [0, {len(supercell)-1}]")
        return center_idx

    cell_center = supercell.lattice.get_cartesian_coords([0.5, 0.5, 0.5])
    best_i = -1
    best_d = 1e18
    for i, site in enumerate(supercell):
        # Prefer non-hydrogen atoms as anchors.
        if site.species_string == "H":
            continue
        d = np.linalg.norm(site.coords - cell_center)
        if d < best_d:
            best_d = d
            best_i = i
    if best_i < 0:
        raise ValueError("Could not find a valid center atom.")
    return best_i


def build_bond_graph(supercell, neighbor_r=2.4):
    all_neigh = supercell.get_all_neighbors(r=neighbor_r)
    graph = [[] for _ in range(len(supercell))]
    for i, neighs in enumerate(all_neigh):
        si = supercell[i].species_string
        ci = supercell[i].coords
        for n in neighs:
            j = n.index
            if j <= i:
                continue
            sj = n.species_string
            d = np.linalg.norm(ci - n.coords)
            if is_valid_bond(si, sj, d):
                graph[i].append(j)
                graph[j].append(i)
    return graph


def unwrap_coordinates(supercell, graph, anchor_idx):
    unwrapped = [None] * len(supercell)
    unwrapped[anchor_idx] = np.array(supercell[anchor_idx].coords)
    q = deque([anchor_idx])

    while q:
        i = q.popleft()
        ui = unwrapped[i]
        fi = np.array(supercell[i].frac_coords)
        for j in graph[i]:
            if unwrapped[j] is not None:
                continue
            fj = np.array(supercell[j].frac_coords)
            dfrac = fj - fi
            dfrac -= np.round(dfrac)
            uc = ui + supercell.lattice.get_cartesian_coords(dfrac)
            unwrapped[j] = uc
            q.append(j)

    # Fill any disconnected nodes with raw cart coords.
    for i in range(len(supercell)):
        if unwrapped[i] is None:
            unwrapped[i] = np.array(supercell[i].coords)
    return unwrapped


def oxygen_or_nitrogen_already_protonated(local_idx, species, coords, cutoff=1.2):
    if species[local_idx] not in {"O", "N"}:
        return False
    p = np.array(coords[local_idx])
    for i, s in enumerate(species):
        if s != "H":
            continue
        if np.linalg.norm(p - np.array(coords[i])) <= cutoff:
            return True
    return False


def place_capping_h(parent_idx, toward_missing, cap_len, species, coords, min_hh=1.5, min_heavy=0.9):
    parent = np.array(coords[parent_idx])
    v = np.array(toward_missing)
    n = np.linalg.norm(v)
    if n < 1e-12:
        return False
    # Place opposite to the missing bond direction.
    u = -v / n
    cand = parent + cap_len * u

    for i, s in enumerate(species):
        d = np.linalg.norm(cand - np.array(coords[i]))
        if s == "H" and d < min_hh:
            return False
        if i != parent_idx and s != "H" and d < min_heavy:
            return False

    species.append("H")
    coords.append(cand)
    return True


def extract_cof_fragment(cif_path, center_idx=-1, radius=6.0, output_path="cof_fragment.xyz"):
    print(f"Loading '{cif_path}'...")
    struct = Structure.from_file(cif_path)

    print("Creating supercell...")
    supercell, _ = make_supercell(struct)

    print("Building bond graph...")
    graph = build_bond_graph(supercell)

    sc_center_idx = pick_center_atom(supercell, center_idx)
    center_pos = np.array(supercell[sc_center_idx].coords)

    print("Unwrapping periodic coordinates...")
    unwrapped = unwrap_coordinates(supercell, graph, sc_center_idx)

    print("Selecting atoms by radius + connectivity...")
    sphere = {
        i for i in range(len(supercell))
        if np.linalg.norm(np.array(unwrapped[i]) - center_pos) <= radius
    }
    sphere.add(sc_center_idx)

    # Keep only connected part touching center inside sphere.
    keep = set([sc_center_idx])
    q = deque([sc_center_idx])
    while q:
        i = q.popleft()
        for j in graph[i]:
            if j in sphere and j not in keep:
                keep.add(j)
                q.append(j)

    # Collect broken bonds for capping.
    broken = []
    for i in keep:
        for j in graph[i]:
            if j not in keep:
                broken.append((i, np.array(unwrapped[j]) - np.array(unwrapped[i])))

    ordered = sorted(keep)
    local = {gi: li for li, gi in enumerate(ordered)}
    species = [supercell[i].species_string for i in ordered]
    coords = [np.array(unwrapped[i]) for i in ordered]

    print("Capping boundary bonds with H...")
    broken_by_parent = {}
    for gi, vec in broken:
        broken_by_parent.setdefault(gi, []).append(vec)

    for gi, vecs in broken_by_parent.items():
        li = local[gi]
        sym = species[li]
        if sym == "H":
            continue
        if sym in {"O", "N"} and oxygen_or_nitrogen_already_protonated(li, species, coords):
            continue
        cap_len = CAP_BOND.get(sym, 1.09)
        # One H per parent anchor in this base COF mode.
        avg = np.zeros(3)
        for v in vecs:
            avg += v
        place_capping_h(li, avg, cap_len, species, coords)

    print(f"Final size: {len(species)} atoms")
    mol = Molecule(species, coords)
    mol.to(filename=output_path, fmt="xyz")
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Base fragmenter for metal-free COFs")
    parser.add_argument("cif_path")
    parser.add_argument("--center", type=int, default=-1, help="Center atom index in generated supercell")
    parser.add_argument("--radius", type=float, default=6.0)
    parser.add_argument("--output", default="cof_fragment.xyz")
    args = parser.parse_args()

    extract_cof_fragment(
        cif_path=args.cif_path,
        center_idx=args.center,
        radius=args.radius,
        output_path=args.output,
    )

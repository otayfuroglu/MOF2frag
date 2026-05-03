# UniFrag

MOF2frag is a Python script designed to cleanly and robustly extract chemically-meaningful symmetric cluster models from Metal-Organic Framework (MOF) `.cif` files. It is specifically optimized for creating stable models for Density Functional Theory (DFT) calculations by using an intelligently-anchored Graph Breadth-First-Search.

## Features
- **Symmetric Cluster Extraction**: Targets a core metal cluster (e.g., an Mg3 or Ni3 node) and strictly completes all organic ligands bonded to those core metals. Disconnected or grazing organic linkers are gracefully ignored to prevent floating geometries.
- **Intact Organic Linkers**: The bond-graph traversal guarantees that your ligands are pulled completely symmetrically without slicing through covalent bonds.
- **Automatic Valence Capping**: Any carboxylate, phenolic, or bound Oxygen atom that loses at least one metal during extraction is automatically "capped" with an artificial Hydrogen atom placed precisely along the averaged bond vector at 0.96 Å. This creates standard and perfectly neutralized closed-shell clusters automatically!

## Installation
Dependencies:
```bash
pip install pymatgen numpy
```

## Usage
Simply run the Python script on your CIF file.
```bash
python fragmenter.py <path_to_cif> --radius 4.0 --output my_fragment.xyz
```

### Arguments
- `cif_path`: Path to your MOF `.cif` file.
- `--center`: The atomic index in the unit cell that acts as the physical origin. Leave it blank (or set to `-1`) and the script will automatically lock onto the most centralized metal atom to avoid periodic boundary wrapping.
- `--radius`: The spatial sphere (in Ångströms) defining the central "core metals". **Radius = 4.0 Å is highly recommended** for MOF-74 analogues to consistently isolate precisely 3 core metal atoms. *Warning: Radii larger than 4.5 Å will trigger a structural warning as they may disrupt molecular symmetries depending on your specific MOF lattice.*
- `--output`: Filepath to save the resulting `.xyz` fragment.

## Example
If you want to extract a 4-linker, 3-Metal core out of an MOF-74 variant:
```bash
python fragmenter.py example/MgMOF74_clean_fromCORE.cif --radius 4.0 --output Mg_cluster.xyz
```
This produces an `.xyz` trajectory containing 3 Mg atoms, 32 C atoms, 24 O atoms, and 23 capping/native Hydrogens!

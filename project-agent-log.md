# Project Agent Log

Chronological handoff log for agents working on UniFrag. Add newest entries at the top. Each entry should include changed files, validation, decisions, and follow-up risks.

## 2026-05-06 - Apply direct coffragmentor combine to normal/min ZnPc
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_cof_zn_pc_series/ZnPc-COF_frag_cof_min.xyz`
  - `test_on_cof_zn_pc_series/ZnPc-DPB_frag_cof_min.xyz`
- Summary:
  - Reverted the failed coffragmentor index/formula helper-selection machinery.
  - Added metallo-PC Path J for both normal and minimized fragments. It directly combines a Zn/N-rich `coffragmentor.py` node molecule with coffragmentor linker molecule image(s), including neighboring-cell images for normal fragments so all four node sides are represented, then duplicates the pair/set along the shortest lattice vector for the ZnPc dimer.
  - Normal ZnPc fragments remain on existing UniFrag Path D.
- Validation:
  - `python -m py_compile fragmentation_oop.py coffragmentor.py` passes.
  - ZnPc-DPB normal Path J: 322 atoms, `Zn2 B16 H80 C192 N16 O16`; minimized Path J: 166 atoms, `Zn2 B4 H32 C96 N16 O16`.
  - ZnPc-COF normal Path J: 210 atoms, `Zn2 B16 H48 C112 N16 O16`; minimized Path J: 138 atoms, `Zn2 B4 H24 C76 N16 O16`.
  - Regression minimized checks: COF-202 70 atoms, COF-300 74 atoms, COF-366 107 atoms.
- Follow-up risks:
  - Path J intentionally trusts coffragmentor node/linker molecules. Visual inspection should decide whether the combined node+linker pair needs an additional neighboring node/image or terminal capping in a later pass.

## 2026-05-04 - ZnPc minimized dimer linker balance
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_cof_zn_pc_series/ZnPc-COF_frag_cof_min.xyz`
- Summary:
  - Updated metallo-PC minimization from one linker globally to one BDBA linker per retained ZnPc layer in the dimer.
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes.
  - ZnPc normal: 242 atoms, `Zn2 B16 H64 C112 N16 O32`.
  - ZnPc minimized: 146 atoms, `Zn2 B4 H40 C76 N16 O8`; each ZnPc layer has one BDBA linker (`B-O = 8`, `O-H = 4`, `B-C = 4`).
  - Regression smoke checks: COF-366 Path D 182 atoms; COF-202 Path B 169 atoms; COF-300 Path C 149 atoms.
- Follow-up risks:
  - Visual inspection should confirm the selected one-linker-per-layer orientation is acceptable for stacked ZnPc variants.

## 2026-05-04 - ZnPc metallo-PC dimer support
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `test_on_cof_zn_pc_series/ZnPc-COF_frag_cof.xyz`
  - `test_on_cof_zn_pc_series/ZnPc-COF_frag_cof_min.xyz`
- Summary:
  - Updated ZnPc metallo-PC Path D to keep the nearest stacked ZnPc core along the short lattice axis, producing a two-layer dimer SBU.
  - Updated final component cleanup so metallo-PC mode preserves the two disconnected stacked principal layers.
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes.
  - ZnPc normal: Path D metallo-PC core dimer, 242 atoms, `Zn2 B16 H64 C112 N16 O32`.
  - ZnPc minimized: Path D metallo-PC core dimer, 146 atoms, `Zn2 B4 H40 C76 N16 O8`.
  - Regression smoke checks: COF-366 Path D 182 atoms; COF-202 Path B 169 atoms; COF-300 Path C 149 atoms.
- Follow-up risks:
  - Dimer selection assumes the stacked partner lies along the shortest lattice axis with about 2.5-4.5 A axial spacing and low perpendicular offset.

## 2026-05-04 - COF fragmentation path tests and agent memory split
- Changed files:
  - `fragmentation_oop.py`
  - `project-memory.md`
  - `project-decisions.md`
  - `project-agent-log.md`
  - `AGENTS.md`
  - COF test outputs under `test_on_cof_2xx_series/`, `test_on_cof_3xx_series/`, `test_on_cof_Por_series/`, and `test_on_cof_zn_pc_series/`
- Summary:
  - Added/validated COF-202 Path B tied-layer handling.
  - Validated COF-300 and COF-320 Path C behavior.
  - Validated COF-366 Path D porphyrin-core behavior.
  - Added ZnPc metallo-PC Path D priority and Zn radius.
  - Tuned ZnPc normal and minimized fragments: normal keeps all BDBA linkers; minimized keeps one full BDBA linker plus ZnPc fused benzene perimeter.
  - Split project coordination docs into `project-memory.md`, `project-decisions.md`, and `project-agent-log.md`.
- Validation:
  - `python -m py_compile fragmentation_oop.py` passes.
  - COF-202 normal: Path B layered set, 169 atoms.
  - COF-300 normal: Path C, 149 atoms.
  - COF-366 normal: Path D porphyrin core, 182 atoms.
  - ZnPc normal: Path D metallo-PC core, 121 atoms, `Zn1 B8 H32 C56 N8 O16`.
  - ZnPc minimized: Path D metallo-PC core, 73 atoms, `Zn1 B2 H20 C38 N8 O4`.
- Decisions made:
  - See `project-decisions.md` for Path B tied-neighbor and metallo-PC priority decisions.
- Follow-up risks:
  - RDKit UFF warns about Zn atom typing during H refinement; generation still succeeds.
  - Visual inspection remains important for new COF families because linker/SBU chemistry varies.

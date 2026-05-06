# Project Decisions

Durable implementation and architecture decisions for UniFrag. This file is the source of truth for decisions; keep entries concise, dated, and actionable.

## Decision 2026-05-04: Keep tied Path B layered neighbors
- Context: COF-202 visually missed one linker when Path B selected only one of two equally close B/O node components at 2.85 A.
- Decision: Path B keeps all same-spacing layered neighbor components within 0.05 A of the best spacing.
- Consequences: COF-202 normal fragment is 169 atoms (`Si2 B3 H75 C83 O6`); minimized fragment is 70 atoms (`Si2 B3 H33 C26 O6`). COF-102 and COF-300 smoke checks remain on Path A/Path C respectively.
- Alternatives considered: Keep only the first nearest layered component; rejected because tie ordering dropped a symmetry-equivalent linker.

## Decision 2026-05-04: Prioritize metallo-PC cores before B/O layered nodes
- Context: ZnPc-COF contains B/O linker nodes, so generic Path B ran before the porphyrin-like N-rich SBU detector and selected a layered B/O fragment.
- Decision: Add Zn bonding radius and detect N-rich metal macrocycles before generic B/O Path B when Zn/Cu/Fe/Co/Ni/Mn are present.
- Consequences: ZnPc-COF uses Path D (Metallo-PC core dimer). Normal keeps the two stacked ZnPc layers plus all connected benzene-1,4-diboronic acid linkers and cuts only far-side O-C bonds so far B-connected O atoms are H-capped: 242 atoms (`Zn2 B16 H64 C112 N16 O32`). Minimized keeps the two stacked ZnPc layers plus one full BDBA linker per layer, retains the ZnPc-core fused benzene perimeter, caps discarded linker attachment points with H, and uses the same far-side O-H termination: 146 atoms (`Zn2 B4 H40 C76 N16 O8`). COF-366 remains Path D porphyrin; COF-202 remains Path B layered set.
- Alternatives considered: Leave ZnPc as Path B; rejected because the intended SBU is the metallo-phthalocyanine core.

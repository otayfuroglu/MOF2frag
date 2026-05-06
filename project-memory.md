# Project Memory Template (Reusable)

Use this file as a persistent engineering memory for this project. Keep entries concise, date-stamped, and actionable.

> Coordination note: agents should also read `project-decisions.md` for durable decisions and `project-agent-log.md` for chronological handoffs.

## 0) Metadata

- Project name:
- Repository URL:
- Default branch:
- Primary language(s):
- Last updated:
- Maintainer(s):

## 1) Project Purpose

### 1.1 Goals
- 
- 

### 1.2 Non-goals
- 
- 

### 1.3 Current scope
- 

## 2) Architecture

### 2.1 Top-level structure
- `path/or/module`: purpose
- `path/or/module`: purpose

### 2.2 Core components
- Component/Class/Service:
  - Responsibility:
  - Inputs:
  - Outputs:
  - Key invariants:

### 2.3 Data flow (high-level)
1. 
2. 
3. 

## 3) Coding Conventions

### 3.1 Style and organization
- 
- 

### 3.2 Error handling
- 

### 3.3 Logging/observability
- 

### 3.4 Performance conventions
- 

## 4) Dependencies and Environment

### 4.1 Runtime dependencies
- Package: why needed
- Package: why needed

### 4.2 Tooling
- Formatter/linter:
- Test framework:
- Build system:

### 4.3 Environment assumptions
- Python/Node/Compiler version:
- OS notes:
- GPU/CPU notes:

## 5) Commands (Copy/Paste)

### 5.1 Setup
```bash
# install deps
```

### 5.2 Run
```bash
# run app/script
```

### 5.3 Test
```bash
# run tests
```

### 5.4 Lint/format
```bash
# lint/format
```

### 5.5 Release/deploy
```bash
# release/deploy
```

## 6) Development Workflow

### 6.1 Branching and PR flow
1. 
2. 
3. 

### 6.2 Change checklist
- [ ] Reproduce baseline behavior
- [ ] Add/adjust tests
- [ ] Validate key scenarios
- [ ] Update docs/memory

### 6.3 Review focus areas
- 
- 

## 7) Testing Strategy

### 7.1 Test tiers
- Unit:
- Integration:
- End-to-end/manual:

### 7.2 Critical regression cases
- Case name: expected outcome
- Case name: expected outcome

### 7.3 Test datasets/fixtures
- `path/to/fixtures`: purpose

## 8) Known Issues and Risks

For each item:
- ID:
- Symptom:
- Root cause (known/suspected):
- Affected area:
- Workaround:
- Status:

## 9) Decisions Log (ADR-lite)

### Decision YYYY-MM-DD: Title
- Context:
- Decision:
- Consequences:
- Alternatives considered:

### Decision 2026-05-04: Keep tied Path B layered neighbors
- Context: COF-202 visually missed one linker when Path B selected only one of two equally close B/O node components at 2.85 A.
- Decision: Path B now keeps all same-spacing layered neighbor components within 0.05 A of the best spacing.
- Consequences: COF-202 normal fragment grows from 124 to 169 atoms; minimized fragment grows from 58 to 70 atoms. COF-102 and COF-300 smoke checks remain on Path A/Path C respectively.
- Alternatives considered: Keep only the first nearest layered component; rejected because tie ordering dropped a symmetry-equivalent linker.

### Decision 2026-05-04: Prioritize metallo-PC cores before B/O layered nodes
- Context: ZnPc-COF contains B/O linker nodes, so the generic Path B detector ran before the porphyrin-like N-rich SBU detector and selected a layered B/O fragment.
- Decision: Add Zn bonding radius and detect N-rich metal macrocycles before generic B/O Path B when Zn/Cu/Fe/Co/Ni/Mn are present.
- Consequences: ZnPc-COF now uses Path D (Metallo-PC core dimer). Normal keeps the two stacked ZnPc layers plus all connected benzene-1,4-diboronic acid linkers and cuts only far-side O-C bonds so far B-connected O atoms are H-capped: 242 atoms, Zn2 B16 H64 C112 N16 O32. Minimized keeps the two stacked ZnPc layers plus one full BDBA linker per layer, retains the ZnPc-core fused benzene perimeter, caps discarded linker attachment points with H, and uses the same far-side O-H termination: 146 atoms, Zn2 B4 H40 C76 N16 O8. COF-366 remains Path D porphyrin; COF-202 remains Path B layered set.
- Alternatives considered: Leave ZnPc as Path B; rejected because the intended SBU is the metallo-phthalocyanine core.

### Decision 2026-05-06: Directly combine coffragmentor node and linker for metallo-PC minimum fragments
- Context: The failed helper-selection approach was undone. `coffragmentor.py` should provide the actual node/linker molecular fragments for ZnPc minimum fragments.
- Decision: Metallo-PC COFs now try Path J: combine one Zn/N-rich coffragmentor node molecule with coffragmentor linker molecule image(s), then duplicate the pair/set along the shortest lattice vector for the ZnPc dimer. Minimized mode keeps one nearest attached linker image; normal mode keeps all neighboring-cell linker images that chemically attach to the node.
- Consequences: ZnPc-DPB normal/min are 322/166 atoms; ZnPc-COF normal/min are 210/138 atoms. Non-metallo-PC COFs fall back to existing paths.
- Alternatives considered: coffragmentor-assisted UniFrag linker selection; reverted.

## 10) Operational Notes

### 10.1 Common failure modes
- 

### 10.2 Recovery steps
1. 
2. 

### 10.3 Permissions/secrets notes
- 

## 11) Migration / Reuse Notes

### 11.1 What to copy to a new project
- 

### 11.2 What must be customized
- Paths/directories
- Dependency versions
- Environment-specific commands
- Domain-specific heuristics

## 12) Quick Start for New Contributor

1. Read sections: 1, 2, 5, 7, 8.
2. Run setup and one smoke test.
3. Validate one known critical case.
4. Make small change and run regression checks.

## Appendix A) Project-Specific Filled Example (Optional)

Use this section only if you want this template to include a concrete reference snapshot.


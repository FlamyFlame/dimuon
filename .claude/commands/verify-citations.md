---
description: >
  Verify every \cite in the internal note resolves in a .bib AND that each cited
  reference actually supports the adjacent claim (ATLAS biblatex/INSPIRE). Use when
  the user says "check citations", "verify references", "are these citations real",
  or before freezing/compiling the note. Read-only. Implements Gate G2 of
  Analysis/docs/academic_writing_workflow.md.
---

Verify citations in the ATLAS internal note. **Read-only — modify nothing.**

Task / scope: $ARGUMENTS  (a section path, a file, or empty = whole note)

## Setup
Read:
- `Analysis/docs/academic_writing_workflow.md` (Gate G2).
- `IntNotes/ANA-HION-2023-07-INT1.tex` — find the `\addbibresource{...}` lines to
  get the bib resource list (currently `bib/ATLAS.bib`, `bib/CMS.bib`,
  `bib/ConfNotes.bib`, `bib/PubNotes.bib`, `ANA-HION-2023-07-INT1.bib`).

## Procedure
Spawn the **`citation-verifier`** agent (Agent tool, `subagent_type:
"general-purpose"`, instructed per `.claude/agents/citation-verifier.md`) with:
- `tex_files`: the in-scope `.tex` files (default: master + every `IntNotes/tex/*.tex`).
- `bib_files`: the resolved `\addbibresource` targets.

The agent does:
1. **Existence/hygiene (offline):** every cited key is defined in some bib;
   report UNDEFINED (CRITICAL), UNCITED user-bib keys (INFO), DUPLICATE/incomplete
   entries (WARNING).
2. **Support (claim↔reference):** for each cited reference adjacent to a claim,
   resolve identity (title/arXiv/DOI) and check the claim is supported — verdict
   VERIFIED / MINOR_DISTORTION / MAJOR_DISTORTION / UNVERIFIABLE / UNVERIFIED_ACCESS.
   Use arXiv / INSPIRE-HEP / DOI lookups where egress allows (`arxiv.org`,
   `inspirehep.net`, `doi.org`); otherwise mark `UNVERIFIED — needs manual`.
3. Flag suspected fabricated references (CRITICAL).

## Output
Relay the agent's Existence + Support tables and Summary. End with:
- **VERDICT: PASS** iff zero UNDEFINED, zero MAJOR_DISTORTION, zero UNVERIFIABLE,
  zero suspected-fabricated.
- **VERDICT: FAIL** otherwise — list exactly what to fix.
- Surface every `UNVERIFIED — needs manual` line for the user (these need a human
  to confirm against the source; they neither pass nor fail automatically).

## Notes
- The note uses biblatex+biber; do not assume bibtex `\bibliography{}`.
- `bib/*.bib` are ATLAS-curated INSPIRE exports — assume they exist; still
  support-check what they are cited for. Fabrication risk is concentrated in the
  user's own `ANA-HION-2023-07-INT1.bib`.

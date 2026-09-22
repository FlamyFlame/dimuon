---
name: academic-writing-gates
description: Apply the dimuon project's publication-grade G1-G7 workflow to ATLAS internal-note and paper writing, editing, citations, numbers, figures, placeholders, and compilation. Use for any IntNotes LaTeX section, paper prose, results description, citation verification, figure synchronization, or note build/freeze task.
---

# Apply academic writing gates

1. Read `Analysis/docs/academic_writing_workflow.md` and the active
   `Analysis/docs/tracking/academic_writing_workflow.md` fully. These override
   generic writing preferences.
2. Read the relevant read-only Claude workflow references:
   `.claude/commands/review-note.md`, `review-paper.md`,
   `verify-citations.md`, `sync-note-figures.md`, `check-note-sync.md`, and
   `compile-note.md`.
3. Load the scientific ground truth, current roadmap/status, placeholder index,
   confirmed metadata, target LaTeX section, master document, and bibliography.
4. Apply all relevant G1-G7 gates:
   - scope and claim traceability
   - citation existence and support
   - independent numerical re-derivation
   - figure provenance and synchronization
   - notation/document consistency
   - reproducible compilation
   - explicit placeholder and readiness disclosure
5. Never convert a roadmap guess, memory, preliminary plot, or placeholder into
   publication fact. Prefer a precise `\todo` or limitation over invented
   completion.
6. Verify citations against primary papers, official ATLAS material, or source
   documents; checking only that a BibTeX key exists is insufficient.
7. Compile when requested and inspect errors, warnings, references, figures, and
   rendered output. Preserve the IntNotes submodule's independent Git state.
8. Record the gate results and remaining blockers in the active tracking record.

Use independent reviewer agents only when the user explicitly requests
subagents; otherwise identify the pass as a main-thread gate review.

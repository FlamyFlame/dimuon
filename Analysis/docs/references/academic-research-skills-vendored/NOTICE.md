# Vendored subset — `academic-research-skills` (ARS)

**Upstream:** https://github.com/imbad0202/academic-research-skills (v3.12.1)
**Author:** Cheng-I Wu. **License:** CC-BY-NC-4.0 (see `LICENSE-CC-BY-NC-4.0`).

This directory holds a **curated, read-only subset** of ARS — only the protocol
reference documents we cherry-pick from (see
`../academic_research_skills_summary.md`). It is **not** the full plugin: the
Python verification scripts, tests, design docs, examples, and the four `SKILL.md`
agent suites are intentionally omitted (they require pip dependencies + live API
egress that this cluster does not provide, and they target generic
markdown/IMRAD/APA7 papers rather than ATLAS LaTeX notes).

We **do not enable ARS as a Claude Code plugin** (rationale: summary doc §Decision).
These files are kept verbatim for attribution and so the ported protocols remain
traceable to their source. Our own implementations live in
`Analysis/docs/academic_writing_workflow.md` and the `.claude/` commands/agents.

Files (all © Cheng-I Wu, CC-BY-NC-4.0):
- `protocols/anti_leakage_protocol.md` — session/materials over model memory.
- `protocols/claim_verification_protocol.md` — claim→source support verdicts.
- `protocols/ai_research_failure_modes.md` — 7-mode hallucination taxonomy (Lu 2026).
- `protocols/integrity_review_protocol.md` — integrity gate structure.
- `protocols/vlm_figure_verification.md` — figure-vs-data fidelity loop.
- `protocols/writing_quality_check.md` — AI-typical phrasing / prose-quality rules.

---
name: review-analysis-code
description: Perform physics-aware review of C++, ROOT, RDataFrame, histogram-filling, correction, selection, and numerical analysis code in the dimuon project. Use when writing, modifying, fixing, or reviewing analysis code, especially weights, cuts, triggers, efficiencies, normalization, year/sample combination, binning, and ROOT object handling.
---

# Review analysis code

1. Use `dimuon-context` and identify the authoritative Physics Procedure.
2. Read `.claude/commands/review-analysis-code.md`,
   `.claude/agents/code-reviewer.md`, `.claude/agents/physics-reviewer.md`, and
   relevant conventions as the latest read-only workflow references. Translate
   Claude tool/subagent syntax to Codex; never edit those files.
3. Inspect the focused diff and all callers/consumers affected by changed
   selections, weights, histograms, object names, or file formats.
4. Check:
   - physics agreement with the procedure and negative constraints
   - correct sample/year/trigger/centrality roles
   - units, normalization, luminosity/TAA factors, and combination formulas
   - OS/SS behavior, truth/reco matching, weights, edge bins, and null handling
   - ROOT ownership/lifetimes, type safety, RDataFrame column/filter ordering,
     histogram registration, and output naming compatibility
   - placeholder propagation and downstream rerun requirements
5. Compile with the repository's ACLiC pattern when ROOT is available. Run the
   smallest meaningful test and inspect output/log evidence.
6. If the user explicitly requested independent agents, use the Codex code and
   physics reviewer agents. Otherwise perform a self-review and label it as such.
7. Fix only when the user requested implementation. Record findings and
   validation in the active tracking document.

Report findings by severity with precise file/line evidence. Do not approve
merely because compilation succeeds.

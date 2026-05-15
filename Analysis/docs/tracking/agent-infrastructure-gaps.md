# Agent Infrastructure — Knowledge Gaps

**Objective**: Document domain knowledge that is needed for fully effective agent review but is not currently available in the repository documentation or auto-memory. These gaps should be filled as the analysis progresses.

## Gap 1: ATLAS publication plotting standards

**What's missing**: Specific ATLAS rules for publication-quality plots beyond what AtlasStyle.C provides:
- Official ATLAS color palette (specific RGB/hex values for data, MC, signal, background)
- Required luminosity label format and placement conventions
- Required sqrt(s) label format
- Ratio panel conventions (when required, Y-axis range, reference line)
- Canvas size requirements for journal submission
- Font size minimums for axis labels, legend text
- When to use "ATLAS Internal" vs "ATLAS Preliminary" vs "ATLAS" label

**Where this would go**: `.claude/conventions/atlas-plotting.md` (currently contains only codebase-derived rules) and `.claude/agents/plot-reviewer.md` (could add ATLAS-specific sub-checklist)

**How to fill**: Consult ATLAS PubCom guidelines or the ATLAS Style Guide TWiki page. Copy the relevant rules into the convention file.

## Gap 2: Systematic uncertainty framework

**What's missing**: This analysis's systematic uncertainty treatment is not documented:
- List of systematic uncertainty sources considered
- How each systematic is evaluated (variation, weight, alternative sample)
- Whether systematics are correlated or uncorrelated across years/centrality bins
- How systematics are propagated to final results
- Treatment in fits (nuisance parameters, priors)

**Where this would go**: `.claude/kb/analysis/` (new article) and `.claude/agents/physics-reviewer.md` (additional checklist items for systematic consistency)

**How to fill**: Document the systematic framework once it is established. The physics-reviewer currently has no systematics checklist items because implementing them without knowing the actual framework would produce incorrect review criteria.

## Gap 3: Background estimation methods

**What's missing**: The specific background estimation strategy for this analysis:
- How combinatorial background is estimated (same-sign subtraction? event mixing? template fit?)
- Signal region definition and any blinding protocol
- Background model validation procedure

**Where this would go**: `.claude/kb/analysis/overview.md` (expand "Signal and backgrounds" section)

**How to fill**: Document once the background estimation method is finalized.

## Gap 4: Fitting framework conventions

**What's missing**: Conventions for RooFit / template fitting in this analysis:
- Fit function choices and their physics motivation
- Fit range conventions
- Goodness-of-fit criteria
- Treatment of nuisance parameters
- When to use binned vs unbinned fits

**Where this would go**: `.claude/kb/procedures/` (new article) and `.claude/conventions/` (new fitting-conventions file)

**How to fill**: Document once the fitting framework is implemented.

## Gap 5: MC correction ordering

**What's missing**: The prescribed order for applying MC corrections:
- Trigger efficiency corrections
- Reconstruction efficiency corrections
- Acceptance corrections
- Resolution corrections / unfolding
- Whether corrections are applied multiplicatively or via reweighting

**Where this would go**: `.claude/kb/procedures/` (new article) and `.claude/agents/physics-reviewer.md` (add ordering checklist item)

**How to fill**: Document the correction chain once it is established. The physics-reviewer currently cannot check "corrections applied in the right order" because the order is not documented.

## Gap 6: Complete variable/branch reference for Analysis output trees

**What's missing**: Full branch listing for the muon-pair trees produced by NTuple processing (Stage 1 output):
- All branch names with types and units
- Which branches exist in which tree variants (data vs truth vs fullsim vs overlay)

**Where this would go**: `.claude/kb/data/variables.md` (expand)

**How to fill**: Generate from code inspection or `TTree::Print()` on a representative file for each tree type.

## Status

| Gap | Blocking? | Priority |
|-----|-----------|----------|
| 1. ATLAS plotting standards | No — current plot reviewer uses practical readability checks | Low (fill when preparing for publication) |
| 2. Systematics framework | No — physics reviewer covers existing known checks | Medium (fill when systematics work begins) |
| 3. Background estimation | No — investigation reviewer can still operate | Medium |
| 4. Fitting conventions | No — no fitting code exists yet | Low (fill when fitting starts) |
| 5. MC correction ordering | No — existing corrections are code-enforced | Medium |
| 6. Branch reference | No — key branches documented, full list is supplementary | Low |

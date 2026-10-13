---
name: review-plot
description: Review dimuon physics figures, histograms, fits, canvases, labels, and plot-producing code for numerical and publication quality. Use whenever a plot is created, modified, regenerated, looks unexpected, changes normalization or labels, or will enter slides, an internal note, or a paper.
---

# Review a physics plot

1. Load the relevant physics procedure and read
   `.claude/commands/review-plot.md`, `.claude/agents/plot-reviewer.md`, and
   `.claude/conventions/atlas-plotting.md` as read-only workflow references.
2. Identify the exact input ROOT file, histogram/object names, selection,
   correction stage, normalization, and output path.
3. Open every target image with visual inspection tooling. Do not review only
   the plotting source.
4. Check numerical sanity, empty/overflow bins, fit status, error bars,
   denominators, ratios, log scales, axis ranges, units, bin widths, and
   expected sample/year/centrality ordering.
5. Check labels against the implemented quantity: collision system, energy,
   luminosity, trigger, working point, correction stage, preliminary/placeholder
   status, ATLAS label, legend, and uncertainty meaning.
6. Compare against source histograms or a known-good prior output where useful.
7. If the user explicitly requested independent agents, use the plot reviewer.
   Otherwise state that the review is a main-thread self-review.
8. Record reviewed files, observed ranges, problems, fixes, and regenerated
   outputs in the active tracking document.

Reject a visually polished figure if its physics quantity, provenance, or
normalization is ambiguous.

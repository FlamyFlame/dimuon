---
name: review-investigation
description: Investigate discrepancies, surprising plots, mismatched yields, failed jobs, normalization shifts, trigger anomalies, sample differences, and unknown root causes in the dimuon analysis. Use when the user asks to debug, understand, compare, diagnose, or explain an unexpected analysis result.
---

# Investigate an unexpected result

1. Load authoritative context and read
   `.claude/commands/review-investigation.md` plus relevant reviewer prompts as
   read-only workflow references.
2. Create or resume an investigation tracking document before testing a
   hypothesis. State the observable symptom quantitatively and preserve the
   original output.
3. List competing hypotheses with discriminating tests. Include physics,
   configuration, sample, selection, normalization, stale-output, and technical
   failure modes.
4. Change one variable at a time. Prefer read-only inspection and minimal
   reproductions before altering production code.
5. Append every result immediately with exact paths, commands, event counts,
   ratios, timestamps, and relevant line references. Move disproven hypotheses
   to `Ruled out` with evidence.
6. Cross-check historical tracking documents before asserting a new cause.
7. Distinguish demonstrated facts, plausible explanations, and open questions.
8. Diagnose only unless the user also requested a fix. If a fix is authorized,
   update the tracking design and invoke the appropriate code/plot review skill.

Do not stop at the first plausible explanation. Require a test that separates it
from reasonable alternatives.

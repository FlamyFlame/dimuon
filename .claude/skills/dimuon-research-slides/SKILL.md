---
name: dimuon-research-slides
description: Create, revise, outline, export, or reference academic physics slides and presentations for the dimuon heavy-ion analysis project at /Users/yuhanguo/Documents/physics/heavy-ion/dimuon. Use whenever the user asks for slides, a slide deck, a presentation, a PPTX, Google Slides, prior-slide plot reuse, or meeting-slide work in this project — including group/analysis meetings, EW/Onia meetings, ATLAS approval or editorial-board meetings, and ZDC/detector talks. Enforces mandatory repository context routing through dimuon/.claude/CLAUDE.md and dimuon/AGENTS.md and canonical docs/tracking/knowledge sources, the project audience briefing, minimal research-slide style, robust plot- and prior-slide resolution under Cowork file mounts, outline-versus-direct-build decisions, dual delivery of deck + reusable source + Google Slides, and Obsidian indexing. This is the Claude port of the Codex skill of the same name and must stay behaviorally identical to it.
---

# Dimuon Research Slides (Claude)

Create concise, plot-led academic research decks for this analysis. This skill is the
Claude/Cowork port of the Codex `dimuon-research-slides` skill and must produce the same
decisions and outputs as Codex; only the tooling and path semantics differ (see
**Environment & paths**). When anything here would diverge from the Codex skill, treat the
Codex behavior as canonical and adapt the mechanism, not the outcome.

For PPTX generation use the **`pptx` skill** and its render-and-verify workflow. For
uploading to and inspecting Google Slides / Drive decks, use the connected Google Drive /
file connector.

## Environment & paths (Cowork specifics)

- `PROJECT_ROOT` (the folder the user opens) = `/Users/yuhanguo/Documents/physics/heavy-ion/dimuon`
- `REPO_ROOT` (the git analysis repo) = `/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/dimuon`
- Style kit lives at `PROJECT_ROOT`: `zdc_style.py`, `DECK_STYLE_GUIDE.md`,
  `zdc_calibration_build_reference.py`, `Guo_ZDC_calibration_BRAN.pptx`. A vendored copy of
  `zdc_style.py` also sits in this skill's `assets/`.
- Plots live at `PROJECT_ROOT/plots/` (and some under `PROJECT_ROOT/datasets/`).
- Obsidian vault = `/Users/yuhanguo/Documents/Obsidian Vault` (usually **not** mounted in a
  Cowork session — see **Obsidian bookkeeping** for how to gain access).

**File-mount path rule (critical).** Read/Write/Edit/Grep/Glob use the host paths above
(`/Users/yuhanguo/...`). The bash sandbox sees the same files under
`/sessions/<id>/mnt/...` (e.g. `PROJECT_ROOT` → `/sessions/<id>/mnt/dimuon`). The session id
changes every session. Therefore:

- When **building in-session**, resolve every figure to an absolute path that exists *in the
  current session mount*, and verify it exists (Read/`ls`) before placing it on a slide.
- Never hardcode a stale `/sessions/<old-id>/...` path (the reference build script's `FIG=`
  line is session-specific and will not resolve in a new session — re-point it).
- In the **reusable source you save for the user**, parameterize the figure directory at the
  top (e.g. `FIG = os.environ.get("FIG", "figs")`) and copy the actual plot PNGs into the
  deck's `figs/` folder so the script regenerates the deck on any machine without hunting for
  session paths.

## Required briefing

Before outlining, editing, or generating any deck, read all of
[references/briefing.md](references/briefing.md). Treat it as the standing project brief.
Apply the section matching the meeting type; do not silently substitute a generic corporate or
conference style. The briefing's editorial rules are authoritative for visual and language
style; `zdc_style.py` is the concrete python-pptx kit that implements that minimal style.

## Mandatory repository context preflight

Perform this preflight **before writing a per-slide outline and before directly generating a
detailed deck**. Dimuon is the default subject. The common exception is a **ZDC / detector
talk** (e.g. ZDC calibration, online monitoring): for those, skip the dimuon-analysis routing,
use the ZDC assets in `PROJECT_ROOT` (`zdc_*`, `Guo_ZDC_*`, `zdc_work_2024_slides/`), and still
apply the briefing's minimal style. If the subject is unclear, ask.

1. Read `REPO_ROOT/.claude/CLAUDE.md` first — it is updated most often and carries Tracking
   Documents, Documentation References, Auto-Dispatch Rules, and Active Tracking Docs. Then
   read `REPO_ROOT/AGENTS.md` for the source-authority precedence. Read the *current* files on
   every slide task, not a cached summary.
2. Read every currently **Active Tracking Doc** listed in `CLAUDE.md` fully when the task
   touches its topic. Do not start an outline or build until this is done.
3. Route context with `rg`/Grep from these sources:
   - `REPO_ROOT/Analysis/docs/analysis_overview.md` — stable physics motivation and methodology.
   - The relevant implementation tracking doc's `Physics Procedure` (or the academic workflow's
     `Workflow Specification`) — the highest-authority procedure.
   - `REPO_ROOT/Analysis/docs/tracking/analysis_roadmap_2026_06.md` and
     `analysis_status_summary.md` — current status, plus topic-specific tracking docs (including
     CLOSED ones, per the "Cross-reference other tracking docs" rule in `CLAUDE.md`).
   - `REPO_ROOT/Analysis/README.md`, pipeline docs, and code — implementation provenance.
   - `REPO_ROOT/.claude/kb/index.md` — select only the relevant KB article(s); do not load the
     whole KB.
   - `REPO_ROOT/IntNotes/analysis_metadata.md` and internal-note material — confirmed inputs and
     note context.
4. Apply the authority precedence in `AGENTS.md`/`CLAUDE.md`. Never treat conversation history,
   slide prose, filenames, or memory as stronger evidence than canonical repository sources.
   Re-read a value's canonical source rather than reciting it from memory.
5. Never modify `REPO_ROOT/.claude/` (Claude-owned context) as part of slide work, except the
   one Obsidian index update, which lives in the vault, not the repo.
6. Treat historical slides as historical evidence. If an old slide conflicts with current
   canonical context, preserve the old plot only when requested, label its historical/preliminary
   status, and do not repeat obsolete claims as current facts.

For each reused plot, establish a compact provenance record before writing its new caption:
source deck / date / slide, embedded image title or filename, original caption/notes, observable
and sample, relevant selections / corrections / status, canonical repository sources consulted,
and any historical-versus-current mismatch. Do not invent missing fields.

## Workflow

1. Identify the meeting type from the prompt: group/analysis, EW/Onia, ATLAS approval, or a
   ZDC/detector talk. If it is genuinely ambiguous and the choice changes content depth or index
   destination, ask the user.
2. **Ask clarifying questions when the prompt is generic or unclear.** If audience, scope, which
   results to feature, length, date, or meeting type are not pinned down, ask focused questions
   before doing work — do not guess. Then decide outline vs. direct build:
   - Generic prompt without per-slide content: return a **numbered per-slide outline first** and
     wait for explicit approval before generating files.
   - Detailed prompt with per-slide instructions (`Slide 1: ...`, `Slide 2: ...`): generate the
     deck immediately without a separate outline.
3. Stop and ask a focused question when a plot's meaning, caption, methodology, procedure, result
   interpretation, or requested slide role is unclear. Never infer scientific meaning from visual
   appearance or filenames alone.
4. Resolve plots and prior-slide context using the procedures below.
5. Build the deck. Generate **both** the `.pptx` and the reusable source code (see Deliverables),
   keeping inputs and paths legible enough for later regeneration.
6. **Render and verify.** Render the complete deck to per-slide images and inspect every slide for
   clipping, overlap, illegible labels, distorted plots, bad aspect ratios, and style drift.
   Iterate until clean. For a thorough pass, hand the rendered images to a sub-agent (Task tool)
   for inspection.
7. **Deliver** the `.pptx` + source locally, then upload/convert to Google Slides (see below).
8. After the final artifacts exist and pass inspection, update **exactly one** relevant Obsidian
   index. Do not index an outline or a failed/incomplete build.

## Build tooling & style kit

- Build with python-pptx via `from zdc_style import *` (helpers: `new_deck`/`title_slide`,
  `section_divider`, `header`, `bullets`, `stat`, `card`, `pic_fit`, `caption`, `notes`,
  `finish`). `zdc_calibration_build_reference.py` is the worked example; `DECK_STYLE_GUIDE.md`
  states the conventions. Copy `zdc_style.py` into the deck folder so the build is self-contained.
- Conventions: 16:9 (13.333×7.5"); title/divider slides on navy, content on white; one-line
  Cambria 34pt titles; Calibri 16–18pt body; full-width plot-led layout; captions 13pt italic
  muted; speaker notes on every slide. No accent bars/stripes/underlines. These agree with the
  briefing's minimal style; the briefing wins on any conflict.
- For routine meeting decks favor the briefing's plainest form (one dominant plot + short factual
  caption, minimal chrome). Reserve the richer presets (navy title, three-dot motif, stat cards)
  for polished standalone talks (interview, seminar).
- Output the deck into a dated deck folder under `PROJECT_ROOT`, matching the existing convention
  (`<topic>_<YYYY_MM_DD>_deck/` with `src/`, `figs/`, `output/`). Render QA artifacts as the
  existing decks do.

## Plot discovery and BNL paths

Use `PROJECT_ROOT = /Users/yuhanguo/Documents/physics/heavy-ion/dimuon`.

- Search `plots/` first; also search plot subdirectories under `datasets/` when relevant.
- Verify each chosen plot file exists in the current session mount before placing it (Read/`ls`).
- When the prompt gives a BNL/SDCC remote path: get the remote→local mapping from
  `~/bin/mount_sdcc_plots.sh` (read it via the host file tools, requesting access with
  `request_cowork_directory` if it is not mounted) and map the remote path to its configured
  local mount path — **do not guess by basename**. Verify the mapped local path is present and
  readable in the current session. If it is not mounted, ask the user to run
  `~/bin/mount_sdcc_plots.sh` (the sandbox cannot run it for them) or to mount/sync the
  directory, then re-check. Never fabricate the figure.
- Preserve plot aspect ratio (use `pic_fit`). Prefer the original vector or highest-resolution
  artifact. Do not crop away axis titles, legends, ATLAS labels, qualifiers, or uncertainty
  annotations.
- If multiple similarly named plots remain and their scientific roles are not unambiguous, stop
  and ask which one to use.

## Template and previous-deck lookup

- For group/analysis decks, inspect the most recent applicable analysis-meeting deck listed in
  the group/analysis index and use its layout as the template. Reproduce its restrained visual
  language; do not copy stale scientific content.
- When the user references previous slides, read the relevant Obsidian index **before** searching
  elsewhere, then open the linked deck through the Google Drive/Slides connector.
- If Drive folder listing or metadata lookup fails but direct deck links in the index remain
  readable, treat this as an access/discovery failure: use the direct indexed links plus
  title/date search; do not conclude that no deck exists from an empty folder listing.
- Report access failures distinctly from a clean "no matching deck found" result.

## Resolving a plot referenced by deck, date, and slide number

When the user says, e.g., "use the plot on slide 11 in my EW/Onia presentation on April 25, 2025":

1. Read the relevant Obsidian slide index and look for an exact meeting/date match and direct
   deck link.
2. If no exact match, search Drive/Slides with short title/date variants and inspect nearby
   indexed meetings. Treat nearby-date matches only as candidates, never as silent corrections.
3. Open the candidate deck and confirm its actual title, displayed meeting date, slide count, and
   meeting type.
4. Resolve the user's 1-based slide number to the actual slide, then read that exact slide.
5. Inspect the slide's text, speaker notes, and page elements. Record the embedded image's
   title/filename and source metadata when exposed. Use the exact embedded image or a verified
   source plot from the repo; do not recreate it from the caption.
6. Search the repository with distinctive identifiers from the slide (image filename, observable,
   sample, selection values, method names, caption phrases, plotter names). Route through
   `CLAUDE.md`/`AGENTS.md`, the relevant tracking docs, `analysis_overview.md`, pipeline docs/code,
   and `.claude/kb/index.md`.
7. Use repository context to explain the plot's meaning and purpose, preserving the distinction
   between what the historical slide establishes and what current repository sources establish.
8. Stop for confirmation if the requested date differs from the deck's actual date, more than one
   candidate remains, the slide number is absent, the image cannot be extracted reliably, or the
   slide's claim conflicts with canonical context.

## Deliverables

Always provide **both**:

- The final `.pptx` (saved under `PROJECT_ROOT` in the dated deck folder).
- The source code used to generate it, including any small local config/data file required for
  reproducibility, plus the `figs/` it references. Keep paths regenerable (see the file-mount
  path rule).

Start with a title slide unless the user explicitly requests insertion into an existing deck.
Use a white-background, plot-led layout from the briefing.

### Google Slides export

After the `.pptx` passes QA, upload it to Google Drive via the connected Drive/file connector and
convert it to a Google Slides presentation (Drive converts `.pptx` on import, or set the target
type to `application/vnd.google-apps.presentation`). Capture the resulting Google Slides URL for
the index. If no Drive connector is connected, save the `.pptx` locally, tell the user to import
it, and record the local path in the index instead.

## Obsidian bookkeeping

The vault is `/Users/yuhanguo/Documents/Obsidian Vault`. It is usually not mounted in Cowork —
request access with `request_cowork_directory` (path: the vault) before reading or writing an
index. Use these indexes inside the `Dimuon analysis/` folder:

- Group or analysis meeting: `Dimuon analysis/Group and analysis meeting slides index.md`
- EW/Onia meeting: `Dimuon analysis/EW Onia meeting slides index.md`
- ATLAS approval, editorial board, or group approval: `Dimuon analysis/ATLAS Approval Meeting Slides.md`

Add additional category indexes only when a genuinely new meeting type appears.

Before writing, read the target index and avoid duplicate dated entries. Update its `Updated:`
date when that field exists. Add the newest entry in the same ordering and local style as the
document, including:

- Date and meeting label/title.
- A working link to the final deck; prefer the Google Slides link if uploaded, otherwise the
  local PPTX path.
- A source-code link/path when the deck is local.
- One brief content/status summary bullet and, when useful, one methodology/result/versioning
  bullet.
- For reused plots, the provenance back-link (which dimuon source/tracking doc establishes the
  plot's meaning), so a later "reuse slide N from this deck" request can trace context.

When merely **referencing** prior slides, do not modify an index. When creating a genuinely new
final deck, indexing is part of completion.

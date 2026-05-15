# Plotting Conventions

## Style

- ATLAS style is available at `Analysis/AtlasStyle.C` — load with `SetAtlasStyle()` from `AtlasStyle.h`
- Shared plot helpers: `Analysis/Utilities/PlotUtils.h` (SetStyle, LogAx, adjustLogXRange)
- Variable-to-title and log-scale maps: `Analysis/Utilities/PlotCommonConfig.h`

## Output format

- Save plots as PNG only (`c->SaveAs("....png")`)
- Do NOT also save PDF unless explicitly requested

## Cross-section plots

- PbPb cross-section plots are always combined (all years together), never per individual year

## Differential scaling

- For differential plots, use `h->Scale(NORM_FACTOR, "width")` — never manual bin-width division

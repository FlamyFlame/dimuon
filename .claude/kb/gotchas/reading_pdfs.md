# Gotcha: reading PDFs on this cluster

**Verified 2026-06-14 (BNL SDCC).** How an agent must read a paper PDF.

| Approach | Works? |
|---|---|
| `WebFetch` on arXiv **abstract** page | ✗ returns the abstract only — no tables/cuts/numbers |
| `WebFetch` on arXiv **PDF URL** | ✗ undecoded binary — content not retrievable |
| **Read** tool on a local `.pdf` | ✗ needs `pdftoppm` (poppler) — not installed, can't install |
| `gs -sDEVICE=txtwrite` (ghostscript 9.54) | ✓ extracts text (minor ligature/spacing artifacts; tables extract poorly) |
| `curl` a PDF | ✓ outbound HTTPS works (no proxy) — but you still need `gs` to read it |

**Do this to read a local PDF:**
```bash
gs -q -dNOPAUSE -dBATCH -sDEVICE=txtwrite \
   -dFirstPage=N -dLastPage=M -sOutputFile=out.txt input.pdf
```
For a remote PDF with no local copy: `curl -sSL <url> -o /tmp/p.pdf` then `gs`.

**Implications for the KB** (see `KB_BUILDING_GUIDE.md` §3):
- A URL is **not** runtime-retrievable for deep content → KB summaries must be
  self-sufficient, and PRIMARY references must commit the PDF in-repo.
- Capture key tables/numbers **as text in the summary** — `gs` is weakest exactly
  on tables, which is where the numbers are.
</content>

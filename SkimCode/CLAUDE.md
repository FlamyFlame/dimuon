# Notes for Claude sessions in SkimCode/

See `README.md` for project docs. Assistant-specific hints:

- **Local test first.** `athena.py TrigRates_CA.py --evtMax=20` on a downloaded AOD before any grid submission.
- **Ask before state-changing ops.** Grid submissions, Rucio rule edits, `pbook` actions: confirm with user first.
- **Bump `.v<N>.`** in `grid_sub.sh` before every resubmit — PanDA task names are immutable.
- **Don't pipe `source`.** `source foo.sh | tail` runs source in a subshell; env vars vanish. Redirect instead: `source foo.sh > /tmp/log 2>&1`. Same for `asetup`.

# Setting up the dimuon repo on a new machine

This repo is hosted publicly at `github.com:FlamyFlame/dimuon.git`. Two
internal-note directories are **git submodules** pointing to private CERN
GitLab repositories:

| Path | Submodule URL | Visibility |
|---|---|---|
| `IntNotes/` | `ssh://git@gitlab.cern.ch:7999/atlas-physics-office/HION/ANA-HION-2023-07/ANA-HION-2023-07-INT1.git` | ATLAS internal (CERN GitLab access) |
| `IntNotesRun2DimuonReference/` | `ssh://git@gitlab.cern.ch:7999/yuhang/intnotesrun2dimuonreference.git` | Private (CERN GitLab) |

Cloning the parent repo without CERN credentials leaves those directories
empty — that is the intended security boundary, not an error.

---

## 1. Prerequisites

You need:
- A CERN account with access to `gitlab.cern.ch` (https://gitlab.cern.ch).
- An SSH key registered with both `github.com` and `gitlab.cern.ch`.

### Generate / register SSH key (skip if already done on this machine)

```bash
# Generate a key if you don't have one
ssh-keygen -t ed25519 -C "your-email@cern.ch"

# Print the public key and paste it into BOTH:
#   https://github.com/settings/keys
#   https://gitlab.cern.ch/-/profile/keys
cat ~/.ssh/id_ed25519.pub
```

### CERN GitLab uses a non-standard SSH port (7999)

The submodule URLs already include `:7999`. If you're behind a firewall that
blocks 7999 (some sites do), add this to `~/.ssh/config`:

```
Host gitlab.cern.ch
    Hostname gitlab.cern.ch
    Port 7999
    User git
    IdentityFile ~/.ssh/id_ed25519
```

Test it:

```bash
ssh -T -p 7999 git@gitlab.cern.ch
# Should print:  Welcome to GitLab, @yourusername!
```

---

## 2. Clone the repo with submodules

Two equivalent ways:

```bash
# One-shot
git clone --recurse-submodules git@github.com:FlamyFlame/dimuon.git
cd dimuon

# Or after a plain clone
git clone git@github.com:FlamyFlame/dimuon.git
cd dimuon
git submodule update --init --recursive
```

After this you should see content under `IntNotes/` and
`IntNotesRun2DimuonReference/`.

---

## 3. Verify

```bash
git submodule status
```

Expected output (SHAs will differ):

```
 a1b2c3d... IntNotes (heads/master)
 e4f5g6h... IntNotesRun2DimuonReference (heads/master)
```

A leading `-` means the submodule isn't initialized; `+` means the working
tree has moved past the recorded SHA.

---

## 4. BNL SDCC cluster notes

On `pc.sdcc.bnl.gov` / `usatlas` machines, work in
`/usatlas/u/<user>/workarea/` or wherever you keep your code area.

Standard setup before cloning:

```bash
# (existing) ROOT environment
source /usatlas/u/yuhanguo/setup.sh
```

The CERN GitLab SSH key needs to exist on the cluster too. Either copy it
from your laptop (`scp ~/.ssh/id_ed25519* user@pc.sdcc.bnl.gov:~/.ssh/`) or
generate a new one on the cluster and register it on both
`gitlab.cern.ch/-/profile/keys` and `github.com/settings/keys`.

---

## 5. Working with the submodules day-to-day

A submodule is a separate git repo nested inside the parent. Two things to
remember:

1. **Edits to files under `IntNotes/` or `IntNotesRun2DimuonReference/`**
   are committed and pushed *inside that submodule* — to CERN GitLab, not
   to GitHub:

   ```bash
   cd IntNotes
   git add -p
   git commit -m "Update Sec. 3 motivation"
   git push                # pushes to gitlab.cern.ch
   cd ..
   ```

2. **The parent repo tracks which submodule commit to use.** After
   pushing a new submodule commit, update the pointer in the parent:

   ```bash
   git add IntNotes        # records the new submodule SHA
   git commit -m "Bump IntNotes to vN"
   git push                # pushes to GitHub
   ```

When pulling new changes from GitHub:

```bash
git pull
git submodule update --init --recursive    # sync submodules to recorded SHAs
```

If you want submodules to track the latest `master` of each submodule
instead of the pinned SHA:

```bash
git submodule update --remote --merge
```

---

## 6. Troubleshooting

**`Permission denied (publickey)` when fetching a submodule.**
Your SSH key isn't registered with `gitlab.cern.ch` (or the wrong key is
being offered). Run `ssh -T -p 7999 git@gitlab.cern.ch` to test, and check
`~/.ssh/config`.

**`Connection timed out` on port 22 to gitlab.cern.ch.**
Default SSH port is 22; CERN GitLab uses 7999. Use the `~/.ssh/config`
snippet in §1 or include `:7999` in the URL.

**Submodule directory is empty after clone.**
You forgot `--recurse-submodules`. Run `git submodule update --init
--recursive`. If that fails with permission errors, you don't have CERN
GitLab access — the submodules are intentionally private.

**Want to skip submodules entirely on a machine that has no CERN access.**
That's fine — the parent repo works without them. Just don't run
`submodule update`. The `IntNotes/` and `IntNotesRun2DimuonReference/`
directories stay empty.

---

## 7. Why this layout

- Parent repo is **public on GitHub** so that recruiters / AI screening
  tools can read the analysis code.
- Internal notes live in **private CERN GitLab** repos as submodules so
  that:
  - ATLAS-internal content is never exposed in public history.
  - The ANA-HION-2023-07-INT1 submodule keeps its existing ATLAS Physics
    Office CI (`.gitlab-ci.yml`) on the gitlab side.
  - PIs / collaborators with CERN accounts edit the notes directly on
    gitlab.cern.ch without needing GitHub access.
- The parent's `.claude/` skills, `kb/`, and code context are visible
  locally alongside the note text, so Claude can write with full context
  while keeping the public-vs-private boundary clean.

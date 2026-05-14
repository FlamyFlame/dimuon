# BigPanda/JEDI PbPb Skimming Jobs — April 2025

## Overview

Final set of JEDI tasks for skimming PbPb2023 and PbPb2024 data with FCal + Centrality info, submitted in April 2025.
All tasks used `pathena TrigRates.py` with centrality weighting (`w_centr`).

An earlier attempt in late March 2025 (TaskID **43793693**, created 2025-03-21) failed with only 63/70752 jobs completing (~0.09%) and was not pursued further.

---

## ✅ PbPb2024 Data — 3 Tasks (all 100% succeeded)

| TaskID | Part | Succeeded/Total | Output Dataset | Date Ended |
|--------|------|-----------------|----------------|------------|
| **44157029** | part1 | 23617 / 23617 ✅ | `user.yuhang.TrigRates.dimuon.PbPb2024data.Apr2025.w_centr.v1_part1._MYSTREAM/` | 2025-04-16 |
| **44169186** | part2 (replacement — v1_part2 job exhausted) | 23356 / 23356 ✅ | `user.yuhang.TrigRates.dimuon.PbPb2024data.Apr2025.w_centr.v2_part2._MYSTREAM/` | 2025-04-17 |
| **44157036** | part3 | 23779 / 23779 ✅ | `user.yuhang.TrigRates.dimuon.PbPb2024data.Apr2025.w_centr.v1_part3._MYSTREAM/` | 2025-04-15 |

### Input Datasets (PbPb2024)

- `data24_hi` HardProbes AOD files, multiple run numbers, merge tags `f1537–f1550`, `m2259/m2267`
- Part 1 & 2 used `InDstxt_PbPb2024_HP_part1.txt` / `InDstxt_PbPb2024_HP_part2.txt`; Part 3 used `InDstxt_PbPb2024_HP_part3.txt`

---

## ✅ PbPb2023 Data — 1 Final Task (~98.7% succeeded)

| TaskID | Succeeded/Total | Final Status | Output Dataset | Date Ended |
|--------|-----------------|--------------|----------------|------------|
| **44157121** | 86778 / 87924 | finished | `user.yuhang.TrigRates.dimuon.PbPb2023data.Apr2025.w_centr.v1._MYSTREAM/` | 2025-04-23 |

### Input Dataset (PbPb2023)

- `data23_hi.periodAllYear2.physics_HardProbes.PhysCont.AOD.repro35_v01`

### Notes on Follow-up Task

A follow-up task **44169187** was submitted (to handle throttled output from 44157121), but was **aborted** on 2025-05-18 with 0/87924 jobs completed (killed by prodsys). No output was produced. **Task 44157121 is the final usable 2023 PbPb output.**

---

## Summary Table

| Dataset | TaskID(s) | Success Rate | Output Dataset(s) |
|---------|-----------|--------------|-------------------|
| PbPb2024 part1 | 44157029 | 100% | `user.yuhang.TrigRates.dimuon.PbPb2024data.Apr2025.w_centr.v1_part1._MYSTREAM/` |
| PbPb2024 part2 | 44169186 | 100% | `user.yuhang.TrigRates.dimuon.PbPb2024data.Apr2025.w_centr.v2_part2._MYSTREAM/` |
| PbPb2024 part3 | 44157036 | 100% | `user.yuhang.TrigRates.dimuon.PbPb2024data.Apr2025.w_centr.v1_part3._MYSTREAM/` |
| PbPb2023 | 44157121 | ~98.7% | `user.yuhang.TrigRates.dimuon.PbPb2023data.Apr2025.w_centr.v1._MYSTREAM/` |

---

## Failed / Aborted Tasks (for reference)

| TaskID | Created | Status | Succeeded/Total | Notes |
|--------|---------|--------|-----------------|-------|
| 43793693 | 2025-03-21 | finished | 63 / 70752 | Earlier March attempt — essentially failed |
| 44169187 | 2025-04-15 | aborted | 0 / 87924 | Follow-up to 44157121; killed by prodsys |

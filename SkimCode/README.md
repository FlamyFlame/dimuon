# SkimCode quick README

## What this code does

This package runs the `TrigRates` Athena algorithm to produce a skim ROOT ntuple (`myfile.root`) for heavy-ion / pp dimuon studies.

Main features in `TrigRates`:

- Event cleaning (GRL, detector-error checks, vertex cut)
- Trigger decisions and prescales (including L1 info)
- Muon trigger matching (single- and dimuon paths)
- Muon MCP tools:
	- selection
	- momentum calibration
	- efficiency corrections (Medium + Tight)
	- trigger matching
- Optional storage of tracks, vertex, MET, truth, and MC event info


## Config scripts (JO vs CA)

Under `scripts/`:

- `TrigRates_JO.py`:
	- Legacy job options configuration
	- Use for **R21 / pre-R24** workflows
- `TrigRates_CA.py`:
	- ComponentAccumulator-based configuration
	- Use for **R24/R25 (recommended for modern setup)**
- `TrigRates.py`:
	- Original legacy script (kept for compatibility with existing run folders)
- `TruthTrigRates.py`:
	- Truth-only skimming configuration


## How to run

### R21 (legacy JO)

Important: run `setupATLAS -c centos7` first, then run setup inside the container shell.

```bash
cd /afs/cern.ch/user/y/yuhang/eos/dimuon/SkimCode
setupATLAS -c centos7
```

Inside the container:

```bash
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
cd build_21
acmSetup --sourcedir=../source AthAnalysis,21.2.200
cd ../run_powheg_pp17fullsim
cp ../scripts/TrigRates_JO.py ./TrigRates_JO.py
athena TrigRates_JO.py --evtMax=100
```

If you still use an old tested local script in the run folder, you can run:

```bash
athena TrigRates.py --evtMax=100
```


### R25 (CA)

```bash
cd /afs/cern.ch/user/y/yuhang/eos/dimuon/SkimCode
source setup_25.sh
cd build_25
acm compile
cd ..
```

Run in a target folder (example `run_24pp`):

```bash
cd run_24pp
cp ../scripts/TrigRates_CA.py ./TrigRates_CA.py
athena TrigRates_CA.py --evtMax=100
```

If `athena` is not accepted by your environment wrapper, fallback:

```bash
python TrigRates_CA.py --evtMax=100
```


## Output

- Default output file: `myfile.root`
- Tree name: `HeavyIonD3PD`


## Notes

- Keep R21 and R25 environments separate (do not mix setup scripts in one shell).
- If C++ in `source/HFtrigValidation/` changes, always rebuild (`acm compile`) in the matching release build directory.
- Recommended workflow:
	- R21 -> `TrigRates_JO.py`
	- R25 -> `TrigRates_CA.py`


## Dataset toggle cheat-sheet

Use exactly one dataset toggle as `True` in the chosen script, keep all others `False`.

### Common run folders

- `run_powheg_pp17fullsim` (R21 JO):
	- script: `TrigRates_JO.py` (or legacy local `TrigRates.py`)
	- toggle: `do_pp_MC_fullsim_17 = True`

- `run_23hi` (R25 CA):
	- script: `TrigRates_CA.py`
	- toggle: `do_hi2023 = True`

- `run_24hi` (R25 CA):
	- script: `TrigRates_CA.py`
	- toggle: `do_hi2024 = True`

- `run_24pp` (R25 CA):
	- script: `TrigRates_CA.py`
	- toggle: `do_pp2024 = True`


### Available toggles

- Data:
	- `do_hi2015`, `do_hi2018`, `do_hi2023`, `do_hi2024`
	- `do_pp2015`, `do_pp2017`, `do_pp2024`
- MC fullsim:
	- `do_pp_MC_fullsim_17`, `do_pp_MC_fullsim_24`


### Quick sanity checks before running

- Exactly one `do_*` toggle is `True`
- Input file path for that toggle exists
- Correct script for release:
	- R21: JO (`TrigRates_JO.py`)
	- R25: CA (`TrigRates_CA.py`)



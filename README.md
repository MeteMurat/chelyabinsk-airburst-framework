# Chelyabinsk-class airburst framework

MATLAB source and numerical replay harness for **An Auditable Shape-Informed Reduced-Order Framework for Chelyabinsk-Class Meteoroid Entry, Fragmentation, and Airburst Consequence Screening** (manuscript AISR-D-26-02161).

**Read [validation/REPLAY_STATUS.md](validation/REPLAY_STATUS.md) before interpreting this repository as a completed reproduction.** The initial repository is explicitly a replay candidate, not a claim that a new MATLAB run has passed. The status is updated only from actual MATLAB result files.

This repository separates:

1. Byte-identical historical scientific kernels and the original full-1000 MATLAB driver.
2. A new portable test harness, using the frozen input design and comparing newly computed outputs with archived reference CSVs.
3. Third-party/raw binary inputs that must be supplied locally and are **not uploaded to GitHub**.

## What the replay tests

- Binary STL processing, 19.8 m volume-equivalent scaling, all 500 orientation samples with the MATLAB `twister` random stream, and the local-radius proxy.
- Nominal G4 and G0-G4 at capacity 128, base horizon 360 s and time step 0.02 s.
- Independent reference-vs-diagnostic kernel agreement for nominal scalar values, time series, deposition profile, and RNG state.
- Eight explicitly frozen Monte Carlo rows as a preliminary test.
- The original 1,000-member campaign, including its original horizon, closure, semantic, pilot, and parity checks.
- Pressure radii computed from the **new** source terms with the frozen global coefficient.
- 4,000 exposure records recomputed from the exact native population raster and the **new** pressure radii.

Expected CSVs are regression targets. They are not substituted for newly executed simulations. The replay does not tune parameters to make a failed comparison pass. Numerical reproduction does not independently validate the physical approximations against observations.

## Required software

The targeted local configuration is **MATLAB R2025b and Mapping Toolbox** (`readgeoraster`, `intrinsicToGeographic`). Other releases have not been claimed as tested. The default replay is sequential CPU computation and does not open a parallel pool or request GPU arrays. This does not reconstruct the hardware scheduling of every historical run.

Windows PowerShell is used by the convenience launcher. MATLAB can also be called directly on other supported systems. Python is not required for the MATLAB scientific calculation. Git and GitHub CLI are only needed for the optional publication workflow.

## Required inputs -- not in the public repository

Place the exact files at:

```
04_data/geometry/Meteor.stl
04_data/population/turkiye.tif
```

See [docs/DATA_REQUIREMENTS.md](docs/DATA_REQUIREMENTS.md) and `inputs/required_binary_inputs.json` for SHA256 identities and access limitations. A different grid, year or geometry must not silently be substituted.

## Run in a clean folder

Do not add old FOSTRAD/meteor folders with `genpath`. The harness limits the path to MathWorks folders and the explicit final-core/reproduction folders; it verifies which kernel is actually resolved. It does not call `savepath` or change installed MATLAB settings.

From a cloned repository with the required local binary inputs:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\reproduce.ps1 -Mode smoke -MatlabExe "C:\Mathworks.Matlab.R2025b.x64\bin\matlab.exe"
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\reproduce.ps1 -Mode full1000 -MatlabExe "C:\Mathworks.Matlab.R2025b.x64\bin\matlab.exe"
```

Or from MATLAB, starting in this repository:

```matlab
addpath(fullfile(pwd,'reproduction'));
run_clean_replay(pwd,'smoke');
run_clean_replay(pwd,'full1000');
```

The full campaign's historical driver supports checkpoints. Resume only in an unchanged workspace with the same source and input hashes. The supplied outer launcher enforces a source-manifest identity. Never copy archived *output* matrices into checkpoint directories.

## Outputs and success criteria

Local outputs are stored in `replay_outputs/` and `05_outputs/`; these are ignored by Git. Public-safe status files are written under `validation/`.

- `PASS_SMOKE_REPLAY_FULL1000_PENDING`: geometry, nominal/G0-G4, eight MC cases and their downstream transformations passed; **not** a full-1000 pass.
- `PASS_CLEAN_MATLAB_FORWARD_REPLAY_FULL1000_AND_POPULATION`: all numerical replay checks, 1,000 MC rows and 4,000 exposure rows passed.
- `FAIL_CLEAN_MATLAB_REPLAY`: a real failed check, dependency issue or execution error; inspect local `checks.csv` and console logs. Do not relabel it PASS.

A historical terminal-incomplete flag is not automatically a failed airburst result. The original residual, semantic and conservation gates remain in force; terminal active mass is not promoted to a published physical endpoint.

The old nominal exposure wrapper rounded energy to 331.6302 kt, giving a tiny radius difference relative to the unrounded nominal pressure calculation. The new chained replay uses the unrounded result, records that difference and requires the archived cell/population counts to match. No coefficient is recalibrated. Physical circular footprint and raster support area have separate output columns.

## Scope and provenance

`01_final_core/` contains the executable historical core. `historical_drivers/` preserves three original orchestration files for review; these are **not automatically added to MATLAB's path** and may require historical audit contracts. The portable harness does not fabricate those contracts. See [docs/REPLAY_DESIGN.md](docs/REPLAY_DESIGN.md) and `docs/source_origins.json`.

The full replay designation covers the numerical chain and generated summary outputs, not a byte-for-byte reconstruction of every typeset manuscript page or graphic. Geometry and original Monte Carlo scripts also regenerate their own diagnostic figures.

No FOSTRAD panel-aerodynamics code is required by the default final solver. Exploratory FOSTRAD code, the obsolete 40-ID/120-s runs and the unseeded historical kernel are deliberately not mixed into the executable path.

See [RIGHTS_AND_THIRD_PARTY_NOTICES.md](RIGHTS_AND_THIRD_PARTY_NOTICES.md). Public availability does not imply a newly granted blanket software/data licence.

## Lightweight archive checks without MATLAB

`python tests/verify_reference_fixtures.py` verifies source/input/output fixture identities and 4,000 archived pressure transformations. This check does not run the entry solver and does not mark the MATLAB replay as complete.

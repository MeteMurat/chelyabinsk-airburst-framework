# Replay design and identity

## Unmodified scientific sources

- `meteor_single_run.m`: reference entry/fragmentation core; SHA256 `BE3DAC129E6A9DC9D0CCF1D8F95D72E7E0A11EA19E041942E71A20CD47697ABC`.
- `meteor_single_run_r3b3_diagnostic.m`: diagnostic core; SHA256 `58EF0D1C6F45476C600356961B3BCABAF9C7A3D2FE0569E20401D9436E68C59A`.
- `run_r3c2c_r2_corrected_full1000.m`: historical full campaign driver; SHA256 `D3140A7D002FC33F074298FC47868FF73948A0D5D9390F1EBB40B90FD04D3082`.
- Both geometry scripts are copied byte-for-byte. Their local helper functions, rather than FOSTRAD `STLcurvature.m`, define the geometry/radius outputs.

All exact file identities are in `SOURCE_MANIFEST.json` and `docs/source_origins.json`. Two expected-summary CSVs omit only their private `stl_file` path column; numeric strings are preserved.

## New portable orchestration

`reproduction/run_clean_replay.m` is a new harness, not a newly claimed historical execution file. It uses the historical nominal parameter contract and the frozen MC input design. Nominal/G0-G4 are direct core calls. Full-1000 invokes the unmodified original driver with genuine archived design/pilot/horizon inputs. The downstream pressure and raster-centroid transformations implement the archived algorithms without requiring historical folder names or manufactured PASS files.

The named historical horizon-gate CSV is a preserved upstream input, not a result of the new run. The new run's acceptance comes from actual core calls, original campaign gates and comparisons with the held reference outputs.

## Random streams and computation

The current batch stream is explicitly set to `twister`; the historical kernels then set their frozen seeds. Input design generation is not re-randomized: the exact 1,000 physical-input rows and fragmentation seeds are supplied. This freezes the experiment rather than constructing another MC sample.

Default execution is CPU/sequential. GPU availability, historical author attestations about acceleration and the replay's execution backend are distinct claims. No past worker count or GPU capacity is inferred by this test.

## Checks

- Source hashes and resolved kernel paths, before and after computation.
- Both binary-input hashes.
- Full geometry summary, including 500 stochastic orientation samples; local-radius summary.
- G0-G4 and nominal metrics; nominal reference/diagnostic time series, deposition profiles and random states.
- Eight predefined MC cases before the full campaign.
- Every archived publishable metric across 1,000 newly computed runs; capacity-loss, terminal flag, and original full-campaign gates.
- Pressure radii and circular areas computed from newly obtained source terms.
- Native-grid exposure, counts and support areas from the newly obtained radii.

General numerical regression tolerance is max(1e-8 absolute, 1e-10 relative), geometry max(1e-9 absolute, 1e-9 relative), and integer/count comparisons are exact. Population comparisons allow 1e-4 person absolute / 1e-11 relative for summation; valid/included cell counts are exact. Conservation and residual gates independently retain the historical thresholds. Tolerances are explicit in source, not adjusted after a mismatch.

## Not claimed

This is not independent physical validation, a new hazard model, a new uncertainty design, proof of historical acceleration use, or a statement that every third-party datum is publicly redistributable. New MATLAB results remain pending until the actual local MATLAB process passes. No manuscript/ledger file is changed by this package.

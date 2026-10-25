# Rights and third-party boundary

The author requested public source availability for peer review. No new blanket licence is asserted for all historical files or third-party inputs by this packaging operation; existing notices and rights remain with their respective holders. A specific open-source licence for the author-owned implementation can be selected separately by the rights holders.

The default executable source does not redistribute the FOSTRAD suite, its `STLcurvature.m`, or other exploratory third-party modules. This avoids incorrectly relicensing those modules as original project code. FOSTRAD development history is distinct from the frozen reduced-order scientific implementation.

MATLAB and the Mapping Toolbox are external proprietary dependencies and are not included. Raw STL/raster inputs are required locally but excluded from Git publication; see `docs/DATA_REQUIREMENTS.md`. Presence in a historical archive and hash equality establish identity, not a licence grant.

The reference result tables are supplied to identify the paper's numerical outputs and to test numerical replay. Do not conflate these output tables with redistribution of the raw population raster.

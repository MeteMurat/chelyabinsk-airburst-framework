# Exact inputs and public-data boundary

## Geometry

`04_data/geometry/Meteor.stl`

SHA256: `30DC0618D7FEF9F2C4A9AE4D3830EC85EA1B2B9FD8EEAC64260461E223303F76`

Size: 1,344,084 bytes. After duplicate-vertex consolidation: 13,442 vertices, 26,880 faces. This is the input found in the supplied historical development archive and matched to the archived manifest. This identity does not establish that the shape is the uniquely reconstructed historical Chelyabinsk object or establish its redistribution licence.

## Population analysis raster

`04_data/population/turkiye.tif`

SHA256: `084349E49339606D30B95B6AF60EB06DF012297BEF2040F72E6D3C836253BE16`

Size: 3,756,736 bytes. Grid: 755 by 2,299; geographic WGS84/EPSG:4326; 684,872 admissible cells and 684,119 positive cells. Frozen domain population: 87,794,916.9867744.

The manuscript identifies the source as WorldPop Global2 2026 constrained population-count product R2025A v1 alpha, DOI `10.5258/SOTON/WP00845`. The analysis raster's identity was bound to the archived local manifest. **A complete public download-to-exact-analysis-raster construction recipe has not been verified by this packaging step.** Obtaining a generic country raster or a newer release is not a verified substitute for this exact input.

To reproduce the precise paper outputs, obtain the exact analysis input through the corresponding author's documented review/data-access mechanism or establish a separately documented reconstruction from the provider. Respect the provider's licence and the paper's stated non-redistribution boundary. No raw raster is embedded in the GitHub repository.

## Why two delivery surfaces exist

The local replay package provided to the author contains the two already supplied inputs under `LOCAL_ONLY_INPUTS/`. The launcher copies them to the isolated local workspace after hashing. The GitHub publisher uses an explicit text-file allowlist; the two binaries, MATLAB MAT files, raw logs, editor/reviewer correspondence, manuscript, and private evidence ZIPs are excluded.

Thus the author can perform the full local replay without a new archive scan. Public code availability is **not** advertised as one-click complete data availability. The remaining public raw-input acquisition limitation must be stated honestly to reviewers.

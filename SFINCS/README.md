# SFINCS workflows

## SCS-CN infiltration

New global P1/P2/P3 and regional P1/P2 builds enable classic SCS-CN by default
through `scs_infiltration.py`. Supply the external **GCN250 ARCII** raster:

```powershell
$env:SFINCS_GCN250_PATH = 'D:\external\soil\GCN250_ARCII.tif'
```

On Linux: `export SFINCS_GCN250_PATH=/data/soil/GCN250_ARCII.tif`.
For the TOML launcher, set `[infiltration].cn_raster` in the configuration.
To explicitly retain a no-additional-infiltration build, set
`SFINCS_INFILTRATION=none` (or TOML `method="none"`). This does not remove
infiltration already present in a supplied base model.

The module median-reprojects CN to the original rotated/projected grid and
writes `sfincs.scs` in active-cell `sfincs.ind` order: S (inches) =
max(1000/CN - 10, 0), rounded to 3 decimals, little-endian float32. It sets
`sfacinf=0.2`, removes zero `qinf`, and sets `storecumprcp=1`, required to avoid
the known classic-CN cumulative-precipitation issue in the solver version used
for this study. Missing CN is explicitly counted and assigned CN100/S0
(no infiltration), not silently interpolated. Soil parameters are fixed
across scenarios; no dry-period recovery or future land-cover change is added.

Every build includes `infiltration_metadata.json` with the grid/source signature,
binary checksum and missing-CN count. Immutable files may remain hardlinked;
the module replaces config/binary files atomically instead of editing hardlinks
in place. Existing simulated directories are rejected. New FULL and all
removal experiments must use consistent infiltration; do not mix old FULL
outputs with new infiltration removal outputs. This is **surface-water
infiltration**, not a groundwater seepage model. Soil data are external and
are not redistributed with the repository.

Bundled precomputed examples retain their original provenance; enabling the
module affects newly built cases, not archived example outputs.

Use `global_build/P1_build_sfincs_models_from_partition.py` to construct the 97
reusable `GTC_####` grids. Then use the historical/future compound-forcing
builders in `global_build/` to assemble runnable events.

For a first test, return to the repository root and run:

```bash
python example/PRD_single_TC/run_prd_single_tc.py --dry-run
python example/PRD_single_TC/run_prd_single_tc.py
```

The regional `P1_PRD_building.py`, `P1_YRD_building.py`, `P1_BoB_building.py` and
`P1_Misp_building.py` scripts are retained for provenance. They expect to run in
a HydroMT-SFINCS examples directory containing the documented `data/` tree;
set `CAMA_MAP_ROOT` for their CaMa-Flood map input. These four scripts are
alternative P1 regional builders. After building the selected region, run
`P2_sfincs_couple_adcirc_cama_tc_rain_future_batch.py` to assemble its coupled
forcing cases.

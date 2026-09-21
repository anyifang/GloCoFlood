# SFINCS workflows

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

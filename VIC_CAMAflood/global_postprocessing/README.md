# Global CaMa-Flood to SFINCS inflow processing

The original `VIC_CAMAflood/auto_cama_vic_coupled.py` workflow prepares and
runs VIC/CaMa-Flood. The scripts here convert the routed global discharge into
the GTC-domain products consumed by `SFINCS/global_build`.

`P1_plot_sfincs_inflow_diagnostics.py` matches SFINCS source points to independent
upstream CaMa-Flood cells, reads historical and CMIP6 `outflw*.bin` files,
constructs daily caches, derives monthly climatologies and TC-month weights,
and writes diagnostic tables and figures.

Example:

```bash
python P1_plot_sfincs_inflow_diagnostics.py \
  --input_dir external/camaflood/output \
  --map_dir external/camaflood/glb_15min \
  --models_dir external/sfincs/gtc_models \
  --out_dir work/gtc_cama_inflow_figures \
  --ibtracs_path external/tracks/IBTrACS.ALL.v04r01.nc \
  --figs9_basin_events_csv external/tracks/05_basin_landfall_events.csv
```

Use `--ids GTC_0009 --max_files_per_case 1 --max_steps_per_file 10` for a
small PRD smoke test. The same defaults can be supplied with
`GLOCOFLOOD_CAMA_ROOT`, `GLOCOFLOOD_SFINCS_MODEL_ROOT`,
`GLOCOFLOOD_IBTRACS_NC` and `GLOCOFLOOD_TC_EVENTS_CSV`.

After daily caches and the 12-month weight table exist, create the per-inlet
P50/P90 table:

```bash
python P2_build_gtc_inlet_tc_flow_summary_from_cache.py \
  --root work/gtc_cama_inflow_figures \
  --cache-key <key-printed-by-the-diagnostic-run>
```

The SFINCS builder consumes the generated `tables/` and `cache/` directories;
these products are intentionally ignored by Git.

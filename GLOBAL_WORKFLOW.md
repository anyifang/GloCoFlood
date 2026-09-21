# Reproducible global compound-flood workflow

The public workflow is organized by data dependency rather than by individual
case folder:

```text
downscaled TC tracks
        |
        v
ADCIRC/global_build  ------> fort.63 water levels and fort.22 metadata
        |                                      |
        |                                      v
VIC_CAMAflood ----------------------> SFINCS/global_build
  routed discharge + GTC P50/P90       boundary + flow + C15/TCR rainfall
                                               |
                                               v
                                 runnable historical/SSP SFINCS cases
```

## Execution order

1. Build the 97 reusable SFINCS domains with
   `SFINCS/global_build/P1_build_sfincs_models_from_partition.py` (or supply
   existing domains), and build the global ADCIRC blocks.
2. Run VIC and CaMa-Flood using `VIC_CAMAflood/auto_cama_vic_coupled.py` and
   the appropriate HPC wrapper.
3. Run `VIC_CAMAflood/global_postprocessing/P1_plot_sfincs_inflow_diagnostics.py`
   and
   `VIC_CAMAflood/global_postprocessing/P2_build_gtc_inlet_tc_flow_summary_from_cache.py`.
4. Generate historical ADCIRC TC cases with
   `ADCIRC/global_build/P5_build_adcirc_era5_tc_run_files_from_selected_tracks.m`.
5. Select/stage return-period cases, run ADCIRC, and collect `fort.63` using
   ADCIRC stages P6, P7 and P8; use P9 to rebuild the future forcing.
6. Copy `SFINCS/global_build/config/paths.example.toml` to `paths.toml`, point
   it to the products above, and validate it.
7. Run the fixed PRD one-TC example before launching the four global scenario
   builders.

```bash
cd SFINCS/global_build
python run_global_build_from_config.py --config config/paths.toml --check-only
python run_global_build_from_config.py --config config/paths.toml --workers 8
```

The self-contained PRD example is run separately from the repository root:

```bash
python example/PRD_single_TC/run_prd_single_tc.py --dry-run
python example/PRD_single_TC/run_prd_single_tc.py
```

All large global inputs and outputs stay outside the repository; only the
compact PRD example inputs are included. Configuration files
contain paths only; credentials must never be added to them.

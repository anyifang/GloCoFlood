# Global ADCIRC tropical-cyclone workflow

This directory extends the original regional ADCIRC scripts to the global
ADCIRC block inventory used by the SFINCS production builder. Large meshes,
tracks, lookup tables, model binaries and generated `fort.*` cases are not
included.

## Main stages

1. `P1_build_adcirc_model_blocks_from_partitions.m` groups the global coastal
   SFINCS partitions into the 41 regional ADCIRC blocks and writes the block
   catalog and inner/outer boundary rings.
2. `P2_run_global_tc_adcirc_block_meshes.m` generates the corresponding
   OceanMesh2D `fort.14` meshes.
3. `P3_build_tide_adcirc_fort13151922_fix_stable.m` converts each mesh into a
   stable static library containing `fort.13`, `fort.14`, `fort.15`,
   tide-only `fort.19`, background `fort.22` and a batch template.
4. `P4_count_and_plot_era5_tc_for_adcirc_inner_domains.m` selects the downscaled
   historical TC events that cross each exact inner refined domain.
5. `P5_build_adcirc_era5_tc_run_files_from_selected_tracks.m` creates historical
   event directories and generates C15 `fort.22`, tidal `fort.19`, control and
   static mesh files. For a multi-process build, run
   `P5_parallel_build_adcirc_era5_tc_run_files_block_workers.m` instead; it
   invokes the same P5 builder by block and must not be run after the serial
   P5 build.
6. `P6_prepare_return_period_tc_event_reruns.m` reads the selected controlling
   TC tables, removes duplicate events and stages 100-, 200- and 500-year
   reruns.
7. `P7_submit_selected_return_period_tc_events.sh` submits the staged ADCIRC jobs.
8. `P8_collect_return_out_fort63_maxele63.sh` collects their water-level outputs.
9. `P9_build_future_return_period_tc_event_reruns.m` rebuilds `fort.19` and
   `fort.22` for SSP1-2.6, SSP2-4.5 and SSP3-7.0 using the precomputed
   block/TC intensity and local sea-level-rise tables.

The controlling-TC selection tables and future scaling tables are scientific
intermediate products, not bundled data. Point the staging and future scripts
to those products through the variables below.

The versioned `catalog/` directory contains the lightweight 41-block metadata,
all boundary-ring CSV files and the 41-block-to-97-GTC membership table used by
the published workflow. It does not contain meshes or simulation output.

## Configuration

Copy `configure_paths.example.m` to `configure_paths.m`, adjust the
repository-relative `external/` and `work/` layout when needed, and run it
before the workflow scripts:

```matlab
run('configure_paths.m');
run('P5_build_adcirc_era5_tc_run_files_from_selected_tracks.m');
```

For a clean rebuild from the upstream coastal partitions, run the first four
stages before the historical case builder. If the published 41-block layout is
being reused, keep `ADCIRC_BLOCK_CATALOG_ROOT` pointed at `catalog/` and begin
with mesh/static-input generation.

For eight external workers:

```matlab
run('configure_paths.m');
setenv('ADCIRC_TC_PARALLEL_LAUNCH', '0');
run('P5_parallel_build_adcirc_era5_tc_run_files_block_workers.m');
```

Inspect the generated plan and then launch its generated worker script. Each
worker writes distinct blocks and has a private TMD directory.

## Wind convention

`vmax_trks` is read as the total 1-min near-surface maximum wind. Translation
speed is removed to obtain the symmetric surface core used by the Cv/Rmax and
Holland-B calculations. C15 supplies the radial profile of axisymmetric 1-min
surface tangential wind; the translation vector is restored once, and the
complete field is converted to the 10-min surface wind written to `fort.22`
with the configured `0.893` averaging-period factor. No additional `0.85`
multiplier is applied. The external lookup-field name `C15.vg` is not a
gradient, upper-tropospheric or environmental wind variable; local variables
use the explicit name `c15_surface_wind_1min_ms`.

## Required external inputs

- global ADCIRC block meshes and block/event tables;
- historical and future downscaled TC tracks;
- C15 lookup profiles, calibrated Cv table and wind-pressure fit table;
- TMD tidal model files;
- controlling-TC return-period tables;
- future TC-intensity and local-SLR scaling tables;
- an ADCIRC executable and site-specific batch template.

Copy `sub_intel.example.sh` to a local `sub_intel.sh` and set `ADCIRC_BIN` for
the target cluster. No institution-specific executable or queue path is stored
in the public template.

# Global historical and future SFINCS case builder

## Scope

This directory contains the production-scale builder for the 97 global SFINCS
domains and their selected tropical-cyclone events. It creates four scenarios:

- `historical`
- `ssp126`
- `ssp245`
- `ssp370`

Each runnable case combines:

1. an existing SFINCS domain template;
2. ADCIRC boundary water levels;
3. CaMa-Flood river discharge;
4. C15-driven, physics-based TCR rainfall;
5. scenario-specific terrain/VLM information where applicable.

Large source datasets and generated cases are deliberately excluded. The sole
exception is the compact PRD one-TC input case used by the direct-run example.

## Scripts

| Script | Purpose |
|---|---|
| `P1_build_sfincs_models_from_partition.py` | Build the 97 reusable SFINCS grids from the published coastal partitions. |
| `plot_sfincs_partition_models.py` | Plot and inspect generated grids, boundaries and matched CaMa-Flood inflows. |
| `P2_build_sfincs_historical_return_period_compound_forcing.py` | Historical event builder and shared implementation. |
| `P3_build_sfincs_future_return_period_compound_forcing.py` | SSP1-2.6, SSP2-4.5 and SSP3-7.0 extension. |
| `run_four_compound_builds_parallel.py` | Four- or eight-process orchestration, cache preparation and logging. |
| `run_global_build_from_config.py` | Portable TOML-configured entry point. |
| `../../example/PRD_single_TC/run_prd_single_tc.py` | Directly runnable one-TC Pearl River Delta example, with optional full reconstruction. |
| `utilities/build_sfincs_forcing_decomposition_cases.py` | Optional leave-one-driver-out cases. |
| `utilities/check_sfincs_map_result.py` | Optional result integrity check. |

The numbered production sequence is P1 (reusable domains), P2 (historical
forcing) and P3 (future forcing). `run_four_compound_builds_parallel.py` and
`run_global_build_from_config.py` are alternative launchers for P2/P3, not
additional stages to run afterwards. The plotting and `utilities/` scripts
are optional diagnostics.

## C15/TCR wind convention

The track variable `vmax_trks` is the 1-min near-surface wind. Rmax and
Holland-B retain the surface-wind definitions used by the ADCIRC builder.
Physics-based TCR requires gradient-level wind, so the builder first computes
`gradient_core = max(vmax_trks - translation_speed, 0) / 0.9`. This gradient
peak selects the C15 lookup profile, which is then used directly as the TCR
gradient-wind profile; `C15.vg` is not divided by `0.9` a second time. The
same `0.9` relationship is used by the TCR humidity diagnosis. ADCIRC's
independent `0.893` conversion from 1-min to 10-min surface wind is not
applied to TCR rainfall.


## Installation

Create an environment with Python 3.11 or newer and install the packages in
`requirements.txt`. CLIMADA and CLIMADA Petals must be installed with their
TCR data dependencies.

## Configuration

Copy:

```text
config/paths.example.toml -> config/paths.toml
```

The supplied values are repository-relative (`external/...` for inputs and
`work/...` for outputs). Relative values are resolved from the repository
root, independently of the current working directory. Edit this local file if
your data layout differs; it is intentionally ignored by Git.

Validate the configuration without launching builders:

```bash
python run_global_build_from_config.py --config config/paths.toml --check-only
```

## Build the 97 reusable SFINCS domains

The published partition polygons and metadata are included in
`global_tc_coastal_model_partitions/`. Check the selection without accessing
the large terrain and river datasets:

```bash
python P1_build_sfincs_models_from_partition.py --dry-run --limit 1 --no-plot
```

For a real build, supply a HydroMT-SFINCS `examples` directory containing
`data/FABDEM_Coastal_200km`, `data/coastaline`, `data/river`, `data/Landcover`
and the MDT NetCDF used by the workflow, plus a CaMa-Flood map root containing
`glb_15min`:

```bash
python P1_build_sfincs_models_from_partition.py \
  --hydromt-examples external/HydroMT-SFINCS/examples \
  --cama-base-dir external/camaflood/map_v420 \
  --cama-res 15 \
  --res 200 \
  --folder-naming model-id \
  --output-dir work/global_sfincs_partition_models_cama_15min
```

The same two external paths may be provided through
`HYDROMT_SFINCS_EXAMPLES` and `CAMA_MAP_ROOT`. Model folders default to stable
`GTC_####` identifiers so that the historical and future forcing builders can
consume them directly.

## Pearl River Delta quick example

Immediately after cloning, prepare and verify the bundled historical PRD TC
case with:

```bash
python example/PRD_single_TC/run_prd_single_tc.py
```

Use `--dry-run` first when checking a new installation. To reconstruct the
same case from the complete source datasets, configure `config/paths.toml` and
add `--full-build`. See `example/PRD_single_TC/README.md` for the fixed domain,
block, event and output layout.

## Recommended smoke test

Build one event without writing output:

```bash
python run_global_build_from_config.py --config config/paths.toml --dry-run --limit-events 1 --workers 4
```

Then generate one real event:

```bash
python run_global_build_from_config.py --config config/paths.toml --limit-events 1 --workers 4
```

## Full build

```bash
python run_global_build_from_config.py --config config/paths.toml --workers 8 --continue-on-error
```

Use `--overwrite` only when existing complete cases must be regenerated.

## Output layout

```text
<output_root>/
|-- historical/
|-- ssp126/
|-- ssp245/
|-- ssp370/
|-- _cache/
`-- _parallel_build_logs/
```

The runnable scenario directories use relative references and hard links where
possible. Preserve the directory tree when transferring cases to another
filesystem or HPC system.

## Reproducibility warning

The scripts are open, but the results also depend on external ADCIRC outputs,
CaMa-Flood products, downscaled TC tracks, monthly environmental fields, C15
lookup tables, SFINCS domain templates, VLM-adjusted terrain and the exact
CLIMADA versions. Record checksums or dataset versions for a reproducible
publication archive.

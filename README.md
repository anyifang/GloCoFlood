<h1 align="center">GloCoFlood</h1>

<p align="center">
  An open workflow for global-to-regional compound-flood modelling with
  ADCIRC, VIC, CaMa-Flood and SFINCS.
</p>

<p align="center">
  <img alt="MATLAB" src="https://img.shields.io/badge/MATLAB-ADCIRC-orange">
  <img alt="Python" src="https://img.shields.io/badge/Python-SFINCS%20%7C%20VIC-blue">
  <img alt="Shell" src="https://img.shields.io/badge/Shell-HPC%20workflow-lightgrey">
  <img alt="Status" src="https://img.shields.io/badge/status-research%20workflow-informational">
</p>

GloCoFlood couples storm surge and tides, river discharge and
tropical-cyclone rainfall to construct high-resolution coastal inundation
cases. The public workflow contains the production scripts used for 41 ADCIRC
blocks and 97 reusable global SFINCS coastal domains, together with a compact
Pearl River Delta (PRD) example that can be prepared immediately after the
repository is downloaded.

> **Scope.** This is a reproducible research workflow rather than a packaged
> modelling service. Large licensed datasets, model executables and global
> simulation outputs are not distributed in the repository.


## Coupled modelling framework
<img width="5223" height="4130" alt="fig1" src="https://github.com/user-attachments/assets/f5b953df-eb7f-4ba9-a5b4-312dd4c88542" />

| Component | Function | Principal products |
| --- | --- | --- |
| [ADCIRC](ADCIRC/README.md) | Regional meshes, tides and TC wind-pressure forcing | `fort.13`, `fort.14`, `fort.15`, `fort.19`, `fort.22`, `fort.63` |
| [VIC + CaMa-Flood](VIC_CAMAflood/README.md) | Rainfall-runoff modelling and river routing | Runoff fields and routed inlet discharge |
| [SFINCS](SFINCS/README.md) | Coupling and high-resolution inundation modelling | Runnable historical and future compound-flood cases |

The complete dependency order is documented in
[GLOBAL_WORKFLOW.md](GLOBAL_WORKFLOW.md).

## Repository structure

```text
GloCoFlood/
|-- ADCIRC/
|   |-- Global_autofunction/          MATLAB geometry and mesh helpers
|   `-- global_build/                 P1-P9 global ADCIRC workflow
|-- SFINCS/
|   |-- global_build/                 97-domain and compound-case builders
|   `-- P1_* / P2_*                   original regional reference scripts
|-- VIC_CAMAflood/
|   |-- global_postprocessing/        routed-flow products for SFINCS
|   `-- auto_cama_vic_coupled.py      regional VIC/CaMa-Flood workflow
|-- data_download/                    supported data-download utilities
|-- example/PRD_single_TC/            bundled, directly verifiable example
|-- check_repository.py               dependency-free repository checker
|-- QUICKSTART.md                     shortest first-use route
|-- GLOBAL_WORKFLOW.md                end-to-end production sequence
|-- VALIDATION.md                     checks completed for this release
`-- scripts_manifest.csv              public script inventory
```

Files prefixed with `P1_`, `P2_`, and so on are ordered by dependency within
their directory. Files sharing a number are independent tasks or alternative
implementations at the same stage. In particular, the serial and parallel
ADCIRC P5 builders are alternatives and must not both be used for the same
build.

## Global production workflow

### 1. Acquire and register input data

The utilities in [`data_download/`](data_download/README.md) support official
IBTrACS and ERA5/CDS requests and check the proposed external-data layout:

```bash
python data_download/P1_download_ibtracs.py
python data_download/P1_download_cds_era5.py tcr-environment \
  --start-year 1975 --end-year 2014
python data_download/P2_check_data_layout.py external
```

Downloads remain outside version control. FABDEM, CaMa-Flood maps, tidal data,
C15/CLIMADA support data and project-specific downscaled historical and future
TC tracks must be obtained separately from their official providers under the
applicable licences.

### 2. Build ADCIRC domains and event forcing

[`ADCIRC/global_build/`](ADCIRC/global_build/README.md) provides the ordered
P1-P9 workflow to:

1. group the 97 coastal domains into 41 ADCIRC blocks;
2. generate OceanMesh2D meshes and static/tidal `fort.*` inputs;
3. select historical TC events and build event-specific forcing;
4. stage return-period simulations and collect `fort.63`; and
5. reconstruct future SSP1-2.6, SSP2-4.5 and SSP3-7.0 cases.

Copy `ADCIRC/global_build/configure_paths.example.m` to the ignored local file
`configure_paths.m`, set the external-data and work locations, and follow the
numbered scripts in the component README.
<img width="6004" height="5223" alt="adcirc_fort14_oceanmesh_blog_panels_6col_sharedcb_v6" src="https://github.com/user-attachments/assets/04b7d8a9-c08e-42c3-bfdc-a51781d59c64" />

### 3. Generate river discharge

[`VIC_CAMAflood/`](VIC_CAMAflood/README.md) prepares VIC meteorological
forcing, converts runoff for CaMa-Flood and routes river discharge. The global
post-processing stages then match independent upstream CaMa-Flood cells to
SFINCS inlet points and construct the TC-weighted P50/P90 inflow products:

```bash
python VIC_CAMAflood/global_postprocessing/P1_plot_sfincs_inflow_diagnostics.py --help
python VIC_CAMAflood/global_postprocessing/P2_build_gtc_inlet_tc_flow_summary_from_cache.py --help
```

### 4. Build the SFINCS domains

<img width="3285" height="1849" alt="image" src="https://github.com/user-attachments/assets/b94de808-9dc2-445b-9b09-f8165dd4dd9a" />

[`SFINCS/global_build/P1_build_sfincs_models_from_partition.py`](SFINCS/global_build/P1_build_sfincs_models_from_partition.py)
builds reusable `GTC_####` domains from the included partition inventory. A
non-writing selection check can be run before the large external terrain and
river datasets are connected:

```bash
python SFINCS/global_build/P1_build_sfincs_models_from_partition.py \
  --dry-run --limit 1 --no-plot
```

### 5. Assemble historical and future compound cases

Copy
`SFINCS/global_build/config/paths.example.toml` to the ignored local file
`paths.toml`, edit the external paths, and validate the configuration:

```bash
python SFINCS/global_build/run_global_build_from_config.py \
  --config SFINCS/global_build/config/paths.toml --check-only
```

Start with one event:

```bash
python SFINCS/global_build/run_global_build_from_config.py \
  --config SFINCS/global_build/config/paths.toml \
  --limit-events 1 --workers 4
```

After that smoke test passes, launch the full historical and three-scenario
build with eight workers:

```bash
python SFINCS/global_build/run_global_build_from_config.py \
  --config SFINCS/global_build/config/paths.toml \
  --workers 8 --continue-on-error
```

The builders cache the invariant ADCIRC-to-SFINCS boundary mapping and reuse
static domain files where possible. Preserve the generated directory tree
when transferring cases because `sfincs.inp` may contain relative references.


## Software and data requirements

The full production workflow requires local installations of MATLAB,
OceanMesh2D, ADCIRC, Python, VIC, CaMa-Flood, HydroMT-SFINCS/SFINCS, CLIMADA
and CLIMADA Petals, together with their model-specific datasets. Exact package
and executable versions should be recorded for each scientific release.

The repository intentionally excludes:

- model executables and cluster-specific modules;
- global forcing archives and licence-controlled datasets;
- generated ADCIRC and SFINCS case trees;
- bulk NetCDF, raster and binary results; and
- credentials or machine-specific path configuration.

Do not commit `paths.toml`, `configure_paths.m`, credentials or generated
outputs. The supplied `.gitignore` excludes these common local products.

## Validation and reproducibility

The release checker uses only the Python standard library:

```bash
python check_repository.py
```

It checks the manifest, 97 SFINCS partitions, 41 ADCIRC blocks, all 97
block-domain memberships and every checksum in the bundled PRD case. See
[VALIDATION.md](VALIDATION.md) for the completed code, path, syntax and
fresh-copy tests. These checks validate repository integration; they do not
replace scientific validation of external datasets or full ADCIRC,
CaMa-Flood and SFINCS simulations.

### Application in four representative estuarine and deltaic regions
To demonstrate the transferability of the modelling framework, GloCoFlood was
applied to four contrasting, TC-prone coastal systems: the Mississippi River
Delta (MRD), Yangtze River Delta (YRD), Bay of Bengal (BoB) and Pearl River
Delta (PRD). The same coupled workflow is used in each region, while the model
domains, mesh resolution, river network and forcing boundaries are adapted to
the local coastal setting.

<img width="3163" height="3499" alt="global_domain_examples" src="https://github.com/user-attachments/assets/8ba6a4c2-b89d-447f-8bdf-a46396db4c75" />


## Quick start


Run these commands from the repository root:

```bash
python check_repository.py
python -m pip install -r example/PRD_single_TC/requirements.txt
python example/PRD_single_TC/run_prd_single_tc.py
```

The first command checks the published inventories and bundled-data
checksums. The remaining commands prepare and validate one complete historical
PRD compound-forcing case under:

```text
example/PRD_single_TC/output/
```

This default example does not require the global source archives and does not
launch SFINCS. If a SFINCS executable is available, run the hydrodynamic model
as well with:

```bash
python example/PRD_single_TC/run_prd_single_tc.py \
  --sfincs-executable /path/to/sfincs

<img width="9990" height="3456" alt="prd_sfincs_domain" src="https://github.com/user-attachments/assets/2de7e5b9-a751-4987-b16c-a2be41a4f42c" />
```

Windows users may run `run_prd_example.bat`; Linux and macOS users may run
`./run_prd_example.sh`. See [QUICKSTART.md](QUICKSTART.md) and the
[PRD example guide](example/PRD_single_TC/README.md) for details.

## Author

Anyifang Zhang  
Southern University of Science and Technology  
Contact: `zhangayf@sustech.edu.cn`

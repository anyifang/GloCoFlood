<h1 align="center">Global Compound Flood Modeling Workflow</h1>

<p align="center">
  A script-level workflow for coupling ADCIRC, SFINCS, VIC, and CaMa-Flood to simulate compound coastal flooding driven by storm surge, tides, river discharge, and tropical-cyclone rainfall.
</p>

<p align="center">
  <img alt="MATLAB" src="https://img.shields.io/badge/MATLAB-ADCIRC-orange">
  <img alt="Python" src="https://img.shields.io/badge/Python-SFINCS%20%7C%20VIC-blue">
  <img alt="Shell" src="https://img.shields.io/badge/Shell-HPC%20workflow-lightgrey">
  <img alt="Status" src="https://img.shields.io/badge/status-research%20workflow-informational">
</p>

## Overview

This repository documents the main modeling scripts used to build a global-to-regional compound flood workflow. The workflow links:

- `ADCIRC` for storm surge and tidal water-level modeling.
- `VIC + CaMa-Flood` for rainfall-runoff generation and river discharge routing.
- `SFINCS` for high-resolution coastal and deltaic inundation simulation.

Only the core modeling scripts are included. Large forcing datasets, model executables, generated case folders, NetCDF outputs, rasters, and figures are intentionally excluded.

## Workflow

<img width="5222" height="4132" alt="fig1" src="https://github.com/user-attachments/assets/f97b2529-758b-4898-ad56-21fba8fb3f1f" />

## Example
<img width="3163" height="3499" alt="fig2" src="https://github.com/user-attachments/assets/d1c17bcf-6f36-4ab9-8ef6-5cc7dcc7aa49" />

## Repository Layout

```text
.
|-- ADCIRC/
|   |-- Global_auto.m
|   |-- Global_autofunction/
|   |-- make_fort13_open_boundary_boost.m
|   |-- make_fort14_to_cases.m
|   |-- make_fort15_from_template.m
|   |-- make_fort19_from_tmd.m
|   `-- make_fort22_IBTrACS_C15.m
|-- SFINCS/
|   |-- BoB_building.py
|   |-- Misp_building.py
|   |-- PRD_building.py
|   |-- YRD_building.py
|   `-- sfincs_couple_adcirc_cama_tc_rain_future_batch.py
`-- VIC_CAMAflood/
    |-- auto_cama_vic_coupled.py
    |-- auto_cama_vic_coupled_module_BoB.sh
    |-- auto_cama_vic_coupled_module_Misp.sh
    |-- auto_cama_vic_coupled_module_PRD.sh
    |-- auto_cama_vic_coupled_module_YRD.sh
    |-- Autorun_Camaflood_exp.sh
    `-- plot_vic_runoff.py
```

## Model Components

| Folder | Role | Main outputs |
| --- | --- | --- |
| `ADCIRC/` | Build regional ADCIRC domains, meshes, open boundaries, and storm/tide forcing files. | ADCIRC mesh and `fort.*` input files. |
| `VIC_CAMAflood/` | Prepare VIC domains, meteorological forcing, VIC runoff, and CaMa-Flood routing setup. | VIC runoff and CaMa-Flood river discharge forcing. |
| `SFINCS/` | Build regional SFINCS domains and assemble coupled surge, river, and rainfall forcing. | SFINCS model folders and coupled inundation cases. |

## ADCIRC

`ADCIRC/Global_auto.m` is the main ADCIRC domain and mesh driver. It selects the target estuary or delta, constructs the computational domain, filters coastline data, selects DEM inputs, builds the OceanMesh2D mesh, interpolates bathymetry/topography, and exports the ADCIRC mesh.

Helper functions in `ADCIRC/Global_autofunction/` support automatic domain generation, coastline filtering, fine DEM selection, morphological smoothing, and open-boundary construction.

The `make_fort*.m` scripts prepare ADCIRC case inputs:

| Script | Purpose |
| --- | --- |
| `make_fort13_open_boundary_boost.m` | Generate spatial nodal attributes such as friction and open-boundary adjustments. |
| `make_fort14_to_cases.m` | Copy the generated mesh into case directories. |
| `make_fort15_from_template.m` | Create ADCIRC control files from a template. |
| `make_fort19_from_tmd.m` | Generate tidal boundary forcing from tide model data. |
| `make_fort22_IBTrACS_C15.m` | Generate tropical-cyclone wind and pressure forcing from TC tracks. |

## VIC + CaMa-Flood

`VIC_CAMAflood/auto_cama_vic_coupled.py` is the main Python workflow for building regional VIC and CaMa-Flood inputs. It prepares the regional CaMa-Flood map, creates VIC domain and parameter files, converts ERA5 or prepared meteorological data into VIC forcing, runs VIC, converts VIC runoff into CaMa-Flood forcing, and prepares CaMa-Flood execution scripts.

The region-specific shell workflows wrap the coupled setup for HPC environments:

| Script | Region |
| --- | --- |
| `auto_cama_vic_coupled_module_BoB.sh` | Bay of Bengal |
| `auto_cama_vic_coupled_module_Misp.sh` | Mississippi |
| `auto_cama_vic_coupled_module_PRD.sh` | Pearl River Delta |
| `auto_cama_vic_coupled_module_YRD.sh` | Yangtze River Delta |

`Autorun_Camaflood_exp.sh` is the CaMa-Flood run template used after runoff forcing and regional map files are prepared.

`plot_vic_runoff.py` is kept as an optional quick-look utility for checking VIC runoff fields before routing.

## SFINCS

The SFINCS build scripts create regional high-resolution flood models:

| Script | Region |
| --- | --- |
| `BoB_building.py` | Bay of Bengal |
| `Misp_building.py` | Mississippi |
| `PRD_building.py` | Pearl River Delta |
| `YRD_building.py` | Yangtze River Delta |

These scripts create SFINCS grids, load and adjust terrain, generate model masks and boundary cells, connect river inflow points from CaMa-Flood, prepare land-cover and roughness data, and link external forcing files.

`sfincs_couple_adcirc_cama_tc_rain_future_batch.py` assembles coupled SFINCS cases by combining ADCIRC water levels, CaMa-Flood discharge, and tropical-cyclone rainfall.

## Suggested Execution Order

1. Update local paths, model installation paths, data paths, and HPC module settings.
2. Run `ADCIRC/Global_auto.m` to generate the ADCIRC regional mesh.
3. Run the `ADCIRC/make_fort*.m` scripts to prepare ADCIRC input files.
4. Run the appropriate `VIC_CAMAflood/auto_cama_vic_coupled_module_*.sh` script, or call `auto_cama_vic_coupled.py` directly.
5. Run `VIC_CAMAflood/Autorun_Camaflood_exp.sh` to route runoff with CaMa-Flood.
6. Run the target `SFINCS/*_building.py` script to build the regional SFINCS domain.
7. Run `SFINCS/sfincs_couple_adcirc_cama_tc_rain_future_batch.py` to generate coupled SFINCS forcing and cases.

## External Requirements

The workflow expects local installations or datasets for:

- ADCIRC
- OceanMesh2D
- SFINCS and HydroMT-SFINCS
- VIC Image Driver
- CaMa-Flood
- ERA5 meteorological forcing
- IBTrACS tropical cyclone tracks
- TMD tidal data
- DEM, coastline, land-cover, and river-network datasets

## Notes for Reuse

Many scripts contain hard-coded paths from the original computing environment. Before reuse, check:

- Data root directories.
- Model executable paths.
- HPC module names and queue settings.
- Region names, bounding boxes, and estuary-specific parameters.
- Input/output folder conventions.

## Citation

If you use or adapt this workflow, please cite the associated study or contact the repository author for citation details.

## Author

Anyifang Zhang

# VIC and CaMa-Flood workflows

The portable postprocessing used by the 97-domain production workflow is in
`global_postprocessing/`. It converts routed CaMa-Flood results into the
per-GTC inlet statistics consumed by the SFINCS forcing builders.

`auto_cama_vic_coupled.py` and the four regional shell wrappers are retained
as regional setup/run workflows. Configure their data and executables through
environment variables rather than editing user-specific absolute paths. The
main variables are:

- `CAMA_DIR`;
- `HYDROMT_SFINCS_EXAMPLES`;
- `VIC_EXEC`;
- `VIC_DATA_ROOT`;
- `ERA5_DAILY_ROOT`;
- `VIC_READY_FORCING_PREFIX`;
- `LAND_DEM_NC`.
- `GLOBAL_SOIL_TXT`.

The regional legacy driver can fall back to uniform placeholder soil
properties when `GLOBAL_SOIL_TXT` is absent. That fallback is retained only
to reproduce the original setup workflow and is not suitable for production
or scientific simulations. Supply a validated global soil-parameter table and
check the generated VIC parameter file before running VIC/CaMa-Flood.

HPC module names differ between systems and may also need to be overridden in
the shell environment through `MODULE_INIT_SH`, `VIC_RUNTIME_MODULE`,
`MPI_RUNTIME_MODULE` and `NETCDF_LIB_DIR`. The scripts contain no
institution-specific module or library path. See `data_download/README.md` for
ERA5 request scripts.

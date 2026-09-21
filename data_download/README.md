# Downloading external data

Run these scripts from the repository root. Downloads default to `external/`,
which is ignored by Git.

## Automated downloads

The two `P1_` scripts are independent acquisition tasks and may be run in
either order. Run `P2_check_data_layout.py` after the required downloads and
manual datasets have been placed under `external/`.

Download the official NOAA global IBTrACS archive:

```bash
python data_download/P1_download_ibtracs.py
```

This is the raw best-track archive. The production historical forcing builder
uses the project-specific downscaled track NetCDF containing `vmax_trks`,
`u250_trks`, `v250_trks`, `u850_trks` and `v850_trks`; raw IBTrACS is an input
to that preprocessing and is not a drop-in replacement.

Official source: [NOAA NCEI IBTrACS v04r01 archive](https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/).

Preview monthly ERA5 pressure-level requests used by TCR:

```bash
python data_download/P1_download_cds_era5.py tcr-environment \
  --start-year 1975 --end-year 2014
```

These requests write the annual `era5_t600_monthly_<year>.nc` filenames read
directly by the historical forcing builder. They include temperature at
600 hPa and 250/850-hPa winds used by preprocessing and diagnostics.

Preview six-hourly ERA5 single-level forcing for VIC:

```bash
python data_download/P1_download_cds_era5.py vic-forcing \
  --start-year 2000 --end-year 2000 --area 26 110 20 116
```

ERA5-Land rainfall can be requested with `era5-land-rainfall`. CDS requests
are previews unless `--execute` is supplied. Before execution, install
`data_download/requirements.txt`, create a Copernicus Climate Data Store
account, accept the relevant dataset terms and configure the current CDS API
credentials. Credentials are never accepted as command-line arguments and
must not be committed.

Official CDS datasets:

- [ERA5 monthly averaged pressure-level data](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels-monthly-means)
- [ERA5 hourly single-level data](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels)
- [ERA5-Land hourly data](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land)

## Manual or licence-controlled inputs

The following inputs are intentionally not downloaded automatically:

- FABDEM and other elevation products;
- HydroMT-SFINCS coastline, river and land-cover inputs;
- CaMa-Flood map packages;
- ADCIRC-compatible tidal model data;
- C15 lookup profiles and CLIMADA tropical-cyclone rainfall support data;
- downscaled historical/future synthetic tracks and CMIP6 environmental data.

Relevant official project pages include [FABDEM](https://www.fabdem.org/),
[CaMa-Flood](https://hydro.iis.u-tokyo.ac.jp/~yamadai/cama-flood/) and
[HydroMT-SFINCS](https://deltares.github.io/hydromt_sfincs/latest/). Always
verify the current access terms and required citation at the provider.

Obtain these from their official providers under the applicable terms, then
set their locations in `SFINCS/global_build/config/paths.toml` or through the
environment variables documented in the component READMEs. Do not substitute
lower-resolution look-alike products: grid definitions and variable semantics
must match the production scripts.

Check a proposed external-data layout with:

```bash
python data_download/P2_check_data_layout.py external
```

The checker reports missing inputs but does not modify them.

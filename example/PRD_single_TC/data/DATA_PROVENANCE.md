# Bundled PRD example data

The bundled files form one complete SFINCS input case for `GTC_0009` in the
Pearl River Delta. The forcing period is 4 July 1975 19:00 UTC to 9 July 1975
19:00 UTC. The event identifier is
`track_000261_tracks_GL_era5_197501_201412_000261` and its ADCIRC block is
`ADC_WNP_04`.

The files were extracted from the production GloCoFlood workflow rather than
invented for the example:

- `sfincs.dep`, `sfincs.msk`, `sfincs.ind`, `sfincs.man` and `sfincs.bnd`
  describe the reusable PRD grid and boundary;
- `sfincs.bzs` contains the event water-level boundary derived from ADCIRC
  `fort.63`;
- `sfincs.src` and `sfincs.dis` contain the CaMa-Flood-based river inflow;
- `sfincs_wind.nc`, `sfincs_pressure.nc` and `sfincs_precipitation.nc` contain
  the TC wind, pressure and TCR precipitation forcing generated with the same
  settings as the global historical build.

`MANIFEST.csv` records byte sizes and SHA-256 checksums. Users publishing or
redistributing the example remain responsible for citing and complying with
the licences of SFINCS, ADCIRC, CaMa-Flood, IBTrACS and the atmospheric source
datasets used by their own full reconstruction.

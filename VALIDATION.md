# Integration validation

Validated locally on 2026-09-21 without launching ADCIRC, CaMa-Flood or SFINCS
simulations.

- The repository is based on the downloaded `main` branch of
  `anyifang/GloCoFlood`.
- All 22 Python files pass AST parsing.
- All 14 MATLAB files parse successfully with MATLAB R2024b (`mtree`).
- All 10 shell scripts pass `bash -n` with Git Bash.
- All local Markdown links and all published-manifest paths resolve.
- The versioned ADCIRC catalog contains 41 blocks, the 97-domain membership
  mapping, inner/outer GeoJSON and all block boundary-ring CSV files.
- The fixed PRD example resolves `GTC_0009` to `ADC_WNP_04` and one historical
  event. Its bundled mode was executed end to end and produced a verified
  12-file, 15.83-MiB SFINCS input case without launching the solver.
- A fresh-copy test was run from an unrelated temporary directory. Both the
  documented Python command and `run_prd_example.bat` generated the PRD case;
  every file reference in the generated `sfincs.inp` is case-relative.
- The documented two-package requirements file was also installed into a new
  clean virtual environment; the PRD preparation and NetCDF validation then
  completed successfully with no pre-existing project environment.
- No user-specific absolute path or supplied credential occurs in the release
  directory. Default production paths are repository-relative, and generated
  caches/logs have been removed from the release tree.
- `check_repository.py` validates 97 SFINCS partitions, 41 ADCIRC blocks, 97
  memberships and every checksum in the 12-file bundled PRD case.
- IBTrACS, ERA5 pressure-level, ERA5 single-level and ERA5-Land download
  scripts pass non-mutating request previews.

These checks validate code structure, path resolution and input discovery.
They do not replace scientific validation of the external datasets or a full
production simulation.

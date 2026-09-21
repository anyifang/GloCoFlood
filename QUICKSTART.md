# GloCoFlood quick start

This page is the shortest route from a fresh download to a verified example.
Run all commands from the repository root.

## 1. Check the download

```bash
python check_repository.py
```

The checker uses only the Python standard library. It verifies the published
manifest, the 97 SFINCS partitions, 41 ADCIRC blocks, 97 block memberships and
all SHA-256 checksums in the bundled PRD example. Missing optional scientific
packages are reported as information.

## 2. Prepare the bundled Pearl River Delta case

```bash
python -m pip install -r example/PRD_single_TC/requirements.txt
python example/PRD_single_TC/run_prd_single_tc.py
```

The output is written below `example/PRD_single_TC/output/`. This step does not
need the large global input archives and does not run the SFINCS executable.

To launch the hydrodynamic simulation as well:

```bash
python example/PRD_single_TC/run_prd_single_tc.py --sfincs-executable external/bin/sfincs
```

Windows users can alternatively run `run_prd_example.bat`; Linux and macOS
users can run `./run_prd_example.sh` after making it executable.

## 3. Choose the full workflow you need

- Build the 97 reusable model grids: see `SFINCS/global_build/README.md`.
- Build ADCIRC blocks and TC forcing: see `ADCIRC/global_build/README.md`.
- Prepare CaMa-Flood inflows: see
  `VIC_CAMAflood/global_postprocessing/README.md`.
- Assemble historical and future compound-forcing cases: see
  `GLOBAL_WORKFLOW.md` and `SFINCS/global_build/config/paths.example.toml`.

The full workflows require external model executables and datasets. They are
not downloaded automatically because of their size and separate licences.
See `data_download/README.md` for the supported IBTrACS and ERA5 download
scripts and the list of licence-controlled inputs that must be obtained
manually.

## Repository conventions

- `global_build` contains the portable production workflows.
- Older top-level regional scripts are retained for provenance and may still
  require site-specific path edits.
- Never commit `paths.toml`, credentials, model outputs or downloaded external
  datasets; these are excluded by `.gitignore`.

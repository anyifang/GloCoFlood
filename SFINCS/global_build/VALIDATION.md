# Release validation

Validated on 2026-09-20 with the production data inventory and the
`sfincs_tcr` Python environment.

Checks completed:

- the reusable-domain builder and plotting helper compile successfully;
- the included partition GeoJSON/CSV load all 97 coastal domains;
- a one-domain `--dry-run --limit 1 --no-plot` resolves `GTC_0064` without
  requiring external terrain or CaMa-Flood data;
- all Python files compile successfully;
- the TOML configuration loader resolves and validates every required input;
- `--dry-run --limit-events 1 --workers 4` starts all four scenario builders;
- historical, SSP1-2.6, SSP2-4.5 and SSP3-7.0 all exit successfully;
- the inventory contains 97 mapped/runnable SFINCS models in 41 ADCIRC blocks;
- the tested future ADCIRC inventories each contain all 1,425 expected non-empty
  `fort.63` files;
- the C15 library contains 19,186 lookup profiles;
- the PRD one-TC wrapper resolves exactly one mapped domain (`GTC_0009`), one
  ADCIRC block (`ADC_WNP_04`) and one fixed historical event in dry-run mode;
- its output verifier successfully opens an existing production instance of
  that PRD case and confirms all required SFINCS and meteorological files;
- the smoke test writes no runnable SFINCS case.

This validates the code path and input discovery, not the scientific contents
of external datasets and not a complete 3,849-case build.

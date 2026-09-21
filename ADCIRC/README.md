# ADCIRC workflows

The supported public workflow is under [`global_build/`](global_build/). It
contains the scripts used to define the 41 global ADCIRC blocks, build meshes
and static `fort.*` inputs, select historical tropical cyclones, generate
event forcing, stage return-period events, and rebuild future scenarios.

Start with [`global_build/README.md`](global_build/README.md), then copy
`global_build/configure_paths.example.m` to an untracked local
`configure_paths.m` and set the external data and executable paths.

`Global_autofunction/` contains reusable geometry and coastline helper
functions retained from the regional domain-building workflow. The original
regional driver and its site-specific `fort.*` preparation scripts are not
part of this portable release; their production replacements are the scripts
in `global_build/`.

Large meshes, TMD tidal data, C15 lookup data, TC-track archives, ADCIRC
executables and generated event directories are external inputs and are not
stored in this repository.

%% Portable path configuration for the global ADCIRC workflow.
% Copy this file to configure_paths.m. All defaults are constructed from the
% repository root; edit only the relative layout if your external/ directory
% differs. Run it in the same MATLAB session as the global-build script.

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
external_root = fullfile(repo_root, 'external');
work_root = fullfile(repo_root, 'work');

setenv('ADCIRC_TC_MESH_ROOT', fullfile(external_root, 'adcirc', 'adcirc_fort14_meshes_stable'));
setenv('ADCIRC_PARTITION_INPUT_ROOT', fullfile(repo_root, 'SFINCS', 'global_build', 'global_tc_coastal_model_partitions'));
setenv('ADCIRC_BLOCK_OUTPUT_ROOT', fullfile(work_root, 'global_tc_adcirc_model_blocks'));
setenv('ADCIRC_BLOCK_CATALOG_ROOT', fullfile(repo_root, 'ADCIRC', 'global_build', 'catalog'));
setenv('OCEANMESH2D_ROOT', fullfile(external_root, 'OceanMesh2D-Projection'));
setenv('ADCIRC_MESH_OUTPUT_ROOT', fullfile(work_root, 'adcirc_fort14_meshes'));
setenv('ADCIRC_MESH_FIGURE_ROOT', fullfile(work_root, 'adcirc_mesh_figures'));
setenv('ADCIRC_STATIC_SOURCE_MESH_ROOT', fullfile(work_root, 'adcirc_fort14_meshes'));
setenv('ADCIRC_STATIC_OUTPUT_MESH_ROOT', fullfile(work_root, 'adcirc_fort14_meshes_stable'));
setenv('ADCIRC_FORT15_TEMPLATE', fullfile(external_root, 'adcirc', 'templates', 'fort.15'));
setenv('ADCIRC_SUBMIT_TEMPLATE', fullfile(repo_root, 'ADCIRC', 'global_build', 'sub_intel.example.sh'));
setenv('ADCIRC_TIDELEVEL_ROOT', fullfile(external_root, 'tide_gauge_validation_data'));
setenv('ADCIRC_FORT13_CF_UPRANGE_SHP', fullfile(external_root, 'adcirc', 'CF_uprange.shp'));
setenv('ADCIRC_TC_BLOCKS_CSV', fullfile(repo_root, 'ADCIRC', 'global_build', 'catalog', 'global_tc_adcirc_blocks.csv'));
setenv('ADCIRC_TC_EVENT_CSV', fullfile(external_root, 'tracks', 'era5_tc_events_by_adcirc_inner_domain_vmax33_matlab.csv'));
setenv('ADCIRC_TC_TRACK_NC', fullfile(external_root, 'tracks', 'tracks_GL_era5_197501_201412.nc'));
setenv('ADCIRC_TC_BASIN_CV_CSV', fullfile(external_root, 'C15_predata', 'accepted_Cv_distribution_summary.csv'));
setenv('ADCIRC_TC_WIND_PRESSURE_FIT_CSV', fullfile(external_root, 'C15_predata', 'TC_wind_pressure_tailfit_1980_2020.csv'));
setenv('ADCIRC_TC_C15_PREDATA_DIR', fullfile(external_root, 'C15_predata'));
setenv('ADCIRC_TC_TMD_ROOT', fullfile(external_root, 'tides', 'TMD2.5', 'TMD'));
setenv('ADCIRC_TC_CASE_OUTPUT_ROOT', fullfile(work_root, 'adcirc_era5_tc_selected_run_files'));

% External MATLAB worker preparation.
setenv('ADCIRC_TC_PARALLEL_WORKERS', '8');
setenv('ADCIRC_TC_PARALLEL_TMD_TEMPLATE', fullfile(external_root, 'tides', 'TMD2.5'));
setenv('ADCIRC_TC_PARALLEL_WORKER_ROOT', fullfile(work_root, 'adcirc_workers'));
setenv('ADCIRC_TC_PARALLEL_LAUNCH', '0');

% Return-period event staging.
setenv('PREPARE_RP_RERUN_RUNFILE_ROOT', fullfile(work_root, 'adcirc_era5_tc_selected_run_files'));
setenv('PREPARE_RP_RERUN_SELECTION_DIR', fullfile(external_root, 'adcirc', 'block_return_period_tc_selection'));
setenv('PREPARE_RP_RERUN_OUTPUT_ROOT', fullfile(work_root, 'return_period_tc_event_reruns'));
setenv('PREPARE_RP_RERUN_SUBMIT_TEMPLATE', fullfile(repo_root, 'ADCIRC', 'global_build', 'sub_intel.example.sh'));

% Future forcing rebuild.
setenv('ADCIRC_FUTURE_SOURCE_ROOT', fullfile(work_root, 'return_period_tc_event_reruns'));
setenv('ADCIRC_FUTURE_OUTPUT_ROOT', fullfile(work_root, 'return_period_tc_event_reruns'));
setenv('ADCIRC_FUTURE_EXACT_INNER_FORCING_DIR', fullfile(external_root, 'adcirc', 'exact_inner_tc_return_period_forcing'));
setenv('ADCIRC_FUTURE_NUM_WORKERS', '8');

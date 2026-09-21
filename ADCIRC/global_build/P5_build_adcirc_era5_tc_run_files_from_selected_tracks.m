%% P5_build_adcirc_era5_tc_run_files_from_selected_tracks.m
% Build ADCIRC run directories for ERA5 downscaled TC events selected by
% P4_count_and_plot_era5_tc_for_adcirc_inner_domains.m.
%
% Each output case contains fort.13, fort.14, fort.15, fort.19, and fort.22.
% Large static files are hard-linked by default; fort.15 is copied and
% relabeled for each case; fort.22 is regenerated with C15 TC winds.
%
% Useful environment overrides:
%   ADCIRC_TC_EVENT_CSV        full path to event CSV
%   ADCIRC_TC_CASE_OUTPUT_ROOT full path to output directory
%   ADCIRC_TC_CASE_BLOCKS      comma-separated block IDs, e.g. ADC_WNP_05,ADC_NATL_01
%   ADCIRC_TC_CASE_MAX_CASES   maximum number of cases for testing
%   ADCIRC_TC_CASE_MAX_PER_BLOCK maximum cases per ADCIRC block for testing
%   ADCIRC_TC_CASE_OVERWRITE   1/true to rebuild existing complete cases
%   ADCIRC_TC_TMD_ROOT         worker-specific TMD root, e.g. external/tides/TMD2.50/TMD
%   ADCIRC_TC_GENERATE_FORT19_FROM_TMD 0/false to skip fort.19 generation
%   ADCIRC_TC_SUMMARY_CSV_NAME worker-specific summary CSV name
%   ADCIRC_TC_CONFIG_MAT_NAME  worker-specific config MAT name
%   ADCIRC_TC_WORKER_ID        optional worker label written to logs/config
%   ADCIRC_TC_RANDOM_DAY_SEED   common deterministic TC calendar seed

clearvars;
clc;

SCRIPT_DIR = fileparts(mfilename('fullpath'));
if isempty(SCRIPT_DIR)
    SCRIPT_DIR = pwd;
end

%% =========================== USER SETTINGS =============================
% Change the variables in this block for normal use. The functions below
% should not need manual edits unless the file format itself changes.

P = struct();

% Input and output paths.
P.mesh_root = fullfile(SCRIPT_DIR, 'external', 'adcirc_fort14_meshes_stable');
P.blocks_csv = fullfile(SCRIPT_DIR, 'global_tc_adcirc_blocks.csv');
P.event_csv = fullfile(SCRIPT_DIR, 'era5_tc_events_by_adcirc_inner_domain_vmax33_matlab.csv');
P.track_nc = fullfile(SCRIPT_DIR, 'external', 'tracks_GL_era5_197501_201412.nc');
P.output_root = fullfile(SCRIPT_DIR, 'output', 'adcirc_era5_tc_selected_run_files');
P.worker_id = "";

% Wind model and auxiliary data.
% Cv is taken from the accepted mean produced by
% valid_U10_timevaryCvfix_C15_global_all_stat_min5.m.
P.basin_cv_csv = fullfile(SCRIPT_DIR, 'external', 'accepted_Cv_distribution_summary.csv');
P.basin_cv_mat = '';
P.wind_pressure_fit_csv = fullfile(SCRIPT_DIR, 'external', 'TC_wind_pressure_tailfit_1980_2020.csv');
P.c15_predata_dir = fullfile(SCRIPT_DIR, 'external', 'C15_predata');

% Tide generation. fort.19 is rebuilt from TMD for each TC calendar window.
P.generate_fort19_from_tmd = true;
P.tmd_root = fullfile(SCRIPT_DIR, 'external', 'TMD2.5', 'TMD');
P.tmd_latlon_file = fullfile(P.tmd_root, 'LAT_LON', 'lat_lon');
P.tmd_out_dir = fullfile(P.tmd_root, 'OUT');
P.clean_tmd_before_each_case = true;
P.allow_tmd_trim_or_pad = true;

% Case selection controls.
P.selected_blocks = strings(0, 1);  % e.g. ["ADC_WNP_05"; "ADC_NATL_01"]; empty = all blocks.
P.max_cases = inf;                  % global maximum case count; inf = no global limit.
P.max_cases_per_block = inf;         % first N selected TCs per block; inf = no per-block limit.
P.overwrite = true;                 % false = skip complete existing cases.
P.fast_block_resume = true;          % true = resume each block from its last existing case.
P.use_environment_overrides = true; % true = allow ADCIRC_TC_* environment variables to override these settings.

% Case directory and file behavior.
P.static_file_mode = 'hardlink';     % hardlink | copy
P.required_static_files = {'fort.13', 'fort.14'};
P.optional_helper_files = {'sub_intel.sh'};
P.case_complete_files = {'fort.13', 'fort.14', 'fort.15', 'fort.19', 'fort.22', 'fort22_meta.txt', 'storm_track_forcing_meta.mat'};
P.fort15_name = 'fort.15';
P.fort22_name = 'fort.22';
P.meta_file_name = 'fort22_meta.txt';
P.track_meta_mat_name = 'storm_track_forcing_meta.mat';
P.completion_marker_name = 'case_complete.ok';
P.min_fort19_record_bytes = 8;
P.min_fort22_record_bytes = 10;
P.summary_csv_name = 'adcirc_era5_tc_case_build_summary.csv';
P.config_mat_name = 'adcirc_era5_tc_case_build_config.mat';
P.case_name_prefix = 'track';
P.case_name_max_chars = 100;
P.fort15_run_desc_max_chars = 32;

% Time window. Each track's global LMI is aligned once to a deterministic
% calendar date. Every block keeps that same absolute TC timeline, then uses
% its own max_time_index to define the local impact time and +/- run window.
P.pre_impact_days = 3;
P.post_impact_days = 2;
P.random_day_seed = 20260705;
P.time_alignment_scheme = 'global_lmi_anchor_block_impact_v2';
P.default_track_dt_seconds = 3600;

% Meteorological grid used by fort.22. When force_rebuild_met_grid is true,
% dlon/dlat are applied even if the mesh block already has fort22_meta.txt.
P.force_rebuild_met_grid = true;
P.fallback_met_margin_deg = 0.20;
P.fallback_dlon = 0.10;
P.fallback_dlat = 0.10;
P.fallback_window_start = datetime(2021, 1, 1, 0, 0, 0);
P.fallback_window_days = 28;
P.sync_fort15_output_end_to_rnday = true;
P.force_fort15_fort63_output = true;
P.fort15_fort63_output_start_mode = 'fixed'; % fixed | impact | source
P.fort15_fort63_output_start_day = 2.5;
P.fort15_fort63_output_end_day = 5.0;
P.auto_update_fort15_dramp = true;

% Wind-pressure and C15 model controls.
P.wind_pressure_fit_mode = 'track_basin'; % track_basin | global
P.cv_weight_field = 'n_accepted_tc';
P.Pn_hPa = 1013;
P.wpr_pc_min_hpa = 850;
P.wpr_pc_max_hpa = P.Pn_hPa - 1;
P.Re = 6371000;
P.Vmax_list = 15:1:120;
P.Rmax_list = 20:1:200;
P.standard_pressure_pa = 101300;
P.background_fort22_line = '0 0 101300';
P.fort22_zero_uv_tolerance = 1e-8;
P.fort22_zero_pressure_tolerance_pa = 1e-4;

% Tunable constants inside the C15 implementation.
P.c15_min_lookup_vm_ms = 15;
P.c15_max_lookup_vm_ms = 120;
P.c15_min_rmax_m = 30e3;
P.c15_holland_b_min = 1.0;
P.c15_holland_b_max = 2.5;
P.c15_one_min_to_ten_min_factor = 0.893;
P.c15_blend_inner_radius_m = 500e3;
P.c15_blend_outer_radius_m = 700e3;
P.c15_beta_inner_base_deg = 10;
P.c15_beta_inner_slope_deg = 10;
P.c15_beta_mid_base_deg = 20;
P.c15_beta_mid_slope_deg = 25;
P.c15_beta_outer_deg = 25;
P.c15_beta_mid_radius_factor = 1.2;

% Derived output paths. Usually do not edit these names; edit output_root
% and the *_name variables above instead.
P.summary_csv = fullfile(P.output_root, P.summary_csv_name);
P.config_mat = fullfile(P.output_root, P.config_mat_name);

%% ========================= END USER SETTINGS ============================

if P.use_environment_overrides
    P = apply_environment_overrides(P);
end

assert(exist(P.event_csv, 'file') == 2, 'Missing event CSV: %s', P.event_csv);
assert(exist(P.track_nc, 'file') == 2, 'Missing track NetCDF: %s', P.track_nc);
assert(exist(P.mesh_root, 'dir') == 7, 'Missing mesh root: %s', P.mesh_root);
assert(exist(P.blocks_csv, 'file') == 2, 'Missing block CSV: %s', P.blocks_csv);
assert(exist(P.wind_pressure_fit_csv, 'file') == 2, 'Missing wind-pressure fit CSV: %s', P.wind_pressure_fit_csv);
assert(exist(P.c15_predata_dir, 'dir') == 7, 'Missing C15 predata directory: %s', P.c15_predata_dir);
if P.generate_fort19_from_tmd
    assert_tmd_ready(P.tmd_root, P.tmd_latlon_file);
end

ensure_dir(P.output_root);
save(P.config_mat, 'P', '-v7.3');

Events = readtable(P.event_csv, 'TextType', 'string');
Events = validate_and_filter_events(Events, P);
EventsBeforeResume = height(Events);
ResumePlan = table();
if ~P.overwrite && P.fast_block_resume
    [Events, ResumePlan] = fast_block_resume_events(Events, P);
end
fprintf('\n============================================================\n');
fprintf('Build ADCIRC ERA5 TC run files\n');
fprintf('Event CSV : %s\n', P.event_csv);
fprintf('Track NC  : %s\n', P.track_nc);
fprintf('Mesh root : %s\n', P.mesh_root);
fprintf('Output    : %s\n', P.output_root);
fprintf('Cases     : %d\n', height(Events));
if ~P.overwrite && P.fast_block_resume
    fprintf('FastResume: kept %d/%d cases after timing-aware block resume checks\n', ...
        height(Events), EventsBeforeResume);
    if height(ResumePlan) > 0
        disp(ResumePlan(:, {'block_id', 'n_events', 'last_existing_position', ...
            'resume_position', 'n_cases_to_build', 'status'}));
    end
end
if isfinite(P.max_cases_per_block)
    fprintf('Per block : first %d events per block\n', P.max_cases_per_block);
end
fprintf('Static    : %s for fort.13/fort.14; copy fort.15\n', P.static_file_mode);
if P.generate_fort19_from_tmd
    fprintf('Tide      : global LMI calendar anchor uses deterministic seed %d; each block uses its local impact window\n', ...
        P.random_day_seed);
    fprintf('TMD root  : %s\n', P.tmd_root);
else
    fprintf('Tide      : fort.19 must already exist in the case directory or be added manually\n');
end
if strlength(string(P.worker_id)) > 0
    fprintf('Worker    : %s\n', string(P.worker_id));
end
fprintf('Window    : impact - %.2f days to impact + %.2f days (RNDAY=%.2f)\n', ...
    P.pre_impact_days, P.post_impact_days, P.pre_impact_days + P.post_impact_days);
fprintf('Time rule : %s\n', P.time_alignment_scheme);
fprintf('Overwrite : %d\n', P.overwrite);
fprintf('============================================================\n');

if height(Events) == 0
    SummaryRows = {};
    Summary = summary_rows_to_table(SummaryRows);
    writetable(Summary, P.summary_csv);
    save(P.config_mat, 'P', 'Events', 'Summary', 'ResumePlan', '-v7.3');
    fprintf('\nDone. No cases need to be built after fast block resume.\nSummary: %s\n', P.summary_csv);
    return;
end

fprintf('\nLoading C15 constants and lookup tables...\n');
[Cv_global, CvTable] = load_global_weighted_mean_Cv(P.basin_cv_csv, P.basin_cv_mat, P.cv_weight_field);
WindPressureFits = load_wind_pressure_fits(P.wind_pressure_fit_csv, P.wpr_pc_min_hpa, P.wpr_pc_max_hpa);
WindC15Lib = load_wind_c15_library(P.c15_predata_dir, P.Vmax_list, P.Rmax_list);
fprintf('  Cv_global = %.4f (source: %s; field: %s)\n', ...
    Cv_global, P.basin_cv_csv, P.cv_weight_field);
fprintf('  wind-pressure fit mode: %s\n', P.wind_pressure_fit_mode);
fprintf('  C15 lookup tables loaded: %d\n', WindC15Lib.Count);

trackInfo = prepare_track_nc_info(P.track_nc, P.default_track_dt_seconds);

SummaryRows = {};
for i = 1:height(Events)
    event = table2struct(Events(i, :));
    blockId = char(event.block_id);
    tcId = char(event.track_id);
    trackIndex = double(event.track_index);
    caseName = make_case_name(trackIndex, tcId, P);
    caseDir = fullfile(P.output_root, blockId, caseName);
    blockDir = fullfile(P.mesh_root, blockId);

    fprintf('\n[%d/%d] %s | %s\n', i, height(Events), blockId, tcId);
    row = make_empty_summary_row(event, caseDir);

    try
        if exist(blockDir, 'dir') ~= 7
            error('Missing mesh block directory: %s', blockDir);
        end
        ensure_dir(caseDir);

        if ~P.overwrite && ~P.fast_block_resume
            [isComplete, completeMessage, completeInfo] = case_is_complete(caseDir, blockDir, event, P);
            if isComplete
                fprintf('  -> complete case exists, skipped (%s)\n', char(completeMessage));
                row = fill_summary_row_from_complete_info(row, completeInfo, caseDir, P);
                row.status = "skipped_existing";
                row.message = "";
                SummaryRows(end + 1, :) = summary_struct_to_cell(row); %#ok<SAGROW>
                close_case_file_handles(P, char(row.status));
                continue;
            elseif strlength(completeMessage) > 0
                fprintf('  -> existing case incomplete: %s; rebuilding\n', char(completeMessage));
            end
        end

        clear_case_completion_marker(caseDir, P);

        [Grid, Window, sourceMeta] = read_case_grid_and_window(blockDir, P);

        WindPressureFit = select_wind_pressure_fit_for_track(event, trackInfo, WindPressureFits, P);
        fprintf('  -> wind-pressure fit: track basin %s -> %s (A=%.4g, B=%.4g, Pref=%.1f hPa)\n', ...
            WindPressureFit.track_basin_code, WindPressureFit.basin, ...
            WindPressureFit.A, WindPressureFit.B, WindPressureFit.prefHpa);

        tc = read_era5_track_as_tc(P.track_nc, trackInfo, event, WindPressureFit);
        Window = set_tc_calendar_window(Window, event, trackInfo, tc, P);
        tc = assign_tc_absolute_time(tc, Window);
        fprintf('  -> fort.22 grid: NWLON=%d, NWLAT=%d, WTIMINC=%d s, nSteps=%d\n', ...
            Grid.NWLON, Grid.NWLAT, Window.WTIMINC, numel(Window.time_data));
        fprintf('  -> TC global LMI anchor: index=%d, time=%s\n', ...
            Window.global_anchor_index, datestr(Window.global_anchor_time, 31));
        fprintf('  -> block impact: index=%d, time=%s; local window=%s to %s\n', ...
            tc.event_time_index, datestr(Window.impact_time, 31), ...
            datestr(Window.start, 31), datestr(Window.end, 31));

        link_or_copy_static_files(blockDir, caseDir, P);
        copy_and_update_fort15(blockDir, caseDir, blockId, trackIndex, tcId, event, Window, Grid, P);
        copy_optional_helper_files(blockDir, caseDir, P);
        [fort22Path, maxWindDomain, maxWindTime, activeSteps] = write_fort22_for_case( ...
            fullfile(caseDir, P.fort22_name), tc, Grid, Window, Cv_global, WindC15Lib, WindPressureFit, P);
        if P.generate_fort19_from_tmd
            fort14Path = fullfile(caseDir, 'fort.14');
            G = read_fort14_nodes_and_open_boundary(fort14Path);
            generate_fort19_from_tmd_for_window(G, fullfile(caseDir, 'fort.19'), Window, P);
        end

        write_case_meta(fullfile(caseDir, P.meta_file_name), blockId, event, tc, ...
            Grid, Window, sourceMeta, Cv_global, CvTable, WindPressureFit, ...
            fort22Path, maxWindDomain, maxWindTime, activeSteps, P);

        save(fullfile(caseDir, P.track_meta_mat_name), ...
            'P', 'event', 'tc', 'Grid', 'Window', 'sourceMeta', ...
            'Cv_global', 'WindPressureFit', 'maxWindDomain', 'maxWindTime', 'activeSteps', '-v7.3');
        write_case_completion_marker(caseDir, event, Grid, Window, fort22Path, activeSteps, maxWindDomain, P);

        row.status = "built";
        row.message = "";
        row.case_dir = string(caseDir);
        row.fort22 = string(fort22Path);
        row.window_start = string(datestr(Window.start, 31));
        row.window_end = string(datestr(Window.end, 31));
        row.nsteps = numel(Window.time_data);
        row.nwlon = Grid.NWLON;
        row.nwlat = Grid.NWLAT;
        row.active_steps = activeSteps;
        row.domain_max_wind_ms = maxWindDomain;
        row.domain_max_wind_time = string(fmt_time_safe(maxWindTime));
        fprintf('  -> built, active wind steps=%d, domain max wind=%.2f m/s\n', activeSteps, maxWindDomain);
    catch ME
        row.status = "failed";
        row.message = string(ME.message);
        fprintf('  !! failed: %s\n', ME.message);
    end

    close_case_file_handles(P, char(row.status));
    SummaryRows(end + 1, :) = summary_struct_to_cell(row); %#ok<SAGROW>
end

Summary = summary_rows_to_table(SummaryRows);
writetable(Summary, P.summary_csv);
save(P.config_mat, 'P', 'Events', 'Summary', 'ResumePlan', '-v7.3');

fprintf('\nDone.\nSummary: %s\n', P.summary_csv);

%% ========================================================================
function P = apply_environment_overrides(P)

v = strtrim(string(getenv('ADCIRC_TC_WORKER_ID')));
if strlength(v) > 0
    P.worker_id = v;
end

v = strtrim(string(getenv('ADCIRC_TC_EVENT_CSV')));
if strlength(v) > 0
    P.event_csv = char(v);
end

v = strtrim(string(getenv('ADCIRC_TC_BLOCKS_CSV')));
if strlength(v) > 0
    P.blocks_csv = char(v);
end

v = strtrim(string(getenv('ADCIRC_TC_TRACK_NC')));
if strlength(v) > 0
    P.track_nc = char(v);
end

v = strtrim(string(getenv('ADCIRC_TC_BASIN_CV_CSV')));
if strlength(v) > 0
    P.basin_cv_csv = char(v);
end

v = strtrim(string(getenv('ADCIRC_TC_BASIN_CV_MAT')));
if strlength(v) > 0
    P.basin_cv_mat = char(v);
end

v = strtrim(string(getenv('ADCIRC_TC_WIND_PRESSURE_FIT_CSV')));
if strlength(v) > 0
    P.wind_pressure_fit_csv = char(v);
end

v = strtrim(string(getenv('ADCIRC_TC_C15_PREDATA_DIR')));
if strlength(v) > 0
    P.c15_predata_dir = char(v);
end

v = strtrim(string(getenv('ADCIRC_TC_MESH_ROOT')));
if strlength(v) > 0
    P.mesh_root = char(v);
end

v = strtrim(string(getenv('ADCIRC_TC_CASE_OUTPUT_ROOT')));
if strlength(v) > 0
    P.output_root = char(v);
end

v = strtrim(string(getenv('ADCIRC_TC_CASE_BLOCKS')));
if strlength(v) > 0
    parts = split(v, ',');
    parts = strip(parts);
    P.selected_blocks = parts(strlength(parts) > 0);
end

v = str2double(strtrim(string(getenv('ADCIRC_TC_CASE_MAX_CASES'))));
if isfinite(v) && v >= 0
    P.max_cases = floor(v);
end

v = str2double(strtrim(string(getenv('ADCIRC_TC_CASE_MAX_PER_BLOCK'))));
if isfinite(v) && v >= 0
    P.max_cases_per_block = floor(v);
end

v = strtrim(string(getenv('ADCIRC_TC_CASE_OVERWRITE')));
if any(strcmpi(v, ["1", "true", "yes", "y"]))
    P.overwrite = true;
elseif any(strcmpi(v, ["0", "false", "no", "n"]))
    P.overwrite = false;
end

v = strtrim(string(getenv('ADCIRC_TC_FAST_BLOCK_RESUME')));
if any(strcmpi(v, ["1", "true", "yes", "y", "on"]))
    P.fast_block_resume = true;
elseif any(strcmpi(v, ["0", "false", "no", "n", "off"]))
    P.fast_block_resume = false;
end

v = strtrim(string(getenv('ADCIRC_TC_GENERATE_FORT19_FROM_TMD')));
if any(strcmpi(v, ["1", "true", "yes", "y", "on"]))
    P.generate_fort19_from_tmd = true;
elseif any(strcmpi(v, ["0", "false", "no", "n", "off"]))
    P.generate_fort19_from_tmd = false;
end

v = strtrim(string(getenv('ADCIRC_TC_TMD_ROOT')));
if strlength(v) > 0
    P.tmd_root = char(v);
    P.tmd_latlon_file = fullfile(P.tmd_root, 'LAT_LON', 'lat_lon');
    P.tmd_out_dir = fullfile(P.tmd_root, 'OUT');
end

v = strtrim(string(getenv('ADCIRC_TC_SUMMARY_CSV_NAME')));
if strlength(v) > 0
    P.summary_csv_name = char(v);
end

v = strtrim(string(getenv('ADCIRC_TC_CONFIG_MAT_NAME')));
if strlength(v) > 0
    P.config_mat_name = char(v);
end

v = str2double(strtrim(string(getenv('ADCIRC_TC_FORT63_START_DAY'))));
if isfinite(v)
    P.fort15_fort63_output_start_day = v;
    P.fort15_fort63_output_start_mode = 'fixed';
end

v = str2double(strtrim(string(getenv('ADCIRC_TC_FORT63_END_DAY'))));
if isfinite(v)
    P.fort15_fort63_output_end_day = v;
end

v = str2double(strtrim(string(getenv('ADCIRC_TC_RANDOM_DAY_SEED'))));
if isfinite(v) && v >= 0
    P.random_day_seed = floor(v);
end

P.summary_csv = fullfile(P.output_root, P.summary_csv_name);
P.config_mat = fullfile(P.output_root, P.config_mat_name);
end

%% ========================================================================
function Events = validate_and_filter_events(Events, P)

need = ["block_id", "track_id", "track_index", "year", "max_inner_vmax_ms", ...
    "lmi_ms", "max_lon", "max_lat", "max_time_index"];
missing = setdiff(need, string(Events.Properties.VariableNames));
if ~isempty(missing)
    error('Event CSV missing columns: %s', strjoin(missing, ', '));
end

Events.block_id = string(Events.block_id);
Events.track_id = string(Events.track_id);
Events.track_index = double(Events.track_index);
Events.max_time_index = double(Events.max_time_index);

ok = isfinite(Events.track_index) & Events.track_index >= 1 & isfinite(Events.max_time_index);
Events = Events(ok, :);

if ~isempty(P.selected_blocks)
    Events = Events(ismember(Events.block_id, P.selected_blocks), :);
end

if height(Events) > 0
    key = Events.block_id + "|" + string(round(Events.track_index));
    [~, ia] = unique(key, 'stable');
    Events = Events(ia, :);
end

Events = order_events_by_selected_blocks(Events, P.selected_blocks);

if isfinite(P.max_cases_per_block)
    Events = limit_events_per_block(Events, P.max_cases_per_block);
end

Events = group_events_by_block(Events, P.selected_blocks);

if isfinite(P.max_cases) && height(Events) > P.max_cases
    Events = Events(1:P.max_cases, :);
end
end

%% ========================================================================
function Events = limit_events_per_block(Events, maxPerBlock)

if isempty(Events) || ~isfinite(maxPerBlock)
    return;
end
maxPerBlock = floor(maxPerBlock);
if maxPerBlock <= 0
    Events = Events([], :);
    return;
end

keep = false(height(Events), 1);
blocks = unique(Events.block_id, 'stable');
for i = 1:numel(blocks)
    idx = find(Events.block_id == blocks(i));
    idx = idx(1:min(maxPerBlock, numel(idx)));
    keep(idx) = true;
end
Events = Events(keep, :);
end

%% ========================================================================
function Events = group_events_by_block(Events, selectedBlocks)

if isempty(Events)
    return;
end

if nargin >= 2 && ~isempty(selectedBlocks)
    blocks = string(selectedBlocks(:));
    blocks = blocks(ismember(blocks, Events.block_id));
else
    blocks = unique(Events.block_id, 'stable');
end
order = [];
for i = 1:numel(blocks)
    order = [order; find(Events.block_id == blocks(i))]; %#ok<AGROW>
end
Events = Events(order, :);
end

%% ========================================================================
function Events = order_events_by_selected_blocks(Events, selectedBlocks)

if isempty(Events) || isempty(selectedBlocks)
    return;
end

selectedBlocks = string(selectedBlocks(:));
selectedBlocks = selectedBlocks(strlength(selectedBlocks) > 0);
if isempty(selectedBlocks)
    return;
end

order = [];
for i = 1:numel(selectedBlocks)
    order = [order; find(Events.block_id == selectedBlocks(i))]; %#ok<AGROW>
end
if numel(order) < height(Events)
    used = false(height(Events), 1);
    used(order) = true;
    order = [order; find(~used)]; %#ok<AGROW>
end
Events = Events(order, :);
end

%% ========================================================================
function [EventsOut, ResumePlan] = fast_block_resume_events(Events, P)

ResumePlan = empty_resume_plan();
EventsOut = Events;
if isempty(Events)
    return;
end

blocks = unique(Events.block_id, 'stable');
keep = false(height(Events), 1);
Rows = {};

for ib = 1:numel(blocks)
    blockId = blocks(ib);
    idx = find(Events.block_id == blockId);
    nBlockEvents = numel(idx);
    blockDir = fullfile(P.mesh_root, char(blockId));

    lastExistingLocal = 0;
    lastExistingCaseDir = "";
    for j = 1:nBlockEvents
        event = table2struct(Events(idx(j), :));
        caseDir = event_case_dir(event, P);
        if exist(caseDir, 'dir') == 7
            lastExistingLocal = j;
            lastExistingCaseDir = string(caseDir);
        end
    end

    if lastExistingLocal == 0
        resumeLocal = 1;
        status = "no_existing_cases";
        msg = "build from first event";
        lastComplete = false;
    else
        event = table2struct(Events(idx(lastExistingLocal), :));
        [lastComplete, msg] = case_is_complete(char(lastExistingCaseDir), blockDir, event, P);
        if lastComplete
            resumeLocal = lastExistingLocal + 1;
            if resumeLocal > nBlockEvents
                status = "block_complete";
            else
                status = "last_existing_complete";
            end
        elseif contains(string(msg), "time alignment scheme mismatch")
            % Old directories may exist for the entire block. During the
            % migration, find the first case that does not carry a completed
            % marker for the current timing scheme. This preserves restart
            % progress without accepting any legacy-aligned case.
            resumeLocal = first_case_without_current_timing_marker(Events, idx, P);
            status = "time_alignment_upgrade";
            if resumeLocal > nBlockEvents
                status = "block_complete";
                lastComplete = true;
                msg = "all cases have current timing completion markers";
            else
                msg = sprintf('resume timing upgrade at block position %d', resumeLocal);
            end
        else
            resumeLocal = lastExistingLocal;
            status = "last_existing_incomplete";
        end
    end

    if resumeLocal <= nBlockEvents
        keep(idx(resumeLocal:end)) = true;
        resumeGlobal = idx(resumeLocal);
        nCasesToBuild = nBlockEvents - resumeLocal + 1;
        resumeCaseDir = string(event_case_dir(table2struct(Events(resumeGlobal, :)), P));
    else
        resumeGlobal = NaN;
        nCasesToBuild = 0;
        resumeCaseDir = "";
    end

    Rows(end + 1, :) = {blockId, nBlockEvents, lastExistingLocal, lastExistingCaseDir, ...
        logical(lastComplete), resumeLocal, resumeGlobal, resumeCaseDir, ...
        nCasesToBuild, status, string(msg)}; %#ok<AGROW>
end

EventsOut = Events(keep, :);
if ~isempty(Rows)
    ResumePlan = cell2table(Rows, 'VariableNames', resume_plan_var_names());
end
end

%% ========================================================================
function resumeLocal = first_case_without_current_timing_marker(Events, idx, P)

resumeLocal = 1;
for j = 1:numel(idx)
    event = table2struct(Events(idx(j), :));
    caseDir = event_case_dir(event, P);
    markerFile = fullfile(caseDir, P.completion_marker_name);
    if exist(markerFile, 'file') ~= 2
        resumeLocal = j;
        return;
    end
    requiredScheme = sprintf('time_alignment_scheme = %s', P.time_alignment_scheme);
    if ~text_file_contains(markerFile, 'status = complete') || ...
            ~text_file_contains(markerFile, requiredScheme)
        resumeLocal = j;
        return;
    end
    resumeLocal = j + 1;
end
end

%% ========================================================================
function T = empty_resume_plan()

T = cell2table(cell(0, numel(resume_plan_var_names())), 'VariableNames', resume_plan_var_names());
end

%% ========================================================================
function names = resume_plan_var_names()

names = {'block_id', 'n_events', 'last_existing_position', 'last_existing_case_dir', ...
    'last_existing_complete', 'resume_position', 'resume_global_index', ...
    'resume_case_dir', 'n_cases_to_build', 'status', 'message'};
end

%% ========================================================================
function caseDir = event_case_dir(event, P)

blockId = char(event.block_id);
tcId = char(event.track_id);
trackIndex = double(event.track_index);
caseName = make_case_name(trackIndex, tcId, P);
caseDir = fullfile(P.output_root, blockId, caseName);
end

%% ========================================================================
function trackInfo = prepare_track_nc_info(ncFile, defaultDtSeconds)

years = clean_fill_to_nan(double(ncread(ncFile, 'tc_years')));
years = years(:);
trackInfo.nTrack = numel(years);
trackInfo.years = years;
trackInfo.months = read_optional_nc_vector(ncFile, 'tc_month', trackInfo.nTrack);
trackInfo.basins = read_optional_nc_string_vector(ncFile, 'tc_basins', trackInfo.nTrack);
[trackInfo.trackDim, trackInfo.nTime] = infer_track_dim_from_nc(ncFile, 'lon_trks', trackInfo.nTrack);

if any(strcmp({ncinfo(ncFile).Variables.Name}, 'time'))
    t = clean_fill_to_nan(double(ncread(ncFile, 'time')));
    t = t(:);
else
    t = [];
end
if numel(t) ~= trackInfo.nTime || ~any(isfinite(t))
    t = (0:trackInfo.nTime-1).' * defaultDtSeconds;
end
if numel(t) >= 2
    dt = median(diff(t(isfinite(t))));
    if ~isfinite(dt) || dt <= 0
        dt = defaultDtSeconds;
    end
else
    dt = defaultDtSeconds;
end
trackInfo.time_seconds = t(:);
trackInfo.dt_seconds = dt;
end

%% ========================================================================
function v = read_optional_nc_vector(ncFile, varName, nTrack)

names = {ncinfo(ncFile).Variables.Name};
if ~any(strcmp(names, varName))
    v = nan(nTrack, 1);
    return;
end
v = clean_fill_to_nan(double(ncread(ncFile, varName)));
v = v(:);
if numel(v) < nTrack
    v(end + 1:nTrack, 1) = NaN;
elseif numel(v) > nTrack
    v = v(1:nTrack);
end
end

%% ========================================================================
function v = read_optional_nc_string_vector(ncFile, varName, nTrack)

names = {ncinfo(ncFile).Variables.Name};
if ~any(strcmp(names, varName))
    v = strings(nTrack, 1);
    return;
end

raw = ncread(ncFile, varName);
if isstring(raw)
    v = raw(:);
elseif iscell(raw)
    v = string(raw(:));
elseif ischar(raw)
    if size(raw, 1) == nTrack
        v = string(cellstr(raw));
    elseif size(raw, 2) == nTrack
        v = string(cellstr(raw.'));
    else
        v = string(raw(:));
    end
else
    v = string(raw(:));
end
v = upper(strtrim(v));
if numel(v) < nTrack
    v(end + 1:nTrack, 1) = "";
elseif numel(v) > nTrack
    v = v(1:nTrack);
end
end

%% ========================================================================
function [trackDim, nTime] = infer_track_dim_from_nc(ncFile, varName, nTrack)

info = ncinfo(ncFile, varName);
sz = info.Size;
if numel(sz) ~= 2
    error('%s must be 2-D in %s.', varName, ncFile);
end
if sz(1) == nTrack
    trackDim = 1;
    nTime = sz(2);
elseif sz(2) == nTrack
    trackDim = 2;
    nTime = sz(1);
else
    error('Cannot infer track dimension for %s in %s.', varName, ncFile);
end
end

%% ========================================================================
function [Grid, Window, sourceMeta] = read_case_grid_and_window(blockDir, P)

metaFile = fullfile(blockDir, P.meta_file_name);
sourceMeta = struct('meta_file', metaFile, 'source', "fort22_meta");
if exist(metaFile, 'file') == 2
    meta = read_key_value_meta(metaFile);

    Window = struct();
    Window.start = parse_meta_datetime(get_meta_str(meta, 'WindowStart'));
    Window.end = parse_meta_datetime(get_meta_str(meta, 'WindowEnd'));
    Window.WTIMINC = round(get_meta_num(meta, 'WTIMINC'));
    nSteps = get_meta_num(meta, 'nSteps');
    if ~isfinite(nSteps) || nSteps <= 0
        nSteps = round(seconds(Window.end - Window.start) / Window.WTIMINC) + 1;
    end
    Window.time_data = Window.start + seconds((0:round(nSteps)-1).' .* Window.WTIMINC);

    if isfield(P, 'force_rebuild_met_grid') && P.force_rebuild_met_grid
        fort14 = fullfile(blockDir, 'fort.14');
        Mesh = read_fort14_nodes_only(fort14);
        Grid = build_met_grid_from_mesh(Mesh, P.fallback_met_margin_deg, P.fallback_dlon, P.fallback_dlat);
        sourceMeta.source = "fort14_grid_meta_window";
    else
        Grid = struct();
        Grid.NWLON = round(get_meta_num(meta, 'NWLON'));
        Grid.NWLAT = round(get_meta_num(meta, 'NWLAT'));
        Grid.WLONMIN = get_meta_num(meta, 'WLONMIN');
        Grid.WLATMAX = get_meta_num(meta, 'WLATMAX');
        Grid.DLON = get_meta_num(meta, 'WLONINC');
        Grid.DLAT = get_meta_num(meta, 'WLATINC');
        Grid.met_lon = Grid.WLONMIN + (0:Grid.NWLON-1) .* Grid.DLON;
        Grid.met_lat = Grid.WLATMAX - (0:Grid.NWLAT-1) .* Grid.DLAT;
        [Grid.MET_LON, Grid.MET_LAT] = meshgrid(Grid.met_lon, Grid.met_lat);
    end
    return;
end

fort14 = fullfile(blockDir, 'fort.14');
Mesh = read_fort14_nodes_only(fort14);
Grid = build_met_grid_from_mesh(Mesh, P.fallback_met_margin_deg, P.fallback_dlon, P.fallback_dlat);
Window = struct();
Window.start = P.fallback_window_start;
Window.end = Window.start + days(P.fallback_window_days);
Window.WTIMINC = P.default_track_dt_seconds;
Window.time_data = Window.start + seconds((0:round(seconds(Window.end - Window.start) / Window.WTIMINC)).' .* Window.WTIMINC);
sourceMeta.source = "fallback_from_fort14";
end

%% ========================================================================
function Window = set_tc_calendar_window(Window, event, trackInfo, tc, P)

trackIndex = round(double(event.track_index));
tcYear = safe_index(trackInfo.years, trackIndex);
tcMonth = safe_index(trackInfo.months, trackIndex);

if ~isfinite(tcYear) && isfield(event, 'year')
    tcYear = double(event.year);
end
if ~isfinite(tcMonth) && isfield(event, 'month')
    tcMonth = double(event.month);
end

tcYear = round(tcYear);
tcMonth = round(tcMonth);
if ~isfinite(tcYear) || tcYear < 1800 || tcYear > 2300
    error('Cannot determine valid TC year for track_index=%d.', trackIndex);
end
if ~isfinite(tcMonth) || tcMonth < 1 || tcMonth > 12
    error('Cannot determine valid TC month for track_index=%d.', trackIndex);
end

if ~isfinite(Window.WTIMINC) || Window.WTIMINC <= 0
    Window.WTIMINC = P.default_track_dt_seconds;
end

randomDay = deterministic_random_day(tcYear, tcMonth, trackIndex, P.random_day_seed);
globalAnchorTime = datetime(tcYear, tcMonth, randomDay, 0, 0, 0);
impactTime = globalAnchorTime + seconds(tc.event_relative_seconds);
windowStart = impactTime - days(P.pre_impact_days);
windowEnd = impactTime + days(P.post_impact_days);
durationSeconds = seconds(windowEnd - windowStart);
nSteps = round(durationSeconds / Window.WTIMINC) + 1;
if abs((nSteps - 1) * Window.WTIMINC - durationSeconds) > 1e-6
    warning('Run duration %.3f s is not exactly divisible by WTIMINC=%d s; final time will follow WindowEnd.', ...
        durationSeconds, Window.WTIMINC);
end

Window.source_start = Window.start;
Window.source_end = Window.end;
Window.tc_year = tcYear;
Window.tc_month = tcMonth;
Window.random_day = randomDay;
Window.time_alignment_scheme = P.time_alignment_scheme;
Window.global_anchor_index = tc.global_anchor_index;
Window.global_anchor_time = globalAnchorTime;
Window.block_impact_index = tc.event_time_index;
Window.impact_time = impactTime;
Window.pre_impact_days = P.pre_impact_days;
Window.post_impact_days = P.post_impact_days;
Window.start = windowStart;
Window.end = windowEnd;
Window.time_data = Window.start + seconds((0:nSteps-1).' .* Window.WTIMINC);
end

%% ========================================================================
function tc = assign_tc_absolute_time(tc, Window)

tc.time = Window.global_anchor_time + seconds(tc.relative_time_seconds);
tc.event_aligned_time = Window.impact_time;
tc.global_anchor_time = Window.global_anchor_time;
tc.time_alignment_scheme = Window.time_alignment_scheme;
end

%% ========================================================================
function day = deterministic_random_day(tcYear, tcMonth, trackIndex, seed0)

maxDay = eomday(tcYear, tcMonth);
% The calendar day must be intrinsic to the TC. Do not include block_id:
% the same track can force several ADCIRC blocks and must retain one tide phase.
key = sprintf('%d_%04d_%02d_%d', round(trackIndex), tcYear, tcMonth, round(seed0));
bytes = double(uint8(key));
hash = uint32(2166136261);
for i = 1:numel(bytes)
    hash = bitxor(hash, uint32(bytes(i)));
    hash = uint32(mod(double(hash) * 16777619, 2^32));
end
seed = double(mod(hash, uint32(2^31 - 1)));
oldState = rng;
cleanupObj = onCleanup(@() rng(oldState)); %#ok<NASGU>
rng(seed, 'twister');
day = randi(maxDay);
end

%% ========================================================================
function meta = read_key_value_meta(metaFile)

txt = fileread(metaFile);
lines = regexp(txt, '\r\n|\n|\r', 'split');
meta = containers.Map('KeyType', 'char', 'ValueType', 'char');
for i = 1:numel(lines)
    line = strtrim(lines{i});
    if isempty(line) || startsWith(line, '[')
        continue;
    end
    tok = regexp(line, '^([^=]+?)\s*=\s*(.*)$', 'tokens', 'once');
    if isempty(tok)
        continue;
    end
    key = regexprep(strtrim(tok{1}), '\s+', '');
    val = strtrim(tok{2});
    meta(key) = val;
end
end

%% ========================================================================
function s = get_meta_str(meta, key)

if isKey(meta, key)
    s = meta(key);
else
    s = '';
end
end

%% ========================================================================
function x = get_meta_num(meta, key)

x = NaN;
s = get_meta_str(meta, key);
if isempty(s)
    return;
end
x = sscanf(s, '%f', 1);
if isempty(x)
    x = NaN;
end
end

%% ========================================================================
function t = parse_meta_datetime(s)

s = strtrim(char(s));
if isempty(s)
    t = NaT;
    return;
end
try
    t = datetime(s, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
catch
    t = datetime(s);
end
try
    t.TimeZone = '';
catch
end
end

%% ========================================================================
function [tf, msg, info] = case_is_complete(caseDir, blockDir, event, P)

tf = false;
msg = "";
info = struct('fort22', string(fullfile(caseDir, P.fort22_name)), ...
    'window_start', "", 'window_end', "", 'nsteps', NaN, 'nwlon', NaN, ...
    'nwlat', NaN, 'active_steps', NaN, 'domain_max_wind_ms', NaN, ...
    'domain_max_wind_time', "");

if exist(caseDir, 'dir') ~= 7
    msg = "case directory missing";
    return;
end

need = P.case_complete_files;
for i = 1:numel(need)
    name = need{i};
    f = fullfile(caseDir, name);
    bytes = file_size_bytes(f);
    if ~isfinite(bytes)
        msg = sprintf('missing %s', name);
        return;
    end
    minBytes = 1;
    if strcmpi(name, P.track_meta_mat_name)
        minBytes = 1024;
    end
    if bytes < minBytes
        msg = sprintf('%s is too small (%g bytes)', name, bytes);
        return;
    end
end

for i = 1:numel(P.required_static_files)
    name = P.required_static_files{i};
    src = fullfile(blockDir, name);
    dst = fullfile(caseDir, name);
    srcBytes = file_size_bytes(src);
    dstBytes = file_size_bytes(dst);
    if ~isfinite(srcBytes)
        msg = sprintf('source %s is missing', name);
        return;
    end
    if srcBytes ~= dstBytes
        msg = sprintf('%s size mismatch with source block', name);
        return;
    end
end

matFile = fullfile(caseDir, P.track_meta_mat_name);
try
    matVars = string({whos('-file', matFile).name});
catch ME
    msg = sprintf('cannot read %s: %s', P.track_meta_mat_name, ME.message);
    return;
end
requiredMatVars = ["P", "event", "tc", "Grid", "Window"];
missingMatVars = requiredMatVars(~ismember(requiredMatVars, matVars));
if ~isempty(missingMatVars)
    msg = sprintf('%s missing variables: %s', P.track_meta_mat_name, strjoin(missingMatVars, ', '));
    return;
end

metaFile = fullfile(caseDir, P.meta_file_name);
try
    meta = read_key_value_meta(metaFile);
catch ME
    msg = sprintf('cannot read %s: %s', P.meta_file_name, ME.message);
    return;
end

if ~text_file_contains(metaFile, '[ForcingTrackPoints]')
    msg = sprintf('%s is missing final section', P.meta_file_name);
    return;
end

blockMeta = string(get_meta_str(meta, 'block'));
if strlength(blockMeta) == 0 || blockMeta ~= string(event.block_id)
    msg = 'metadata block does not match event';
    return;
end

schemeMeta = string(get_meta_str(meta, 'time_alignment_scheme'));
if strlength(schemeMeta) == 0 || schemeMeta ~= string(P.time_alignment_scheme)
    msg = sprintf('time alignment scheme mismatch: found "%s", require "%s"', ...
        char(schemeMeta), P.time_alignment_scheme);
    return;
end

trackIdMeta = string(get_meta_str(meta, 'track_id'));
if strlength(trackIdMeta) == 0 || trackIdMeta ~= string(event.track_id)
    msg = 'metadata track_id does not match event';
    return;
end

trackIndexMeta = round(get_meta_num(meta, 'track_index'));
if ~isfinite(trackIndexMeta) || trackIndexMeta ~= round(double(event.track_index))
    msg = 'metadata track_index does not match event';
    return;
end

nSteps = round(get_meta_num(meta, 'nSteps'));
nwlon = round(get_meta_num(meta, 'NWLON'));
nwlat = round(get_meta_num(meta, 'NWLAT'));
wtiminc = round(get_meta_num(meta, 'WTIMINC'));
if ~all(isfinite([nSteps, nwlon, nwlat, wtiminc])) || any([nSteps, nwlon, nwlat, wtiminc] <= 0)
    msg = 'metadata has invalid nSteps/NWLON/NWLAT/WTIMINC';
    return;
end

fort15File = fullfile(caseDir, P.fort15_name);
fort15Token = sprintf('track_index=%d', round(double(event.track_index)));
if ~text_file_contains(fort15File, fort15Token)
    msg = sprintf('%s does not contain expected track_index marker', P.fort15_name);
    return;
end

fort22File = fullfile(caseDir, P.fort22_name);
fort22Lines = double(nSteps) .* double(nwlon) .* double(nwlat);
fort22Bytes = file_size_bytes(fort22File);
if fort22Bytes < fort22Lines .* P.min_fort22_record_bytes
    msg = sprintf('%s is shorter than expected for %d steps on %dx%d grid', ...
        P.fort22_name, nSteps, nwlat, nwlon);
    return;
end
if ~text_file_ends_with_newline(fort22File) || ~last_nonempty_line_has_numbers(fort22File, 3)
    msg = sprintf('%s tail is incomplete or malformed', P.fort22_name);
    return;
end

fort19File = fullfile(caseDir, 'fort.19');
try
    nOpen = count_fort14_open_boundary_nodes_cached(fullfile(blockDir, 'fort.14'));
catch ME
    msg = sprintf('cannot count open-boundary nodes: %s', ME.message);
    return;
end
fort19Rows = 1 + double(nSteps) .* double(nOpen);
fort19Bytes = file_size_bytes(fort19File);
if fort19Bytes < fort19Rows .* P.min_fort19_record_bytes
    msg = sprintf('fort.19 is shorter than expected for %d steps and %d open nodes', nSteps, nOpen);
    return;
end
firstFort19 = first_nonempty_line(fort19File);
firstFort19Value = sscanf(firstFort19, '%f', 1);
if isempty(firstFort19Value) || abs(firstFort19Value - wtiminc) > 0.5
    msg = 'fort.19 first line does not match WTIMINC';
    return;
end
if ~text_file_ends_with_newline(fort19File) || ~last_nonempty_line_has_numbers(fort19File, 1)
    msg = 'fort.19 tail is incomplete or malformed';
    return;
end

info.fort22 = string(fort22File);
info.window_start = string(get_meta_str(meta, 'WindowStart'));
info.window_end = string(get_meta_str(meta, 'WindowEnd'));
info.nsteps = nSteps;
info.nwlon = nwlon;
info.nwlat = nwlat;
info.active_steps = get_meta_num(meta, 'active_wind_steps');
info.domain_max_wind_ms = get_meta_num(meta, 'domain_max_wind_ms');
info.domain_max_wind_time = string(get_meta_str(meta, 'domain_max_wind_time'));

markerFile = fullfile(caseDir, P.completion_marker_name);
if exist(markerFile, 'file') == 2
    msg = "validated marker";
else
    msg = "validated legacy files";
end
tf = true;
end

%% ========================================================================
function row = fill_summary_row_from_complete_info(row, info, caseDir, P)

row.case_dir = string(caseDir);
row.fort22 = string(fullfile(caseDir, P.fort22_name));
if isfield(info, 'fort22') && strlength(info.fort22) > 0
    row.fort22 = info.fort22;
end
row.window_start = info.window_start;
row.window_end = info.window_end;
row.nsteps = info.nsteps;
row.nwlon = info.nwlon;
row.nwlat = info.nwlat;
row.active_steps = info.active_steps;
row.domain_max_wind_ms = info.domain_max_wind_ms;
row.domain_max_wind_time = info.domain_max_wind_time;
end

%% ========================================================================
function clear_case_completion_marker(caseDir, P)

markerFile = fullfile(caseDir, P.completion_marker_name);
if exist(markerFile, 'file') == 2
    try
        delete(markerFile);
    catch
    end
end
end

%% ========================================================================
function write_case_completion_marker(caseDir, event, Grid, Window, fort22Path, activeSteps, maxWindDomain, P)

markerFile = fullfile(caseDir, P.completion_marker_name);
fid = fopen(markerFile, 'wt');
if fid < 0
    error('Cannot write completion marker: %s', markerFile);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'status = complete\n');
fprintf(fid, 'completed_local = %s\n', datestr(now, 31));
fprintf(fid, 'block = %s\n', char(event.block_id));
fprintf(fid, 'track_id = %s\n', char(event.track_id));
fprintf(fid, 'track_index = %d\n', round(double(event.track_index)));
fprintf(fid, 'time_alignment_scheme = %s\n', P.time_alignment_scheme);
fprintf(fid, 'global_anchor_index = %d\n', Window.global_anchor_index);
fprintf(fid, 'global_anchor_time = %s\n', datestr(Window.global_anchor_time, 31));
fprintf(fid, 'block_impact_index = %d\n', Window.block_impact_index);
fprintf(fid, 'block_impact_time = %s\n', datestr(Window.impact_time, 31));
fprintf(fid, 'fort22 = %s\n', fort22Path);
fprintf(fid, 'nSteps = %d\n', numel(Window.time_data));
fprintf(fid, 'NWLON = %d\n', Grid.NWLON);
fprintf(fid, 'NWLAT = %d\n', Grid.NWLAT);
fprintf(fid, 'WTIMINC = %d\n', Window.WTIMINC);
fprintf(fid, 'active_wind_steps = %d\n', activeSteps);
fprintf(fid, 'domain_max_wind_ms = %.6f\n', maxWindDomain);
fprintf(fid, 'generator = P5_build_adcirc_era5_tc_run_files_from_selected_tracks.m\n');
end

%% ========================================================================
function bytes = file_size_bytes(path)

D = dir(path);
if isempty(D) || D(1).isdir
    bytes = NaN;
else
    bytes = double(D(1).bytes);
end
end

%% ========================================================================
function tf = text_file_contains(path, token)

tf = false;
try
    txt = fileread(path);
catch
    return;
end
tf = contains(txt, token);
end

%% ========================================================================
function tf = text_file_ends_with_newline(path)

tf = false;
fid = fopen(path, 'r');
if fid < 0
    return;
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
if fseek(fid, 0, 'eof') ~= 0
    return;
end
n = ftell(fid);
if n <= 0 || fseek(fid, -1, 'eof') ~= 0
    return;
end
b = fread(fid, 1, 'uint8');
tf = ~isempty(b) && any(double(b) == [10, 13]);
end

%% ========================================================================
function tf = last_nonempty_line_has_numbers(path, nRequired)

line = last_nonempty_line(path);
vals = sscanf(line, '%f');
tf = numel(vals) >= nRequired;
end

%% ========================================================================
function line = first_nonempty_line(path)

line = '';
fid = fopen(path, 'rt');
if fid < 0
    return;
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
while true
    s = fgetl(fid);
    if ~ischar(s)
        return;
    end
    s = strtrim(s);
    if ~isempty(s)
        line = s;
        return;
    end
end
end

%% ========================================================================
function line = last_nonempty_line(path)

line = '';
fid = fopen(path, 'r');
if fid < 0
    return;
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
if fseek(fid, 0, 'eof') ~= 0
    return;
end
fileBytes = ftell(fid);
if fileBytes <= 0
    return;
end
chunkBytes = min(fileBytes, 65536);
if fseek(fid, -chunkBytes, 'eof') ~= 0
    return;
end
bytes = fread(fid, chunkBytes, '*uint8').';
txt = char(bytes);
lines = regexp(txt, '\r\n|\n|\r', 'split');
for i = numel(lines):-1:1
    s = strtrim(lines{i});
    if ~isempty(s)
        line = s;
        return;
    end
end
end

%% ========================================================================
function nOpen = count_fort14_open_boundary_nodes_cached(fort14_file)

persistent cache
if isempty(cache)
    cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
end
key = char(fort14_file);
if isKey(cache, key)
    nOpen = cache(key);
    return;
end
nOpen = count_fort14_open_boundary_nodes(fort14_file);
cache(key) = nOpen;
end

%% ========================================================================
function nOpen = count_fort14_open_boundary_nodes(fort14_file)

fid = fopen(fort14_file, 'rt');
if fid < 0
    error('Cannot open fort.14: %s', fort14_file);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fgetl(fid);
line2 = fgetl(fid);
v = sscanf(strtrim(line2), '%f');
if numel(v) < 2
    error('Cannot parse fort.14 line 2.');
end
NE = round(v(1));
NP = round(v(2));
for i = 1:NP
    if ~ischar(fgetl(fid))
        error('Unexpected EOF in fort.14 node table.');
    end
end
for i = 1:NE
    if ~ischar(fgetl(fid))
        error('Unexpected EOF in fort.14 element table.');
    end
end

s = fgetl(fid);
if ~ischar(s)
    error('Cannot read open-boundary count from fort.14.');
end
NOPE = sscanf(strtrim(s), '%d', 1);
if isempty(NOPE)
    NOPE = 0;
end

s = fgetl(fid);
if ~ischar(s)
    error('Cannot read open-boundary node total from fort.14.');
end
NETA = sscanf(strtrim(s), '%d', 1);
if isempty(NETA)
    NETA = 0;
end

nOpen = 0;
for ib = 1:NOPE
    header = fgetl(fid);
    if ~ischar(header)
        error('Unexpected EOF before open-boundary header %d.', ib);
    end
    hv = sscanf(strtrim(header), '%d', 1);
    if isempty(hv)
        error('Cannot parse open-boundary header %d.', ib);
    end
    nvdll = round(hv);
    nOpen = nOpen + nvdll;
    for k = 1:nvdll
        if ~ischar(fgetl(fid))
            error('Unexpected EOF in open-boundary %d.', ib);
        end
    end
end

if NETA > 0 && nOpen ~= NETA
    error('fort.14 open-boundary count mismatch: headers sum to %d, NETA=%d.', nOpen, NETA);
end
end

%% ========================================================================
function link_or_copy_static_files(blockDir, caseDir, P)

for i = 1:numel(P.required_static_files)
    name = P.required_static_files{i};
    src = fullfile(blockDir, name);
    dst = fullfile(caseDir, name);
    if exist(src, 'file') ~= 2
        error('Missing required static file: %s', src);
    end
    if exist(dst, 'file') == 2
        if P.overwrite
            delete(dst);
        else
            continue;
        end
    end
    if strcmpi(P.static_file_mode, 'copy')
        copyfile(src, dst);
    else
        ok = create_hardlink(dst, src);
        if ~ok
            warning('Hardlink failed for %s; copying instead.', name);
            copyfile(src, dst);
        end
    end
end
end

%% ========================================================================
function ok = create_hardlink(linkPath, targetPath)

ok = false;
if ispc
    cmd = sprintf('cmd /c mklink /H "%s" "%s"', linkPath, targetPath);
else
    cmd = sprintf('ln "%s" "%s"', targetPath, linkPath);
end
[status, ~] = system(cmd);
ok = (status == 0) && exist(linkPath, 'file') == 2;
end

%% ========================================================================
function copy_optional_helper_files(blockDir, caseDir, P)

for i = 1:numel(P.optional_helper_files)
    name = P.optional_helper_files{i};
    src = fullfile(blockDir, name);
    dst = fullfile(caseDir, name);
    if exist(src, 'file') ~= 2
        continue;
    end
    copyfile(src, dst);
end
end

%% ========================================================================
function copy_and_update_fort15(blockDir, caseDir, blockId, trackIndex, tcId, event, Window, Grid, P)

src = fullfile(blockDir, P.fort15_name);
dst = fullfile(caseDir, P.fort15_name);
if exist(src, 'file') ~= 2
    error('Missing fort.15: %s', src);
end
txt = fileread(src);
lines = regexp(txt, '\r\n|\n|\r', 'split');
if ~isempty(lines) && isempty(lines{end})
    lines(end) = [];
end

desc = sprintf('%s TC%06d', blockId, round(trackIndex));
if numel(desc) > P.fort15_run_desc_max_chars
    desc = desc(1:P.fort15_run_desc_max_chars);
end
if ~isempty(lines)
    lines{1} = sprintf(' %-32s ! 32 CHARACTER ALPHANUMERIC RUN DESCRIPTION', desc);
end

idx_nws = find_line_contains(lines, '! NWS - WIND STRESS AND BAROMETRIC PRESSURE OPTION PARAMETER');
if ~isempty(idx_nws)
    lines{idx_nws} = ' 6                                   ! NWS - WIND STRESS AND BAROMETRIC PRESSURE OPTION PARAMETER';
end

idx_met = find_line_any(lines, {'WTIMINC', 'WITMINC', 'STIMINC'});
if ~isempty(idx_met)
    lines{idx_met} = sprintf( ...
        ' %d %d %.2f %.2f %.2f %.2f %d    ! NWLAT NWLON WLATMAX WLONMIN WLATINC WLONINC WTIMINC', ...
        Grid.NWLAT, Grid.NWLON, Grid.WLATMAX, Grid.WLONMIN, Grid.DLAT, Grid.DLON, Window.WTIMINC);
end

RNDAY = days(Window.end - Window.start);
idx_rnday = find_line_contains(lines, '! RNDAY - TOTAL LENGTH OF SIMULATION (IN DAYS)');
if ~isempty(idx_rnday)
    lines{idx_rnday} = sprintf(' %.8f                                 ! RNDAY - TOTAL LENGTH OF SIMULATION (IN DAYS)', RNDAY);
end

if P.auto_update_fort15_dramp
    idx_dramp = find_line_contains(lines, '! DRAMP - DURATION OF RAMP FUNCTION (IN DAYS)');
    if ~isempty(idx_dramp)
        DRAMP = min(1.0, RNDAY);
        lines{idx_dramp} = sprintf(' %.8f                                 ! DRAMP - DURATION OF RAMP FUNCTION (IN DAYS)', DRAMP);
    end
end

if P.sync_fort15_output_end_to_rnday
    output_tokens = { ...
        'NOUTE,TOUTSE,TOUTFE', ...
        'NOUTV,TOUTSV,TOUTFV', ...
        'NOUTM,TOUTSM,TOUTFM', ...
        'NOUTGV,TOUTSGV,TOUTFGV', ...
        'NOUTGW,TOUTSGW,TOUTFGW'};
    for k = 1:numel(output_tokens)
        idx_out = find_line_contains(lines, output_tokens{k});
        if ~isempty(idx_out)
            lines{idx_out} = update_output_schedule_line(lines{idx_out}, RNDAY);
        end
    end
end

idx_ge = find_line_contains(lines, 'NOUTGE,TOUTSGE,TOUTFGE');
if ~isempty(idx_ge)
    lines{idx_ge} = update_fort63_output_schedule_line(lines{idx_ge}, Window, RNDAY, P);
end

lines{end + 1} = sprintf('! ERA5 TC case generated by P5_build_adcirc_era5_tc_run_files_from_selected_tracks.m');
lines{end + 1} = sprintf('! block=%s ; track_id=%s ; track_index=%d', blockId, tcId, round(trackIndex));
lines{end + 1} = sprintf('! selected_inner_vmax_ms=%.3f ; selected_lon=%.6f ; selected_lat=%.6f', ...
    double(event.max_inner_vmax_ms), double(event.max_lon), double(event.max_lat));
lines{end + 1} = sprintf('! TimeAlignment=%s', Window.time_alignment_scheme);
lines{end + 1} = sprintf('! TCYear=%d ; TCMonth=%02d ; GlobalAnchorDay=%02d ; GlobalLMIAnchor=%s', ...
    Window.tc_year, Window.tc_month, Window.random_day, datestr(Window.global_anchor_time, 31));
lines{end + 1} = sprintf('! GlobalLMIIndex=%d ; BlockImpactIndex=%d ; BlockImpact=%s', ...
    Window.global_anchor_index, Window.block_impact_index, datestr(Window.impact_time, 31));
lines{end + 1} = sprintf('! WindowStart=%s ; WindowEnd=%s ; WTIMINC=%d', ...
    datestr(Window.start, 31), datestr(Window.end, 31), Window.WTIMINC);
fort63StartDay = fort63_output_start_day(Window, RNDAY, P);
fort63EndDay = fort63_output_end_day(Window, RNDAY, P);
lines{end + 1} = sprintf('! fort63_output_window_days=%.8f to %.8f', ...
    fort63StartDay, fort63EndDay);
lines{end + 1} = sprintf('! NWLAT=%d ; NWLON=%d ; WLATMAX=%.6f ; WLONMIN=%.6f ; WLATINC=%.6f ; WLONINC=%.6f', ...
    Grid.NWLAT, Grid.NWLON, Grid.WLATMAX, Grid.WLONMIN, Grid.DLAT, Grid.DLON);

fid = fopen(dst, 'wt');
if fid < 0
    error('Cannot write %s', dst);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
for i = 1:numel(lines)
    fprintf(fid, '%s\n', lines{i});
end
end

%% ========================================================================
function tc = read_era5_track_as_tc(ncFile, trackInfo, event, WindPressureFit)

trackIndex = round(double(event.track_index));
if trackIndex < 1 || trackIndex > trackInfo.nTrack
    error('Track index %d outside NetCDF track range 1..%d.', trackIndex, trackInfo.nTrack);
end
nTime = trackInfo.nTime;

if trackInfo.trackDim == 1
    start = [trackIndex 1];
    count = [1 nTime];
    lon = clean_fill_to_nan(double(ncread(ncFile, 'lon_trks', start, count)));
    lat = clean_fill_to_nan(double(ncread(ncFile, 'lat_trks', start, count)));
    wind = clean_fill_to_nan(double(ncread(ncFile, 'vmax_trks', start, count)));
else
    start = [1 trackIndex];
    count = [nTime 1];
    lon = clean_fill_to_nan(double(ncread(ncFile, 'lon_trks', start, count))).';
    lat = clean_fill_to_nan(double(ncread(ncFile, 'lat_trks', start, count))).';
    wind = clean_fill_to_nan(double(ncread(ncFile, 'vmax_trks', start, count))).';
end

lon = wrap_to_180_local(lon(:));
lat = lat(:);
wind = wind(:);
timeSeconds = trackInfo.time_seconds(:);
if numel(timeSeconds) ~= numel(lon)
    timeSeconds = (0:numel(lon)-1).' .* trackInfo.dt_seconds;
end

valid = isfinite(lon) & isfinite(lat) & isfinite(wind) & isfinite(timeSeconds) & ...
    lat >= -90 & lat <= 90;
if ~any(valid)
    error('Track %d has no finite lon/lat/wind samples.', trackIndex);
end

origIndex = find(valid);
lon = lon(valid);
lat = lat(valid);
wind = wind(valid);
timeSeconds = timeSeconds(valid);

eventIndex = round(double(event.max_time_index));
eventIndex = choose_valid_event_index(eventIndex, origIndex, lon, lat, wind, event);
eventPos = find(origIndex == eventIndex, 1, 'first');
if isempty(eventPos)
    [~, eventPos] = min(abs(origIndex - eventIndex));
end

[timeSeconds, order] = sort(timeSeconds);
lon = lon(order);
lat = lat(order);
wind = wind(order);
origIndex = origIndex(order);
eventPos = find(origIndex == eventIndex, 1, 'first');
if isempty(eventPos)
    [~, eventPos] = min(abs(origIndex - eventIndex));
    eventIndex = origIndex(eventPos);
end

% A track has one block-independent absolute timeline. Anchor the global
% maximum intensity (LMI) to the deterministic calendar date, then retain
% every original time offset. A block's max_time_index only selects that
% block's local impact time; it never shifts the track itself.
[~, globalAnchorPos] = max(wind);
globalAnchorIndex = origIndex(globalAnchorPos);
relativeTimeSeconds = timeSeconds - timeSeconds(globalAnchorPos);

lon_unwrapped = rad2deg(unwrap(deg2rad(lon)));
pressure_hpa = wind_to_pressure_hpa(wind, WindPressureFit);

tc = struct();
tc.time = NaT(numel(timeSeconds), 1);
tc.relative_time_seconds = relativeTimeSeconds(:);
tc.lon = wrap_to_180_local(lon_unwrapped(:));
tc.lon_unwrapped = lon_unwrapped(:);
tc.lat = lat(:);
tc.wind = wind(:);
tc.p = pressure_hpa(:);
tc.name = char(event.track_id);
tc.sid = char(event.track_id);
tc.year = safe_index(trackInfo.years, trackIndex);
tc.month = safe_index(trackInfo.months, trackIndex);
tc.basin = safe_index(trackInfo.basins, trackIndex);
tc.pressure_fit_basin = char(WindPressureFit.basin);
tc.pressure_fit_A = WindPressureFit.A;
tc.pressure_fit_B = WindPressureFit.B;
tc.pressure_fit_pref_hPa = WindPressureFit.prefHpa;
tc.track_index = trackIndex;
tc.original_time_index = origIndex(:);
tc.event_time_index = eventIndex;
tc.event_relative_seconds = relativeTimeSeconds(eventPos);
tc.event_aligned_time = NaT;
tc.global_anchor_index = globalAnchorIndex;
tc.global_anchor_position = globalAnchorPos;
tc.global_anchor_vmax_ms = wind(globalAnchorPos);
tc.global_anchor_time = NaT;
tc.event_lon = double(event.max_lon);
tc.event_lat = double(event.max_lat);
tc.event_vmax_ms = double(event.max_inner_vmax_ms);
tc.track_dt_seconds = trackInfo.dt_seconds;
end

%% ========================================================================
function eventIndex = choose_valid_event_index(eventIndex, origIndex, lon, lat, wind, event)

if isfinite(eventIndex) && any(origIndex == eventIndex)
    return;
end

eventLon = double(event.max_lon);
eventLat = double(event.max_lat);
if isfinite(eventLon) && isfinite(eventLat)
    dlon = wrap_to_180_local(lon - eventLon);
    d = hypot(dlon .* cosd(0.5 .* (lat + eventLat)), lat - eventLat);
    [~, k] = min(d);
    eventIndex = origIndex(k);
    return;
end

[~, k] = max(wind);
eventIndex = origIndex(k);
end

%% ========================================================================
function x = safe_index(v, idx)

if isempty(v) || idx < 1 || idx > numel(v)
    x = NaN;
else
    x = v(idx);
end
end

%% ========================================================================
function [fort22Path, maxWindDomain, maxWindTime, activeSteps] = write_fort22_for_case( ...
    fort22Path, tc, Grid, Window, Cv_tc, WindC15Lib, WindPressureFit, P)

if exist(fort22Path, 'file') == 2
    delete(fort22Path);
end

fid = fopen(fort22Path, 'wt');
if fid < 0
    error('Cannot write fort.22: %s', fort22Path);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

timeData = Window.time_data(:);
dtStep = seconds(Window.WTIMINC);
dtHr = Window.WTIMINC / 3600;
maxWindDomain = 0;
maxWindTime = NaT;
activeSteps = 0;

for it = 1:numel(timeData)
    t = timeData(it);
    [hasTc, stateNow, stateNext] = get_tc_state_for_step(tc, t, dtStep, WindPressureFit);
    if hasTc
        [U, V, Pres] = calc_c15_uvp_field_grid( ...
            stateNow.lat, stateNow.lon, stateNext.lat, stateNext.lon, dtHr, ...
            stateNow.wind, stateNow.p, Grid.MET_LAT, Grid.MET_LON, ...
            Cv_tc, P.Pn_hPa, P.Re, WindC15Lib, WindPressureFit, P);
        windMag = hypot(U, V);
        stepMax = max(windMag(:), [], 'omitnan');
        if isfinite(stepMax) && stepMax > 0
            activeSteps = activeSteps + 1;
        end
        if isfinite(stepMax) && stepMax > maxWindDomain
            maxWindDomain = stepMax;
            maxWindTime = t;
        end
        write_fort22_snapshot_nws6(fid, U, V, Pres, P);
    else
        write_background_step_nws6(fid, Grid.NWLAT, Grid.NWLON, P);
    end
end
end

%% ========================================================================
function write_case_meta(metaFile, blockId, event, tc, Grid, Window, sourceMeta, ...
    Cv_global, CvTable, WindPressureFit, fort22Path, maxWindDomain, maxWindTime, activeSteps, P) %#ok<INUSD>

fid = fopen(metaFile, 'wt');
if fid < 0
    error('Cannot write %s', metaFile);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'ADCIRC fort.22 (NWS=6) ERA5 TC metadata\n');
fprintf(fid, '----------------------------------------\n');
fprintf(fid, 'block = %s\n', blockId);
fprintf(fid, 'case_source = P4_count_and_plot_era5_tc_for_adcirc_inner_domains.m\n');
fprintf(fid, 'event_csv = %s\n', P.event_csv);
fprintf(fid, 'track_nc = %s\n', P.track_nc);
fprintf(fid, 'fort22 = %s\n', fort22Path);
fprintf(fid, 'static_file_mode = %s\n', P.static_file_mode);
fprintf(fid, 'source_block_meta = %s\n', sourceMeta.meta_file);
fprintf(fid, 'source_block_meta_mode = %s\n', sourceMeta.source);
fprintf(fid, 'time_alignment_scheme = %s\n', P.time_alignment_scheme);
fprintf(fid, 'alignment = global track LMI aligned once; block max_time_index defines local impact time without shifting track\n');
fprintf(fid, 'simulation_window_rule = impact_time - %.3f days to impact_time + %.3f days\n', ...
    P.pre_impact_days, P.post_impact_days);
fprintf(fid, 'calendar_rule = global LMI uses track year/month plus deterministic day keyed only by track_index; random_day_seed=%d\n', ...
    P.random_day_seed);

fprintf(fid, '\n[Event]\n');
fprintf(fid, 'track_id = %s\n', char(event.track_id));
fprintf(fid, 'track_index = %d\n', round(double(event.track_index)));
fprintf(fid, 'year = %.0f\n', double(event.year));
fprintf(fid, 'selected_inner_vmax_ms = %.6f\n', double(event.max_inner_vmax_ms));
fprintf(fid, 'track_lmi_ms = %.6f\n', double(event.lmi_ms));
fprintf(fid, 'selected_lon = %.8f\n', double(event.max_lon));
fprintf(fid, 'selected_lat = %.8f\n', double(event.max_lat));
fprintf(fid, 'selected_max_time_index = %.0f\n', double(event.max_time_index));

fprintf(fid, '\n[Grid]\n');
fprintf(fid, 'WindowStart = %s\n', datestr(Window.start, 31));
fprintf(fid, 'WindowEnd = %s\n', datestr(Window.end, 31));
fprintf(fid, 'SourceWindowStart = %s\n', datestr(Window.source_start, 31));
fprintf(fid, 'SourceWindowEnd = %s\n', datestr(Window.source_end, 31));
fprintf(fid, 'WTIMINC = %d\n', Window.WTIMINC);
fprintf(fid, 'nSteps = %d\n', numel(Window.time_data));
fprintf(fid, 'TC_calendar_year = %d\n', Window.tc_year);
fprintf(fid, 'TC_calendar_month = %d\n', Window.tc_month);
fprintf(fid, 'TC_global_anchor_day = %d\n', Window.random_day);
fprintf(fid, 'TC_random_day = %d\n', Window.random_day);
fprintf(fid, 'TC_global_anchor_index = %d\n', Window.global_anchor_index);
fprintf(fid, 'TC_global_anchor_time = %s\n', datestr(Window.global_anchor_time, 31));
fprintf(fid, 'TC_block_impact_index = %d\n', Window.block_impact_index);
fprintf(fid, 'TC_block_impact_time = %s\n', datestr(Window.impact_time, 31));
% Retained as a compatibility alias; it is now the block-specific impact.
fprintf(fid, 'TC_aligned_impact_time = %s\n', datestr(Window.impact_time, 31));
fprintf(fid, 'PreImpactDays = %.8f\n', Window.pre_impact_days);
fprintf(fid, 'PostImpactDays = %.8f\n', Window.post_impact_days);
fprintf(fid, 'NWLON = %d\n', Grid.NWLON);
fprintf(fid, 'NWLAT = %d\n', Grid.NWLAT);
fprintf(fid, 'WLONMIN = %.8f\n', Grid.WLONMIN);
fprintf(fid, 'WLATMAX = %.8f\n', Grid.WLATMAX);
fprintf(fid, 'WLONINC = %.8f\n', Grid.DLON);
fprintf(fid, 'WLATINC = %.8f\n', Grid.DLAT);

fprintf(fid, '\n[C15]\n');
fprintf(fid, 'Cv_used = %.8f\n', Cv_global);
fprintf(fid, 'Cv_source_csv = %s\n', P.basin_cv_csv);
fprintf(fid, 'Cv_weight = %s\n', P.cv_weight_field);
fprintf(fid, 'PressureFit_mode = %s\n', P.wind_pressure_fit_mode);
fprintf(fid, 'PressureFit_selection = %s\n', WindPressureFit.selection_method);
fprintf(fid, 'PressureFit_track_basin_code = %s\n', WindPressureFit.track_basin_code);
fprintf(fid, 'PressureFit_basin = %s\n', WindPressureFit.basin);
fprintf(fid, 'PressureFit_source = %s\n', WindPressureFit.source_csv);
fprintf(fid, 'PressureFit_model = %s\n', WindPressureFit.model);
fprintf(fid, 'PressureFit_pref_hPa = %.8f\n', WindPressureFit.prefHpa);
fprintf(fid, 'PressureFit_A = %.12g\n', WindPressureFit.A);
fprintf(fid, 'PressureFit_B = %.12g\n', WindPressureFit.B);
fprintf(fid, 'PressureFit_pressure_clamp_hPa = %.2f, %.2f\n', WindPressureFit.pcMinHpa, WindPressureFit.pcMaxHpa);

fprintf(fid, '\n[Tide]\n');
fprintf(fid, 'fort19_generation = %s\n', bool_to_string(P.generate_fort19_from_tmd));
fprintf(fid, 'TMD_root = %s\n', P.tmd_root);
fprintf(fid, 'TMD_start_time = %s\n', datestr(Window.start, 31));
fprintf(fid, 'TMD_end_time = %s\n', datestr(Window.end, 31));
fprintf(fid, 'TMD_time_step_seconds = %d\n', Window.WTIMINC);
fprintf(fid, 'TMD_nSteps = %d\n', numel(Window.time_data));

fprintf(fid, '\n[AlignedTrack]\n');
fprintf(fid, 'event_aligned_time = %s\n', datestr(tc.event_aligned_time, 31));
fprintf(fid, 'global_anchor_time = %s\n', datestr(tc.global_anchor_time, 31));
fprintf(fid, 'global_anchor_index = %d\n', tc.global_anchor_index);
fprintf(fid, 'global_anchor_vmax_ms = %.6f\n', tc.global_anchor_vmax_ms);
fprintf(fid, 'event_relative_to_global_anchor_seconds = %.6f\n', tc.event_relative_seconds);
fprintf(fid, 'track_year = %.0f\n', tc.year);
fprintf(fid, 'track_month = %.0f\n', tc.month);
fprintf(fid, 'track_basin = %s\n', char(string(tc.basin)));
fprintf(fid, 'pressure_fit_basin = %s\n', tc.pressure_fit_basin);
fprintf(fid, 'track_dt_seconds = %.0f\n', tc.track_dt_seconds);
fprintf(fid, 'track_start = %s\n', datestr(tc.time(1), 31));
fprintf(fid, 'track_end = %s\n', datestr(tc.time(end), 31));
fprintf(fid, 'track_points = %d\n', numel(tc.time));
fprintf(fid, 'track_max_vmax_ms = %.6f\n', max(tc.wind, [], 'omitnan'));
fprintf(fid, 'track_min_pressure_hPa = %.6f\n', min(tc.p, [], 'omitnan'));

fprintf(fid, '\n[ForcingResult]\n');
fprintf(fid, 'active_wind_steps = %d\n', activeSteps);
fprintf(fid, 'domain_max_wind_ms = %.6f\n', maxWindDomain);
fprintf(fid, 'domain_max_wind_time = %s\n', fmt_time_safe(maxWindTime));

fprintf(fid, '\n[ForcingTrackPoints]\n');
fprintf(fid, 'index,original_time_index,time,lon180,lat,vmax_ms,pressure_hPa\n');
for i = 1:numel(tc.time)
    fprintf(fid, '%d,%d,%s,%.8f,%.8f,%.6f,%.6f\n', i, tc.original_time_index(i), ...
        datestr(tc.time(i), 31), tc.lon(i), tc.lat(i), tc.wind(i), tc.p(i));
end
end

%% ========================================================================
function assert_tmd_ready(tmd_root, latlon_file)

if exist(tmd_root, 'dir') ~= 7
    error('TMD root does not exist: %s', tmd_root);
end
latlon_dir = fileparts(latlon_file);
if exist(latlon_dir, 'dir') ~= 7
    error('TMD LAT_LON directory does not exist: %s', latlon_dir);
end
end

%% ========================================================================
function generate_fort19_from_tmd_for_window(G, fort19_file, Window, P)

WTIMINC = Window.WTIMINC;
dt_min = round(WTIMINC / 60);
if abs(WTIMINC / 60 - dt_min) > 1e-8
    warning('WTIMINC is not an integer minute; rounded dt_min to %d.', dt_min);
end

nSteps = numel(Window.time_data);
if nSteps <= 0
    nSteps = round(seconds(Window.end - Window.start) / WTIMINC) + 1;
end
if nSteps <= 0
    error('Invalid time window for fort.19.');
end
if isempty(G.open_nodes_ordered)
    error('Cannot generate fort.19 because fort.14 has no open-boundary nodes.');
end

write_tmd_latlon_file(P.tmd_latlon_file, G.open_lat_ordered, G.open_lon_ordered, ...
    Window.start, dt_min, nSteps);

close_open_file_handles();
if P.clean_tmd_before_each_case
    cleanup_tmd_outputs(P.tmd_root, P.tmd_out_dir);
end
ensure_blank_data_out(fullfile(P.tmd_root, 'data.out'));
run_tmd_ato(P.tmd_root);
close_open_file_handles();

tide = read_tmd_series_matrix(P.tmd_root, numel(G.open_nodes_ordered), nSteps, P.allow_tmd_trim_or_pad);
write_fort19_file(fort19_file, tide, WTIMINC);
end

%% ========================================================================
function close_case_file_handles(P, status)

nClosed = close_open_file_handles();
if nClosed > 0
    fprintf('  -> closed %d lingering file handle(s) after %s\n', nClosed, status);
end

if isfield(P, 'tmd_root') && exist(P.tmd_root, 'dir') == 7
    close_open_file_handles();
end
end

%% ========================================================================
function nClosed = close_open_file_handles()

nClosed = 0;
fids = fopen('all');
for i = 1:numel(fids)
    fid = fids(i);
    if fid <= 2
        continue;
    end
    try
        fclose(fid);
        nClosed = nClosed + 1;
    catch
    end
end
end

%% ========================================================================
function write_tmd_latlon_file(latlon_file, open_lat, open_lon, start_time, dt_min, nSteps)

fid = fopen(latlon_file, 'wt');
if fid < 0
    error('Cannot write TMD input file: %s', latlon_file);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '   lat   lon            yy   mm  dd hh mi sec dt(min) TSLength\n');

[yy, mm, dd] = ymd(start_time);
hh = hour(start_time);
mi = minute(start_time);
ss = round(second(start_time));

for i = 1:numel(open_lat)
    fprintf(fid, '  %.10f  %.10f', open_lat(i), open_lon(i));
    fprintf(fid, '    %04d %02d %02d %02d %02d %02d %d %d\n', ...
        yy, mm, dd, hh, mi, ss, dt_min, nSteps);
end
end

%% ========================================================================
function cleanup_tmd_outputs(tmd_root, tmd_out_dir)

delete_if_exists(fullfile(tmd_root, 'data_*.mat'));
delete_if_exists(fullfile(tmd_root, 'data.out'));
delete_if_exists(fullfile(tmd_root, 'data*'));

if exist(tmd_out_dir, 'dir') == 7
    delete_if_exists(fullfile(tmd_out_dir, 'data.out'));
    delete_if_exists(fullfile(tmd_out_dir, 'data*'));
end
end

%% ========================================================================
function delete_if_exists(pattern)

D = dir(pattern);
for i = 1:numel(D)
    if D(i).isdir
        continue;
    end
    try
        delete(fullfile(D(i).folder, D(i).name));
    catch
    end
end
end

%% ========================================================================
function ensure_blank_data_out(data_out_file)

fid = fopen(data_out_file, 'wt');
if fid < 0
    error('Cannot write %s', data_out_file);
end
fprintf(fid, '%s\n', ' ');
fclose(fid);
end

%% ========================================================================
function run_tmd_ato(tmd_root)

old_dir = pwd;
cleanupObj = onCleanup(@() cd(old_dir)); %#ok<NASGU>
cd(tmd_root);

if exist('TMD_ato', 'file') ~= 2 && exist('TMD_ato', 'builtin') ~= 5
    error('TMD_ato is not on the MATLAB path.');
end

TMD_ato;
end

%% ========================================================================
function tide = read_tmd_series_matrix(tmd_root, nOpen, nSteps, allow_trim_or_pad)

tide = nan(nOpen, nSteps);
for i = 1:nOpen
    f = fullfile(tmd_root, sprintf('data_%d.mat', i));
    if exist(f, 'file') ~= 2
        error('Missing TMD output: %s', f);
    end

    S = load(f);
    if ~isfield(S, 'TimeSeries')
        error('%s does not contain TimeSeries.', f);
    end

    ts = real(S.TimeSeries(:));
    if numel(ts) == nSteps
        tide(i, :) = ts.';
    elseif allow_trim_or_pad
        tmp = nan(nSteps, 1);
        if numel(ts) > nSteps
            tmp(:) = ts(1:nSteps);
        elseif isempty(ts)
            tmp(:) = 0;
        else
            tmp(1:numel(ts)) = ts;
            tmp(numel(ts)+1:end) = ts(end);
        end
        tide(i, :) = tmp.';
    else
        error('TMD length mismatch for data_%d.mat: got %d, need %d.', i, numel(ts), nSteps);
    end
end
end

%% ========================================================================
function write_fort19_file(fort19_file, tide, ETIMINC)

[nOpen, nSteps] = size(tide);
fid = fopen(fort19_file, 'wt');
if fid < 0
    error('Cannot write fort.19: %s', fort19_file);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '%d\n', ETIMINC);
for k = 1:nSteps
    for j = 1:nOpen
        fprintf(fid, '%.10f\n', tide(j, k));
    end
end
end

%% ========================================================================
function [Cv_global, T] = load_global_weighted_mean_Cv(csv_file, mat_file, weight_field)

if exist(csv_file, 'file') == 2
    T = readtable(csv_file, 'TextType', 'string');
elseif exist(mat_file, 'file') == 2
    S = load(mat_file);
    if isfield(S, 'T_basinCv')
        T = S.T_basinCv;
    elseif isfield(S, 'basinCvStats')
        T = struct2table(S.basinCvStats);
    else
        error('Cannot find T_basinCv or basinCvStats in %s', mat_file);
    end
else
    error('Cannot find Cv summary: %s or %s', csv_file, mat_file);
end

if ~ismember('mean_Cv', T.Properties.VariableNames)
    error('Cv table must contain mean_Cv.');
end

if height(T) == 1 && ismember('n_accepted_tc', T.Properties.VariableNames) && ...
        ismember('n_Cv_samples', T.Properties.VariableNames)
    Cv_global = T.mean_Cv(1);
    if ~isfinite(Cv_global)
        error('Accepted Cv summary contains a non-finite mean_Cv.');
    end
    return
end

if ~ismember(weight_field, T.Properties.VariableNames)
    w = ones(height(T), 1);
else
    w = T.(weight_field);
end

cv = T.mean_Cv;
good = isfinite(cv) & isfinite(w) & w > 0;
if ~any(good)
    error('No finite Cv values with positive weights.');
end

Cv_global = sum(cv(good) .* w(good)) ./ sum(w(good));
end

%% ========================================================================
function Fits = load_wind_pressure_fits(csv_file, pc_min_hpa, pc_max_hpa)

T = readtable(csv_file, 'TextType', 'string');
need_vars = {'basin', 'prefHpa', 'A', 'B'};
missing = setdiff(need_vars, T.Properties.VariableNames);
if ~isempty(missing)
    error('Wind-pressure fit CSV missing columns: %s', strjoin(missing, ', '));
end

Fits = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 1:height(T)
    basinName = string(T.basin(i));
    if strlength(basinName) == 0
        continue;
    end

    fit = struct();
    fit.source_csv = csv_file;
    fit.basin = basinName;
    fit.prefHpa = T.prefHpa(i);
    fit.A = T.A(i);
    fit.B = T.B(i);
    fit.pcMinHpa = pc_min_hpa;
    fit.pcMaxHpa = pc_max_hpa;
    if ismember('tailQuantile', T.Properties.VariableNames)
        fit.tailQuantile = T.tailQuantile(i);
    else
        fit.tailQuantile = NaN;
    end
    if ismember('pressureBinWidthHpa', T.Properties.VariableNames)
        fit.pressureBinWidthHpa = T.pressureBinWidthHpa(i);
    else
        fit.pressureBinWidthHpa = NaN;
    end
    if ismember('model', T.Properties.VariableNames)
        fit.model = char(string(T.model(i)));
    else
        fit.model = sprintf('V_ms = %.6g * (%.6g - Pc_hPa) ^ %.6g', fit.A, fit.prefHpa, fit.B);
    end
    Fits(char(basinName)) = fit;
end
if ~isKey(Fits, 'GLOBAL')
    error('Cannot find GLOBAL row in wind-pressure fit CSV: %s', csv_file);
end
end

%% ========================================================================
function fit = select_wind_pressure_fit_for_track(event, trackInfo, Fits, P)

trackIndex = round(double(event.track_index));
trackBasin = "";
if trackIndex >= 1 && trackIndex <= numel(trackInfo.basins)
    trackBasin = upper(strtrim(string(trackInfo.basins(trackIndex))));
end
basinId = pressure_basin_id_from_track_code(trackBasin);

if strcmpi(P.wind_pressure_fit_mode, 'global')
    fit = Fits('GLOBAL');
    selectionMethod = "forced_global";
elseif strlength(basinId) > 0 && isKey(Fits, char(basinId))
    fit = Fits(char(basinId));
    selectionMethod = "track_nc.tc_basins";
else
    warning('Track %d has unsupported TC basin code "%s"; using GLOBAL pressure fit.', ...
        trackIndex, trackBasin);
    fit = Fits('GLOBAL');
    selectionMethod = "track_basin_fallback_global";
end
fit.track_basin_code = char(trackBasin);
fit.selection_method = char(selectionMethod);
end

%% ========================================================================
function basinId = pressure_basin_id_from_track_code(trackBasin)

switch upper(strtrim(string(trackBasin)))
    case "NA"
        basinId = "BASIN_NATL";
    case "EP"
        basinId = "BASIN_ENP";
    case "WP"
        basinId = "BASIN_WNP";
    case "NI"
        basinId = "BASIN_NIO";
    case "SI"
        basinId = "BASIN_SIO";
    case {"SP", "AU"}
        basinId = "BASIN_AUSSP";
    otherwise
        basinId = "";
end
end

%% ========================================================================
function p_hpa = wind_to_pressure_hpa(wind_ms, fit)

p_hpa = nan(size(wind_ms));
wind_ms = double(wind_ms);
good = isfinite(wind_ms) & wind_ms >= 0;
if ~any(good)
    return;
end
delta_p = (wind_ms(good) ./ fit.A) .^ (1 ./ fit.B);
p = fit.prefHpa - delta_p;
p = max(fit.pcMinHpa, min(fit.pcMaxHpa, p));
p_hpa(good) = p;
end

%% ========================================================================
function WindC15Lib = load_wind_c15_library(predata_dir, Vmax_list, Rmax_list)

WindC15Lib = containers.Map;
for Vmax = Vmax_list
    for Rmaxkm = Rmax_list
        fname = fullfile(predata_dir, ...
            ['Wind_C15_data_Vmax', num2str(Vmax), '_Rmax', num2str(Rmaxkm), '.mat']);
        if exist(fname, 'file') == 2
            S = load(fname);
            key = sprintf('%d_%d', Vmax, Rmaxkm);
            WindC15Lib(key) = S.Wind_C15_data;
        end
    end
end
if WindC15Lib.Count == 0
    error('No C15 lookup tables found under %s', predata_dir);
end
end

%% ========================================================================
function [has_tc, state_now, state_next] = get_tc_state_for_step(tc, t, dt_step, WindPressureFit)

if isempty(tc.time) || t < tc.time(1) || t > tc.time(end)
    has_tc = false;
    state_now = empty_tc_state();
    state_next = empty_tc_state();
    return;
end

state_now = interp_tc_state(tc, t, WindPressureFit);
t2 = t + dt_step;
if t2 > tc.time(end)
    state_next = state_now;
else
    state_next = interp_tc_state(tc, t2, WindPressureFit);
    if ~isfinite(state_next.lat) || ~isfinite(state_next.lon) || ~isfinite(state_next.wind)
        state_next = state_now;
    end
end

if ~isfinite(state_now.lat) || ~isfinite(state_now.lon) || ~isfinite(state_now.wind)
    has_tc = false;
else
    has_tc = true;
end
end

%% ========================================================================
function s = interp_tc_state(tc, t_query, WindPressureFit)

t_src = datenum(tc.time(:));
tq = datenum(t_query);

[t_src_u, IA] = unique(t_src, 'stable');
lat_u = tc.lat(IA);
if isfield(tc, 'lon_unwrapped') && numel(tc.lon_unwrapped) == numel(tc.lon)
    lon_u = tc.lon_unwrapped(IA);
else
    lon_u = rad2deg(unwrap(deg2rad(tc.lon(IA))));
end
wind_u = tc.wind(IA);
p_u = tc.p(IA);

lat = interp1(t_src_u, lat_u, tq, 'linear', NaN);
lon = wrap_to_180_local(interp1(t_src_u, lon_u, tq, 'linear', NaN));
wind = interp1(t_src_u, wind_u, tq, 'linear', NaN);
p = interp1(t_src_u, p_u, tq, 'linear', NaN);

if ~isfinite(p) || p <= 0
    p = wind_to_pressure_hpa(wind, WindPressureFit);
end
s = struct('lat', lat, 'lon', lon, 'wind', wind, 'p', p);
end

%% ========================================================================
function s = empty_tc_state()

s = struct('lat', NaN, 'lon', NaN, 'wind', NaN, 'p', NaN);
end

%% ========================================================================
function [U, V, Pres, Rmax_m] = calc_c15_uvp_field_grid( ...
    lat_c, lon_c, lat_n, lon_n, dt_hr, v_c, p_c, ...
    LAT, LON, Cv, Pn_hPa, Re, WindC15Lib, WindPressureFit, Params) %#ok<INUSD>

U = zeros(size(LAT));
V = zeros(size(LAT));
Pres = ones(size(LAT)) * Params.standard_pressure_pa;
Rmax_m = NaN;

if ~isfinite(lat_c) || ~isfinite(lon_c) || ~isfinite(v_c)
    return;
end
if ~isfinite(p_c) || p_c <= 0
    p_c = wind_to_pressure_hpa(v_c, WindPressureFit);
end
if ~isfinite(p_c) || p_c >= Pn_hPa
    return;
end

if dt_hr <= 0 || ~isfinite(dt_hr)
    vmc = 0;
    fai = 0;
else
    [dis, alpha] = fast_track_step(lat_c, lon_c, lat_n, lon_n);
    vmc = dis / (dt_hr * 3600);
    fai = alpha;
end

vm_raw = v_c - vmc;
min_lookup_vm = Params.c15_min_lookup_vm_ms;
weak_scale = min(max(vm_raw, 0) ./ min_lookup_vm, 1);
if weak_scale <= 0
    return;
end
vm = max(vm_raw, min_lookup_vm);
vm = min(vm, Params.c15_max_lookup_vm_ms);
p_model_hPa = Pn_hPa - weak_scale .* (Pn_hPa - p_c);

Rmax_m = Cv * 51.6 * exp(-0.0223 * vm + 0.0281 * abs(lat_c)) * 1000;
Rmax_m = max(Rmax_m, Params.c15_min_rmax_m);

B = (v_c^2) * 1.15 * exp(1) / max(Pn_hPa - p_model_hPa, 1) / 100;
B = max(Params.c15_holland_b_min, min(B, Params.c15_holland_b_max));

dlon_grid = wrap_to_180_local(LON - lon_c);
dx = dlon_grid .* cosd(0.5 * (LAT + lat_c)) * 111320;
dy = (LAT - lat_c) * 110540;
r = hypot(dx, dy);
r = max(r, 1.0);

cta = atan2d(dy, dx);
cta(cta < 0) = cta(cta < 0) + 360;

Pg_TC = (p_model_hPa + (Pn_hPa - p_model_hPa) .* exp(-(Rmax_m ./ r).^B)) * 100;

Vmax = round(vm);
Rmaxkm = round(Rmax_m / 1000);
key = sprintf('%d_%d', Vmax, Rmaxkm);

if isKey(WindC15Lib, key)
    C15 = WindC15Lib(key);
    rr = double(C15.rr(:));
    % C15.vg is the lookup file's axisymmetric 1-min surface tangential
    % wind speed, not a gradient, upper-tropospheric or environmental wind.
    c15_surface_wind_1min_table_ms = double(C15.vg(:));
    c15_surface_wind_1min_ms = interp1( ...
        rr, c15_surface_wind_1min_table_ms, r, 'linear', 0);
else
    c15_surface_wind_1min_ms = zeros(size(r));
end

c15_surface_wind_1min_ms(~isfinite(c15_surface_wind_1min_ms) | ...
    c15_surface_wind_1min_ms < 0) = 0;
c15_surface_wind_1min_ms = weak_scale .* c15_surface_wind_1min_ms;

beta = zeros(size(r));
I1 = r < Rmax_m;
I2 = r >= Rmax_m & r < Params.c15_beta_mid_radius_factor * Rmax_m;
I3 = r >= Params.c15_beta_mid_radius_factor * Rmax_m;
beta(I1) = Params.c15_beta_inner_base_deg + Params.c15_beta_inner_slope_deg .* (r(I1) ./ Rmax_m);
beta(I2) = Params.c15_beta_mid_base_deg + Params.c15_beta_mid_slope_deg .* (r(I2) ./ Rmax_m - 1);
beta(I3) = Params.c15_beta_outer_deg;

vmoc = weak_scale .* vmc .* r .* Rmax_m ./ (r.^2 + Rmax_m^2);
hemisphere_sign = sign(lat_c);
if hemisphere_sign == 0
    hemisphere_sign = 1;
end
rotation_angle = hemisphere_sign .* (90 + beta);
Vx_TC = c15_surface_wind_1min_ms .* cosd(cta + rotation_angle) + ...
    vmoc .* cosd(fai);
Vy_TC = c15_surface_wind_1min_ms .* sind(cta + rotation_angle) + ...
    vmoc .* sind(fai);
% Convert the complete 1-min surface wind vector to the 10-min wind used by
% ADCIRC fort.22. This is an averaging-period conversion, not a
% gradient-to-surface wind reduction.
Vx_TC = Params.c15_one_min_to_ten_min_factor .* Vx_TC;
Vy_TC = Params.c15_one_min_to_ten_min_factor .* Vy_TC;

R1 = Params.c15_blend_inner_radius_m;
R2 = Params.c15_blend_outer_radius_m;
lamda = zeros(size(r));
I_mid = r >= R1 & r <= R2;
I_far = r > R2;
lamda(I_mid) = (r(I_mid) - R1) ./ (R2 - R1);
lamda(I_far) = 1.0;

U = (1 - lamda) .* Vx_TC;
V = (1 - lamda) .* Vy_TC;
Pres = (1 - lamda) .* Pg_TC + lamda .* Params.standard_pressure_pa;
end

%% ========================================================================
function [dis_m, alpha_deg] = fast_track_step(lat1, lon1, lat2, lon2)

dlon = wrap_to_180_local(lon2 - lon1);
dx = dlon * cosd(0.5 * (lat1 + lat2)) * 111320;
dy = (lat2 - lat1) * 110540;
dis_m = hypot(dx, dy);
alpha_deg = atan2d(dy, dx);
if alpha_deg < 0
    alpha_deg = alpha_deg + 360;
end
end

%% ========================================================================
function write_background_step_nws6(fid, NWLAT, NWLON, P)

for k = 1:NWLAT
    for j = 1:NWLON
        fprintf(fid, '%s\n', P.background_fort22_line);
    end
end
end

%% ========================================================================
function write_fort22_snapshot_nws6(fid, U, V, Pres, P)

[nlat, nlon] = size(U);
tol_uv = P.fort22_zero_uv_tolerance;
tol_p = P.fort22_zero_pressure_tolerance_pa;
for k = 1:nlat
    for j = 1:nlon
        u = U(k, j);
        v = V(k, j);
        p = Pres(k, j);
        if abs(u) < tol_uv && abs(v) < tol_uv && abs(p - P.standard_pressure_pa) < tol_p
            fprintf(fid, '%s\n', P.background_fort22_line);
        else
            fprintf(fid, '%.1f %.1f %.0f\n', u, v, p);
        end
    end
end
end

%% ========================================================================
function G = read_fort14_nodes_and_open_boundary(fort14_file)

fid = fopen(fort14_file, 'rt');
if fid < 0
    error('Cannot open fort.14: %s', fort14_file);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

G = struct();
G.title = strtrim(fgetl(fid));

line2 = strtrim(fgetl(fid));
a = sscanf(line2, '%d %d');
if numel(a) < 2
    error('Cannot parse fort.14 line 2.');
end
G.NE = a(1);
G.NP = a(2);

node_id = zeros(G.NP, 1);
lon = zeros(G.NP, 1);
lat = zeros(G.NP, 1);
depth = zeros(G.NP, 1);

for i = 1:G.NP
    s = fgetl(fid);
    v = sscanf(s, '%f');
    if numel(v) < 4
        error('Cannot parse node line %d in fort.14.', i);
    end
    node_id(i) = round(v(1));
    lon(i) = v(2);
    lat(i) = v(3);
    depth(i) = v(4);
end

max_node_id = max(node_id);
id_to_idx = zeros(max_node_id, 1);
id_to_idx(node_id) = (1:G.NP).';

elem_id = zeros(G.NE, 1);
elem_node_idx = nan(G.NE, 3);
for i = 1:G.NE
    s = fgetl(fid);
    v = sscanf(s, '%f');
    if numel(v) < 5
        error('Cannot parse element line %d in fort.14.', i);
    end
    elem_id(i) = round(v(1));
    nverts = round(v(2));
    if nverts == 3
        ids = round(v(3:5));
        if all(ids >= 1) && all(ids <= max_node_id) && all(id_to_idx(ids) > 0)
            elem_node_idx(i, :) = id_to_idx(ids);
        end
    end
end

s = strtrim(fgetl(fid));
G.NOPE = sscanf(s, '%d', 1);
if isempty(G.NOPE)
    G.NOPE = 0;
end

s = strtrim(fgetl(fid));
G.NETA = sscanf(s, '%d', 1);
if isempty(G.NETA)
    G.NETA = 0;
end

open_node_ids = [];
for ib = 1:G.NOPE
    header = strtrim(fgetl(fid));
    hv = sscanf(header, '%d');
    if isempty(hv)
        error('Cannot parse open-boundary header %d.', ib);
    end
    nvdll = hv(1);

    tmp = zeros(nvdll, 1);
    for k = 1:nvdll
        line = strtrim(fgetl(fid));
        vv = sscanf(line, '%d');
        if isempty(vv)
            error('Cannot parse open-boundary node %d in boundary %d.', k, ib);
        end
        tmp(k) = vv(1);
    end
    open_node_ids = [open_node_ids; tmp]; %#ok<AGROW>
end

if any(open_node_ids < 1 | open_node_ids > max_node_id | id_to_idx(open_node_ids) <= 0)
    error('Open-boundary node IDs in fort.14 are outside the node table.');
end
open_nodes = id_to_idx(open_node_ids);

G.node_id = node_id;
G.lon = lon;
G.lat = lat;
G.depth = depth;
G.elem_id = elem_id;
G.elem_node_idx = elem_node_idx;
G.open_node_ids_ordered = open_node_ids(:);
G.open_nodes_ordered = open_nodes(:);
G.open_lon_ordered = lon(open_nodes);
G.open_lat_ordered = lat(open_nodes);
end

%% ========================================================================
function Mesh = read_fort14_nodes_only(fort14_file)

fid = fopen(fort14_file, 'rt');
if fid < 0
    error('Cannot open fort.14: %s', fort14_file);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fgetl(fid);
tmp = sscanf(strtrim(fgetl(fid)), '%f');
NE = tmp(1);
NP = tmp(2);

node_id = zeros(NP, 1);
lon = zeros(NP, 1);
lat = zeros(NP, 1);
dep = zeros(NP, 1);

for i = 1:NP
    a = sscanf(fgetl(fid), '%f');
    node_id(i) = a(1);
    lon(i) = a(2);
    lat(i) = a(3);
    dep(i) = a(4);
end

Mesh = struct('NE', NE, 'NP', NP, 'node_id', node_id, 'lon', lon, 'lat', lat, 'dep', dep);
end

%% ========================================================================
function Grid = build_met_grid_from_mesh(Mesh, margin_deg, dlon, dlat)

Grid = struct();
Grid.lon_min = min(Mesh.lon);
Grid.lon_max = max(Mesh.lon);
Grid.lat_min = min(Mesh.lat);
Grid.lat_max = max(Mesh.lat);

Grid.WLONMIN = Grid.lon_min - margin_deg;
Grid.WLONMAX = Grid.lon_max + margin_deg;
Grid.WLATMIN = Grid.lat_min - margin_deg;
Grid.WLATMAX = Grid.lat_max + margin_deg;
Grid.DLON = dlon;
Grid.DLAT = dlat;

Grid.met_lon = Grid.WLONMIN : dlon : Grid.WLONMAX;
Grid.met_lat = Grid.WLATMAX : -dlat : Grid.WLATMIN;

Grid.NWLON = numel(Grid.met_lon);
Grid.NWLAT = numel(Grid.met_lat);

[Grid.MET_LON, Grid.MET_LAT] = meshgrid(Grid.met_lon, Grid.met_lat);
end

%% ========================================================================
function x = clean_fill_to_nan(x)

x = double(x);
x(abs(x) > 1e20) = NaN;
x(x <= -9990) = NaN;
end

%% ========================================================================
function lon = wrap_to_180_local(lon)

lon = mod(double(lon) + 180, 360) - 180;
end

%% ========================================================================
function name = make_case_name(trackIndex, tcId, P)

name = sprintf('%s_%06d_%s', P.case_name_prefix, round(trackIndex), sanitize_name(tcId));
if numel(name) > P.case_name_max_chars
    name = name(1:P.case_name_max_chars);
end
end

%% ========================================================================
function s = sanitize_name(s)

s = char(s);
s = strtrim(s);
if isempty(s)
    s = 'NONAME';
end
s = regexprep(s, '[^\w\-]+', '_');
end

%% ========================================================================
function row = make_empty_summary_row(event, caseDir)

row = struct();
row.block_id = string(event.block_id);
row.track_id = string(event.track_id);
row.track_index = double(event.track_index);
row.year = double(event.year);
row.max_inner_vmax_ms = double(event.max_inner_vmax_ms);
row.lmi_ms = double(event.lmi_ms);
row.max_lon = double(event.max_lon);
row.max_lat = double(event.max_lat);
row.max_time_index = double(event.max_time_index);
row.case_dir = string(caseDir);
row.fort22 = "";
row.window_start = "";
row.window_end = "";
row.nsteps = NaN;
row.nwlon = NaN;
row.nwlat = NaN;
row.active_steps = NaN;
row.domain_max_wind_ms = NaN;
row.domain_max_wind_time = "";
row.status = "pending";
row.message = "";
end

%% ========================================================================
function C = summary_struct_to_cell(row)

C = {row.block_id, row.track_id, row.track_index, row.year, ...
    row.max_inner_vmax_ms, row.lmi_ms, row.max_lon, row.max_lat, row.max_time_index, ...
    row.case_dir, row.fort22, row.window_start, row.window_end, row.nsteps, ...
    row.nwlon, row.nwlat, row.active_steps, row.domain_max_wind_ms, ...
    row.domain_max_wind_time, row.status, row.message};
end

%% ========================================================================
function T = summary_rows_to_table(C)

varNames = {'block_id', 'track_id', 'track_index', 'year', ...
    'max_inner_vmax_ms', 'lmi_ms', 'max_lon', 'max_lat', 'max_time_index', ...
    'case_dir', 'fort22', 'window_start', 'window_end', 'nsteps', ...
    'nwlon', 'nwlat', 'active_steps', 'domain_max_wind_ms', ...
    'domain_max_wind_time', 'status', 'message'};

if isempty(C)
    T = cell2table(cell(0, numel(varNames)), 'VariableNames', varNames);
else
    T = cell2table(C, 'VariableNames', varNames);
end
end

%% ========================================================================
function idx = find_line_contains(lines, token)

idx = [];
for i = 1:numel(lines)
    if contains(lines{i}, token)
        idx = i;
        return;
    end
end
end

%% ========================================================================
function idx = find_line_any(lines, tokens)

idx = [];
for i = 1:numel(lines)
    for j = 1:numel(tokens)
        if contains(lines{i}, tokens{j})
            idx = i;
            return;
        end
    end
end
end

%% ========================================================================
function new_line = update_output_schedule_line(line, rnday)

excl = strfind(line, '!');
if isempty(excl)
    data_part = strtrim(line);
    comment_part = '';
else
    data_part = strtrim(line(1:excl(1)-1));
    comment_part = line(excl(1):end);
end

nums = sscanf(data_part, '%f');
if numel(nums) < 4
    new_line = line;
    return;
end

NOUT = nums(1);
TOUTS = nums(2);
NSPOOL = nums(4);
if TOUTS > rnday
    TOUTS = max(0, rnday - 1/24);
end

new_line = sprintf(' %g %.2f %.2f %g %s', NOUT, TOUTS, rnday, NSPOOL, comment_part);
end

%% ========================================================================
function new_line = update_fort63_output_schedule_line(line, Window, rnday, P)

excl = strfind(line, '!');
if isempty(excl)
    data_part = strtrim(line);
    comment_part = '';
else
    data_part = strtrim(line(1:excl(1)-1));
    comment_part = line(excl(1):end);
end

nums = sscanf(data_part, '%f');
if numel(nums) < 4
    new_line = line;
    return;
end

NOUTGE = nums(1);
if isfield(P, 'force_fort15_fort63_output') && P.force_fort15_fort63_output
    NOUTGE = max(1, NOUTGE);
end

TOUTSGE = fort63_output_start_day(Window, rnday, P);
if ~isfinite(TOUTSGE)
    TOUTSGE = nums(2);
end

TOUTFGE = fort63_output_end_day(Window, rnday, P);
if ~isfinite(TOUTFGE)
    TOUTFGE = nums(3);
end
TOUTFGE = max(0, min(TOUTFGE, rnday));
TOUTSGE = max(0, min(TOUTSGE, max(0, TOUTFGE - 1/24)));

NSPOOLGE = nums(4);
new_line = sprintf(' %g %.2f %.2f %g %s', NOUTGE, TOUTSGE, TOUTFGE, NSPOOLGE, comment_part);
end

%% ========================================================================
function startDay = fort63_output_start_day(Window, rnday, P) %#ok<INUSD>

startDay = NaN;
mode = "impact";
if isfield(P, 'fort15_fort63_output_start_mode')
    mode = lower(string(P.fort15_fort63_output_start_mode));
end

if mode == "fixed" && ...
        isfield(P, 'fort15_fort63_output_start_day') && ...
        isfinite(P.fort15_fort63_output_start_day)
    startDay = double(P.fort15_fort63_output_start_day);
elseif mode == "impact" && isfield(Window, 'impact_time') && isdatetime(Window.impact_time) && ~isnat(Window.impact_time)
    startDay = days(Window.impact_time - Window.start);
elseif mode == "impact" && isfield(Window, 'pre_impact_days')
    startDay = Window.pre_impact_days;
elseif mode == "impact" && isfield(P, 'pre_impact_days')
    startDay = P.pre_impact_days;
end

if isfinite(startDay)
    startDay = max(0, min(startDay, rnday));
end
end

%% ========================================================================
function endDay = fort63_output_end_day(Window, rnday, P) %#ok<INUSD>

endDay = NaN;
if isfield(P, 'fort15_fort63_output_end_day') && isfinite(P.fort15_fort63_output_end_day)
    endDay = double(P.fort15_fort63_output_end_day);
end

if ~isfinite(endDay) && isfinite(rnday)
    endDay = rnday;
end

if isfinite(endDay)
    endDay = max(0, min(endDay, rnday));
end
end

%% ========================================================================
function s = bool_to_string(tf)

if tf
    s = 'true';
else
    s = 'false';
end
end

%% ========================================================================
function s = fmt_time_safe(t)

if isempty(t) || (isdatetime(t) && isnat(t))
    s = 'NaT';
else
    s = datestr(t, 31);
end
end

%% ========================================================================
function ensure_dir(d)

if exist(d, 'dir') ~= 7
    mkdir(d);
end
end

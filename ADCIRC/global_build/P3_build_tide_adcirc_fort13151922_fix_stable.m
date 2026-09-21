%% P3_build_tide_adcirc_fort13151922_fix_stable.m
% Build a stable ADCIRC static mesh library from adcirc_fort14_meshes.
%
% For each mesh block directory:
%   1) copy fort.14 into adcirc_fort14_meshes_steable
%   2) generate strengthened open-boundary fort.13
%   3) select tide-gauge stations inside inner_ring.csv
%   4) choose an off-season one-month window covered by the selected gauge data
%   5) write a tide-only fort22_meta.txt and selected_tide_stations.csv
%   6) generate fort.15 from Examples/fort.15 with elevation output at gauges
%   7) generate zero-wind/standard-pressure fort.22 (NWS=6 background)
%   8) generate matching-window fort.19 open-boundary tide from TMD

clearvars;
close all;
clc;

%% ------------------------- user controls -------------------------
SCRIPT_DIR = fileparts(mfilename('fullpath'));
if isempty(SCRIPT_DIR)
    SCRIPT_DIR = pwd;
end

PROJECT_ROOT = fileparts(fileparts(SCRIPT_DIR));
SOURCE_MESH_ROOT = env_path('ADCIRC_STATIC_SOURCE_MESH_ROOT', ...
    fullfile(SCRIPT_DIR, 'output', 'adcirc_fort14_meshes'));
OUTPUT_MESH_ROOT = env_path('ADCIRC_STATIC_OUTPUT_MESH_ROOT', ...
    fullfile(SCRIPT_DIR, 'output', 'adcirc_fort14_meshes_stable'));
MESH_ROOT = SOURCE_MESH_ROOT;
FORT15_TEMPLATE = env_path('ADCIRC_FORT15_TEMPLATE', ...
    fullfile(SCRIPT_DIR, 'external', 'fort.15'));
SUB_INTEL_TEMPLATE = env_path('ADCIRC_SUBMIT_TEMPLATE', ...
    fullfile(SCRIPT_DIR, 'sub_intel.sh'));
if exist(SUB_INTEL_TEMPLATE, 'file') ~= 2
    SUB_INTEL_TEMPLATE = fullfile(PROJECT_ROOT, 'Examples', 'sub_intel.sh');
end
TIDELEVEL_ROOT = env_path('ADCIRC_TIDELEVEL_ROOT', ...
    fullfile(SCRIPT_DIR, 'external', 'Tidelevel'));

TMD_ROOT = env_path('ADCIRC_TC_TMD_ROOT', ...
    fullfile(SCRIPT_DIR, 'external', 'TMD2.5', 'TMD'));
TMD_LATLON_FILE = fullfile(TMD_ROOT, 'LAT_LON', 'lat_lon');
TMD_OUT_DIR = fullfile(TMD_ROOT, 'OUT');

MAKE_FORT22 = true;
MAKE_FORT13 = true;
MAKE_FORT15 = true;
MAKE_FORT19 = true;
COPY_SUB_INTEL = true;

% Stable-library controls.
% Default behavior strengthens open-boundary friction for every block and
% writes a complete server-ready input set into OUTPUT_MESH_ROOT.
STRONG_OPEN_BOUNDARY_FORT13_ONLY = false;
USE_STRONG_OPEN_BOUNDARY_FOR_ALL_BLOCKS = true;
COPY_EXISTING_FORT13_BLOCKS = {'ADC_AUSSP_02','ADC_WNP_05'};

OVERWRITE_EXISTING = true;
MAX_MESH_BLOCKS = Inf;   % set to a small number for a pilot run
DRY_RUN = false;          % true: select stations/windows only; write no fort files

FORCE_NWS6 = true;
SYNC_OUTPUT_END_TO_RNDAY = true;
AUTO_UPDATE_DRAMP = true;
ZERO_VELOCITY_MET_STATIONS = true;

CLEAN_TMD_BEFORE_EACH_CASE = true;
ALLOW_TMD_TRIM_OR_PAD = true;

PREFERRED_VALIDATION_YEAR = 2021;
USE_UHSLC_STATIONS = true;
USE_NOAA_STATIONS = true;
MAX_VALIDATION_STATIONS_PER_BLOCK = 10;
ALLOW_NO_STATION_CONVERGENCE_RUN = true;
STATION_SNAP_MAX_DISTANCE_M = 5000.0;

DT_HOUR = 1;
WTIMINC = DT_HOUR * 3600;
MET_GRID_MARGIN_DEG = 1.0;
DLON = 1.0;
DLAT = 1.0;
NO_TC_FORT22_LINE = '0 0 101300';

SUMMARY_CSV = fullfile(OUTPUT_MESH_ROOT, 'tide_only_condition_build_summary.csv');

MAX_MESH_BLOCKS = env_double('ADCIRC_TIDE_ONLY_MAX_BLOCKS', MAX_MESH_BLOCKS);
DRY_RUN = env_bool('ADCIRC_TIDE_ONLY_DRY_RUN', DRY_RUN);
STRONG_OPEN_BOUNDARY_FORT13_ONLY = env_bool( ...
    'ADCIRC_STRONG_OPEN_BOUNDARY_FORT13_ONLY', STRONG_OPEN_BOUNDARY_FORT13_ONLY);
ALLOW_NO_STATION_CONVERGENCE_RUN = env_bool( ...
    'ADCIRC_ALLOW_NO_STATION_CONVERGENCE_RUN', ALLOW_NO_STATION_CONVERGENCE_RUN);
STATION_SNAP_MAX_DISTANCE_M = env_double( ...
    'ADCIRC_STATION_SNAP_MAX_DISTANCE_M', STATION_SNAP_MAX_DISTANCE_M);

if STRONG_OPEN_BOUNDARY_FORT13_ONLY
    MAKE_FORT13 = true;
    MAKE_FORT15 = false;
    MAKE_FORT19 = false;
    MAKE_FORT22 = false;
    COPY_SUB_INTEL = false;
    OVERWRITE_EXISTING = true;
end

NEED_TIDE_SETUP = MAKE_FORT15 || MAKE_FORT19 || MAKE_FORT22;

Fort13 = default_fort13_options();
Fort13.use_strong_open_boundary = USE_STRONG_OPEN_BOUNDARY_FOR_ALL_BLOCKS;
Fort13 = apply_fort13_env_overrides(Fort13);

%% ------------------------- checks -------------------------
if ~exist(MESH_ROOT, 'dir')
    error('Missing mesh root: %s', MESH_ROOT);
end
if ~DRY_RUN
    ensure_dir(OUTPUT_MESH_ROOT);
end
if MAKE_FORT15 && ~exist(FORT15_TEMPLATE, 'file')
    error('Missing fort.15 template: %s', FORT15_TEMPLATE);
end
if COPY_SUB_INTEL && ~exist(SUB_INTEL_TEMPLATE, 'file')
    error('Missing sub_intel.sh template: %s', SUB_INTEL_TEMPLATE);
end
if MAKE_FORT19
    assert_tmd_ready(TMD_ROOT, TMD_LATLON_FILE);
end
if NEED_TIDE_SETUP && ~exist(TIDELEVEL_ROOT, 'dir')
    error('Missing tide-level validation root: %s', TIDELEVEL_ROOT);
end

if MAKE_FORT15
    fort15_template_lines = read_text_lines(FORT15_TEMPLATE);
else
    fort15_template_lines = {};
end

if NEED_TIDE_SETUP
    TideStations = load_tide_station_catalog(TIDELEVEL_ROOT, ...
        USE_UHSLC_STATIONS, USE_NOAA_STATIONS);
    if isempty(TideStations)
        error('No tide-level validation stations were loaded from: %s', TIDELEVEL_ROOT);
    end
else
    TideStations = table();
end

MeshDirs = list_mesh_block_dirs(MESH_ROOT);
if isempty(MeshDirs)
    error('No mesh block directories with fort.14 were found under: %s', MESH_ROOT);
end
ONLY_BLOCK_NAME = strtrim(string(getenv('ADCIRC_TIDE_ONLY_BLOCK_NAME')));
if strlength(ONLY_BLOCK_NAME) > 0
    keep_block = strcmpi(string({MeshDirs.name}), ONLY_BLOCK_NAME);
    MeshDirs = MeshDirs(keep_block);
    if isempty(MeshDirs)
        error('No mesh block named %s was found under: %s', ONLY_BLOCK_NAME, MESH_ROOT);
    end
end
MeshDirs = MeshDirs(1:min(numel(MeshDirs), MAX_MESH_BLOCKS));

fprintf('Tide-only ADCIRC input build\n');
fprintf('  source mesh root: %s\n', MESH_ROOT);
fprintf('  output mesh root: %s\n', OUTPUT_MESH_ROOT);
fprintf('  mesh blocks: %d\n', numel(MeshDirs));
fprintf('  strong open-boundary fort.13-only mode: %d\n', STRONG_OPEN_BOUNDARY_FORT13_ONLY);
fprintf('  strong open-boundary for all blocks: %d\n', USE_STRONG_OPEN_BOUNDARY_FOR_ALL_BLOCKS);
fprintf('  copy existing fort.13 blocks: %s\n', strjoin(string(COPY_EXISTING_FORT13_BLOCKS), ', '));
fprintf('  fort.15 template: %s\n', FORT15_TEMPLATE);
fprintf('  sub_intel.sh template: %s\n', SUB_INTEL_TEMPLATE);
fprintf('  tide stations loaded: %d\n', height(TideStations));
fprintf('  no-station convergence fallback: %d\n', ALLOW_NO_STATION_CONVERGENCE_RUN);
fprintf('  TMD root: %s\n', TMD_ROOT);
fprintf('  fort.22 mode: zero wind, standard pressure (%s)\n', NO_TC_FORT22_LINE);
fprintf('  fort.13 open-boundary Cf target/cap: %.4f / %.4f, band %.1f-%.1f km\n', ...
    Fort13.cf_default * Fort13.max_boost, Fort13.cf_cap, Fort13.inner_km, Fort13.outer_km);
fprintf('  fort.13 near-land open-boundary Cf target/cap: %.4f / %.4f, band %.1f-%.1f km\n', ...
    Fort13.near_land_target_cf, Fort13.near_land_cf_cap, ...
    Fort13.near_land_inner_km, Fort13.near_land_outer_km);
fprintf('  dry run: %d\n', DRY_RUN);

AllRows = {};

%% ------------------------- process mesh blocks -------------------------
for iblock = 1:numel(MeshDirs)
    block_name = MeshDirs(iblock).name;
    source_block_dir = fullfile(MESH_ROOT, block_name);
    block_dir = fullfile(OUTPUT_MESH_ROOT, block_name);

    fprintf('\n====================================================\n');
    fprintf('Block %d/%d: %s\n', iblock, numel(MeshDirs), block_name);
    fprintf('Source directory: %s\n', source_block_dir);
    fprintf('Output directory: %s\n', block_dir);
    fprintf('====================================================\n');

    row = init_summary_row(block_name, block_dir);

    try
        if ~DRY_RUN
            ensure_dir(block_dir);
        end

        fort14_file = resolve_block_fort14(source_block_dir, block_name);
        row.fort14 = true;
        fprintf('  fort.14: %s\n', fort14_file);

        out14 = fullfile(block_dir, 'fort.14');
        if ~DRY_RUN && (OVERWRITE_EXISTING || exist(out14, 'file') ~= 2)
            copyfile(fort14_file, out14, 'f');
            fprintf('  copied fort.14: %s\n', out14);
        end

        G = read_fort14_nodes_and_open_boundary(fort14_file);
        row.open_boundary_nodes = numel(G.open_nodes_ordered);
        fprintf('  mesh NP=%d, NE=%d, open-boundary nodes=%d\n', ...
            G.NP, G.NE, row.open_boundary_nodes);

        if NEED_TIDE_SETUP
            [~, ~, Grid] = build_met_grid_from_mesh(G, ...
                MET_GRID_MARGIN_DEG, DLON, DLAT);
            fprintf('  ADCIRC bbox lon [%.3f, %.3f], lat [%.3f, %.3f]\n', ...
                Grid.lon_min, Grid.lon_max, Grid.lat_min, Grid.lat_max);
            fprintf('  fort.22 grid NWLON=%d, NWLAT=%d, WTIMINC=%d s\n', ...
                Grid.NWLON, Grid.NWLAT, WTIMINC);
            row.NWLON = Grid.NWLON;
            row.NWLAT = Grid.NWLAT;

            [offseason_month, season_label] = offseason_month_for_block(block_name);
            [SelectedStations, WindowStart, WindowEnd, station_note] = select_tide_stations_for_block( ...
                source_block_dir, block_name, TideStations, offseason_month, ...
                PREFERRED_VALIDATION_YEAR, MAX_VALIDATION_STATIONS_PER_BLOCK, ...
                ALLOW_NO_STATION_CONVERGENCE_RUN);
            [SelectedStations, station_mesh_note] = enforce_selected_stations_on_mesh( ...
                SelectedStations, G, STATION_SNAP_MAX_DISTANCE_M);
            if isempty(SelectedStations) && ALLOW_NO_STATION_CONVERGENCE_RUN
                [SelectedStations, WindowStart, WindowEnd, fallback_note] = no_station_convergence_selection( ...
                    TideStations, offseason_month, PREFERRED_VALIDATION_YEAR, station_mesh_note);
                station_note = sprintf('%s; %s', station_note, fallback_note);
            else
                station_note = sprintf('%s; %s', station_note, station_mesh_note);
            end

            selected_station_csv = fullfile(block_dir, 'selected_tide_stations.csv');
            if ~DRY_RUN
                writetable(SelectedStations, selected_station_csv);
            end

            nSteps = round(seconds(WindowEnd - WindowStart) / WTIMINC) + 1;
            row.window_start = datestr(WindowStart, 31);
            row.window_end = datestr(WindowEnd, 31);
            row.selected_station_count = height(SelectedStations);
            if height(SelectedStations) > 0
                row.primary_station = char(SelectedStations.station_label(1));
                row.primary_station_source = char(SelectedStations.source(1));
            else
                row.primary_station = 'none';
                row.primary_station_source = 'none';
            end
            row.selected_station_csv = selected_station_csv;
            fprintf('  selected tide stations: %d (%s)\n', ...
                row.selected_station_count, station_note);
            fprintf('  primary station: %s [%s]\n', ...
                row.primary_station, row.primary_station_source);
            fprintf('  tide window: %s -> %s (%s, nSteps=%d)\n', ...
                row.window_start, row.window_end, season_label, nSteps);

            meta_file = fullfile(block_dir, 'fort22_meta.txt');
            if ~DRY_RUN && (OVERWRITE_EXISTING || ~exist(meta_file, 'file'))
                write_tide_only_fort22_meta(meta_file, fort14_file, block_name, Grid, ...
                    WTIMINC, WindowStart, WindowEnd, season_label, nSteps, ...
                    NO_TC_FORT22_LINE, SelectedStations, selected_station_csv);
            end
            row.fort22_meta = ~DRY_RUN;
        else
            Grid = [];
            SelectedStations = table();
            meta_file = fullfile(block_dir, 'fort22_meta.txt');
            fprintf('  tide setup skipped; fort.13-only mode\n');
        end

        if MAKE_FORT22 && ~DRY_RUN
            out22 = fullfile(block_dir, 'fort.22');
            if OVERWRITE_EXISTING || ~exist(out22, 'file')
                nLines = write_no_tc_fort22_file(out22, Grid.NWLAT, Grid.NWLON, ...
                    nSteps, NO_TC_FORT22_LINE);
            else
                nLines = count_text_lines(out22);
            end
            row.fort22 = true;
            row.fort22_lines = nLines;
            fprintf('  wrote fort.22 lines: %d\n', nLines);
        end

        if MAKE_FORT13 && ~DRY_RUN
            out13 = fullfile(block_dir, 'fort.13');
            if OVERWRITE_EXISTING || ~exist(out13, 'file')
                if any(strcmpi(block_name, COPY_EXISTING_FORT13_BLOCKS))
                    src13 = fullfile(source_block_dir, 'fort.13');
                    if exist(src13, 'file') ~= 2
                        error('Missing source fort.13 for copy-existing block %s: %s', block_name, src13);
                    end
                    copyfile(src13, out13, 'f');
                    fprintf('  copied existing locally strengthened fort.13: %s\n', src13);
                else
                    BlockFort13 = fort13_options_for_block(Fort13, block_name);
                    fprintf('  fort.13 strengthened open-boundary Cf target/cap: %.4f / %.4f, band %.1f-%.1f km\n', ...
                        BlockFort13.cf_default * BlockFort13.max_boost, BlockFort13.cf_cap, ...
                        BlockFort13.inner_km, BlockFort13.outer_km);
                    fprintf('  fort.13 strengthened near-land Cf target/cap: %.4f / %.4f, band %.1f-%.1f km\n', ...
                        BlockFort13.near_land_target_cf, BlockFort13.near_land_cf_cap, ...
                        BlockFort13.near_land_inner_km, BlockFort13.near_land_outer_km);
                    backup13 = fullfile(block_dir, 'fort.13.before_strong_open_boundary');
                    if exist(out13, 'file') == 2 && exist(backup13, 'file') ~= 2
                        copyfile(out13, backup13);
                        fprintf('  backed up original fort.13: %s\n', backup13);
                    end
                    make_fort13_for_mesh(G, out13, BlockFort13);
                end
            end
            row.fort13 = true;
        end

        if MAKE_FORT15 && ~DRY_RUN
            out15 = fullfile(block_dir, 'fort.15');
            if OVERWRITE_EXISTING || ~exist(out15, 'file')
                generate_fort15_from_meta(fort15_template_lines, meta_file, out15, ...
                    block_name, FORCE_NWS6, SYNC_OUTPUT_END_TO_RNDAY, ...
                    AUTO_UPDATE_DRAMP, SelectedStations, ZERO_VELOCITY_MET_STATIONS);
            end
            row.fort15 = true;
        end

        if MAKE_FORT19 && ~DRY_RUN
            out19 = fullfile(block_dir, 'fort.19');
            if OVERWRITE_EXISTING || ~exist(out19, 'file')
                generate_fort19_from_tmd(G, meta_file, out19, TMD_ROOT, ...
                    TMD_LATLON_FILE, TMD_OUT_DIR, CLEAN_TMD_BEFORE_EACH_CASE, ...
                    ALLOW_TMD_TRIM_OR_PAD);
            end
            row.fort19 = true;
        end

        if COPY_SUB_INTEL && ~DRY_RUN
            out_sub = fullfile(block_dir, 'sub_intel.sh');
            if OVERWRITE_EXISTING || ~exist(out_sub, 'file')
                copyfile(SUB_INTEL_TEMPLATE, out_sub, 'f');
            end
            row.sub_intel = true;
            fprintf('  copied sub_intel.sh\n');
        end

        row.status = 'ok';
        row.message = '';
        fprintf('  -> done\n');

    catch ME
        row.status = 'fail';
        row.message = ME.message;
        fprintf('  -> failed: %s\n', ME.message);
    end

    AllRows(end+1,:) = summary_row_to_cell(row); %#ok<AGROW>
    fclose all;
end

Summary = cell2table(AllRows, 'VariableNames', summary_var_names());
if ~DRY_RUN
    writetable(Summary, SUMMARY_CSV);
end

fprintf('\n====================================================\n');
fprintf('ADCIRC condition build finished.\n');
if DRY_RUN
    fprintf('Dry-run summary was not written to disk.\n');
else
    fprintf('Summary: %s\n', SUMMARY_CSV);
end
fprintf('====================================================\n');

%% ========================================================================
% local functions
%% ========================================================================

function val = env_bool(name, default_value)
s = getenv(name);
if isempty(s)
    val = default_value;
    return;
end
s = lower(strtrim(s));
val = any(strcmp(s, {'1','true','yes','on'}));
end

function val = env_path(name, default_value)
raw = strtrim(string(getenv(name)));
if strlength(raw) > 0
    val = char(raw);
else
    val = char(default_value);
end
end

function val = env_double(name, default_value)
s = getenv(name);
if isempty(s)
    val = default_value;
    return;
end
x = str2double(s);
if isfinite(x)
    val = x;
else
    val = default_value;
end
end

function ensure_dir(path_name)
if exist(path_name, 'dir') ~= 7
    mkdir(path_name);
end
end

function Opt = default_fort13_options()
Opt.slope_limiter_default = 0.01;
Opt.cf_default = 0.0025;
Opt.use_shallow_manning = true;
Opt.manning_n = 0.01;
Opt.shallow_depth_max = 20.0;
Opt.depth_floor = 5.0;
Opt.cf_cap = 0.050;
Opt.inner_km = 30.0;
Opt.outer_km = 50.0;
Opt.max_boost = 0.020 / Opt.cf_default;
Opt.use_near_land_open_boundary_boost = true;
Opt.open_land_distance_threshold_deg = 1.0;
Opt.near_land_inner_km = 30.0;
Opt.near_land_outer_km = 50.0;
Opt.near_land_target_cf = 0.04;  % 0.03 
Opt.near_land_cf_cap = 0.050;
Opt.use_cf_uprange = true;
Opt.cf_uprange_shp = fullfile(fileparts(mfilename('fullpath')), 'external', 'CF_uprange.shp');
Opt.local_boost = 10.0;
Opt.local_cf_cap = 0.050;
Opt.write_only_nondefault = true;
Opt.chunk_size = 20000;
Opt.use_strong_open_boundary = false;
end

function Opt = apply_fort13_env_overrides(Opt)
Opt.cf_uprange_shp = env_path('ADCIRC_FORT13_CF_UPRANGE_SHP', Opt.cf_uprange_shp);
Opt.cf_cap = env_double('ADCIRC_FORT13_CF_CAP', Opt.cf_cap);
Opt.inner_km = env_double('ADCIRC_FORT13_OPEN_INNER_KM', Opt.inner_km);
Opt.outer_km = env_double('ADCIRC_FORT13_OPEN_OUTER_KM', Opt.outer_km);
open_target_cf = env_double('ADCIRC_FORT13_OPEN_TARGET_CF', Opt.cf_default * Opt.max_boost);
Opt.max_boost = max(1.0, open_target_cf / Opt.cf_default);

Opt.near_land_inner_km = env_double('ADCIRC_FORT13_NEAR_LAND_INNER_KM', Opt.near_land_inner_km);
Opt.near_land_outer_km = env_double('ADCIRC_FORT13_NEAR_LAND_OUTER_KM', Opt.near_land_outer_km);
Opt.near_land_target_cf = env_double('ADCIRC_FORT13_NEAR_LAND_TARGET_CF', Opt.near_land_target_cf);
Opt.near_land_cf_cap = env_double('ADCIRC_FORT13_NEAR_LAND_CF_CAP', Opt.near_land_cf_cap);

Opt.local_boost = env_double('ADCIRC_FORT13_LOCAL_BOOST', Opt.local_boost);
Opt.local_cf_cap = env_double('ADCIRC_FORT13_LOCAL_CF_CAP', Opt.local_cf_cap);

if Opt.outer_km <= Opt.inner_km
    error('ADCIRC_FORT13_OPEN_OUTER_KM must be larger than ADCIRC_FORT13_OPEN_INNER_KM.');
end
if Opt.near_land_outer_km <= Opt.near_land_inner_km
    error('ADCIRC_FORT13_NEAR_LAND_OUTER_KM must be larger than ADCIRC_FORT13_NEAR_LAND_INNER_KM.');
end
end

function Opt = fort13_options_for_block(Opt, block_name)
if ~Opt.use_strong_open_boundary
    return;
end

Opt.cf_cap = max(Opt.cf_cap, 0.080);
Opt.max_boost = max(Opt.max_boost, 0.050 / Opt.cf_default);

Opt.near_land_target_cf = max(Opt.near_land_target_cf, 0.070);
Opt.near_land_cf_cap = max(Opt.near_land_cf_cap, 0.080);
Opt.local_cf_cap = max(Opt.local_cf_cap, 0.080);
end

%% ------------------------- mesh block and season setup -------------------------

function assert_tmd_ready(tmd_root, latlon_file)
if ~exist(tmd_root, 'dir')
    error('TMD root does not exist: %s', tmd_root);
end
latlon_dir = fileparts(latlon_file);
if ~exist(latlon_dir, 'dir')
    error('TMD LAT_LON directory does not exist: %s', latlon_dir);
end
end

function D = list_mesh_block_dirs(mesh_root)
D0 = dir(mesh_root);
D0 = D0([D0.isdir]);
D0 = D0(~ismember({D0.name}, {'.', '..'}));

keep = false(numel(D0), 1);
for i = 1:numel(D0)
    d = fullfile(mesh_root, D0(i).name);
    keep(i) = exist(fullfile(d, 'fort.14'), 'file') == 2 || ...
        exist(fullfile(d, ['ADCIRC_Capsule_Mesh_' D0(i).name '.14']), 'file') == 2;
end
D0 = D0(keep);
[~, order] = sort({D0.name});
D0 = D0(order);

D = D0;
end

function fort14_file = resolve_block_fort14(block_dir, block_name)
cand = {
    fullfile(block_dir, 'fort.14')
    fullfile(block_dir, ['ADCIRC_Capsule_Mesh_' block_name '.14'])
    };

for i = 1:numel(cand)
    if exist(cand{i}, 'file')
        fort14_file = cand{i};
        return;
    end
end

D = dir(fullfile(block_dir, '*.14'));
if numel(D) == 1
    fort14_file = fullfile(block_dir, D(1).name);
else
    names = string({D.name});
    error('Cannot resolve fort.14 for block %s. Found: %s', ...
        block_name, strjoin(names, ', '));
end
end

function [month_num, season_label] = offseason_month_for_block(block_name)
name = upper(strtrim(block_name));

if startsWith(name, 'ADC_NATL')
    season_label = 'NATL off-season January';
    month_num = 1;
elseif startsWith(name, 'ADC_ENP')
    season_label = 'ENP off-season February';
    month_num = 2;
elseif startsWith(name, 'ADC_WNP')
    season_label = 'WNP low-activity February';
    month_num = 2;
elseif startsWith(name, 'ADC_NIO')
    season_label = 'NIO low-activity February';
    month_num = 2;
elseif startsWith(name, 'ADC_SIO')
    season_label = 'SIO off-season July';
    month_num = 7;
elseif startsWith(name, 'ADC_AUSSP')
    season_label = 'AUSSP off-season July';
    month_num = 7;
else
    season_label = 'default off-season January';
    month_num = 1;
end
end

function [MET_LAT, MET_LON, Grid] = build_met_grid_from_mesh(Mesh, margin_deg, dlon, dlat)
Grid = struct();
Grid.lon_min = min(Mesh.lon);
Grid.lon_max = max(Mesh.lon);
Grid.lat_min = min(Mesh.lat);
Grid.lat_max = max(Mesh.lat);

Grid.WLONMIN = floor((Grid.lon_min - margin_deg) / dlon) * dlon;
Grid.WLONMAX = ceil((Grid.lon_max + margin_deg) / dlon) * dlon;
Grid.WLATMIN = floor((Grid.lat_min - margin_deg) / dlat) * dlat;
Grid.WLATMAX = ceil((Grid.lat_max + margin_deg) / dlat) * dlat;
Grid.DLON = dlon;
Grid.DLAT = dlat;

met_lon = Grid.WLONMIN : dlon : Grid.WLONMAX;
met_lat = Grid.WLATMAX : -dlat : Grid.WLATMIN;

Grid.NWLON = numel(met_lon);
Grid.NWLAT = numel(met_lat);
Grid.met_lon = met_lon;
Grid.met_lat = met_lat;

[MET_LON, MET_LAT] = meshgrid(met_lon, met_lat);
end

%% ------------------------- tide station catalog and selection -------------------------

function T = load_tide_station_catalog(tide_root, use_uhslc, use_noaa)
T = empty_station_catalog();

if use_uhslc
    uhslc_summary = fullfile(tide_root, 'UHSLC_hourly_fast_lat_-60_60', ...
        'all_station_time_space_plot', 'station_time_space_summary.csv');
    if exist(uhslc_summary, 'file')
        T = [T; load_uhslc_station_catalog(uhslc_summary)]; %#ok<AGROW>
    else
        warning('UHSLC station summary not found: %s', uhslc_summary);
    end
end

if use_noaa
    noaa_log = fullfile(tide_root, 'NOAA_COOPS_hourly_IBTrACS_3deg_2000_2025', ...
        'download_log.csv');
    if exist(noaa_log, 'file')
        T = [T; load_noaa_station_catalog(noaa_log)]; %#ok<AGROW>
    else
        warning('NOAA download log not found: %s', noaa_log);
    end
end

good = isfinite(T.station_lat) & isfinite(T.station_lon) & strlength(T.file_path) > 0;
T = T(good, :);
end

function T = empty_station_catalog()
T = table('Size', [0, 13], ...
    'VariableTypes', {'string','string','string','double','double', ...
    'string','string','datetime','datetime','double','double','string','string'}, ...
    'VariableNames', {'source','station_id','station_name','station_lat','station_lon', ...
    'file_name','file_path','obs_start','obs_end','n_valid_lines','duration_days', ...
    'status','station_label'});
end

function T = load_uhslc_station_catalog(summary_csv)
S = readtable(summary_csv, 'TextType', 'string', 'VariableNamingRule', 'preserve');
n = height(S);
if n == 0
    T = empty_station_catalog();
    return;
end

source = repmat("UHSLC", n, 1);
station_id = string(S.uhslc_id);
station_name = string(S.station_name);
station_lat = double(S.station_lat);
station_lon = double(S.station_lon);
file_name = string(S.file_name);
file_path = string(S.file_path);
obs_start = parse_datetime_column(S.obs_start);
obs_end = parse_datetime_column(S.obs_end);
n_valid_lines = double(S.n_valid_lines);
duration_days = double(S.duration_days);
status = string(S.status);
station_label = strings(n, 1);

T = table(source, station_id, station_name, station_lat, station_lon, ...
    file_name, file_path, obs_start, obs_end, n_valid_lines, ...
    duration_days, status, station_label, ...
    'VariableNames', {'source','station_id','station_name','station_lat','station_lon', ...
    'file_name','file_path','obs_start','obs_end','n_valid_lines','duration_days', ...
    'status','station_label'});
T.station_label = make_station_labels(T);

good = strcmpi(T.status, "ok") & ~isnat(T.obs_start) & ~isnat(T.obs_end);
T = T(good, :);
end

function T = load_noaa_station_catalog(download_log_csv)
S = readtable(download_log_csv, 'TextType', 'string', 'VariableNamingRule', 'preserve');
downloaded = logical_from_column(S.downloaded);
S = S(downloaded, :);
n = height(S);
if n == 0
    T = empty_station_catalog();
    return;
end

source = repmat("NOAA", n, 1);
station_id = string(S.station_id);
station_name = string(S.station_name);
station_lat = double(S.lat);
station_lon = double(S.lon);
file_path = string(S.file);
file_name = strings(n, 1);
for i = 1:n
    [~, nm, ext] = fileparts(file_path(i));
    file_name(i) = string([char(nm) char(ext)]);
end
obs_start = NaT(n, 1);
obs_end = NaT(n, 1);
n_valid_lines = nan(n, 1);
duration_days = nan(n, 1);
status = repmat("needs_scan", n, 1);
station_label = strings(n, 1);

T = table(source, station_id, station_name, station_lat, station_lon, ...
    file_name, file_path, obs_start, obs_end, n_valid_lines, ...
    duration_days, status, station_label, ...
    'VariableNames', {'source','station_id','station_name','station_lat','station_lon', ...
    'file_name','file_path','obs_start','obs_end','n_valid_lines','duration_days', ...
    'status','station_label'});
T.station_label = make_station_labels(T);
end

function tf = logical_from_column(x)
if islogical(x)
    tf = x;
elseif isnumeric(x)
    tf = x ~= 0;
else
    s = lower(strtrim(string(x)));
    tf = s == "true" | s == "1" | s == "yes";
end
tf = tf(:);
end

function labels = make_station_labels(T)
n = height(T);
labels = strings(n, 1);
for i = 1:n
    sid = strtrim(T.station_id(i));
    sname = strtrim(T.station_name(i));
    if strlength(sname) == 0
        labels(i) = T.source(i) + "_" + sid;
    else
        labels(i) = T.source(i) + "_" + sid + "_" + sname;
    end
end
end

function [Selected, WindowStart, WindowEnd, note] = select_tide_stations_for_block( ...
    block_dir, block_name, TideStations, offseason_month, preferred_year, max_stations, ...
    allow_no_station_fallback)

if nargin < 7
    allow_no_station_fallback = false;
end

inner_file = fullfile(block_dir, 'inner_ring.csv');
if exist(inner_file, 'file') ~= 2
    error('Missing inner refined-region polygon: %s', inner_file);
end

R = readtable(inner_file, 'TextType', 'string', 'VariableNamingRule', 'preserve');
if ~all(ismember({'lon','lat'}, R.Properties.VariableNames))
    error('inner_ring.csv must contain lon and lat columns: %s', inner_file);
end
ring_lon = double(R.lon);
ring_lat = double(R.lat);
good_ring = isfinite(ring_lon) & isfinite(ring_lat);
ring_lon = ring_lon(good_ring);
ring_lat = ring_lat(good_ring);
if numel(ring_lon) < 3
    error('inner_ring.csv does not contain a valid polygon: %s', inner_file);
end

inside = inpolygon(TideStations.station_lon, TideStations.station_lat, ring_lon, ring_lat);
Candidates = TideStations(inside, :);
if isempty(Candidates)
    if allow_no_station_fallback
        [Selected, WindowStart, WindowEnd, note] = no_station_convergence_selection( ...
            TideStations, offseason_month, preferred_year, ...
            sprintf('no inner-ring tide station for %s', block_name));
        return;
    end
    error('No tide-level station is inside the inner refined region for %s.', block_name);
end

Candidates = refine_station_time_ranges(Candidates);
Candidates = Candidates(~isnat(Candidates.obs_start) & ~isnat(Candidates.obs_end), :);
if isempty(Candidates)
    if allow_no_station_fallback
        [Selected, WindowStart, WindowEnd, note] = no_station_convergence_selection( ...
            TideStations, offseason_month, preferred_year, ...
            sprintf('no valid inner-ring tide station time span for %s', block_name));
        return;
    end
    error('No inner-region tide station with valid observation time span for %s.', block_name);
end

cent_lon = mean(ring_lon, 'omitnan');
cent_lat = mean(ring_lat, 'omitnan');
Candidates.inner_centroid_dist_deg = hypot(Candidates.station_lon - cent_lon, ...
    Candidates.station_lat - cent_lat);

try
    [WindowStart, WindowEnd] = choose_validation_month(Candidates, offseason_month, preferred_year);
catch ME
    if allow_no_station_fallback
        [Selected, WindowStart, WindowEnd, note] = no_station_convergence_selection( ...
            TideStations, offseason_month, preferred_year, ME.message);
        return;
    end
    rethrow(ME);
end
cover = Candidates.obs_start <= WindowStart & Candidates.obs_end >= WindowEnd;
Selected = Candidates(cover, :);
if isempty(Selected)
    if allow_no_station_fallback
        [Selected, WindowStart, WindowEnd, note] = no_station_convergence_selection( ...
            TideStations, offseason_month, preferred_year, ...
            sprintf('no selected station covers chosen window for %s', block_name));
        return;
    end
    error('Internal error: no selected station covers chosen window for %s.', block_name);
end

[~, order] = sortrows([-Selected.duration_days, Selected.inner_centroid_dist_deg]);
Selected = Selected(order, :);
if isfinite(max_stations) && height(Selected) > max_stations
    Selected = Selected(1:max_stations, :);
end

Selected.validation_window_start = repmat(WindowStart, height(Selected), 1);
Selected.validation_window_end = repmat(WindowEnd, height(Selected), 1);

note = sprintf('inner-ring candidates=%d, covering window=%d', ...
    height(Candidates), height(Selected));
end

function [Selected, WindowStart, WindowEnd, note] = no_station_convergence_selection( ...
    TideStations, month_num, preferred_year, reason)
Selected = TideStations([], :);
[WindowStart, WindowEnd] = default_convergence_month(month_num, preferred_year);
Selected.validation_window_start = repmat(WindowStart, 0, 1);
Selected.validation_window_end = repmat(WindowEnd, 0, 1);
note = sprintf('no measured station; convergence-only one-month tide run (%s)', reason);
end

function [Selected, note] = enforce_selected_stations_on_mesh(Selected, Mesh, max_snap_distance_m)
Selected = ensure_station_mesh_columns(Selected);
if isempty(Selected)
    note = 'no selected station to mesh-check';
    return;
end

if ~isfield(Mesh, 'elem_node_idx') || isempty(Mesh.elem_node_idx)
    error('Mesh triangle connectivity is unavailable; cannot validate station locations.');
end

keep = true(height(Selected), 1);
snapped_count = 0;
dropped_count = 0;
for i = 1:height(Selected)
    lon0 = Selected.original_station_lon(i);
    lat0 = Selected.original_station_lat(i);
    [lon_use, lat_use, inside, snap_distance_m, elem_id, snap_mode] = ...
        snap_point_to_adc_mesh(lon0, lat0, Mesh);

    if inside || snap_distance_m <= max_snap_distance_m
        Selected.station_lon(i) = lon_use;
        Selected.station_lat(i) = lat_use;
        Selected.adcirc_station_lon(i) = lon_use;
        Selected.adcirc_station_lat(i) = lat_use;
        Selected.adcirc_station_inside_mesh(i) = true;
        Selected.adcirc_station_snap_distance_m(i) = snap_distance_m;
        Selected.adcirc_station_element(i) = elem_id;
        Selected.adcirc_station_snap_mode(i) = string(snap_mode);
        if ~inside
            snapped_count = snapped_count + 1;
        end
    else
        keep(i) = false;
        dropped_count = dropped_count + 1;
        Selected.adcirc_station_inside_mesh(i) = false;
        Selected.adcirc_station_snap_distance_m(i) = snap_distance_m;
        Selected.adcirc_station_snap_mode(i) = "outside_mesh_too_far";
    end
end

Selected = Selected(keep, :);
note = sprintf('mesh station check: kept=%d, snapped=%d, dropped=%d, max_snap=%.0f m', ...
    height(Selected), snapped_count, dropped_count, max_snap_distance_m);
end

function T = ensure_station_mesh_columns(T)
n = height(T);
if ~ismember('original_station_lon', T.Properties.VariableNames)
    T.original_station_lon = T.station_lon;
end
if ~ismember('original_station_lat', T.Properties.VariableNames)
    T.original_station_lat = T.station_lat;
end
if ~ismember('adcirc_station_lon', T.Properties.VariableNames)
    T.adcirc_station_lon = nan(n, 1);
end
if ~ismember('adcirc_station_lat', T.Properties.VariableNames)
    T.adcirc_station_lat = nan(n, 1);
end
if ~ismember('adcirc_station_inside_mesh', T.Properties.VariableNames)
    T.adcirc_station_inside_mesh = false(n, 1);
end
if ~ismember('adcirc_station_snap_distance_m', T.Properties.VariableNames)
    T.adcirc_station_snap_distance_m = nan(n, 1);
end
if ~ismember('adcirc_station_element', T.Properties.VariableNames)
    T.adcirc_station_element = nan(n, 1);
end
if ~ismember('adcirc_station_snap_mode', T.Properties.VariableNames)
    T.adcirc_station_snap_mode = strings(n, 1);
end
end

function [lon_use, lat_use, inside, snap_distance_m, elem_id, snap_mode] = ...
    snap_point_to_adc_mesh(lon0, lat0, Mesh)
tri = Mesh.elem_node_idx;
valid = all(isfinite(tri), 2);
tri = tri(valid, :);
elem_ids = Mesh.elem_id(valid);
if isempty(tri)
    error('No valid triangular elements found in mesh.');
end

lon = Mesh.lon;
lat = Mesh.lat;
scale_x = cosd(lat0) * 111320;
scale_y = 110540;

x1 = (lon(tri(:,1)) - lon0) .* scale_x;
y1 = (lat(tri(:,1)) - lat0) .* scale_y;
x2 = (lon(tri(:,2)) - lon0) .* scale_x;
y2 = (lat(tri(:,2)) - lat0) .* scale_y;
x3 = (lon(tri(:,3)) - lon0) .* scale_x;
y3 = (lat(tri(:,3)) - lat0) .* scale_y;

den = (y2 - y3) .* (x1 - x3) + (x3 - x2) .* (y1 - y3);
good_den = abs(den) > eps;
a = nan(size(den));
b = nan(size(den));
a(good_den) = ((y2(good_den) - y3(good_den)) .* (-x3(good_den)) + ...
    (x3(good_den) - x2(good_den)) .* (-y3(good_den))) ./ den(good_den);
b(good_den) = ((y3(good_den) - y1(good_den)) .* (-x3(good_den)) + ...
    (x1(good_den) - x3(good_den)) .* (-y3(good_den))) ./ den(good_den);
c = 1 - a - b;
inside_idx = find(good_den & a >= -1e-10 & b >= -1e-10 & c >= -1e-10, 1, 'first');
if ~isempty(inside_idx)
    lon_use = lon0;
    lat_use = lat0;
    inside = true;
    snap_distance_m = 0.0;
    elem_id = elem_ids(inside_idx);
    snap_mode = 'inside_mesh';
    return;
end

[d2_12, ~, ~] = point_to_segment_distance2_origin(x1, y1, x2, y2);
[d2_23, ~, ~] = point_to_segment_distance2_origin(x2, y2, x3, y3);
[d2_31, ~, ~] = point_to_segment_distance2_origin(x3, y3, x1, y1);

[d2_edge, ~] = min([d2_12, d2_23, d2_31], [], 2);
[d2_min, tri_idx] = min(d2_edge);
snap_distance_m = sqrt(d2_min);

cx = (x1(tri_idx) + x2(tri_idx) + x3(tri_idx)) / 3;
cy = (y1(tri_idx) + y2(tri_idx) + y3(tri_idx)) / 3;
x_use = cx;
y_use = cy;

lon_use = lon0 + x_use / scale_x;
lat_use = lat0 + y_use / scale_y;
inside = false;
elem_id = elem_ids(tri_idx);
snap_mode = 'snapped_to_nearest_mesh_element_centroid';
end

function [d2, qx, qy] = point_to_segment_distance2_origin(ax, ay, bx, by)
vx = bx - ax;
vy = by - ay;
l2 = vx.^2 + vy.^2;
t = zeros(size(l2));
good = l2 > eps;
t(good) = -(ax(good) .* vx(good) + ay(good) .* vy(good)) ./ l2(good);
t = max(0, min(1, t));
qx = ax + t .* vx;
qy = ay + t .* vy;
d2 = qx.^2 + qy.^2;
end

function [WindowStart, WindowEnd] = default_convergence_month(month_num, preferred_year)
if ~isfinite(preferred_year)
    preferred_year = 2021;
end
WindowStart = datetime(preferred_year, month_num, 1, 0, 0, 0);
WindowEnd = WindowStart + calmonths(1);
end

function Candidates = refine_station_time_ranges(Candidates)
for i = 1:height(Candidates)
    needs_scan = isnat(Candidates.obs_start(i)) || isnat(Candidates.obs_end(i));
    if ~needs_scan
        continue;
    end

    f = char(Candidates.file_path(i));
    if exist(f, 'file') ~= 2
        Candidates.status(i) = "missing_file";
        continue;
    end

    try
        Info = quick_scan_station_csv(f);
        Candidates.obs_start(i) = Info.obs_start;
        Candidates.obs_end(i) = Info.obs_end;
        Candidates.n_valid_lines(i) = Info.n_valid_lines;
        Candidates.duration_days(i) = Info.duration_days;
        Candidates.status(i) = "ok";
    catch ME
        Candidates.status(i) = "scan_fail: " + string(ME.message);
    end
end
end

function [WindowStart, WindowEnd] = choose_validation_month(Candidates, month_num, preferred_year)
year_min = min(year(Candidates.obs_start));
year_max = max(year(Candidates.obs_end));
years_all = year_min:year_max;
years_all = years_all(isfinite(years_all));
years_all = years_all(:).';
years_try = unique([preferred_year, sort(years_all, 'descend')], 'stable');

for y = years_try
    if ~isfinite(y)
        continue;
    end
    t0 = datetime(y, month_num, 1, 0, 0, 0);
    t1 = t0 + calmonths(1);
    cover = Candidates.obs_start <= t0 & Candidates.obs_end >= t1;
    if any(cover)
        WindowStart = t0;
        WindowEnd = t1;
        return;
    end
end

error('No inner-region tide station covers a complete off-season month.');
end

function Info = quick_scan_station_csv(file_path)
T = readtable(file_path, 'TextType', 'string', 'VariableNamingRule', 'preserve');
if isempty(T) || width(T) == 0
    error('empty station table');
end

raw_names = string(T.Properties.VariableNames);
norm_names = normalize_header_names(raw_names);
itime = find(contains(norm_names, "time"), 1, 'first');
iwl = find(contains(norm_names, "sea_level") | contains(norm_names, "water_level") | ...
    norm_names == "v", 1, 'first');
if isempty(itime)
    error('missing time column');
end
if isempty(iwl)
    error('missing sea-level column');
end

time_vec = parse_datetime_column(T{:, itime});
wl_vec = numeric_from_column(T{:, iwl});
Ivalid = ~isnat(time_vec) & isfinite(wl_vec);
if ~any(Ivalid)
    error('no valid time-water-level records');
end

t_valid = time_vec(Ivalid);
Info = struct();
Info.obs_start = min(t_valid);
Info.obs_end = max(t_valid);
Info.n_valid_lines = nnz(Ivalid);
Info.duration_days = days(Info.obs_end - Info.obs_start);
end

function t = parse_datetime_column(x)
if isdatetime(x)
    t = x;
    t.TimeZone = '';
    return;
end

xs = string(x);
xs = strtrim(xs);
t = NaT(size(xs));
I = strlength(xs) > 0 & xs ~= "<missing>";
if ~any(I)
    return;
end

xs2 = xs(I);
xs2 = replace(xs2, "Z", "");
tmp = regexprep(cellstr(xs2), '([+-]\d{2}:\d{2})$', '');
xs2 = string(tmp);
xs2 = replace(xs2, "T", " ");

formats = ["yyyy-MM-dd HH:mm:ss", "yyyy-MM-dd HH:mm", ...
    "yyyy/MM/dd HH:mm:ss", "yyyy/MM/dd HH:mm"];
for ifmt = 1:numel(formats)
    try
        t(I) = datetime(xs2, 'InputFormat', formats(ifmt));
        return;
    catch
    end
end

idx = find(I(:)).';
for k = idx
    try
        s = xs(k);
        s = replace(s, "Z", "");
        s = regexprep(char(s), '([+-]\d{2}:\d{2})$', '');
        s = strrep(s, 'T', ' ');
        t(k) = datetime(s);
    catch
        t(k) = NaT;
    end
end
end

function y = numeric_from_column(x)
if isnumeric(x)
    y = double(x);
else
    y = str2double(string(x));
end
y = y(:);
end

function norm_names = normalize_header_names(raw_names)
norm_names = string(raw_names);
norm_names = replace(norm_names, char(65279), "");
norm_names = strtrim(norm_names);
tmp = regexprep(cellstr(norm_names), '\s*\([^)]*\)\s*', '');
norm_names = string(tmp);
norm_names = lower(norm_names);
norm_names = replace(norm_names, " ", "_");
norm_names = replace(norm_names, "-", "_");
tmp = regexprep(cellstr(norm_names), '_+', '_');
tmp = regexprep(tmp, '^_+|_+$', '');
norm_names = string(tmp);
end

%% ------------------------- fort.22 background forcing -------------------------

function write_tide_only_fort22_meta(meta_txt, fort14_file, block_name, Grid, ...
    WTIMINC, WindowStart, WindowEnd, season_label, nSteps, no_tc_line, ...
    SelectedStations, selected_station_csv)
fid = fopen(meta_txt, 'wt');
if fid < 0
    error('Cannot write %s', meta_txt);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'ADCIRC fort.22 (NWS=6) tide-only metadata\n');
fprintf(fid, '------------------------------------------\n');
fprintf(fid, 'block = %s\n', block_name);
fprintf(fid, 'site = %s\n', block_name);
fprintf(fid, 'mesh_tag = %s\n', block_name);
fprintf(fid, 'fort14 = %s\n', fort14_file);
fprintf(fid, 'forcing_mode = tide_only_no_tc\n');
fprintf(fid, 'season = %s\n', season_label);
fprintf(fid, 'tc_id = tide_only_%s\n', block_name);
fprintf(fid, 'source = tide_only_builder\n');
fprintf(fid, 'scenario = no_tc\n');
fprintf(fid, 'validation_station_constraint = inside_inner_refined_region\n');
fprintf(fid, 'selected_station_csv = %s\n', selected_station_csv);
fprintf(fid, 'selected_station_count = %d\n', height(SelectedStations));
if height(SelectedStations) > 0
    fprintf(fid, 'primary_station = %s\n', char(SelectedStations.station_label(1)));
    fprintf(fid, 'primary_station_source = %s\n', char(SelectedStations.source(1)));
    fprintf(fid, 'primary_station_lon = %.8f\n', SelectedStations.station_lon(1));
    fprintf(fid, 'primary_station_lat = %.8f\n', SelectedStations.station_lat(1));
    if ismember('original_station_lon', SelectedStations.Properties.VariableNames)
        fprintf(fid, 'primary_station_original_lon = %.8f\n', SelectedStations.original_station_lon(1));
        fprintf(fid, 'primary_station_original_lat = %.8f\n', SelectedStations.original_station_lat(1));
    end
    if ismember('adcirc_station_snap_distance_m', SelectedStations.Properties.VariableNames)
        fprintf(fid, 'primary_station_snap_distance_m = %.3f\n', ...
            SelectedStations.adcirc_station_snap_distance_m(1));
        fprintf(fid, 'primary_station_snap_mode = %s\n', ...
            char(SelectedStations.adcirc_station_snap_mode(1)));
        fprintf(fid, 'primary_station_element = %.0f\n', ...
            SelectedStations.adcirc_station_element(1));
    end
    fprintf(fid, 'primary_station_obs_start = %s\n', datestr(SelectedStations.obs_start(1), 31));
    fprintf(fid, 'primary_station_obs_end = %s\n', datestr(SelectedStations.obs_end(1), 31));
else
    fprintf(fid, 'primary_station = none\n');
    fprintf(fid, 'primary_station_source = none\n');
    fprintf(fid, 'primary_station_lon = NaN\n');
    fprintf(fid, 'primary_station_lat = NaN\n');
    fprintf(fid, 'primary_station_obs_start = NaT\n');
    fprintf(fid, 'primary_station_obs_end = NaT\n');
end
fprintf(fid, 'WindowStart = %s\n', datestr(WindowStart, 31));
fprintf(fid, 'WindowEnd = %s\n', datestr(WindowEnd, 31));
fprintf(fid, 'WTIMINC = %d\n', WTIMINC);
fprintf(fid, 'NWLON = %d\n', Grid.NWLON);
fprintf(fid, 'NWLAT = %d\n', Grid.NWLAT);
fprintf(fid, 'WLONMIN = %.2f\n', Grid.WLONMIN);
fprintf(fid, 'WLATMAX = %.2f\n', Grid.WLATMAX);
fprintf(fid, 'WLONINC = %.2f\n', Grid.DLON);
fprintf(fid, 'WLATINC = %.2f\n', Grid.DLAT);
fprintf(fid, 'nSteps = %d\n', nSteps);
fprintf(fid, 'fort22_line = %s\n', no_tc_line);
fprintf(fid, 'note = fort.22 contains zero wind and standard pressure only\n');

fprintf(fid, '\n[SelectedStations]\n');
for i = 1:height(SelectedStations)
    orig_lon = SelectedStations.station_lon(i);
    orig_lat = SelectedStations.station_lat(i);
    snap_dist = 0.0;
    snap_mode = "inside_mesh";
    elem_id = NaN;
    if ismember('original_station_lon', SelectedStations.Properties.VariableNames)
        orig_lon = SelectedStations.original_station_lon(i);
        orig_lat = SelectedStations.original_station_lat(i);
    end
    if ismember('adcirc_station_snap_distance_m', SelectedStations.Properties.VariableNames)
        snap_dist = SelectedStations.adcirc_station_snap_distance_m(i);
        snap_mode = SelectedStations.adcirc_station_snap_mode(i);
        elem_id = SelectedStations.adcirc_station_element(i);
    end
    fprintf(fid, '%d source=%s id=%s lon=%.8f lat=%.8f original_lon=%.8f original_lat=%.8f snap_distance_m=%.3f snap_mode=%s elem=%.0f obs_start=%s obs_end=%s label=%s\n', ...
        i, char(SelectedStations.source(i)), char(SelectedStations.station_id(i)), ...
        SelectedStations.station_lon(i), SelectedStations.station_lat(i), ...
        orig_lon, orig_lat, snap_dist, char(snap_mode), elem_id, ...
        datestr(SelectedStations.obs_start(i), 31), datestr(SelectedStations.obs_end(i), 31), ...
        char(SelectedStations.station_label(i)));
end
end

function nLines = write_no_tc_fort22_file(fort22_file, NWLAT, NWLON, nSteps, no_tc_line)
fid = fopen(fort22_file, 'wt');
if fid < 0
    error('Cannot write fort.22: %s', fort22_file);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

nGrid = NWLAT * NWLON;
oneLine = sprintf('%s\n', no_tc_line);
oneStep = repmat(oneLine, 1, nGrid);
for it = 1:nSteps
    fwrite(fid, oneStep, 'char');
end
nLines = nGrid * nSteps;
end

function nLines = count_text_lines(file_name)
fid = fopen(file_name, 'rt');
if fid < 0
    error('Cannot open file: %s', file_name);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

nLines = 0;
while ischar(fgetl(fid))
    nLines = nLines + 1;
end
end

%% ------------------------- summary -------------------------

function row = init_summary_row(block_name, block_dir)
row = struct();
row.block_name = block_name;
row.block_dir = block_dir;
row.fort14 = false;
row.fort13 = false;
row.fort15 = false;
row.fort19 = false;
row.fort22 = false;
row.fort22_meta = false;
row.sub_intel = false;
row.fort22_lines = NaN;
row.window_start = '';
row.window_end = '';
row.NWLON = NaN;
row.NWLAT = NaN;
row.open_boundary_nodes = NaN;
row.selected_station_count = NaN;
row.primary_station = '';
row.primary_station_source = '';
row.selected_station_csv = '';
row.status = 'pending';
row.message = '';
end

function names = summary_var_names()
names = {'BlockName','BlockDir','Fort14','Fort13','Fort15','Fort19', ...
    'Fort22','Fort22Meta','SubIntel','Fort22Lines','WindowStart','WindowEnd', ...
    'NWLON','NWLAT','OpenBoundaryNodes','SelectedStationCount', ...
    'PrimaryStation','PrimaryStationSource','SelectedStationCsv','Status','Message'};
end

function C = summary_row_to_cell(row)
C = {row.block_name, row.block_dir, row.fort14, row.fort13, row.fort15, ...
    row.fort19, row.fort22, row.fort22_meta, row.sub_intel, row.fort22_lines, ...
    row.window_start, row.window_end, row.NWLON, row.NWLAT, ...
    row.open_boundary_nodes, row.selected_station_count, row.primary_station, ...
    row.primary_station_source, row.selected_station_csv, row.status, row.message};
end

%% ------------------------- fort.15 -------------------------

function generate_fort15_from_meta(template_lines, meta_file, out_file, site_name, ...
    force_NWS6, sync_output_end_to_rnday, auto_update_dramp, ...
    SelectedStations, zero_velocity_met_stations)

if nargin < 9
    zero_velocity_met_stations = true;
end

meta = read_fort22_meta(meta_file);
lines = template_lines;

NWLAT = get_required_numeric(meta, 'NWLAT');
NWLON = get_required_numeric(meta, 'NWLON');
WLATMAX = get_required_numeric(meta, 'WLATMAX');
WLONMIN = get_required_numeric(meta, 'WLONMIN');
WLATINC = get_required_numeric(meta, 'WLATINC');
WLONINC = get_required_numeric(meta, 'WLONINC');
WTIMINC = get_required_numeric(meta, 'WTIMINC');

WindowStart = get_required_datetime(meta, 'WindowStart');
WindowEnd = get_required_datetime(meta, 'WindowEnd');

RNDAY = days(WindowEnd - WindowStart);
if RNDAY <= 0
    error('RNDAY <= 0, check %s', meta_file);
end

tc_id = get_meta_string(meta, 'tc_id', get_file_part(fileparts(meta_file)));
cv_used = get_meta_string(meta, 'Cv_used', '');

run_desc = crop_or_pad(sprintf('%s %s', site_name, tc_id), 32);
lines{1} = sprintf(' %-32s ! 32 CHARACTER ALPHANUMERIC RUN DESCRIPTION', run_desc);

idx_nws = find_line_contains(lines, '! NWS - WIND STRESS AND BAROMETRIC PRESSURE OPTION PARAMETER');
if ~isempty(idx_nws) && force_NWS6
    lines{idx_nws} = ' 6                                   ! NWS - WIND STRESS AND BAROMETRIC PRESSURE OPTION PARAMETER';
end

idx_met = find_line_any(lines, {'WTIMINC', 'WITMINC', 'STIMINC'});
if isempty(idx_met)
    error('Cannot find meteorological grid line in fort.15 template.');
end
lines{idx_met} = sprintf( ...
    ' %d %d %.2f %.2f %.2f %.2f %d    ! NWLAT NWLON WLATMAX WLONMIN WLATINC WLONINC WTIMINC', ...
    NWLAT, NWLON, WLATMAX, WLONMIN, WLATINC, WLONINC, WTIMINC);

idx_rnday = find_line_contains(lines, '! RNDAY - TOTAL LENGTH OF SIMULATION (IN DAYS)');
if isempty(idx_rnday)
    error('Cannot find RNDAY line in fort.15 template.');
end
lines{idx_rnday} = sprintf(' %.8f                                 ! RNDAY - TOTAL LENGTH OF SIMULATION (IN DAYS)', RNDAY);

if auto_update_dramp
    idx_dramp = find_line_contains(lines, '! DRAMP - DURATION OF RAMP FUNCTION (IN DAYS)');
    if ~isempty(idx_dramp)
        DRAMP = min(1.0, RNDAY);
        lines{idx_dramp} = sprintf(' %.8f                                 ! DRAMP - DURATION OF RAMP FUNCTION (IN DAYS)', DRAMP);
    end
end

if sync_output_end_to_rnday
    output_tokens = { ...
        'NOUTE,TOUTSE,TOUTFE', ...
        'NOUTV,TOUTSV,TOUTFV', ...
        'NOUTM,TOUTSM,TOUTFM', ...
        'NOUTGE,TOUTSGE,TOUTFGE', ...
        'NOUTGV,TOUTSGV,TOUTFGV', ...
        'NOUTGW,TOUTSGW,TOUTFGW'};

    for k = 1:numel(output_tokens)
        idx_out = find_line_contains(lines, output_tokens{k});
        if ~isempty(idx_out)
            lines{idx_out} = update_output_schedule_line(lines{idx_out}, RNDAY);
        end
    end
end

lines = set_elevation_station_output_block(lines, SelectedStations);
if zero_velocity_met_stations
    lines = zero_station_block(lines, ...
        'NOUTV,TOUTSV,TOUTFV', ...
        'TOTAL NUMBER OF VELOCITY RECORDING STATIONS');
    lines = zero_station_block(lines, ...
        'NOUTM,TOUTSM,TOUTFM', ...
        'TOTAL NUMBER OF meteorological recording stations');
end

source = get_meta_string(meta, 'source', '');
scenario = get_meta_string(meta, 'scenario', '');
nc_file = get_meta_string(meta, 'nc_file', '');
lines{end+1} = sprintf('! auto-generated for tide-only validation ; block=%s ; tc_id=%s ; Cv_used=%s', ...
    site_name, tc_id, cv_used);
if ~isempty(source) || ~isempty(scenario) || ~isempty(nc_file)
    lines{end+1} = sprintf('! source=%s ; scenario=%s ; nc_file=%s', source, scenario, nc_file);
end

write_text_lines(out_file, lines);
end

function s = get_file_part(p)
[~, s] = fileparts(p);
end

function lines = set_elevation_station_output_block(lines, SelectedStations)
idx_sched = find_line_contains(lines, 'NOUTE,TOUTSE,TOUTFE');
if isempty(idx_sched)
    error('Cannot find elevation station output schedule line in fort.15 template.');
end

idx_count = [];
for i = idx_sched+1:min(idx_sched+5, numel(lines))
    if contains(lines{i}, 'TOTAL NUMBER OF ELEVATION RECORDING STATIONS', 'IgnoreCase', true)
        idx_count = i;
        break;
    end
end
if isempty(idx_count)
    idx_count = idx_sched + 1;
end
if idx_count > numel(lines)
    error('Cannot locate elevation station count line in fort.15 template.');
end

old_count = first_integer_from_line(lines{idx_count});
if ~isfinite(old_count) || old_count < 0
    old_count = 0;
end

if height(SelectedStations) > 0
    lines{idx_sched} = set_first_numeric_token(lines{idx_sched}, 1);
else
    lines{idx_sched} = set_first_numeric_token(lines{idx_sched}, 0);
end
lines{idx_count} = sprintf(' %d                                  ! TOTAL NUMBER OF ELEVATION RECORDING STATIONS', ...
    height(SelectedStations));

station_lines = cell(height(SelectedStations), 1);
for i = 1:height(SelectedStations)
    label = sanitize_station_comment(SelectedStations.station_label(i));
    station_lines{i} = sprintf(' %.8f %.8f ! %s', ...
        SelectedStations.station_lon(i), SelectedStations.station_lat(i), label);
end

remove_first = idx_count + 1;
remove_last = min(idx_count + old_count, numel(lines));
if remove_last >= remove_first
    lines(remove_first:remove_last) = [];
end
lines = [lines(1:idx_count); station_lines; lines(idx_count+1:end)];
end

function label = sanitize_station_comment(label)
label = char(string(label));
label = regexprep(label, '[\r\n]', ' ');
label = regexprep(label, '\s+', ' ');
label = strtrim(label);
if numel(label) > 120
    label = label(1:120);
end
end

function lines = zero_station_block(lines, schedule_token, count_token)
idx_sched = find_line_contains(lines, schedule_token);
if isempty(idx_sched)
    return;
end

idx_count = [];
for i = idx_sched+1:min(idx_sched+5, numel(lines))
    if contains(lines{i}, count_token, 'IgnoreCase', true)
        idx_count = i;
        break;
    end
end
if isempty(idx_count)
    idx_count = idx_sched + 1;
end
if idx_count > numel(lines)
    return;
end

old_count = first_integer_from_line(lines{idx_count});
if ~isfinite(old_count) || old_count < 0
    old_count = 0;
end

lines{idx_sched} = set_first_numeric_token(lines{idx_sched}, 0);
lines{idx_count} = set_count_line_to_zero(lines{idx_count});

remove_first = idx_count + 1;
remove_last = min(idx_count + old_count, numel(lines));
if remove_last >= remove_first
    lines(remove_first:remove_last) = [];
end
end

function n = first_integer_from_line(line)
vals = sscanf(strtrim(line), '%f');
if isempty(vals)
    n = NaN;
else
    n = round(vals(1));
end
end

function line = set_count_line_to_zero(line)
excl = strfind(line, '!');
if isempty(excl)
    comment_part = '';
else
    comment_part = strtrim(line(excl(1):end));
end

if isempty(comment_part)
    line = ' 0';
else
    line = sprintf(' 0                                   %s', comment_part);
end
end

function new_line = set_first_numeric_token(line, first_value)
excl = strfind(line, '!');
if isempty(excl)
    data_part = strtrim(line);
    comment_part = '';
else
    data_part = strtrim(line(1:excl(1)-1));
    comment_part = strtrim(line(excl(1):end));
end

nums = sscanf(data_part, '%f');
if isempty(nums)
    new_line = line;
    return;
end
nums(1) = first_value;

data_new = sprintf(' %g', nums);
if isempty(comment_part)
    new_line = data_new;
else
    new_line = sprintf('%s                        %s', data_new, comment_part);
end
end

%% ------------------------- fort.19 -------------------------

function generate_fort19_from_tmd(G, meta_file, fort19_file, tmd_root, latlon_file, ...
    tmd_out_dir, clean_tmd_before_each_case, allow_trim_or_pad)

meta = read_fort22_meta(meta_file);

WindowStart = get_required_datetime(meta, 'WindowStart');
WindowEnd = get_required_datetime(meta, 'WindowEnd');
WTIMINC = get_required_numeric(meta, 'WTIMINC');

dt_min = round(WTIMINC / 60);
if abs(WTIMINC / 60 - dt_min) > 1e-8
    warning('WTIMINC is not an integer minute; rounded dt_min to %d.', dt_min);
end

nSteps = round(seconds(WindowEnd - WindowStart) / WTIMINC) + 1;
if nSteps <= 0
    error('Invalid time window for fort.19.');
end

write_tmd_latlon_file(latlon_file, G.open_lat_ordered, G.open_lon_ordered, ...
    WindowStart, dt_min, nSteps);

if clean_tmd_before_each_case
    cleanup_tmd_outputs(tmd_root, tmd_out_dir);
end
ensure_blank_data_out(fullfile(tmd_root, 'data.out'));
run_tmd_ato(tmd_root);

tide = read_tmd_series_matrix(tmd_root, numel(G.open_nodes_ordered), nSteps, allow_trim_or_pad);
write_fort19_file(fort19_file, tide, WTIMINC);
end

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

function cleanup_tmd_outputs(tmd_root, tmd_out_dir)
delete_if_exists(fullfile(tmd_root, 'data_*.mat'));
delete_if_exists(fullfile(tmd_root, 'data.out'));
delete_if_exists(fullfile(tmd_root, 'data*'));

if exist(tmd_out_dir, 'dir')
    delete_if_exists(fullfile(tmd_out_dir, 'data.out'));
    delete_if_exists(fullfile(tmd_out_dir, 'data*'));
end
end

function delete_if_exists(pattern)
D = dir(pattern);
for i = 1:numel(D)
    try
        delete(fullfile(D(i).folder, D(i).name));
    catch
    end
end
end

function ensure_blank_data_out(data_out_file)
fid = fopen(data_out_file, 'wt');
if fid < 0
    error('Cannot write %s', data_out_file);
end
fprintf(fid, '%s\n', ' ');
fclose(fid);
end

function run_tmd_ato(tmd_root)
old_dir = pwd;
cleanupObj = onCleanup(@() cd(old_dir)); %#ok<NASGU>
cd(tmd_root);

if exist('TMD_ato', 'file') ~= 2 && exist('TMD_ato', 'builtin') ~= 5
    error('TMD_ato is not on the MATLAB path.');
end

TMD_ato;
end

function tide = read_tmd_series_matrix(tmd_root, nOpen, nSteps, allow_trim_or_pad)
tide = nan(nOpen, nSteps);

for i = 1:nOpen
    f = fullfile(tmd_root, sprintf('data_%d.mat', i));
    if ~exist(f, 'file')
        error('Missing TMD output: %s', f);
    end

    S = load(f);
    if ~isfield(S, 'TimeSeries')
        error('%s does not contain TimeSeries.', f);
    end

    ts = real(S.TimeSeries(:));
    if numel(ts) == nSteps
        tide(i,:) = ts.';
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
        tide(i,:) = tmp.';
    else
        error('TMD length mismatch for data_%d.mat: got %d, need %d.', i, numel(ts), nSteps);
    end
end
end

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
        fprintf(fid, '%.10f\n', tide(j,k));
    end
end
end

%% ------------------------- fort.13 -------------------------

function make_fort13_for_mesh(G, out_fort13, Opt)
open_nodes_unique = unique(G.open_nodes_ordered(:));
if isempty(open_nodes_unique)
    error('No open-boundary nodes in fort.14.');
end

lon0 = mean(G.lon, 'omitnan');
lat0 = mean(G.lat, 'omitnan');

[x_all, y_all] = ll2xy_local(G.lon, G.lat, lon0, lat0);
[x_ob, y_ob] = ll2xy_local(G.lon(open_nodes_unique), G.lat(open_nodes_unique), lon0, lat0);

dist_to_open_m = min_distance_to_points_chunked(x_all, y_all, x_ob, y_ob, Opt.chunk_size);
dist_to_open_km = dist_to_open_m / 1000;

depth = G.depth(:);
depth_use = max(depth, Opt.depth_floor);
cf_base = Opt.cf_default * ones(G.NP, 1);

if Opt.use_shallow_manning
    I_shallow = depth > 0 & depth <= Opt.shallow_depth_max;
    cf_manning = 9.81 * Opt.manning_n^2 ./ (depth_use .^ (1/3));
    cf_manning = min(cf_manning, Opt.cf_cap);
    cf_base(I_shallow) = max(cf_base(I_shallow), cf_manning(I_shallow));
end

boost = ones(G.NP, 1);
I1 = dist_to_open_km <= Opt.inner_km;
I2 = dist_to_open_km > Opt.inner_km & dist_to_open_km < Opt.outer_km;

boost(I1) = Opt.max_boost;
boost(I2) = 1 + (Opt.max_boost - 1) * ...
    (Opt.outer_km - dist_to_open_km(I2)) / (Opt.outer_km - Opt.inner_km);

cf_final = cf_base .* boost;
cf_final = max(cf_final, Opt.cf_default);
cf_final = min(cf_final, Opt.cf_cap);

if Opt.use_near_land_open_boundary_boost
    cf_final = apply_near_land_open_boundary_boost(G, cf_final, Opt);
end

local_mask = false(G.NP, 1);
if Opt.use_cf_uprange && exist(Opt.cf_uprange_shp, 'file')
    try
        local_mask = points_in_shapefile(G.lon, G.lat, Opt.cf_uprange_shp);
        if any(local_mask)
            cf_final(local_mask) = cf_final(local_mask) .* Opt.local_boost;
            cf_final(local_mask) = min(cf_final(local_mask), Opt.local_cf_cap);
        end
    catch ME
        warning('Skipped CF_uprange local boost: %s', ME.message);
    end
end

if Opt.write_only_nondefault
    idx_write = find(abs(cf_final - Opt.cf_default) > 1e-12);
else
    idx_write = (1:G.NP).';
end

write_fort13_file(out_fort13, G.NP, Opt.slope_limiter_default, ...
    Opt.cf_default, idx_write, cf_final(idx_write));
end

function cf_final = apply_near_land_open_boundary_boost(G, cf_final, Opt)
if ~isfield(G, 'land_nodes_ordered') || isempty(G.land_nodes_ordered)
    warning('Near-land open-boundary boost skipped: no land-boundary nodes found in fort.14.');
    return;
end

open_nodes_unique = unique(G.open_nodes_ordered(:));
land_nodes_unique = unique(G.land_nodes_ordered(:));

open_lon = G.lon(open_nodes_unique);
open_lat = G.lat(open_nodes_unique);
land_lon = G.lon(land_nodes_unique);
land_lat = G.lat(land_nodes_unique);

dist_open_to_land_deg = min_angular_distance_deg_chunked( ...
    open_lat, open_lon, land_lat, land_lon, Opt.chunk_size);

near_open_nodes = open_nodes_unique(dist_open_to_land_deg <= Opt.open_land_distance_threshold_deg);
fprintf('  near-land open-boundary nodes: %d / %d (threshold %.2f deg)\n', ...
    numel(near_open_nodes), numel(open_nodes_unique), Opt.open_land_distance_threshold_deg);

if isempty(near_open_nodes)
    return;
end

lon0 = mean(G.lon, 'omitnan');
lat0 = mean(G.lat, 'omitnan');
[x_all, y_all] = ll2xy_local(G.lon, G.lat, lon0, lat0);
[x_near_ob, y_near_ob] = ll2xy_local(G.lon(near_open_nodes), G.lat(near_open_nodes), lon0, lat0);

dist_to_near_open_km = min_distance_to_points_chunked( ...
    x_all, y_all, x_near_ob, y_near_ob, Opt.chunk_size) / 1000;

target_cf = Opt.near_land_target_cf;
I_inner = dist_to_near_open_km <= Opt.near_land_inner_km;
I_outer = dist_to_near_open_km > Opt.near_land_inner_km & ...
    dist_to_near_open_km < Opt.near_land_outer_km;

cf_near = cf_final;
cf_near(I_inner) = max(cf_near(I_inner), target_cf);

if any(I_outer)
    w = (Opt.near_land_outer_km - dist_to_near_open_km(I_outer)) ./ ...
        (Opt.near_land_outer_km - Opt.near_land_inner_km);
    target_outer = Opt.cf_default + w .* (target_cf - Opt.cf_default);
    cf_near(I_outer) = max(cf_near(I_outer), target_outer);
end

affected = I_inner | I_outer;
cf_final(affected) = min(cf_near(affected), Opt.near_land_cf_cap);
fprintf('  near-land open-boundary boost affected nodes: %d, target Cf %.4f, cap %.4f\n', ...
    nnz(affected), Opt.near_land_target_cf, Opt.near_land_cf_cap);
end

function write_fort13_file(fname, np, slope_limiter_default, cf_default, idx_write, cf_vals)
fid = fopen(fname, 'wt');
if fid < 0
    error('Cannot write fort.13: %s', fname);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'fort13\n');
fprintf(fid, '%d\n', np);
fprintf(fid, '2\n');

fprintf(fid, 'elemental_slope_limiter\n');
fprintf(fid, '1\n');
fprintf(fid, '1\n');
fprintf(fid, '%.8f\n', slope_limiter_default);

fprintf(fid, 'quadratic_friction_coefficient_at_sea_floor\n');
fprintf(fid, 'm\n');
fprintf(fid, '1\n');
fprintf(fid, '%.8f\n', cf_default);

fprintf(fid, 'elemental_slope_limiter\n');
fprintf(fid, '0\n');

fprintf(fid, 'quadratic_friction_coefficient_at_sea_floor\n');
fprintf(fid, '%d\n', numel(idx_write));

for i = 1:numel(idx_write)
    fprintf(fid, '%d %.8f\n', idx_write(i), cf_vals(i));
end
end

%% ------------------------- mesh parsing -------------------------

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
    node_id(i) = v(1);
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
    elem_id(i) = v(1);
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

s = strtrim(fgetl(fid));
G.NETA = sscanf(s, '%d', 1);

open_nodes = [];
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
    open_nodes = [open_nodes; tmp]; %#ok<AGROW>
end

land_nodes = [];
land_header = fgetl(fid);
if ischar(land_header)
    hv = sscanf(strtrim(land_header), '%d');
    if ~isempty(hv)
        G.NBOU = hv(1);

        nvel_line = fgetl(fid);
        if ischar(nvel_line)
            vv = sscanf(strtrim(nvel_line), '%d');
            if ~isempty(vv)
                G.NVEL = vv(1);
            else
                G.NVEL = 0;
            end
        else
            G.NVEL = 0;
        end

        for ib = 1:G.NBOU
            header = strtrim(fgetl(fid));
            hv = sscanf(header, '%d');
            if isempty(hv)
                error('Cannot parse land-boundary header %d.', ib);
            end
            nvell = hv(1);

            tmp = zeros(nvell, 1);
            for k = 1:nvell
                line = strtrim(fgetl(fid));
                vv = sscanf(line, '%d');
                if isempty(vv)
                    error('Cannot parse land-boundary node %d in boundary %d.', k, ib);
                end
                tmp(k) = vv(1);
            end
            land_nodes = [land_nodes; tmp]; %#ok<AGROW>
        end
    else
        G.NBOU = 0;
        G.NVEL = 0;
    end
else
    G.NBOU = 0;
    G.NVEL = 0;
end

G.node_id = node_id;
G.lon = lon;
G.lat = lat;
G.depth = depth;
G.elem_id = elem_id;
G.elem_node_idx = elem_node_idx;
G.open_nodes_ordered = open_nodes(:);
G.open_lon_ordered = lon(G.open_nodes_ordered);
G.open_lat_ordered = lat(G.open_nodes_ordered);
G.land_nodes_ordered = land_nodes(:);
if isempty(G.land_nodes_ordered)
    G.land_lon_ordered = [];
    G.land_lat_ordered = [];
else
    G.land_lon_ordered = lon(G.land_nodes_ordered);
    G.land_lat_ordered = lat(G.land_nodes_ordered);
end
end

function [x, y] = ll2xy_local(lon, lat, lon0, lat0)
x = (lon - lon0) .* cosd(lat0) * 111320;
y = (lat - lat0) * 110540;
end

function dmin = min_distance_to_points_chunked(x, y, xb, yb, chunk_size)
n = numel(x);
dmin = inf(n, 1);
if isempty(xb)
    return;
end

for i1 = 1:chunk_size:n
    i2 = min(i1 + chunk_size - 1, n);
    xx = x(i1:i2);
    yy = y(i1:i2);

    dx = xx - xb.';
    dy = yy - yb.';
    d2 = dx.^2 + dy.^2;
    dmin(i1:i2) = sqrt(min(d2, [], 2));
end
end

function dmin = min_angular_distance_deg_chunked(lat, lon, lat_ref, lon_ref, chunk_size)
n = numel(lat);
dmin = inf(n, 1);
if isempty(lat_ref)
    return;
end

lat_ref_rad = deg2rad(lat_ref(:)).';
lon_ref_rad = deg2rad(lon_ref(:)).';

for i1 = 1:chunk_size:n
    i2 = min(i1 + chunk_size - 1, n);

    lat_rad = deg2rad(lat(i1:i2));
    lon_rad = deg2rad(lon(i1:i2));

    dlat = lat_rad - lat_ref_rad;
    dlon = lon_rad - lon_ref_rad;
    a = sin(dlat ./ 2).^2 + cos(lat_rad) .* cos(lat_ref_rad) .* sin(dlon ./ 2).^2;
    c = 2 .* atan2(sqrt(a), sqrt(max(1 - a, 0)));

    dmin(i1:i2) = rad2deg(min(c, [], 2));
end
end

function mask = points_in_shapefile(lon, lat, shpfile)
if exist('shaperead', 'file') ~= 2
    error('shaperead is unavailable.');
end

S = shaperead(shpfile);
mask = false(size(lon));

for i = 1:numel(S)
    x = S(i).X(:);
    y = S(i).Y(:);
    if isempty(x) || isempty(y)
        continue;
    end

    nan_idx = find(isnan(x) | isnan(y));
    seg_start = 1;
    split_idx = [nan_idx; numel(x)+1];

    for k = 1:numel(split_idx)
        seg_end = split_idx(k) - 1;
        if seg_end >= seg_start
            xv = x(seg_start:seg_end);
            yv = y(seg_start:seg_end);
            good = isfinite(xv) & isfinite(yv);
            xv = xv(good);
            yv = yv(good);
            if numel(xv) >= 3
                mask = mask | inpolygon(lon, lat, xv, yv);
            end
        end
        seg_start = split_idx(k) + 1;
    end
end
end

%% ------------------------- shared text/meta utilities -------------------------

function lines = read_text_lines(fname)
fid = fopen(fname, 'rt');
if fid < 0
    error('Cannot open file: %s', fname);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

lines = {};
while true
    tline = fgetl(fid);
    if ~ischar(tline)
        break;
    end
    lines{end+1,1} = tline; %#ok<AGROW>
end
end

function write_text_lines(fname, lines)
fid = fopen(fname, 'wt');
if fid < 0
    error('Cannot write file: %s', fname);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

for i = 1:numel(lines)
    fprintf(fid, '%s\n', lines{i});
end
end

function meta = read_fort22_meta(meta_file)
lines = read_text_lines(meta_file);
meta = struct();

for i = 1:numel(lines)
    s = strtrim(lines{i});
    if isempty(s)
        continue;
    end
    if startsWith(s, 'ADCIRC', 'IgnoreCase', true) || startsWith(s, '---') || startsWith(s, '[')
        continue;
    end

    eq_pos = strfind(s, '=');
    if isempty(eq_pos)
        continue;
    end

    k = strtrim(s(1:eq_pos(1)-1));
    v = strtrim(s(eq_pos(1)+1:end));
    k = matlab.lang.makeValidName(k);
    meta.(k) = v;
end
end

function val = get_required_numeric(meta, key)
k = matlab.lang.makeValidName(key);
if ~isfield(meta, k)
    error('meta is missing numeric field: %s', key);
end
val = str2double(strtrim(meta.(k)));
if ~isfinite(val)
    error('meta field %s is not numeric: %s', key, meta.(k));
end
end

function dt = get_required_datetime(meta, key)
k = matlab.lang.makeValidName(key);
if ~isfield(meta, k)
    error('meta is missing datetime field: %s', key);
end

s = strtrim(meta.(k));
s = strrep(s, 'T', ' ');

try
    dt = datetime(s, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
catch
    try
        dt = datetime(s);
    catch
        error('Cannot parse meta datetime %s: %s', key, s);
    end
end
end

function s = get_meta_string(meta, key, default_value)
k = matlab.lang.makeValidName(key);
if isfield(meta, k)
    s = strtrim(meta.(k));
else
    s = default_value;
end
end

function idx = find_line_contains(lines, token)
idx = [];
for i = 1:numel(lines)
    if contains(lines{i}, token)
        idx = i;
        return;
    end
end
end

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

function s = crop_or_pad(s, n)
s = char(s);
if numel(s) > n
    s = s(1:n);
elseif numel(s) < n
    s = [s repmat(' ', 1, n-numel(s))];
end
end

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

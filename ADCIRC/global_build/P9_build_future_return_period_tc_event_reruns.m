%% P9_build_future_return_period_tc_event_reruns.m
% Build future ADCIRC cases using one intensity ratio per unique TC.
% The control block is the exact P4 ADCIRC-inner block where that TC has its
% largest inner-domain impact intensity. Every rerun block for the same TC
% and scenario uses the control block's single intensity ratio.
% Sea-level rise is different: each TC-block-scenario row uses the AR6 median
% SLR inside the intersection of that block's exact P4 inner polygon and the
% 200 km corridor around the full TC track.
% Only fort.19 and fort.22 are changed. fort.14 is hard-linked by default;
% all other original input/helper files are copied without modification.

clearvars;
clc;

SCRIPT_DIR = fileparts(mfilename('fullpath'));
if isempty(SCRIPT_DIR)
    SCRIPT_DIR = pwd;
end

%% =========================== USER SETTINGS =============================
P = struct();
P.source_root = fullfile(SCRIPT_DIR, 'runfile', 'return_period_tc_event_reruns');
P.output_root = P.source_root;
P.return_period_folders = ["200yr";"500yr"];
P.scenarios = ["ssp126"; "ssp245"; "ssp370"];

P.exact_inner_forcing_dir = fullfile(SCRIPT_DIR, 'exact_inner_tc_return_period_forcing');
P.curve_csv = fullfile(P.exact_inner_forcing_dir, '07_exact_inner_block_return_period_curves.csv');
P.block_tc_scaling_csv = fullfile(P.exact_inner_forcing_dir, ...
    '08b_block_tc_control_intensity_scaling.csv');
P.block_tc_slr_csv = fullfile(P.exact_inner_forcing_dir, ...
    '08c_block_tc_local_slr_2061_2100.csv');
P.slr_quantile = 0.5;
P.future_years = (2061:2100).';
P.slr_search_radius_km = 200;
P.slr_spatial_method = "tc_track_corridor_intersect_exact_p4_inner";

P.overwrite = false;
P.hardlink_fort14 = true;
P.selected_blocks = strings(0, 1);
P.max_cases = inf;
P.progress_interval = 10;
P.summary_flush_interval = 25;
P.c15_cache_size = 768;
P.use_parallel = true;
P.num_workers = 8;
P.use_environment_overrides = true;
P.worker_id = "";

P.unchanged_files = {'fort.13', 'fort.14', 'fort.15', ...
    'storm_track_forcing_meta.mat', 'sub_intel.sh'};
P.fort22_meta_file = 'fort22_meta.txt';
P.fort22_meta_version = 'future_fort22_meta_v1';
P.future_meta_mat = 'future_forcing_meta.mat';
P.future_meta_txt = 'future_forcing_meta.txt';
P.completion_marker = 'future_case_complete_unique_tc_intensity_block_tc_slr_v4.ok';
P.summary_csv = fullfile(P.output_root, ...
    'future_return_period_tc_case_build_summary_unique_tc_intensity_block_slr.csv');
P.slr_summary_csv = fullfile(P.output_root, 'future_block_tc_local_slr_2061_2100.csv');

% Optional environment controls for tests/resume:
% ADCIRC_FUTURE_SCENARIOS=ssp126,ssp245
% ADCIRC_FUTURE_RETURN_PERIODS=100yr
% ADCIRC_FUTURE_BLOCKS=ADC_WNP_08
% ADCIRC_FUTURE_MAX_CASES=1
% ADCIRC_FUTURE_OVERWRITE=1
% ADCIRC_FUTURE_OUTPUT_ROOT=work/test_output
% ADCIRC_FUTURE_WORKER_ID=ssp126
% ADCIRC_FUTURE_USE_PARALLEL=1
% ADCIRC_FUTURE_NUM_WORKERS=4
% ADCIRC_FUTURE_EXACT_INNER_FORCING_DIR=external/adcirc/test_forcing
% ADCIRC_FUTURE_SOURCE_ROOT=work/return_period_tc_event_reruns
if P.use_environment_overrides
    P = apply_environment_overrides(P);
end

%% ========================= END USER SETTINGS ============================
validate_inputs(P);
ensure_dir(P.output_root);

fprintf('\n============================================================\n');
fprintf('Build future ADCIRC cases with one max-impact control block per unique TC\n');
fprintf('Source     : %s\n', P.source_root);
fprintf('Output     : %s\n', P.output_root);
fprintf('Return sets: %s\n', strjoin(cellstr(P.return_period_folders), ', '));
fprintf('Scenarios  : %s\n', strjoin(cellstr(P.scenarios), ', '));
fprintf('Future mean: %d-%d, AR6 median, no VLM\n', P.future_years(1), P.future_years(end));
fprintf('============================================================\n\n');

validate_curve_table(read_csv_table(P.curve_csv));
TcScaling = read_csv_table(P.block_tc_scaling_csv);
validate_block_tc_scaling_table(TcScaling, P.scenarios);

fprintf('[1/3] Discovering source cases...\n');
Tasks = discover_tasks(P);
fprintf('      Found %d source block-cases before scenario expansion.\n', height(Tasks));
fprintf('      Planned future cases: %d\n', height(Tasks) * numel(P.scenarios));
validate_task_scaling_coverage(TcScaling, Tasks, P.scenarios);
fprintf('      Block-TC scaling coverage verified: %d selected TCs.\n', ...
    numel(unique(double(Tasks.track_index))));

fprintf('[2/3] Loading P3 block-specific TC local sea-level rise...\n');
TcSlrTable = read_csv_table(P.block_tc_slr_csv);
validate_tc_slr_table(TcSlrTable, Tasks, P);
writetable(TcSlrTable, P.slr_summary_csv);
fprintf('      Unique TC=%d, block-TC-scenario rows=%d, nearest fallbacks=%d\n', ...
    numel(unique(double(TcSlrTable.track_index))), height(TcSlrTable), ...
    nnz(TcSlrTable.nearest_fallback));
fprintf('      Source: %s\n', P.block_tc_slr_csv);
fprintf('      Build copy: %s\n', P.slr_summary_csv);

summaryNames = {'return_period_set', 'scenario', 'block_id', 'control_block_id', ...
    'case_name', 'source_case_dir', 'future_case_dir', 'status', 'message', ...
    'era5_impact_vmax_ms', 'era5_intensity_return_period_yr', ...
    'cmip6_historical_intensity_ms', 'cmip6_future_intensity_ms', ...
    'intensity_scale_factor', 'intensity_change_percent', 'tc_local_slr_m', ...
    'slr_block_event_lon', 'slr_block_event_lat', 'slr_location_count', ...
    'slr_search_radius_used_km', 'slr_nearest_fallback', ...
    'fort22_max_wind_ms', 'fort22_active_steps'};
fprintf('[3/3] Building future cases...\n');
[JobTaskIndex, JobScenarios] = expand_future_jobs(Tasks, P);
nAttempt = numel(JobTaskIndex);
SummaryRows = cell(nAttempt, numel(summaryNames));
JobReturnPeriods = string(Tasks.return_period_set(JobTaskIndex));

if P.use_parallel && nAttempt > 0
    pool = start_fixed_process_pool(P.num_workers);
    fprintf('      Parallel workers: %d\n', pool.NumWorkers);
    progressQueue = parallel.pool.DataQueue;
    report_parallel_progress("reset", nAttempt, P.progress_interval);
    afterEach(progressQueue, @(status) report_parallel_progress( ...
        status, nAttempt, P.progress_interval));
else
    serialAttempted = 0;
    serialBuilt = 0;
    serialSkipped = 0;
    serialFailed = 0;
end

% Use one parfor batch per return period. The end of each parfor is a hard
% barrier, so no 200yr job starts before every 100yr job has returned, and
% no 500yr job starts before every 200yr job has returned.
for ir = 1:numel(P.return_period_folders)
    rp = P.return_period_folders(ir);
    batchPositions = find(JobReturnPeriods == rp);
    nBatch = numel(batchPositions);
    if nBatch == 0
        continue;
    end
    BatchTasks = Tasks(JobTaskIndex(batchPositions), :);
    BatchScenarios = JobScenarios(batchPositions);
    fprintf('      Starting %s batch: %d cases across %d scenarios\n', ...
        rp, nBatch, numel(P.scenarios));

    if P.use_parallel
        BatchRows = cell(nBatch, numel(summaryNames));
        parfor ib = 1:nBatch
            task = BatchTasks(ib, :);
            row = build_future_job( ...
                task, BatchScenarios(ib), TcScaling, TcSlrTable, P);
            BatchRows(ib, :) = row;
            send(progressQueue, string(row{8}));
        end
        SummaryRows(batchPositions, :) = BatchRows;
    else
        for ib = 1:nBatch
            ij = batchPositions(ib);
            task = BatchTasks(ib, :);
            row = build_future_job( ...
                task, BatchScenarios(ib), TcScaling, TcSlrTable, P);
            SummaryRows(ij, :) = row;
            rowStatus = string(row{8});
            serialAttempted = serialAttempted + 1;
            serialBuilt = serialBuilt + double(rowStatus == "built");
            serialSkipped = serialSkipped + double(rowStatus == "skipped_complete");
            serialFailed = serialFailed + double(rowStatus == "failed");
            if mod(serialAttempted, P.progress_interval) == 0 || ...
                    serialAttempted == nAttempt
                fprintf('      attempted=%d built=%d skipped=%d failed=%d\n', ...
                    serialAttempted, serialBuilt, serialSkipped, serialFailed);
            end
        end
    end

    batchStatus = string(SummaryRows(batchPositions, 8));
    fprintf('      Finished %s: built=%d skipped=%d failed=%d\n', rp, ...
        nnz(batchStatus == "built"), ...
        nnz(batchStatus == "skipped_complete"), ...
        nnz(batchStatus == "failed"));

    completedRows = ~cellfun(@isempty, SummaryRows(:, 8));
    write_summary_rows(SummaryRows(completedRows, :), summaryNames, P.summary_csv);
    fprintf('      Checkpoint summary: %s\n', P.summary_csv);
    nBatchFailed = nnz(batchStatus == "failed");
    if nBatchFailed > 0
        error('%d future cases failed in %s. Stopping before the next return-period batch.', ...
            nBatchFailed, rp);
    end
end

status = string(SummaryRows(:, 8));
nBuilt = nnz(status == "built");
nSkipped = nnz(status == "skipped_complete");
nFailed = nnz(status == "failed");

write_summary_rows(SummaryRows, summaryNames, P.summary_csv);
fprintf('\nDone. attempted=%d built=%d skipped=%d failed=%d\n', ...
    nAttempt, nBuilt, nSkipped, nFailed);
fprintf('Summary: %s\n', P.summary_csv);
if nFailed > 0
    error('%d future cases failed. See the summary CSV.', nFailed);
end

%% ========================================================================
function P = apply_environment_overrides(P)

v = strtrim(string(getenv('ADCIRC_FUTURE_SOURCE_ROOT')));
if strlength(v) > 0
    P.source_root = char(v);
end

v = strtrim(string(getenv('ADCIRC_FUTURE_EXACT_INNER_FORCING_DIR')));
if strlength(v) == 0
    v = strtrim(string(getenv('ADCIRC_FUTURE_SFINCS_FORCING_DIR'))); % backward-compatible test override
end
if strlength(v) > 0
    P.exact_inner_forcing_dir = char(v);
    P.curve_csv = fullfile(P.exact_inner_forcing_dir, '07_exact_inner_block_return_period_curves.csv');
    P.block_tc_scaling_csv = fullfile(P.exact_inner_forcing_dir, ...
        '08b_block_tc_control_intensity_scaling.csv');
    P.block_tc_slr_csv = fullfile(P.exact_inner_forcing_dir, ...
        '08c_block_tc_local_slr_2061_2100.csv');
end

v = strtrim(string(getenv('ADCIRC_FUTURE_SCENARIOS')));
if strlength(v) > 0
    P.scenarios = split_csv_env(v);
end
v = strtrim(string(getenv('ADCIRC_FUTURE_RETURN_PERIODS')));
if strlength(v) > 0
    P.return_period_folders = split_csv_env(v);
end
v = strtrim(string(getenv('ADCIRC_FUTURE_BLOCKS')));
if strlength(v) > 0
    P.selected_blocks = split_csv_env(v);
end
v = strtrim(string(getenv('ADCIRC_FUTURE_MAX_CASES')));
if strlength(v) > 0
    x = str2double(v);
    assert(isfinite(x) && x >= 1, 'ADCIRC_FUTURE_MAX_CASES must be >= 1.');
    P.max_cases = floor(x);
end
v = strtrim(string(getenv('ADCIRC_FUTURE_OVERWRITE')));
if strlength(v) > 0
    P.overwrite = any(strcmpi(v, ["1", "true", "yes", "on"]));
end
v = strtrim(string(getenv('ADCIRC_FUTURE_OUTPUT_ROOT')));
if strlength(v) > 0
    P.output_root = char(v);
    P.summary_csv = fullfile(P.output_root, ...
        'future_return_period_tc_case_build_summary_unique_tc_intensity_block_slr.csv');
    P.slr_summary_csv = fullfile(P.output_root, 'future_block_tc_local_slr_2061_2100.csv');
end
v = strtrim(string(getenv('ADCIRC_FUTURE_WORKER_ID')));
if strlength(v) > 0
    P.worker_id = regexprep(v, '[^A-Za-z0-9_-]', '_');
    P.summary_csv = fullfile(P.output_root, sprintf( ...
        'future_return_period_tc_case_build_summary_unique_tc_intensity_block_slr_%s.csv', char(P.worker_id)));
    P.slr_summary_csv = fullfile(P.output_root, sprintf( ...
        'future_block_tc_local_slr_2061_2100_%s.csv', char(P.worker_id)));
end
v = strtrim(string(getenv('ADCIRC_FUTURE_USE_PARALLEL')));
if strlength(v) > 0
    P.use_parallel = any(strcmpi(v, ["1", "true", "yes", "on"]));
end
v = strtrim(string(getenv('ADCIRC_FUTURE_NUM_WORKERS')));
if strlength(v) > 0
    x = str2double(v);
    assert(isfinite(x) && x >= 1 && x == floor(x), ...
        'ADCIRC_FUTURE_NUM_WORKERS must be a positive integer.');
    P.num_workers = x;
end
end

%% ========================================================================
function values = split_csv_env(v)

values = strip(split(v, ','));
values = values(strlength(values) > 0);
values = values(:);
end

%% ========================================================================
function write_summary_rows(rows, names, path)

T = cell2table(rows, 'VariableNames', names);
writetable(T, path);
end

%% ========================================================================
function [taskIndex, scenarios] = expand_future_jobs(Tasks, P)

nTasks = height(Tasks);
nScenarios = numel(P.scenarios);
taskIndex = zeros(nTasks * nScenarios, 1);
scenarios = strings(nTasks * nScenarios, 1);
cursor = 0;
for ir = 1:numel(P.return_period_folders)
    rpTaskIndex = find(string(Tasks.return_period_set) == ...
        P.return_period_folders(ir));
    for is = 1:nScenarios
        rows = cursor + (1:numel(rpTaskIndex));
        taskIndex(rows) = rpTaskIndex;
        scenarios(rows) = P.scenarios(is);
        cursor = cursor + numel(rpTaskIndex);
    end
end
assert(cursor == nTasks * nScenarios, ...
    'Return-period batching did not cover every discovered source task.');
nKeep = min(numel(taskIndex), P.max_cases);
taskIndex = taskIndex(1:nKeep);
scenarios = scenarios(1:nKeep);

keys = string(Tasks.return_period_set(taskIndex)) + "|" + scenarios + "|" + ...
    string(Tasks.block_id(taskIndex)) + "|" + string(Tasks.case_name(taskIndex));
assert(numel(unique(keys)) == numel(keys), ...
    'Future task expansion produced duplicate output case paths.');
end

%% ========================================================================
function pool = start_fixed_process_pool(numWorkers)

assert(exist('gcp', 'file') == 2 && exist('parpool', 'file') == 2, ...
    ['Parallel Computing Toolbox is not installed in this MATLAB. Install it ' ...
    'to use four workers, or set P.use_parallel=false for serial execution.']);
pool = gcp('nocreate');
if ~isempty(pool) && pool.NumWorkers ~= numWorkers
    fprintf('      Replacing existing %d-worker pool with %d workers...\n', ...
        pool.NumWorkers, numWorkers);
    delete(pool);
    pool = [];
end
if isempty(pool)
    pool = parpool('local', numWorkers);
end
assert(pool.NumWorkers == numWorkers, ...
    'Expected %d parallel workers, but the pool has %d.', ...
    numWorkers, pool.NumWorkers);
end

%% ========================================================================
function report_parallel_progress(status, total, interval)

persistent attempted built skipped failed
if status == "reset"
    attempted = 0;
    built = 0;
    skipped = 0;
    failed = 0;
    return;
end
attempted = attempted + 1;
built = built + double(status == "built");
skipped = skipped + double(status == "skipped_complete");
failed = failed + double(status == "failed");
if mod(attempted, interval) == 0 || attempted == total
    fprintf('      attempted=%d built=%d skipped=%d failed=%d\n', ...
        attempted, built, skipped, failed);
end
end

%% ========================================================================
function row = build_future_job(task, scenario, TcScaling, TcSlrTable, P)

sourceCaseDir = char(task.source_case_dir);
futureSet = task.return_period_set + "_" + scenario;
futureCaseDir = fullfile(P.output_root, char(futureSet), ...
    char(task.block_id), char(task.case_name));
row = make_summary_row(task, scenario, futureCaseDir);

try
    scaleInfo = lookup_tc_scaling( ...
        TcScaling, task.track_index, task.block_id, scenario);
    slrInfo = lookup_tc_slr( ...
        TcSlrTable, task.track_index, task.block_id, scenario);
    if is_future_case_complete(futureCaseDir, P) && ~P.overwrite
        [futureMeta, metaUpdated] = ensure_future_fort22_meta( ...
            sourceCaseDir, futureCaseDir, task, scaleInfo, scenario, slrInfo, P);
        row{4} = char(scaleInfo.control_block_id);
        row{8} = "skipped_complete";
        if metaUpdated
            row{9} = "Existing forcing complete; future fort22_meta.txt refreshed";
        else
            row{9} = "Existing complete future case and current fort22 metadata";
        end
        row{10} = double(scaleInfo.era5_impact_vmax_ms);
        row{11} = double(scaleInfo.era5_intensity_return_period_yr);
        row{12} = double(scaleInfo.cmip6_historical_intensity_ms);
        row{13} = double(scaleInfo.cmip6_future_intensity_ms);
        row{14} = double(scaleInfo.intensity_scale_factor);
        row{15} = 100 * (double(scaleInfo.intensity_scale_factor) - 1);
        row{16} = double(slrInfo.tc_local_slr_m);
        row{17} = double(slrInfo.block_event_lon);
        row{18} = double(slrInfo.block_event_lat);
        row{19} = double(slrInfo.slr_location_count);
        row{20} = double(slrInfo.search_radius_used_km);
        row{21} = logical(slrInfo.nearest_fallback);
        row{22} = double(futureMeta.fort22_max_wind_ms);
        row{23} = double(futureMeta.fort22_active_steps);
        return;
    end
    [buildInfo, futureMeta] = build_one_future_case( ...
        sourceCaseDir, futureCaseDir, task, scaleInfo, scenario, slrInfo, P);

    row{4} = char(scaleInfo.control_block_id);
    row{8} = "built";
    row{9} = "";
    row{10} = buildInfo.era5ImpactVmax;
    row{11} = buildInfo.intensityReturnPeriod;
    row{12} = buildInfo.historicalIntensity;
    row{13} = buildInfo.futureIntensity;
    row{14} = buildInfo.scaleFactor;
    row{15} = 100 * (buildInfo.scaleFactor - 1);
    row{16} = slrInfo.tc_local_slr_m;
    row{17} = slrInfo.block_event_lon;
    row{18} = slrInfo.block_event_lat;
    row{19} = slrInfo.slr_location_count;
    row{20} = slrInfo.search_radius_used_km;
    row{21} = slrInfo.nearest_fallback;
    row{22} = futureMeta.fort22_max_wind_ms;
    row{23} = futureMeta.fort22_active_steps;
catch ME
    row{8} = "failed";
    row{9} = string(getReport(ME, 'basic', 'hyperlinks', 'off'));
    warning('Future case failed: %s | %s', futureCaseDir, ME.message);
end
end

%% ========================================================================
function validate_inputs(P)

assert(exist(P.source_root, 'dir') == 7, 'Missing source root: %s', P.source_root);
assert(exist(P.curve_csv, 'file') == 2, 'Missing curve CSV: %s', P.curve_csv);
assert(exist(P.block_tc_scaling_csv, 'file') == 2, ...
    ['Missing block-TC control-block scaling CSV: %s. Run ' ...
    'plot_future_return_period_forcing_changes_by_sfincs.m first.'], ...
    P.block_tc_scaling_csv);
assert(exist(P.block_tc_slr_csv, 'file') == 2, ...
    ['Missing P3 block-TC SLR CSV: %s. Run ' ...
    'P3_plot_future_return_period_forcing_changes_by_sfincs.m first.'], ...
    P.block_tc_slr_csv);
validScenarios = ["ssp126", "ssp245", "ssp370"];
assert(all(ismember(P.scenarios, validScenarios)), ...
    'Only ssp126, ssp245, and ssp370 are supported.');
for i = 1:numel(P.return_period_folders)
    d = fullfile(P.source_root, char(P.return_period_folders(i)));
    assert(exist(d, 'dir') == 7, 'Missing source return-period folder: %s', d);
end
end

%% ========================================================================
function validate_curve_table(T)

need = {'region_id', 'experiment', 'return_period_yr', 'intensity_ms'};
missing = setdiff(need, T.Properties.VariableNames);
assert(isempty(missing), 'Curve CSV missing columns: %s', strjoin(missing, ', '));
end

%% ========================================================================
function validate_block_tc_scaling_table(T, scenarios)

need = {'scenario', 'track_index', 'block_id', 'control_block_id', 'control_block_label', ...
    'era5_impact_vmax_ms', 'event_lon', 'event_lat', ...
    'era5_intensity_return_period_yr', 'cmip6_historical_intensity_ms', ...
    'cmip6_future_intensity_ms', 'intensity_scale_factor', 'status'};
missing = setdiff(need, T.Properties.VariableNames);
assert(isempty(missing), 'Block-TC scaling CSV missing columns: %s', strjoin(missing, ', '));
key = string(T.track_index) + "|" + string(T.block_id) + "|" + string(T.scenario);
assert(numel(unique(key)) == height(T), ...
    'Block-TC scaling CSV contains duplicate track_index + block_id + scenario rows.');
assert(all(ismember(scenarios, unique(string(T.scenario)))), ...
    'Block-TC scaling CSV does not contain every requested scenario.');

% The table is expanded by block for a common schema with SLR, but the TC
% intensity change must remain identical across all blocks for one TC/SSP.
tcScenarioKey = string(T.track_index) + "|" + string(T.scenario);
[groupIndex, groupKey] = findgroups(tcScenarioKey);
factorRange = splitapply(@(x) max(x) - min(x), ...
    double(T.intensity_scale_factor), groupIndex);
assert(all(isfinite(factorRange)) && all(abs(factorRange) <= 1e-12), ...
    'Block-TC scaling CSV has inconsistent block factors for %s.', ...
    strjoin(cellstr(groupKey(abs(factorRange) > 1e-12)), ', '));
for ig = 1:numel(groupKey)
    I = groupIndex == ig;
    assert(numel(unique(string(T.control_block_id(I)))) == 1, ...
        'Block-TC scaling CSV has inconsistent control blocks for %s.', groupKey(ig));
end
end

%% ========================================================================
function validate_tc_slr_table(T, Tasks, P)

need = {'track_index', 'block_id', 'scenario', 'tc_local_slr_m', ...
    'block_event_lon', 'block_event_lat', 'slr_location_count', ...
    'nominal_search_radius_km', 'search_radius_used_km', ...
    'nearest_fallback', 'spatial_method', 'track_point_count', ...
    'year_start', 'year_end', 'quantile'};
missing = setdiff(need, T.Properties.VariableNames);
assert(isempty(missing), 'P3 block-TC SLR CSV missing columns: %s', ...
    strjoin(missing, ', '));
key = string(T.track_index) + "|" + string(T.block_id) + "|" + string(T.scenario);
assert(numel(unique(key)) == height(T), ...
    'P3 block-TC SLR CSV contains duplicate track_index + block_id + scenario rows.');
taskKey = string(Tasks.track_index) + "|" + string(Tasks.block_id);
[~, taskRows] = unique(taskKey, 'stable');
for it = 1:numel(taskRows)
    task = Tasks(taskRows(it), :);
    for is = 1:numel(P.scenarios)
        I = double(T.track_index) == double(task.track_index) & ...
            string(T.block_id) == string(task.block_id) & ...
            string(T.scenario) == string(P.scenarios(is));
        assert(nnz(I) == 1, ...
            'Missing P3 SLR for track %.0f block %s scenario %s.', ...
            task.track_index, task.block_id, P.scenarios(is));
        assert(isfinite(double(T.tc_local_slr_m(I))), ...
            'Invalid P3 SLR for track %.0f block %s scenario %s.', ...
            task.track_index, task.block_id, P.scenarios(is));
        assert(double(T.year_start(I)) == P.future_years(1) && ...
            double(T.year_end(I)) == P.future_years(end) && ...
            abs(double(T.quantile(I)) - P.slr_quantile) < 1e-12 && ...
            abs(double(T.nominal_search_radius_km(I)) - P.slr_search_radius_km) < 1e-12 && ...
            string(T.spatial_method(I)) == P.slr_spatial_method && ...
            double(T.track_point_count(I)) >= 1, ...
            'P3 SLR configuration mismatch for track %.0f block %s scenario %s.', ...
            task.track_index, task.block_id, P.scenarios(is));
    end
end
end

%% ========================================================================
function Tasks = discover_tasks(P)

rpCol = strings(0, 1);
blockCol = strings(0, 1);
caseCol = strings(0, 1);
sourceCol = strings(0, 1);
trackCol = zeros(0, 1);

for ir = 1:numel(P.return_period_folders)
    rp = P.return_period_folders(ir);
    rpDir = fullfile(P.source_root, char(rp));
    blockDirs = dir(rpDir);
    blockDirs = blockDirs([blockDirs.isdir]);
    blockDirs = blockDirs(~ismember({blockDirs.name}, {'.', '..'}));
    [~, order] = sort({blockDirs.name});
    blockDirs = blockDirs(order);
    for ib = 1:numel(blockDirs)
        blockId = string(blockDirs(ib).name);
        if ~isempty(P.selected_blocks) && ~ismember(blockId, P.selected_blocks)
            continue;
        end
        blockDir = fullfile(rpDir, char(blockId));
        caseDirs = dir(blockDir);
        caseDirs = caseDirs([caseDirs.isdir]);
        caseDirs = caseDirs(~ismember({caseDirs.name}, {'.', '..'}));
        [~, caseOrder] = sort({caseDirs.name});
        caseDirs = caseDirs(caseOrder);
        for ic = 1:numel(caseDirs)
            sourceCaseDir = fullfile(blockDir, caseDirs(ic).name);
            if exist(fullfile(sourceCaseDir, 'storm_track_forcing_meta.mat'), 'file') ~= 2
                continue;
            end
            rpCol(end + 1, 1) = rp; %#ok<AGROW>
            blockCol(end + 1, 1) = blockId; %#ok<AGROW>
            caseCol(end + 1, 1) = string(caseDirs(ic).name); %#ok<AGROW>
            sourceCol(end + 1, 1) = string(sourceCaseDir); %#ok<AGROW>
            trackCol(end + 1, 1) = parse_track_index(caseDirs(ic).name); %#ok<AGROW>
        end
    end
end
Tasks = table(rpCol, blockCol, caseCol, sourceCol, trackCol, ...
    'VariableNames', {'return_period_set', 'block_id', 'case_name', ...
    'source_case_dir', 'track_index'});
end

%% ========================================================================
function trackIndex = parse_track_index(caseName)

token = regexp(char(caseName), '^track_(\d+)_', 'tokens', 'once');
assert(~isempty(token), 'Cannot parse track index from case name: %s', caseName);
trackIndex = str2double(token{1});
assert(isfinite(trackIndex), 'Invalid track index in case name: %s', caseName);
end

function validate_task_scaling_coverage(TcScaling, Tasks, scenarios)

taskKey = string(Tasks.track_index) + "|" + string(Tasks.block_id);
[~, taskRows] = unique(taskKey, 'stable');
for it = 1:numel(taskRows)
    task = Tasks(taskRows(it), :);
    for is = 1:numel(scenarios)
        I = double(TcScaling.track_index) == double(task.track_index) & ...
            string(TcScaling.block_id) == string(task.block_id) & ...
            string(TcScaling.scenario) == string(scenarios(is));
        assert(nnz(I) == 1, ...
            ['Expected one control-block scaling row for track %.0f block %s ' ...
            'scenario %s; found %d.'], ...
            task.track_index, task.block_id, scenarios(is), nnz(I));
        assert(string(TcScaling.status(I)) == "ok", ...
            ['Control-block scaling is not usable for track %.0f block %s scenario %s: %s. ' ...
            'Run the full historical/future curve calculation before building cases.'], ...
            task.track_index, task.block_id, scenarios(is), string(TcScaling.status(I)));
        factor = double(TcScaling.intensity_scale_factor(I));
        assert(isfinite(factor) && factor > 0, ...
            ['Invalid control-block intensity scale factor for track %.0f block %s ' ...
            'scenario %s.'], task.track_index, task.block_id, scenarios(is));
    end
end
end

%% ========================================================================
function row = lookup_tc_scaling(TcScaling, trackIndex, blockId, scenario)

I = double(TcScaling.track_index) == double(trackIndex) & ...
    string(TcScaling.block_id) == string(blockId) & ...
    string(TcScaling.scenario) == string(scenario);
assert(nnz(I) == 1, ...
    ['Cannot uniquely locate control-block scaling for track %.0f block %s ' ...
    'scenario %s.'], trackIndex, blockId, scenario);
row = TcScaling(I, :);
assert(string(row.status) == "ok", ...
    'Control-block scaling is unusable for track %.0f scenario %s: %s', ...
    trackIndex, scenario, string(row.status));
assert(strlength(string(row.control_block_id)) > 0, ...
    'Missing max-impact control block for track %.0f.', trackIndex);
assert(isfinite(double(row.intensity_scale_factor)) && double(row.intensity_scale_factor) > 0, ...
    'Invalid unique-TC intensity scale factor for track %.0f.', trackIndex);
end

%% ========================================================================
function row = lookup_tc_slr(TcSlrTable, trackIndex, blockId, scenario)

I = double(TcSlrTable.track_index) == double(trackIndex) & ...
    string(TcSlrTable.block_id) == string(blockId) & ...
    string(TcSlrTable.scenario) == string(scenario);
assert(nnz(I) == 1, ...
    'Cannot uniquely locate TC-local SLR for track %.0f block %s scenario %s.', ...
    trackIndex, blockId, scenario);
row = TcSlrTable(I, :);
end

%% ========================================================================
function row = make_summary_row(task, scenario, futureCaseDir)

row = {char(task.return_period_set), char(scenario), char(task.block_id), '', ...
    char(task.case_name), char(task.source_case_dir), futureCaseDir, '', '', ...
    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, false, NaN, NaN};
end

%% ========================================================================
function tf = is_future_case_complete(caseDir, P)

tf = exist(fullfile(caseDir, P.completion_marker), 'file') == 2 && ...
    file_has_bytes(fullfile(caseDir, 'fort.19'), 8) && ...
    file_has_bytes(fullfile(caseDir, 'fort.22'), 10) && ...
    exist(fullfile(caseDir, P.future_meta_mat), 'file') == 2;
end

%% ========================================================================
function tf = file_has_bytes(path, minBytes)

d = dir(path);
tf = ~isempty(d) && d(1).bytes >= minBytes;
end

%% ========================================================================
function [info, FutureMeta] = build_one_future_case(sourceCaseDir, futureCaseDir, ...
    task, scaleInfo, scenario, slrInfo, P)

metaPath = fullfile(sourceCaseDir, 'storm_track_forcing_meta.mat');
S = load(metaPath);
need = {'tc', 'Grid', 'Window', 'Cv_global', 'WindPressureFit', 'P', 'event'};
missing = setdiff(need, fieldnames(S));
assert(isempty(missing), 'Track metadata missing fields: %s', strjoin(missing, ', '));

era5ImpactVmax = double(scaleInfo.era5_impact_vmax_ms);
intensityRp = double(scaleInfo.era5_intensity_return_period_yr);
histIntensity = double(scaleInfo.cmip6_historical_intensity_ms);
futureIntensity = double(scaleInfo.cmip6_future_intensity_ms);
scaleFactor = double(scaleInfo.intensity_scale_factor);
assert(isfinite(scaleFactor) && scaleFactor > 0, ...
    'Invalid future intensity scale factor %.6g.', scaleFactor);

futureTc = build_scaled_future_tc(S.tc, scaleFactor, S.WindPressureFit);

ensure_dir(futureCaseDir);
copy_unchanged_case_files(sourceCaseDir, futureCaseDir, P);

fort19Source = fullfile(sourceCaseDir, 'fort.19');
fort19Dest = fullfile(futureCaseDir, 'fort.19');
write_slr_adjusted_fort19(fort19Source, fort19Dest, slrInfo.tc_local_slr_m);

fort22Dest = fullfile(futureCaseDir, 'fort.22');
[maxWind, maxWindTime, activeSteps] = write_future_fort22( ...
    fort22Dest, futureTc, S.Grid, S.Window, S.Cv_global, ...
    S.WindPressureFit, S.P, P);

FutureMeta = struct();
FutureMeta.source_case_dir = sourceCaseDir;
FutureMeta.future_case_dir = futureCaseDir;
FutureMeta.return_period_set = char(task.return_period_set);
FutureMeta.scenario = char(scenario);
FutureMeta.block_id = char(task.block_id);
FutureMeta.control_block_id = char(scaleInfo.control_block_id);
FutureMeta.control_block_label = char(scaleInfo.control_block_label);
FutureMeta.control_event_lon = double(scaleInfo.event_lon);
FutureMeta.control_event_lat = double(scaleInfo.event_lat);
FutureMeta.tc_intensity_spatial_method = ...
    ['one max-impact control block per unique TC, selected from all exact P4 ' ...
    'ADCIRC-inner events; one ratio shared by every rerun block for that TC and scenario'];
FutureMeta.era5_impact_vmax_ms = era5ImpactVmax;
FutureMeta.era5_intensity_return_period_yr = intensityRp;
FutureMeta.cmip6_historical_intensity_ms = histIntensity;
FutureMeta.cmip6_future_intensity_ms = futureIntensity;
FutureMeta.intensity_scale_factor = scaleFactor;
FutureMeta.intensity_change_percent = 100 * (scaleFactor - 1);
FutureMeta.curve_interpolation = ...
    'inverse ERA5 and forward CMIP6 interpolation in log10(return period)';
FutureMeta.tc_local_slr_m = slrInfo.tc_local_slr_m;
FutureMeta.slr_block_id = char(slrInfo.block_id);
FutureMeta.slr_block_event_lon = slrInfo.block_event_lon;
FutureMeta.slr_block_event_lat = slrInfo.block_event_lat;
FutureMeta.slr_block_event_vmax_ms = slrInfo.block_event_vmax_ms;
FutureMeta.slr_nominal_search_radius_km = slrInfo.nominal_search_radius_km;
FutureMeta.slr_search_radius_used_km = slrInfo.search_radius_used_km;
FutureMeta.slr_location_count = slrInfo.slr_location_count;
FutureMeta.slr_nearest_fallback = slrInfo.nearest_fallback;
FutureMeta.slr_spatial_mask = char(slrInfo.spatial_method);
FutureMeta.slr_track_point_count = slrInfo.track_point_count;
FutureMeta.slr_quantile = P.slr_quantile;
FutureMeta.slr_year_start = P.future_years(1);
FutureMeta.slr_year_end = P.future_years(end);
FutureMeta.fort22_max_wind_ms = maxWind;
FutureMeta.fort22_max_wind_time = maxWindTime;
FutureMeta.fort22_active_steps = activeSteps;
FutureMeta = add_track_change_metadata(FutureMeta, S.tc, futureTc);
FutureMeta.created_at = datetime('now');
save(fullfile(futureCaseDir, P.future_meta_mat), 'FutureMeta');
write_future_meta_text(fullfile(futureCaseDir, P.future_meta_txt), FutureMeta);
write_future_fort22_meta(fullfile(futureCaseDir, P.fort22_meta_file), ...
    sourceCaseDir, futureCaseDir, task, scenario, scaleInfo, slrInfo, ...
    S, futureTc, FutureMeta, P);

write_completion_marker(futureCaseDir, P);

info = struct('era5ImpactVmax', era5ImpactVmax, ...
    'intensityReturnPeriod', intensityRp, ...
    'historicalIntensity', histIntensity, ...
    'futureIntensity', futureIntensity, ...
    'scaleFactor', scaleFactor);
end

%% ========================================================================
function [FutureMeta, updated] = ensure_future_fort22_meta( ...
    sourceCaseDir, futureCaseDir, task, scaleInfo, scenario, slrInfo, P)

futureMetaPath = fullfile(futureCaseDir, P.future_meta_mat);
M = load(futureMetaPath, 'FutureMeta');
assert(isfield(M, 'FutureMeta'), 'Missing FutureMeta in %s.', futureMetaPath);
FutureMeta = M.FutureMeta;
scaleFactor = double(scaleInfo.intensity_scale_factor);
metaTextPath = fullfile(futureCaseDir, P.fort22_meta_file);
if future_fort22_meta_is_current(metaTextPath, task, scenario, scaleFactor, P)
    updated = false;
    return;
end

sourceMetaPath = fullfile(sourceCaseDir, 'storm_track_forcing_meta.mat');
S = load(sourceMetaPath);
need = {'tc', 'Grid', 'Window', 'WindPressureFit', 'event'};
missing = setdiff(need, fieldnames(S));
assert(isempty(missing), 'Track metadata missing fields: %s', strjoin(missing, ', '));
futureTc = build_scaled_future_tc(S.tc, scaleFactor, S.WindPressureFit);
FutureMeta = add_track_change_metadata(FutureMeta, S.tc, futureTc);
FutureMeta.fort22_metadata_refreshed_at = datetime('now');
save(futureMetaPath, 'FutureMeta');
write_future_meta_text(fullfile(futureCaseDir, P.future_meta_txt), FutureMeta);
write_future_fort22_meta(metaTextPath, sourceCaseDir, futureCaseDir, ...
    task, scenario, scaleInfo, slrInfo, S, futureTc, FutureMeta, P);
write_completion_marker(futureCaseDir, P);
updated = true;
end

%% ========================================================================
function futureTc = build_scaled_future_tc(sourceTc, scaleFactor, WindPressureFit)

futureTc = sourceTc;
futureTc.wind = double(sourceTc.wind) .* scaleFactor;
futureTc.p = wind_to_pressure_hpa(futureTc.wind, WindPressureFit);
if isfield(futureTc, 'event_vmax_ms')
    futureTc.event_vmax_ms = double(sourceTc.event_vmax_ms) .* scaleFactor;
end
if isfield(futureTc, 'global_anchor_vmax_ms')
    futureTc.global_anchor_vmax_ms = ...
        double(sourceTc.global_anchor_vmax_ms) .* scaleFactor;
end
end

%% ========================================================================
function M = add_track_change_metadata(M, sourceTc, futureTc)

M.historical_track_max_wind_ms = max(double(sourceTc.wind), [], 'omitnan');
M.future_track_max_wind_ms = max(double(futureTc.wind), [], 'omitnan');
M.historical_track_min_pressure_hpa = min(double(sourceTc.p), [], 'omitnan');
M.future_track_min_pressure_hpa = min(double(futureTc.p), [], 'omitnan');
end

%% ========================================================================
function tf = future_fort22_meta_is_current(path, task, scenario, scaleFactor, P)

tf = false;
if exist(path, 'file') ~= 2
    return;
end
try
    content = string(fileread(path));
catch
    return;
end
required = [ ...
    "metadata_version = " + string(P.fort22_meta_version)
    "forcing_period = future"
    "return_period_set = " + string(task.return_period_set)
    "scenario = " + string(scenario)
    "block_id = " + string(task.block_id)
    "intensity_scale_factor = " + string(sprintf('%.10f', scaleFactor))
    ];
tf = true;
for i = 1:numel(required)
    tf = tf && contains(content, required(i));
end
end

%% ========================================================================
function write_future_fort22_meta(path, sourceCaseDir, futureCaseDir, ...
    task, scenario, scaleInfo, slrInfo, S, futureTc, M, P)

tempPath = [path, '.tmp'];
fid = fopen(tempPath, 'wt');
assert(fid >= 0, 'Cannot write %s.', tempPath);
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'ADCIRC fort.22 (NWS=6) future TC metadata\n');
fprintf(fid, 'metadata_version = %s\n', P.fort22_meta_version);
fprintf(fid, 'forcing_period = future\n');
fprintf(fid, 'return_period_set = %s\n', char(task.return_period_set));
fprintf(fid, 'scenario = %s\n', char(scenario));
fprintf(fid, 'block_id = %s\n', char(task.block_id));
fprintf(fid, 'case_name = %s\n', char(task.case_name));
fprintf(fid, 'source_case_dir = %s\n', sourceCaseDir);
fprintf(fid, 'future_case_dir = %s\n', futureCaseDir);
fprintf(fid, 'source_historical_fort22_meta = %s\n', ...
    fullfile(sourceCaseDir, 'fort22_meta.txt'));
fprintf(fid, 'fort22 = %s\n', fullfile(futureCaseDir, 'fort.22'));

fprintf(fid, '\n[FutureIntensityScaling]\n');
fprintf(fid, 'control_block_id = %s\n', char(scaleInfo.control_block_id));
fprintf(fid, 'era5_impact_vmax_ms = %.8f\n', double(scaleInfo.era5_impact_vmax_ms));
fprintf(fid, 'era5_intensity_return_period_yr = %.8f\n', ...
    double(scaleInfo.era5_intensity_return_period_yr));
fprintf(fid, 'cmip6_historical_intensity_ms = %.8f\n', ...
    double(scaleInfo.cmip6_historical_intensity_ms));
fprintf(fid, 'cmip6_future_intensity_ms = %.8f\n', ...
    double(scaleInfo.cmip6_future_intensity_ms));
fprintf(fid, 'intensity_scale_factor = %.10f\n', ...
    double(scaleInfo.intensity_scale_factor));
fprintf(fid, 'intensity_change_percent = %.8f\n', ...
    100 .* (double(scaleInfo.intensity_scale_factor) - 1));
fprintf(fid, 'wind_change = future track wind equals historical track wind times scale factor\n');
fprintf(fid, ['pressure_change = center pressure recomputed from future wind; ' ...
    'stronger wind increases pressure deficit and lowers center pressure\n']);
fprintf(fid, 'pressure_fit_model = %s\n', char(string(S.WindPressureFit.model)));
fprintf(fid, 'pressure_fit_pref_hPa = %.8f\n', double(S.WindPressureFit.prefHpa));
fprintf(fid, 'pressure_fit_A = %.12g\n', double(S.WindPressureFit.A));
fprintf(fid, 'pressure_fit_B = %.12g\n', double(S.WindPressureFit.B));

fprintf(fid, '\n[HistoricalFutureTrackSummary]\n');
fprintf(fid, 'historical_track_max_wind_ms = %.8f\n', M.historical_track_max_wind_ms);
fprintf(fid, 'future_track_max_wind_ms = %.8f\n', M.future_track_max_wind_ms);
fprintf(fid, 'historical_track_min_pressure_hPa = %.8f\n', ...
    M.historical_track_min_pressure_hpa);
fprintf(fid, 'future_track_min_pressure_hPa = %.8f\n', M.future_track_min_pressure_hpa);

fprintf(fid, '\n[GridAndTime]\n');
fprintf(fid, 'WindowStart = %s\n', datestr(S.Window.start, 31));
fprintf(fid, 'WindowEnd = %s\n', datestr(S.Window.end, 31));
fprintf(fid, 'WTIMINC = %d\n', S.Window.WTIMINC);
fprintf(fid, 'nSteps = %d\n', numel(S.Window.time_data));
fprintf(fid, 'NWLON = %d\n', S.Grid.NWLON);
fprintf(fid, 'NWLAT = %d\n', S.Grid.NWLAT);
fprintf(fid, 'WLONMIN = %.8f\n', S.Grid.WLONMIN);
fprintf(fid, 'WLATMAX = %.8f\n', S.Grid.WLATMAX);
fprintf(fid, 'WLONINC = %.8f\n', S.Grid.DLON);
fprintf(fid, 'WLATINC = %.8f\n', S.Grid.DLAT);
fprintf(fid, 'fort22_uv_precision_ms = 0.1\n');
fprintf(fid, 'fort22_pressure_precision_pa = 1\n');

fprintf(fid, '\n[ForcingResult]\n');
fprintf(fid, 'active_wind_steps = %d\n', M.fort22_active_steps);
fprintf(fid, 'domain_max_wind_ms = %.8f\n', M.fort22_max_wind_ms);
fprintf(fid, 'domain_max_wind_time = %s\n', fmt_time_safe(M.fort22_max_wind_time));
fprintf(fid, 'tc_local_slr_m_for_fort19 = %.10f\n', double(slrInfo.tc_local_slr_m));
fprintf(fid, 'slr_spatial_mask = %s\n', char(slrInfo.spatial_method));

fprintf(fid, '\n[FutureForcingTrackPoints]\n');
fprintf(fid, ['index,original_time_index,time,lon180,lat,historical_vmax_ms,' ...
    'future_vmax_ms,historical_pressure_hPa,future_pressure_hPa\n']);
nTrack = numel(futureTc.time);
if isfield(futureTc, 'original_time_index') && ...
        numel(futureTc.original_time_index) == nTrack
    originalIndex = double(futureTc.original_time_index(:));
else
    originalIndex = (1:nTrack).';
end
for i = 1:nTrack
    fprintf(fid, '%d,%d,%s,%.8f,%.8f,%.6f,%.6f,%.6f,%.6f\n', ...
        i, round(originalIndex(i)), datestr(futureTc.time(i), 31), ...
        double(futureTc.lon(i)), double(futureTc.lat(i)), ...
        double(S.tc.wind(i)), double(futureTc.wind(i)), ...
        double(S.tc.p(i)), double(futureTc.p(i)));
end

clear cleanupObj;
movefile(tempPath, path, 'f');
end

%% ========================================================================
function write_completion_marker(caseDir, P)

fid = fopen(fullfile(caseDir, P.completion_marker), 'wt');
assert(fid >= 0, 'Cannot write completion marker in %s.', caseDir);
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'completed=%s\n', datestr(now, 31));
end

%% ========================================================================
function value = fmt_time_safe(t)

if isempty(t) || (isdatetime(t) && all(isnat(t)))
    value = 'NaT';
else
    value = datestr(t(1), 31);
end
end

%% ========================================================================
function copy_unchanged_case_files(sourceDir, destDir, P)

for i = 1:numel(P.unchanged_files)
    name = P.unchanged_files{i};
    src = fullfile(sourceDir, name);
    dst = fullfile(destDir, name);
    if exist(src, 'file') ~= 2
        if strcmp(name, 'sub_intel.sh')
            continue;
        end
        error('Missing unchanged source file: %s', src);
    end
    if exist(dst, 'file') == 2
        continue;
    end
    if strcmp(name, 'fort.14') && P.hardlink_fort14
        ok = create_hardlink(dst, src);
        if ok
            continue;
        end
        warning('Hardlink failed; copying fort.14: %s', src);
    end
    copyfile(src, dst);
end
end

%% ========================================================================
function ok = create_hardlink(linkPath, targetPath)

ok = false;
if ispc
    cmd = sprintf('cmd /c mklink /H "%s" "%s" >nul', linkPath, targetPath);
else
    cmd = sprintf('ln "%s" "%s"', targetPath, linkPath);
end
[status, ~] = system(cmd);
ok = status == 0 && exist(linkPath, 'file') == 2;
end

%% ========================================================================
function write_slr_adjusted_fort19(sourcePath, destPath, slrM)

assert(exist(sourcePath, 'file') == 2, 'Missing source fort.19: %s', sourcePath);
tempPath = [destPath, '.tmp'];
fidIn = fopen(sourcePath, 'rt');
assert(fidIn >= 0, 'Cannot read %s.', sourcePath);
cleanupIn = onCleanup(@() fclose(fidIn)); %#ok<NASGU>
fidOut = fopen(tempPath, 'wt');
assert(fidOut >= 0, 'Cannot write %s.', tempPath);
cleanupOut = onCleanup(@() fclose(fidOut)); %#ok<NASGU>

firstLine = fgetl(fidIn);
assert(ischar(firstLine), 'Empty source fort.19: %s', sourcePath);
fprintf(fidOut, '%s\n', firstLine);
nValues = 0;
while true
    line = fgetl(fidIn);
    if ~ischar(line)
        break;
    end
    value = str2double(strtrim(line));
    assert(isfinite(value), 'Non-numeric fort.19 value in %s.', sourcePath);
    fprintf(fidOut, '%.10f\n', value + slrM);
    nValues = nValues + 1;
end
assert(nValues > 0, 'No water-level records found in %s.', sourcePath);
clear cleanupOut cleanupIn;
movefile(tempPath, destPath, 'f');
end

%% ========================================================================
function [maxWindDomain, maxWindTime, activeSteps] = write_future_fort22( ...
    destPath, tc, Grid, Window, Cv, WindPressureFit, ModelP, P)

ModelP = fill_model_defaults(ModelP);
tempPath = [destPath, '.tmp'];
if exist(tempPath, 'file') == 2
    delete(tempPath);
end
fid = fopen(tempPath, 'wt');
assert(fid >= 0, 'Cannot write future fort.22: %s', tempPath);
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
            Cv, ModelP.Pn_hPa, ModelP.Re, WindPressureFit, ModelP, P);
        windMag = hypot(U, V);
        stepMax = max(windMag(:), [], 'omitnan');
        if isfinite(stepMax) && stepMax > 0
            activeSteps = activeSteps + 1;
        end
        if isfinite(stepMax) && stepMax > maxWindDomain
            maxWindDomain = stepMax;
            maxWindTime = t;
        end
        write_fort22_snapshot_nws6(fid, U, V, Pres, ModelP);
    else
        write_background_step_nws6(fid, Grid.NWLAT, Grid.NWLON, ModelP);
    end
end
clear cleanupObj;
movefile(tempPath, destPath, 'f');
end

%% ========================================================================
function P = fill_model_defaults(P)

defaults = struct('Pn_hPa', 1013, 'Re', 6371000, 'standard_pressure_pa', 101300, ...
    'background_fort22_line', '0 0 101300', 'fort22_zero_uv_tolerance', 1e-8, ...
    'fort22_zero_pressure_tolerance_pa', 1e-4, 'c15_min_lookup_vm_ms', 15, ...
    'c15_max_lookup_vm_ms', 120, 'c15_min_rmax_m', 30e3, ...
    'c15_holland_b_min', 1.0, 'c15_holland_b_max', 2.5, ...
    'c15_one_min_to_ten_min_factor', 0.893, ...
    'c15_blend_inner_radius_m', 500e3, 'c15_blend_outer_radius_m', 700e3, ...
    'c15_beta_inner_base_deg', 10, 'c15_beta_inner_slope_deg', 10, ...
    'c15_beta_mid_base_deg', 20, 'c15_beta_mid_slope_deg', 25, ...
    'c15_beta_outer_deg', 25, 'c15_beta_mid_radius_factor', 1.2);
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(P, names{i})
        P.(names{i}) = defaults.(names{i});
    end
end
end

%% ========================================================================
function [hasTc, stateNow, stateNext] = get_tc_state_for_step(tc, t, dtStep, WindPressureFit)

if isempty(tc.time) || t < tc.time(1) || t > tc.time(end)
    hasTc = false;
    stateNow = empty_tc_state();
    stateNext = empty_tc_state();
    return;
end
stateNow = interp_tc_state(tc, t, WindPressureFit);
t2 = t + dtStep;
if t2 > tc.time(end)
    stateNext = stateNow;
else
    stateNext = interp_tc_state(tc, t2, WindPressureFit);
    if ~isfinite(stateNext.lat) || ~isfinite(stateNext.lon) || ~isfinite(stateNext.wind)
        stateNext = stateNow;
    end
end
hasTc = isfinite(stateNow.lat) && isfinite(stateNow.lon) && isfinite(stateNow.wind);
end

%% ========================================================================
function s = interp_tc_state(tc, tQuery, WindPressureFit)

tSrc = datenum(tc.time(:));
tq = datenum(tQuery);
[tSrc, ia] = unique(tSrc, 'stable');
lat = tc.lat(ia);
if isfield(tc, 'lon_unwrapped') && numel(tc.lon_unwrapped) == numel(tc.lon)
    lon = tc.lon_unwrapped(ia);
else
    lon = rad2deg(unwrap(deg2rad(tc.lon(ia))));
end
wind = tc.wind(ia);
pressure = tc.p(ia);
s.lat = interp1(tSrc, lat, tq, 'linear', NaN);
s.lon = wrap_to_180_local(interp1(tSrc, lon, tq, 'linear', NaN));
s.wind = interp1(tSrc, wind, tq, 'linear', NaN);
s.p = interp1(tSrc, pressure, tq, 'linear', NaN);
if ~isfinite(s.p) || s.p <= 0
    s.p = wind_to_pressure_hpa(s.wind, WindPressureFit);
end
end

%% ========================================================================
function s = empty_tc_state()

s = struct('lat', NaN, 'lon', NaN, 'wind', NaN, 'p', NaN);
end

%% ========================================================================
function [U, V, Pres] = calc_c15_uvp_field_grid( ...
    latC, lonC, latN, lonN, dtHr, vC, pC, LAT, LON, Cv, PnHpa, Re, ...
    WindPressureFit, Params, BuildP) %#ok<INUSD>

U = zeros(size(LAT));
V = zeros(size(LAT));
Pres = ones(size(LAT)) * Params.standard_pressure_pa;
if ~isfinite(latC) || ~isfinite(lonC) || ~isfinite(vC)
    return;
end
if ~isfinite(pC) || pC <= 0
    pC = wind_to_pressure_hpa(vC, WindPressureFit);
end
if ~isfinite(pC) || pC >= PnHpa
    return;
end

if dtHr <= 0 || ~isfinite(dtHr)
    vmc = 0;
    fai = 0;
else
    [distance, fai] = fast_track_step(latC, lonC, latN, lonN);
    vmc = distance / (dtHr * 3600);
end
vmRaw = vC - vmc;
weakScale = min(max(vmRaw, 0) ./ Params.c15_min_lookup_vm_ms, 1);
if weakScale <= 0
    return;
end
vm = max(vmRaw, Params.c15_min_lookup_vm_ms);
vm = min(vm, Params.c15_max_lookup_vm_ms);
pModelHpa = PnHpa - weakScale .* (PnHpa - pC);

rmaxM = Cv * 51.6 * exp(-0.0223 * vm + 0.0281 * abs(latC)) * 1000;
rmaxM = max(rmaxM, Params.c15_min_rmax_m);
B = (vC^2) * 1.15 * exp(1) / max(PnHpa - pModelHpa, 1) / 100;
B = max(Params.c15_holland_b_min, min(B, Params.c15_holland_b_max));

dlon = wrap_to_180_local(LON - lonC);
dx = dlon .* cosd(0.5 * (LAT + latC)) * 111320;
dy = (LAT - latC) * 110540;
r = max(hypot(dx, dy), 1.0);
cta = atan2d(dy, dx);
cta(cta < 0) = cta(cta < 0) + 360;
pgTc = (pModelHpa + (PnHpa - pModelHpa) .* exp(-(rmaxM ./ r).^B)) * 100;

vmaxKey = round(vm);
rmaxKeyKm = round(rmaxM / 1000);
C15 = get_c15_lookup(Params.c15_predata_dir, vmaxKey, rmaxKeyKm, BuildP.c15_cache_size);
if isempty(C15)
    c15_surface_wind_1min_ms = zeros(size(r));
else
    % C15.vg is the lookup file's axisymmetric 1-min surface tangential
    % wind speed, not a gradient, upper-tropospheric or environmental wind.
    c15_surface_wind_1min_ms = interp1( ...
        double(C15.rr(:)), double(C15.vg(:)), r, 'linear', 0);
end
c15_surface_wind_1min_ms(~isfinite(c15_surface_wind_1min_ms) | ...
    c15_surface_wind_1min_ms < 0) = 0;
c15_surface_wind_1min_ms = weakScale .* c15_surface_wind_1min_ms;

beta = zeros(size(r));
i1 = r < rmaxM;
i2 = r >= rmaxM & r < Params.c15_beta_mid_radius_factor * rmaxM;
i3 = r >= Params.c15_beta_mid_radius_factor * rmaxM;
beta(i1) = Params.c15_beta_inner_base_deg + ...
    Params.c15_beta_inner_slope_deg .* (r(i1) ./ rmaxM);
beta(i2) = Params.c15_beta_mid_base_deg + ...
    Params.c15_beta_mid_slope_deg .* (r(i2) ./ rmaxM - 1);
beta(i3) = Params.c15_beta_outer_deg;

vmoc = weakScale .* vmc .* r .* rmaxM ./ (r.^2 + rmaxM^2);
hemisphereSign = sign(latC);
if hemisphereSign == 0
    hemisphereSign = 1;
end
rotationAngle = hemisphereSign .* (90 + beta);
vxTc = c15_surface_wind_1min_ms .* cosd(cta + rotationAngle) + ...
    vmoc .* cosd(fai);
vyTc = c15_surface_wind_1min_ms .* sind(cta + rotationAngle) + ...
    vmoc .* sind(fai);
% Convert the complete 1-min surface wind vector to the 10-min wind used by
% ADCIRC fort.22. This is an averaging-period conversion, not a
% gradient-to-surface wind reduction.
vxTc = Params.c15_one_min_to_ten_min_factor .* vxTc;
vyTc = Params.c15_one_min_to_ten_min_factor .* vyTc;

r1 = Params.c15_blend_inner_radius_m;
r2 = Params.c15_blend_outer_radius_m;
lambda = zeros(size(r));
iMid = r >= r1 & r <= r2;
iFar = r > r2;
lambda(iMid) = (r(iMid) - r1) ./ (r2 - r1);
lambda(iFar) = 1;
U = (1 - lambda) .* vxTc;
V = (1 - lambda) .* vyTc;
Pres = (1 - lambda) .* pgTc + lambda .* Params.standard_pressure_pa;
end

%% ========================================================================
function C15 = get_c15_lookup(predataDir, vmax, rmaxKm, maxCacheSize)

persistent Cache CacheOrder CacheDir
if isempty(Cache) || ~strcmp(CacheDir, predataDir)
    Cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    CacheOrder = cell(0, 1);
    CacheDir = predataDir;
end
key = sprintf('%d_%d', vmax, rmaxKm);
if isKey(Cache, key)
    C15 = Cache(key);
    return;
end
path = fullfile(predataDir, ...
    sprintf('Wind_C15_data_Vmax%d_Rmax%d.mat', vmax, rmaxKm));
if exist(path, 'file') == 2
    S = load(path, 'Wind_C15_data');
    C15 = S.Wind_C15_data;
else
    C15 = [];
end
if Cache.Count >= maxCacheSize && ~isempty(CacheOrder)
    remove(Cache, CacheOrder{1});
    CacheOrder(1) = [];
end
Cache(key) = C15;
CacheOrder{end + 1, 1} = key;
end

%% ========================================================================
function [distanceM, alphaDeg] = fast_track_step(lat1, lon1, lat2, lon2)

dlon = wrap_to_180_local(lon2 - lon1);
dx = dlon * cosd(0.5 * (lat1 + lat2)) * 111320;
dy = (lat2 - lat1) * 110540;
distanceM = hypot(dx, dy);
alphaDeg = atan2d(dy, dx);
if alphaDeg < 0
    alphaDeg = alphaDeg + 360;
end
end

%% ========================================================================
function write_background_step_nws6(fid, nLat, nLon, P)

for k = 1:nLat
    for j = 1:nLon
        fprintf(fid, '%s\n', P.background_fort22_line);
    end
end
end

%% ========================================================================
function write_fort22_snapshot_nws6(fid, U, V, Pres, P)

[nLat, nLon] = size(U);
for k = 1:nLat
    for j = 1:nLon
        u = U(k, j);
        v = V(k, j);
        pressure = Pres(k, j);
        if abs(u) < P.fort22_zero_uv_tolerance && ...
                abs(v) < P.fort22_zero_uv_tolerance && ...
                abs(pressure - P.standard_pressure_pa) < P.fort22_zero_pressure_tolerance_pa
            fprintf(fid, '%s\n', P.background_fort22_line);
        else
            fprintf(fid, '%.1f %.1f %.0f\n', u, v, pressure);
        end
    end
end
end

%% ========================================================================
function pHpa = wind_to_pressure_hpa(windMs, fit)

pHpa = nan(size(windMs));
windMs = double(windMs);
good = isfinite(windMs) & windMs >= 0;
if ~any(good)
    return;
end
deltaP = (windMs(good) ./ fit.A) .^ (1 ./ fit.B);
pressure = fit.prefHpa - deltaP;
pressure = max(fit.pcMinHpa, min(fit.pcMaxHpa, pressure));
pHpa(good) = pressure;
end

%% ========================================================================
function lon = wrap_to_180_local(lon)

lon = mod(double(lon) + 180, 360) - 180;
end

%% ========================================================================
function write_future_meta_text(path, M)

fid = fopen(path, 'wt');
assert(fid >= 0, 'Cannot write %s.', path);
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'Future ADCIRC forcing metadata\n');
fprintf(fid, 'source_case_dir = %s\n', M.source_case_dir);
fprintf(fid, 'return_period_set = %s\n', M.return_period_set);
fprintf(fid, 'scenario = %s\n', M.scenario);
fprintf(fid, 'block_id = %s\n', M.block_id);
fprintf(fid, 'control_block_id = %s\n', M.control_block_id);
fprintf(fid, 'control_block_label = %s\n', M.control_block_label);
fprintf(fid, 'control_event_lon = %.8f\n', M.control_event_lon);
fprintf(fid, 'control_event_lat = %.8f\n', M.control_event_lat);
fprintf(fid, 'tc_intensity_spatial_method = %s\n', M.tc_intensity_spatial_method);
fprintf(fid, 'era5_impact_vmax_ms = %.8f\n', M.era5_impact_vmax_ms);
fprintf(fid, 'era5_intensity_return_period_yr = %.8f\n', M.era5_intensity_return_period_yr);
fprintf(fid, 'cmip6_historical_intensity_ms = %.8f\n', M.cmip6_historical_intensity_ms);
fprintf(fid, 'cmip6_future_intensity_ms = %.8f\n', M.cmip6_future_intensity_ms);
fprintf(fid, 'intensity_scale_factor = %.10f\n', M.intensity_scale_factor);
fprintf(fid, 'intensity_change_percent = %.8f\n', M.intensity_change_percent);
if isfield(M, 'historical_track_max_wind_ms')
    fprintf(fid, 'historical_track_max_wind_ms = %.8f\n', M.historical_track_max_wind_ms);
    fprintf(fid, 'future_track_max_wind_ms = %.8f\n', M.future_track_max_wind_ms);
    fprintf(fid, 'historical_track_min_pressure_hPa = %.8f\n', ...
        M.historical_track_min_pressure_hpa);
    fprintf(fid, 'future_track_min_pressure_hPa = %.8f\n', ...
        M.future_track_min_pressure_hpa);
end
fprintf(fid, 'tc_local_slr_m = %.10f\n', M.tc_local_slr_m);
fprintf(fid, 'slr_block_id = %s\n', M.slr_block_id);
fprintf(fid, 'slr_block_event_lon = %.8f\n', M.slr_block_event_lon);
fprintf(fid, 'slr_block_event_lat = %.8f\n', M.slr_block_event_lat);
fprintf(fid, 'slr_block_event_vmax_ms = %.8f\n', M.slr_block_event_vmax_ms);
fprintf(fid, 'slr_spatial_mask = %s\n', M.slr_spatial_mask);
fprintf(fid, 'slr_nominal_search_radius_km = %.3f\n', M.slr_nominal_search_radius_km);
fprintf(fid, 'slr_search_radius_used_km = %.6f\n', M.slr_search_radius_used_km);
fprintf(fid, 'slr_location_count = %d\n', M.slr_location_count);
fprintf(fid, 'slr_nearest_fallback = %d\n', M.slr_nearest_fallback);
fprintf(fid, 'slr_track_point_count = %d\n', M.slr_track_point_count);
fprintf(fid, 'slr_quantile = %.3f\n', M.slr_quantile);
fprintf(fid, 'slr_years = %d-%d\n', M.slr_year_start, M.slr_year_end);
fprintf(fid, 'fort22_max_wind_ms = %.8f\n', M.fort22_max_wind_ms);
fprintf(fid, 'fort22_active_steps = %d\n', M.fort22_active_steps);
fprintf(fid, 'changed_files = fort.19, fort.22\n');
fprintf(fid, 'unchanged_original_inputs = fort.13, fort.14, fort.15\n');
end

%% ========================================================================
function ensure_dir(path)

if exist(path, 'dir') ~= 7
    mkdir(path);
end
end

%% ========================================================================
function T = read_csv_table(path)

T = readtable(path, 'Delimiter', ',', 'TextType', 'string', ...
    'VariableNamingRule', 'preserve');
end

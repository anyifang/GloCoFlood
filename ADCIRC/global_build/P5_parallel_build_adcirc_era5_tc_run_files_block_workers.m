%% P5_parallel_build_adcirc_era5_tc_run_files_block_workers.m
% Prepare and optionally launch external MATLAB workers for building ERA5 TC
% ADCIRC forcing cases. This uses coarse external parallelism: each worker
% gets fixed ADCIRC blocks and a private TMD copy, so no two workers write to
% the same TMD data.out / LAT_LON / data_*.mat files.

clearvars;
clc;

SCRIPT_DIR = fileparts(mfilename('fullpath'));
if isempty(SCRIPT_DIR)
    SCRIPT_DIR = pwd;
end

P = struct();
P.base_build_script = fullfile(SCRIPT_DIR, 'P5_build_adcirc_era5_tc_run_files_from_selected_tracks.m');
P.event_csv = fullfile(SCRIPT_DIR, 'era5_tc_events_by_adcirc_inner_domain_vmax33_matlab.csv');
P.output_root = fullfile(SCRIPT_DIR, 'output', 'adcirc_era5_tc_selected_run_files');
P.worker_root = fullfile(SCRIPT_DIR, 'adcirc_era5_tc_parallel_workers');
P.worker_script_dir = fullfile(P.worker_root, 'worker_scripts');
P.worker_log_dir = fullfile(P.worker_root, 'logs');
P.worker_plan_csv = fullfile(P.worker_root, 'parallel_worker_plan.csv');
P.worker_plan_mat = fullfile(P.worker_root, 'parallel_worker_plan.mat');
P.master_launch_bat = fullfile(P.worker_root, 'launch_all_parallel_workers.bat');

% Keep this modest by default. Each worker is a full MATLAB process and loads
% the C15 lookup tables. Increase only when RAM and disk I/O are sufficient.
P.worker_count = 8;
P.launch_workers = false;
P.monitor_after_launch = false;

% TMD2.5 is used as the template. The script creates TMD2.50, TMD2.51, ...
% under the same parent folder. Each worker uses one private ...\TMD folder.
P.tmd_template_package = fullfile(SCRIPT_DIR, 'external', 'TMD2.5');
P.tmd_copy_parent = fileparts(P.tmd_template_package);
P.tmd_copy_start_index = 50;
P.create_missing_tmd_copies = true;

P.matlab_exe = 'matlab';
P.generate_fort19_from_tmd = true;
P.overwrite_existing_cases = false;  % false lets interrupted/restarted runs skip complete cases.
% Must be identical for every worker. It anchors each TC's global LMI once;
% block-specific impact indices then select different local run windows.
P.random_day_seed = 20260705;
P.selected_blocks = strings(0, 1);   % empty = all blocks in event CSV.
P.block_priority_prefixes = ["ADC_NA_", "ADC_WNP_"];  % workers process North Atlantic first, then WNP, then other blocks.

% Optional testing controls. Leave as inf for the full 30k-case run.
P.max_cases_per_block = inf;

P = apply_environment_overrides(P);

assert(exist(P.base_build_script, 'file') == 2, 'Missing base build script: %s', P.base_build_script);
assert(exist(P.event_csv, 'file') == 2, 'Missing event CSV: %s', P.event_csv);
assert(exist(P.tmd_template_package, 'dir') == 7, 'Missing TMD template package: %s', P.tmd_template_package);

ensure_dir(P.worker_root);
ensure_dir(P.worker_script_dir);
ensure_dir(P.worker_log_dir);
ensure_dir(P.output_root);

Events = readtable(P.event_csv, 'TextType', 'string');
Events.block_id = string(Events.block_id);
if ~isempty(P.selected_blocks)
    Events = Events(ismember(Events.block_id, P.selected_blocks), :);
end
assert(height(Events) > 0, 'No events remain after block filtering.');

Plan = build_worker_plan(Events, P);
prepare_tmd_copies(Plan, P);
Plan = write_worker_files(Plan, P, SCRIPT_DIR);
writetable(Plan, P.worker_plan_csv);
save(P.worker_plan_mat, 'Plan', 'P', '-v7.3');
write_master_launch_bat(Plan, P);
write_merge_script(Plan, P, SCRIPT_DIR);

fprintf('\n============================================================\n');
fprintf('External worker plan prepared\n');
fprintf('Event CSV      : %s\n', P.event_csv);
fprintf('Output root    : %s\n', P.output_root);
fprintf('Worker count   : %d\n', height(Plan));
fprintf('Plan CSV       : %s\n', P.worker_plan_csv);
fprintf('Plan MAT       : %s\n', P.worker_plan_mat);
fprintf('Launch BAT     : %s\n', P.master_launch_bat);
fprintf('TMD template   : %s\n', P.tmd_template_package);
fprintf('TC date seed   : %d (shared by all workers)\n', P.random_day_seed);
fprintf('Launch workers : %d\n', P.launch_workers);
fprintf('============================================================\n');
disp(Plan(:, {'worker_id', 'n_blocks', 'n_events', 'tmd_root', 'summary_csv_name'}));

if P.launch_workers
    fprintf('\nLaunching workers...\n');
    for i = 1:height(Plan)
        cmd = sprintf('start "ADCIRC_TC_%s" /min cmd /c ""%s""', ...
            char(Plan.worker_id(i)), char(Plan.bat_file(i)));
        [status, msg] = system(cmd);
        if status ~= 0
            warning('Failed to launch %s: %s', char(Plan.worker_id(i)), msg);
        else
            fprintf('  launched %s -> %s\n', char(Plan.worker_id(i)), char(Plan.log_file(i)));
        end
    end
    if P.monitor_after_launch
        fprintf('\nEntering MATLAB progress monitor. Press Ctrl+C to stop monitoring only; workers keep running.\n');
        monitor_adcirc_era5_tc_parallel_workers(P.worker_root);
    end
else
    fprintf('\nWorkers were not launched. To start them, run:\n  %s\n', P.master_launch_bat);
end

%% ========================================================================
function P = apply_environment_overrides(P)

v = str2double(strtrim(string(getenv('ADCIRC_TC_PARALLEL_WORKERS'))));
if isfinite(v) && v >= 1
    P.worker_count = floor(v);
end

v = strtrim(string(getenv('ADCIRC_TC_PARALLEL_LAUNCH')));
if any(strcmpi(v, ["1", "true", "yes", "y", "on"]))
    P.launch_workers = true;
elseif any(strcmpi(v, ["0", "false", "no", "n", "off"]))
    P.launch_workers = false;
end

v = strtrim(string(getenv('ADCIRC_TC_PARALLEL_MONITOR')));
if any(strcmpi(v, ["1", "true", "yes", "y", "on"]))
    P.monitor_after_launch = true;
elseif any(strcmpi(v, ["0", "false", "no", "n", "off"]))
    P.monitor_after_launch = false;
end

v = strtrim(string(getenv('ADCIRC_TC_PARALLEL_CREATE_TMD')));
if any(strcmpi(v, ["1", "true", "yes", "y", "on"]))
    P.create_missing_tmd_copies = true;
elseif any(strcmpi(v, ["0", "false", "no", "n", "off"]))
    P.create_missing_tmd_copies = false;
end

v = strtrim(string(getenv('ADCIRC_TC_CASE_OUTPUT_ROOT')));
if strlength(v) > 0
    P.output_root = char(v);
end

v = strtrim(string(getenv('ADCIRC_TC_PARALLEL_WORKER_ROOT')));
if strlength(v) > 0
    P.worker_root = char(v);
    P.worker_script_dir = fullfile(P.worker_root, 'worker_scripts');
    P.worker_log_dir = fullfile(P.worker_root, 'logs');
    P.worker_plan_csv = fullfile(P.worker_root, 'parallel_worker_plan.csv');
    P.worker_plan_mat = fullfile(P.worker_root, 'parallel_worker_plan.mat');
    P.master_launch_bat = fullfile(P.worker_root, 'launch_all_parallel_workers.bat');
end

v = strtrim(string(getenv('ADCIRC_TC_EVENT_CSV')));
if strlength(v) > 0
    P.event_csv = char(v);
end

v = strtrim(string(getenv('ADCIRC_TC_PARALLEL_BLOCKS')));
if strlength(v) > 0
    parts = split(v, ',');
    P.selected_blocks = strip(parts(strlength(strip(parts)) > 0));
end

v = strtrim(string(getenv('ADCIRC_TC_PARALLEL_BLOCK_PRIORITY')));
if strlength(v) > 0
    parts = split(v, ',');
    P.block_priority_prefixes = strip(parts(strlength(strip(parts)) > 0));
end

v = strtrim(string(getenv('ADCIRC_TC_PARALLEL_TMD_TEMPLATE')));
if strlength(v) > 0
    P.tmd_template_package = char(v);
    P.tmd_copy_parent = fileparts(P.tmd_template_package);
end

v = strtrim(string(getenv('ADCIRC_TC_PARALLEL_TMD_PARENT')));
if strlength(v) > 0
    P.tmd_copy_parent = char(v);
end

v = str2double(strtrim(string(getenv('ADCIRC_TC_PARALLEL_TMD_START_INDEX'))));
if isfinite(v) && v >= 0
    P.tmd_copy_start_index = floor(v);
end

v = strtrim(string(getenv('ADCIRC_TC_PARALLEL_MATLAB_EXE')));
if strlength(v) > 0
    P.matlab_exe = char(v);
end

v = strtrim(string(getenv('ADCIRC_TC_CASE_OVERWRITE')));
if any(strcmpi(v, ["1", "true", "yes", "y", "on"]))
    P.overwrite_existing_cases = true;
elseif any(strcmpi(v, ["0", "false", "no", "n", "off"]))
    P.overwrite_existing_cases = false;
end

v = str2double(strtrim(string(getenv('ADCIRC_TC_CASE_MAX_PER_BLOCK'))));
if isfinite(v) && v >= 0
    P.max_cases_per_block = floor(v);
end

v = str2double(strtrim(string(getenv('ADCIRC_TC_RANDOM_DAY_SEED'))));
if isfinite(v) && v >= 0
    P.random_day_seed = floor(v);
end
end

%% ========================================================================
function Plan = build_worker_plan(Events, P)

[blocks, ~, ib] = unique(Events.block_id, 'stable');
counts = accumarray(ib, 1);
nWorkers = min(P.worker_count, numel(blocks));

blockPriority = block_priority_rank(blocks, P.block_priority_prefixes);
priorityGroups = sort(unique(blockPriority));

workerBlocks = cell(nWorkers, 1);
workerCounts = zeros(nWorkers, 1);
for ig = 1:numel(priorityGroups)
    idx = find(blockPriority == priorityGroups(ig));
    [~, order] = sort(counts(idx), 'descend');
    idx = idx(order);
    for j = 1:numel(idx)
        k = idx(j);
        [~, iw] = min(workerCounts);
        workerBlocks{iw} = [workerBlocks{iw}; blocks(k)]; %#ok<AGROW>
        workerCounts(iw) = workerCounts(iw) + counts(k);
    end
end

worker_id = strings(nWorkers, 1);
tmd_package = strings(nWorkers, 1);
tmd_root = strings(nWorkers, 1);
block_ids = strings(nWorkers, 1);
n_blocks = zeros(nWorkers, 1);
n_events = zeros(nWorkers, 1);
summary_csv_name = strings(nWorkers, 1);
config_mat_name = strings(nWorkers, 1);
summary_csv = strings(nWorkers, 1);
wrapper_m = strings(nWorkers, 1);
bat_file = strings(nWorkers, 1);
log_file = strings(nWorkers, 1);

for iw = 1:nWorkers
    worker_id(iw) = sprintf('worker_%02d', iw);
    copyName = sprintf('TMD2.%02d', P.tmd_copy_start_index + iw - 1);
    tmd_package(iw) = string(fullfile(P.tmd_copy_parent, copyName));
    tmd_root(iw) = string(fullfile(char(tmd_package(iw)), 'TMD'));
    block_ids(iw) = strjoin(workerBlocks{iw}, ',');
    n_blocks(iw) = numel(workerBlocks{iw});
    n_events(iw) = workerCounts(iw);
    summary_csv_name(iw) = sprintf('adcirc_era5_tc_case_build_summary_%s.csv', worker_id(iw));
    config_mat_name(iw) = sprintf('adcirc_era5_tc_case_build_config_%s.mat', worker_id(iw));
    summary_csv(iw) = string(fullfile(P.output_root, char(summary_csv_name(iw))));
    wrapper_m(iw) = string(fullfile(P.worker_script_dir, sprintf('%s_build.m', worker_id(iw))));
    bat_file(iw) = string(fullfile(P.worker_script_dir, sprintf('%s_build.bat', worker_id(iw))));
    log_file(iw) = string(fullfile(P.worker_log_dir, sprintf('%s_build.log', worker_id(iw))));
end

Plan = table(worker_id, n_blocks, n_events, block_ids, tmd_package, tmd_root, ...
    summary_csv_name, config_mat_name, summary_csv, wrapper_m, bat_file, log_file);
end

%% ========================================================================
function rank = block_priority_rank(blocks, priorityPrefixes)

blocks = string(blocks);
priorityPrefixes = string(priorityPrefixes(:));
priorityPrefixes = priorityPrefixes(strlength(priorityPrefixes) > 0);

rank = repmat(numel(priorityPrefixes) + 1, numel(blocks), 1);
for i = 1:numel(priorityPrefixes)
    hit = startsWith(blocks, priorityPrefixes(i));
    rank(hit & rank > i) = i;
end
end

%% ========================================================================
function prepare_tmd_copies(Plan, P)

for i = 1:height(Plan)
    pkg = char(Plan.tmd_package(i));
    root = char(Plan.tmd_root(i));
    if exist(root, 'dir') == 7
        continue;
    end
    if ~P.create_missing_tmd_copies
        warning('TMD copy is missing and creation is disabled: %s', root);
        continue;
    end
    fprintf('Creating TMD copy for %s:\n  %s\n  -> %s\n', ...
        char(Plan.worker_id(i)), P.tmd_template_package, pkg);
    ok = copyfile(P.tmd_template_package, pkg);
    if ~ok || exist(root, 'dir') ~= 7
        error('Failed to create TMD copy: %s', pkg);
    end
end
end

%% ========================================================================
function Plan = write_worker_files(Plan, P, scriptDir)

for i = 1:height(Plan)
    if P.generate_fort19_from_tmd && exist(char(Plan.tmd_root(i)), 'dir') ~= 7
        warning('%s TMD root is missing: %s', char(Plan.worker_id(i)), char(Plan.tmd_root(i)));
    end

    wrapperLines = {
        sprintf('%% Auto-generated by P5_parallel_build_adcirc_era5_tc_run_files_block_workers.m')
        sprintf('setenv(''ADCIRC_TC_WORKER_ID'', ''%s'');', char(Plan.worker_id(i)))
        sprintf('setenv(''ADCIRC_TC_CASE_BLOCKS'', ''%s'');', escape_matlab_string(char(Plan.block_ids(i))))
        sprintf('setenv(''ADCIRC_TC_EVENT_CSV'', ''%s'');', escape_matlab_string(P.event_csv))
        sprintf('setenv(''ADCIRC_TC_CASE_OUTPUT_ROOT'', ''%s'');', escape_matlab_string(P.output_root))
        sprintf('setenv(''ADCIRC_TC_RANDOM_DAY_SEED'', ''%d'');', P.random_day_seed)
        sprintf('setenv(''ADCIRC_TC_TMD_ROOT'', ''%s'');', escape_matlab_string(char(Plan.tmd_root(i))))
        sprintf('setenv(''ADCIRC_TC_SUMMARY_CSV_NAME'', ''%s'');', char(Plan.summary_csv_name(i)))
        sprintf('setenv(''ADCIRC_TC_CONFIG_MAT_NAME'', ''%s'');', char(Plan.config_mat_name(i)))
        sprintf('setenv(''ADCIRC_TC_GENERATE_FORT19_FROM_TMD'', ''%d'');', double(P.generate_fort19_from_tmd))
        sprintf('setenv(''ADCIRC_TC_CASE_OVERWRITE'', ''%d'');', double(P.overwrite_existing_cases))
        sprintf('setenv(''ADCIRC_TC_CASE_MAX_PER_BLOCK'', ''%s'');', numeric_env_value(P.max_cases_per_block))
        sprintf('cd(''%s'');', escape_matlab_string(scriptDir))
        sprintf('run(''%s'');', escape_matlab_string(P.base_build_script))
        };
    write_lines(char(Plan.wrapper_m(i)), wrapperLines);

    matlabCmd = quote_executable(P.matlab_exe);
    % MATLAB on Windows still writes to the inherited console when -logfile
    % is used.  The short-lived launcher cmd can close that stream and cause
    % fprintf to fail with "iostream stream error".  Redirect stdout/stderr
    % at the shell level so every worker keeps a valid private output stream.
    batchLogCmd = sprintf('%s -wait -batch "run(''%s'')" > "%s" 2>&1', ...
        matlabCmd, char(Plan.wrapper_m(i)), char(Plan.log_file(i)));
    batLines = {
        '@echo off'
        'setlocal'
        sprintf('cd /d "%s"', scriptDir)
        batchLogCmd
        'endlocal'
        };
    write_lines(char(Plan.bat_file(i)), batLines);
end
end

%% ========================================================================
function write_master_launch_bat(Plan, P)

lines = {
    '@echo off'
    'setlocal'
    sprintf('cd /d "%s"', P.worker_script_dir)
    };
for i = 1:height(Plan)
    lines{end + 1, 1} = sprintf('start "ADCIRC_TC_%s" /min cmd /c ""%s""', ...
        char(Plan.worker_id(i)), char(Plan.bat_file(i))); %#ok<AGROW>
end
lines{end + 1, 1} = 'endlocal';
write_lines(P.master_launch_bat, lines);
end

%% ========================================================================
function write_merge_script(Plan, P, scriptDir)

mergeFile = fullfile(P.worker_root, 'merge_parallel_worker_summaries.m');
outCsv = fullfile(P.output_root, 'adcirc_era5_tc_case_build_summary_parallel_merged.csv');
lines = {
    '%% Auto-generated summary merge script.'
    'clearvars; clc;'
    sprintf('output_root = ''%s'';', escape_matlab_string(P.output_root))
    sprintf('out_csv = ''%s'';', escape_matlab_string(outCsv))
    'D = dir(fullfile(output_root, ''adcirc_era5_tc_case_build_summary_worker_*.csv''));'
    'Tables = cell(numel(D), 1);'
    'for i = 1:numel(D)'
    '    Tables{i} = readtable(fullfile(D(i).folder, D(i).name), ''TextType'', ''string'');'
    'end'
    'if isempty(Tables)'
    '    warning(''No worker summary CSV files found under %s'', output_root);'
    'else'
    '    Summary = vertcat(Tables{:});'
    '    writetable(Summary, out_csv);'
    '    fprintf(''Merged %d worker summaries, %d rows:\n  %s\n'', numel(D), height(Summary), out_csv);'
    'end'
    };
write_lines(mergeFile, lines);

% Also store a tiny README next to the plan.
readmeFile = fullfile(P.worker_root, 'README_parallel_workers.txt');
readme = {
    'External parallel ERA5 TC forcing build'
    ''
    sprintf('Generated from: %s', scriptDir)
    sprintf('Worker plan: %s', P.worker_plan_csv)
    sprintf('Launch all workers: %s', P.master_launch_bat)
    sprintf('Shared deterministic TC calendar seed: %d', P.random_day_seed)
    'Launch and monitor from the current MATLAB window:'
    '  setenv(''ADCIRC_TC_PARALLEL_LAUNCH'',''1'')'
    '  setenv(''ADCIRC_TC_PARALLEL_MONITOR'',''1'')'
    '  build_adcirc_era5_tc_run_files_parallel_block_workers'
    ''
    sprintf('Monitor existing workers: monitor_adcirc_era5_tc_parallel_workers(''%s'')', P.worker_root)
    sprintf('Merge summaries after workers finish: run %s in MATLAB', mergeFile)
    ''
    'Each worker owns fixed ADCIRC blocks and one private TMD copy.'
    'Do not run two workers with the same TMD root at the same time.'
    };
write_lines(readmeFile, readme);
end

%% ========================================================================
function write_lines(path, lines)

fid = fopen(path, 'wt');
if fid < 0
    error('Cannot write %s', path);
end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
for i = 1:numel(lines)
    fprintf(fid, '%s\n', lines{i});
end
end

%% ========================================================================
function ensure_dir(d)

if exist(d, 'dir') ~= 7
    mkdir(d);
end
end

%% ========================================================================
function s = escape_matlab_string(s)

s = strrep(char(s), '''', '''''');
end

%% ========================================================================
function s = quote_executable(exePath)

exePath = char(exePath);
if contains(exePath, '\') || contains(exePath, '/') || contains(exePath, ' ')
    s = sprintf('"%s"', exePath);
else
    s = exePath;
end
end

%% ========================================================================
function s = numeric_env_value(x)

if isfinite(x)
    s = sprintf('%d', floor(x));
else
    s = 'inf';
end
end

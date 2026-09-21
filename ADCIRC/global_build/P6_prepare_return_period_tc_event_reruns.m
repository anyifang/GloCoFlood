%% P6_prepare_return_period_tc_event_reruns.m
% Build ADCIRC rerun case folders for the events selected by
% plot_nearshore_water_level_return_period.m.
%
% Output layout:
%   return_period_tc_event_reruns/
%     050yr/ADC_WNP_01/track_xxx/fort.13 ...
%     100yr/ADC_WNP_01/track_xxx/fort.13 ...
%
% The script reads:
%   nearshore_water_level_return_period_plots/
%     block_return_period_tc_selection/*_050deg_050yr_controlling_tc_cells.csv
% Case inputs may be stored either as extracted directories under runfile/ or
% as runfile/<block_id>.zip.  For ZIP sources, only the selected track case
% files are extracted; the full archive is never unpacked.
%
% It de-duplicates events within each return-period/block pair because many
% 0.25 degree cells can be controlled by the same typhoon event.

clearvars;
clc;

SCRIPT_DIR = fileparts(mfilename('fullpath'));
if isempty(SCRIPT_DIR)
    SCRIPT_DIR = pwd;
end

P = struct();
P.return_period_years = env_number_list('PREPARE_RP_RERUN_RPS', [100 200 500]); % 100 200 500
P.block_patterns = env_string_list('PREPARE_RP_RERUN_BLOCK_PATTERN', "ADC_*"); 
P.selection_grid_tag = env_grid_tag('PREPARE_RP_RERUN_GRID_TAG', "050deg");
P.runfile_root = env_path('PREPARE_RP_RERUN_RUNFILE_ROOT', ...
    fullfile(SCRIPT_DIR, 'output', 'adcirc_era5_tc_selected_run_files'));
P.selection_dir = first_existing_dir({ ...
    env_path('PREPARE_RP_RERUN_SELECTION_DIR', ''), ...
    fullfile(SCRIPT_DIR, 'nearshore_water_level_return_period_plots', 'block_return_period_tc_selection'), ...
    fullfile(SCRIPT_DIR, 'read_result_nearshore_water_level', 'block_return_period_tc_selection')}, ...
    fullfile(SCRIPT_DIR, 'nearshore_water_level_return_period_plots', 'block_return_period_tc_selection'));
if isfolder(P.runfile_root)
    defaultOutputRoot = fullfile(P.runfile_root, 'return_period_tc_event_reruns');
else
    defaultOutputRoot = fullfile(SCRIPT_DIR, 'return_period_tc_event_reruns');
end
P.output_root = env_path('PREPARE_RP_RERUN_OUTPUT_ROOT', defaultOutputRoot);
P.submit_template = first_existing_file({ ...
    env_path('PREPARE_RP_RERUN_SUBMIT_TEMPLATE', ''), ...
    fullfile(SCRIPT_DIR, 'sub_intel.sh'), ...
    fullfile(P.runfile_root, 'sub_intel.sh')}, ...
    fullfile(SCRIPT_DIR, 'sub_intel.sh'));

% Case input files to copy into every selected rerun directory.
P.required_case_files = {'fort.13', 'fort.14', 'fort.15', 'fort.19', 'fort.22', 'sub_intel.sh'};
P.optional_case_files = {'fort22_meta.txt', 'storm_track_forcing_meta.mat'};

% Keep this false for incremental runs. Missing files are still copied.
P.overwrite_existing_files = true;

% Patch fort.15 so reruns write fort.63 for the requested output window.
% NSPOOLGE is preserved unless P.force_nspoolge is set.
P.patch_fort15_for_fort63_window = true;
P.fort63_start_day = 2.5;
P.fort63_end_day = 5.0;
P.force_nspoolge = [];  % [] means preserve the value already in fort.15.

% Set true to avoid filesystem writes and print planned copies.
P.dry_run = env_flag('PREPARE_RP_RERUN_DRY_RUN', false);
P.overwrite_existing_files = env_flag('PREPARE_RP_RERUN_OVERWRITE', P.overwrite_existing_files);
P.print_dry_run_case_paths = env_flag('PREPARE_RP_RERUN_PRINT_CASES', false);

fprintf('Preparing return-period ADCIRC rerun cases...\n');
fprintf('Selection CSV directory: %s\n', P.selection_dir);
fprintf('Runfile root: %s\n', P.runfile_root);
fprintf('Output root: %s\n', P.output_root);
fprintf('Return periods: %s yr\n', strjoin(string(P.return_period_years), ', '));
fprintf('Block pattern(s): %s\n', strjoin(P.block_patterns, ', '));
fprintf('Selection grid tag: %s\n', P.selection_grid_tag);

if ~isfolder(P.selection_dir)
    error('Selection directory does not exist: %s', P.selection_dir);
end
if ~P.dry_run
    ensure_dir(P.output_root);
end

allRows = empty_manifest_rows();

for rp = P.return_period_years(:).'
    rpTag = return_period_tag(rp);
    pattern = sprintf('ADC_*_%s_%s_controlling_tc_cells.csv', char(P.selection_grid_tag), rpTag);
    selectionFiles = dir(fullfile(P.selection_dir, pattern));
    selectionFiles = filter_selection_files(selectionFiles, P.block_patterns);

    fprintf('\nReturn period %g yr (%s): found %d selection CSV files.\n', ...
        rp, rpTag, numel(selectionFiles));

    rpRows = empty_manifest_rows();
    for iFile = 1:numel(selectionFiles)
        csvFile = fullfile(selectionFiles(iFile).folder, selectionFiles(iFile).name);
        T = readtable(csvFile, 'TextType', 'string');
        if isempty(T)
            fprintf('  %s: empty, skipped.\n', selectionFiles(iFile).name);
            continue;
        end

        EventTable = summarize_selected_events(T, rp, csvFile);
        fprintf('  %s: %d cells -> %d unique events.\n', ...
            selectionFiles(iFile).name, height(T), height(EventTable));

        for iEvent = 1:height(EventTable)
            row = prepare_one_event(EventTable(iEvent, :), rp, rpTag, P, SCRIPT_DIR);
            rpRows(end + 1, 1) = row; %#ok<SAGROW>
            allRows(end + 1, 1) = row; %#ok<SAGROW>
        end
    end

    if ~isempty(rpRows)
        rpTable = struct2table(rpRows);
        rpManifest = fullfile(P.output_root, sprintf('%s_selected_tc_events_manifest.csv', rpTag));
        if ~P.dry_run
            writetable(rpTable, rpManifest);
            fprintf('  Wrote RP manifest: %s\n', rpManifest);
        else
            fprintf('  DRY RUN: would write RP manifest: %s\n', rpManifest);
        end
    end
end

if ~isempty(allRows)
    allTable = struct2table(allRows);
    allManifest = fullfile(P.output_root, 'all_selected_tc_events_manifest.csv');
    if ~P.dry_run
        writetable(allTable, allManifest);
        write_block_manifests(allTable, P.output_root);
        write_submit_helper(P, unique(string(allTable.rp_tag), 'stable'));
        write_readme(P, allTable);
    end

    fprintf('\nDone. Selected %d unique return-period/block events in total.\n', height(allTable));
    print_copy_status_summary(allTable);
    if P.dry_run
        fprintf('DRY RUN: would write main manifest: %s\n', allManifest);
        fprintf('DRY RUN: would write submit helper: %s\n', fullfile(P.output_root, 'P7_submit_selected_return_period_tc_events.sh'));
    else
        fprintf('Main manifest: %s\n', allManifest);
        fprintf('Submit helper: %s\n', fullfile(P.output_root, 'P7_submit_selected_return_period_tc_events.sh'));
    end
else
    fprintf('\nNo events were selected. Check that the selection CSV files exist.\n');
end

function EventTable = summarize_selected_events(T, rp, csvFile)
requiredVars = ["block_id", "selected_track_case"];
for v = requiredVars
    if ~ismember(v, string(T.Properties.VariableNames))
        error('Missing required column "%s" in %s', v, csvFile);
    end
end

if ismember("return_period_yr", string(T.Properties.VariableNames))
    keep = abs(T.return_period_yr - rp) < 1e-6;
    T = T(keep, :);
end

if isempty(T)
    EventTable = table();
    return;
end

blockId = string(T.block_id);
trackCase = string(T.selected_track_case);
eventKey = blockId + "|" + trackCase;
[G, keys] = findgroups(eventKey);

rows = repmat(empty_event_row(), numel(keys), 1);
for i = 1:numel(keys)
    idx = find(G == i);
    firstIdx = idx(1);

    rows(i).block_id = char(blockId(firstIdx));
    rows(i).selected_track_case = char(trackCase(firstIdx));
    rows(i).selected_track_index = get_numeric_value(T, "selected_track_index", firstIdx, NaN);
    rows(i).selected_event_index = get_numeric_value(T, "selected_event_index", firstIdx, NaN);
    rows(i).selected_event_rank = get_numeric_value(T, "selected_event_rank", firstIdx, NaN);
    rows(i).selected_event_return_period_yr = get_numeric_value(T, "selected_event_return_period_yr", firstIdx, NaN);
    rows(i).selected_event_peak_water_level_m = max(get_numeric_vector(T, "selected_event_peak_water_level_m", idx), [], 'omitnan');
    rows(i).controlled_cell_count = numel(idx);
    rows(i).max_cell_return_level_m = max(cell_return_level_values(T, idx), [], 'omitnan');
    rows(i).selection_csv = char(csvFile);
    rows(i).selected_maxele_file = char(get_string_value(T, "selected_maxele_file", firstIdx, ""));
end

EventTable = struct2table(rows);
end

function files = filter_selection_files(files, blockPatterns)
if isempty(files)
    return;
end

keep = false(numel(files), 1);
for i = 1:numel(files)
    blockId = block_id_from_selection_file(files(i).name);
    keep(i) = any(wildcard_match(blockId, blockPatterns));
end
files = files(keep);
end

function blockId = block_id_from_selection_file(fileName)
token = regexp(fileName, '^(ADC_.+)_[0-9]+(?:p[0-9]+)?deg_', 'tokens', 'once');
if isempty(token)
    blockId = string(fileName);
else
    blockId = string(token{1});
end
end

function tf = wildcard_match(textValue, patterns)
textValue = char(textValue);
patterns = string(patterns);
tf = false(size(patterns));
for i = 1:numel(patterns)
    expr = ['^' regexptranslate('wildcard', char(patterns(i))) '$'];
    tf(i) = ~isempty(regexp(textValue, expr, 'once'));
end
end

function row = prepare_one_event(E, rp, rpTag, P, scriptDir)
blockId = char(E.block_id);
trackCase = char(E.selected_track_case);
selectedMaxeleFile = char(E.selected_maxele_file);

[sourcePath, sourceType, zipCasePrefix, sourceStatus] = locate_case_source( ...
    scriptDir, P.runfile_root, blockId, trackCase, selectedMaxeleFile);
destCaseDir = fullfile(P.output_root, rpTag, blockId, trackCase);

copyStatus = "not_started";
missingRequired = strings(0, 1);
copiedFiles = 0;
skippedFiles = 0;
patchedFort15 = false;

if sourceStatus ~= "found"
    copyStatus = "missing_source_case";
else
    if sourceType == "zip"
        [copyStatus, copiedFiles, skippedFiles, missingRequired] = copy_case_files_from_zip( ...
            sourcePath, zipCasePrefix, destCaseDir, P);
    else
        [copyStatus, copiedFiles, skippedFiles, missingRequired] = copy_case_files( ...
            sourcePath, destCaseDir, P);
    end

    if copyStatus ~= "dry_run" && P.patch_fort15_for_fort63_window
        fort15Path = fullfile(destCaseDir, 'fort.15');
        if isfile(fort15Path)
            patchedFort15 = patch_fort15_for_fort63_window(fort15Path, P);
        end
    end
end

row = empty_manifest_row();
row.return_period_yr = rp;
row.rp_tag = char(rpTag);
row.block_id = blockId;
row.selected_track_case = trackCase;
row.selected_track_index = E.selected_track_index;
row.selected_event_index = E.selected_event_index;
row.selected_event_rank = E.selected_event_rank;
row.selected_event_return_period_yr = E.selected_event_return_period_yr;
row.selected_event_peak_water_level_m = E.selected_event_peak_water_level_m;
row.controlled_cell_count = E.controlled_cell_count;
row.max_cell_return_level_m = E.max_cell_return_level_m;
if sourceType == "zip"
    row.source_case_dir = sprintf('%s::%s', sourcePath, zipCasePrefix);
else
    row.source_case_dir = char(sourcePath);
end
row.dest_case_dir = char(destCaseDir);
row.selected_maxele_file = selectedMaxeleFile;
row.selection_csv = char(E.selection_csv);
row.copy_status = char(copyStatus);
row.copied_file_count = copiedFiles;
row.skipped_file_count = skippedFiles;
row.missing_required_files = char(strjoin(missingRequired, ';'));
row.fort15_patched_fort63_window = patchedFort15;
end

function [copyStatus, copiedFiles, skippedFiles, missingRequired] = copy_case_files_from_zip(zipPath, casePrefix, destDir, P)
% Windows tar/libarchive handles some ZIP central directories that MATLAB's
% Java runtime rejects.  It also avoids MATLAB-to-Java byte-array copy
% semantics, which can silently write zero-filled files.
if ispc
    [copyStatus, copiedFiles, skippedFiles, missingRequired] = ...
        copy_case_files_from_zip_with_tar(zipPath, casePrefix, destDir, P);
    return;
end

try
    zipArchive = java.util.zip.ZipFile(java.io.File(char(zipPath)));
catch javaError
    if ispc
        warn_java_zip_fallback_once(zipPath, javaError.message);
        [copyStatus, copiedFiles, skippedFiles, missingRequired] = ...
            copy_case_files_from_zip_with_tar(zipPath, casePrefix, destDir, P);
        return;
    end
    rethrow(javaError);
end

copyStatus = "ok";
copiedFiles = 0;
skippedFiles = 0;
missingRequired = strings(0, 1);

zipCleanup = onCleanup(@() zipArchive.close());
allFiles = [string(P.required_case_files), string(P.optional_case_files)];

for fileName = allFiles
    destFile = fullfile(destDir, fileName);
    isRequired = any(fileName == string(P.required_case_files));
    useSubmitTemplate = fileName == "sub_intel.sh" && isfile(P.submit_template);

    if useSubmitTemplate
        sourceExists = true;
        zipEntry = [];
    else
        entryName = char(string(casePrefix) + "/" + fileName);
        zipEntry = zipArchive.getEntry(entryName);
        sourceExists = ~isempty(zipEntry) && ~zipEntry.isDirectory();
    end

    if ~sourceExists
        if isRequired
            missingRequired(end + 1, 1) = fileName; %#ok<AGROW>
        end
        continue;
    end

    if isfile(destFile) && ~P.overwrite_existing_files
        skippedFiles = skippedFiles + 1;
        continue;
    end

    if P.dry_run
        copiedFiles = copiedFiles + 1;
        continue;
    end

    ensure_dir(destDir);
    if useSubmitTemplate
        [ok, msg] = copyfile(P.submit_template, destFile);
        if ~ok
            error('Failed to copy %s to %s: %s', P.submit_template, destFile, msg);
        end
    else
        extract_one_zip_entry(zipArchive, zipEntry, destFile);
    end
    copiedFiles = copiedFiles + 1;
end

if ~isempty(missingRequired)
    copyStatus = "missing_required_files";
elseif P.dry_run
    copyStatus = "dry_run";
elseif copiedFiles == 0 && skippedFiles > 0
    copyStatus = "already_exists";
end

if P.dry_run && P.print_dry_run_case_paths
    fprintf('    DRY RUN ZIP: %s::%s -> %s (%d files)\n', ...
        zipPath, casePrefix, destDir, copiedFiles);
end
clear zipCleanup;
end

function [copyStatus, copiedFiles, skippedFiles, missingRequired] = copy_case_files_from_zip_with_tar(zipPath, casePrefix, destDir, P)
copyStatus = "ok";
copiedFiles = 0;
skippedFiles = 0;
missingRequired = strings(0, 1);
repairExistingCase = case_files_need_repair(destDir);
if repairExistingCase
    warning('Existing case contains unreadable or zero-filled files and will be repaired: %s', destDir);
end

assert_safe_command_path(zipPath);
assert_safe_command_path(casePrefix);
try
    archiveEntries = cached_tar_archive_entries(zipPath);
catch listError
    warning('Could not list %s using tar: %s', zipPath, listError.message);
    copyStatus = "missing_source_case";
    return;
end
allFiles = [string(P.required_case_files), string(P.optional_case_files)];
filesToExtract = strings(0, 1);
entryNames = strings(0, 1);

for fileName = allFiles
    destFile = fullfile(destDir, fileName);
    isRequired = any(fileName == string(P.required_case_files));
    useSubmitTemplate = fileName == "sub_intel.sh" && isfile(P.submit_template);
    entryName = string(casePrefix) + "/" + fileName;
    sourceExists = useSubmitTemplate || any(archiveEntries == entryName);

    if ~sourceExists
        if isRequired
            missingRequired(end + 1, 1) = fileName; %#ok<AGROW>
        end
        continue;
    end
    if isfile(destFile) && ~P.overwrite_existing_files && ~repairExistingCase
        skippedFiles = skippedFiles + 1;
        continue;
    end
    if useSubmitTemplate
        if ~P.dry_run
            ensure_dir(destDir);
            [ok, msg] = copyfile(P.submit_template, destFile);
            if ~ok
                error('Failed to copy %s to %s: %s', P.submit_template, destFile, msg);
            end
        end
        copiedFiles = copiedFiles + 1;
    else
        filesToExtract(end + 1, 1) = fileName; %#ok<AGROW>
        entryNames(end + 1, 1) = entryName; %#ok<AGROW>
    end
end

if ~P.dry_run && ~isempty(entryNames)
    stagingDir = tempname;
    ensure_dir(stagingDir);
    stagingCleanup = onCleanup(@() remove_temp_dir(stagingDir));
    assert_safe_command_path(stagingDir);
    quotedEntries = compose('"%s"', entryNames);
    extractCommand = sprintf('tar -xf "%s" -C "%s" %s', ...
        char(zipPath), stagingDir, strjoin(quotedEntries, ' '));
    [extractStatus, extractOutput] = system(extractCommand);
    if extractStatus ~= 0
        error('Selective tar extraction failed for %s::%s: %s', ...
            zipPath, casePrefix, strtrim(extractOutput));
    end

    ensure_dir(destDir);
    stagedCaseDir = fullfile(stagingDir, char(replace(string(casePrefix), '/', filesep)));
    for i = 1:numel(filesToExtract)
        sourceFile = fullfile(stagedCaseDir, filesToExtract(i));
        destFile = fullfile(destDir, filesToExtract(i));
        assert_extracted_file_usable(sourceFile, filesToExtract(i));
        [ok, msg] = copyfile(sourceFile, destFile);
        if ~ok
            error('Failed to move selectively extracted file %s to %s: %s', sourceFile, destFile, msg);
        end
    end
    clear stagingCleanup;
end
copiedFiles = copiedFiles + numel(filesToExtract);

if ~isempty(missingRequired)
    copyStatus = "missing_required_files";
elseif P.dry_run
    copyStatus = "dry_run";
elseif copiedFiles == 0 && skippedFiles > 0
    copyStatus = "already_exists";
end

if P.dry_run && P.print_dry_run_case_paths
    fprintf('    DRY RUN ZIP (tar): %s::%s -> %s (%d files)\n', ...
        zipPath, casePrefix, destDir, copiedFiles);
end
end

function needsRepair = case_files_need_repair(caseDir)
needsRepair = false;
if ~isfolder(caseDir)
    return;
end

filesToCheck = ["fort.13", "fort.14", "fort.15", "fort.19", "fort.22", ...
    "fort22_meta.txt", "storm_track_forcing_meta.mat"];
for fileName = filesToCheck
    filePath = fullfile(caseDir, fileName);
    if ~isfile(filePath)
        continue;
    end
    if file_begins_with_zeros(filePath)
        needsRepair = true;
        return;
    end
    if fileName == "storm_track_forcing_meta.mat"
        try
            S = load(filePath, 'tc');
            if ~isfield(S, 'tc')
                needsRepair = true;
                return;
            end
        catch
            needsRepair = true;
            return;
        end
    end
end
end

function assert_extracted_file_usable(filePath, fileName)
if ~isfile(filePath) || dir(filePath).bytes == 0
    error('Selective extraction produced a missing or empty file: %s', filePath);
end
if file_begins_with_zeros(filePath)
    error('Selective extraction produced a zero-filled file: %s', filePath);
end
if string(fileName) == "storm_track_forcing_meta.mat"
    try
        S = load(filePath, 'tc');
    catch loadError
        error('Extracted MAT file is unreadable (%s): %s', filePath, loadError.message);
    end
    if ~isfield(S, 'tc')
        error('Extracted MAT file lacks variable tc: %s', filePath);
    end
end
end

function tf = file_begins_with_zeros(filePath)
fid = fopen(filePath, 'r');
if fid < 0
    tf = true;
    return;
end
fileCleanup = onCleanup(@() fclose(fid));
bytes = fread(fid, 64, '*uint8');
tf = ~isempty(bytes) && all(bytes == 0);
clear fileCleanup;
end

function entries = cached_tar_archive_entries(zipPath)
persistent entryCache
if isempty(entryCache)
    entryCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end
cacheKey = char(zipPath);
if isKey(entryCache, cacheKey)
    entries = entryCache(cacheKey);
    return;
end

assert_safe_command_path(zipPath);
listCommand = sprintf('tar -tf "%s"', char(zipPath));
[listStatus, listOutput] = system(listCommand);
if listStatus ~= 0
    error('tar returned status %d: %s', listStatus, strtrim(listOutput));
end
entries = strip(splitlines(string(listOutput)));
entries(entries == "") = [];
entryCache(cacheKey) = entries;
end

function warn_java_zip_fallback_once(zipPath, errorMessage)
persistent warnedPaths
if isempty(warnedPaths)
    warnedPaths = containers.Map('KeyType', 'char', 'ValueType', 'logical');
end
warningKey = char(zipPath);
if ~isKey(warnedPaths, warningKey)
    warning('MATLAB Java could not open %s (%s). Falling back to selective extraction with Windows tar.', ...
        zipPath, errorMessage);
    warnedPaths(warningKey) = true;
end
end

function assert_safe_command_path(pathValue)
pathValue = char(pathValue);
if contains(pathValue, '"') || contains(pathValue, newline) || contains(pathValue, char(13))
    error('ZIP/tar path contains unsupported command characters: %s', pathValue);
end
end

function remove_temp_dir(dirName)
if isfolder(dirName)
    rmdir(dirName, 's');
end
end

function extract_one_zip_entry(zipArchive, zipEntry, destFile)
inputStream = zipArchive.getInputStream(zipEntry);
outputStream = java.io.FileOutputStream(java.io.File(char(destFile)));
streamCleanup = onCleanup(@() close_zip_streams(inputStream, outputStream));
buffer = zeros(1, 1024 * 1024, 'int8');

while true
    count = inputStream.read(buffer, 0, numel(buffer));
    if count < 0
        break;
    end
    outputStream.write(buffer, 0, count);
end

outputStream.flush();
clear streamCleanup;
end

function close_zip_streams(inputStream, outputStream)
try
    inputStream.close();
catch
end
try
    outputStream.close();
catch
end
end

function [copyStatus, copiedFiles, skippedFiles, missingRequired] = copy_case_files(sourceDir, destDir, P)
copyStatus = "ok";
copiedFiles = 0;
skippedFiles = 0;
missingRequired = strings(0, 1);

if P.dry_run
    if P.print_dry_run_case_paths
        fprintf('    DRY RUN: %s -> %s\n', sourceDir, destDir);
    end
    copyStatus = "dry_run";
    return;
end

ensure_dir(destDir);

allFiles = [string(P.required_case_files), string(P.optional_case_files)];
for fileName = allFiles
    if fileName == "sub_intel.sh" && isfile(P.submit_template)
        sourceFile = P.submit_template;
    else
        sourceFile = fullfile(sourceDir, fileName);
    end
    destFile = fullfile(destDir, fileName);
    isRequired = any(fileName == string(P.required_case_files));

    if ~isfile(sourceFile)
        if isRequired
            missingRequired(end + 1, 1) = fileName; %#ok<AGROW>
        end
        continue;
    end

    if isfile(destFile) && ~P.overwrite_existing_files
        skippedFiles = skippedFiles + 1;
        continue;
    end

    [ok, msg] = copyfile(sourceFile, destFile);
    if ~ok
        error('Failed to copy %s to %s: %s', sourceFile, destFile, msg);
    end
    copiedFiles = copiedFiles + 1;
end

if ~isempty(missingRequired)
    copyStatus = "missing_required_files";
elseif copiedFiles == 0 && skippedFiles > 0
    copyStatus = "already_exists";
end
end

function patched = patch_fort15_for_fort63_window(fort15Path, P)
patched = false;
lines = readlines(fort15Path);

rnday = parse_rnday(lines);
if ~isfinite(rnday) || rnday <= 0
    warning('Could not parse RNDAY from %s. fort.15 was not patched.', fort15Path);
    return;
end

iOut = find(contains(lines, 'NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE'), 1, 'first');
if isempty(iOut)
    warning('Could not find NOUTGE line in %s. fort.15 was not patched.', fort15Path);
    return;
end

oldLine = char(lines(iOut));
numbers = regexp(extract_before_comment(oldLine), number_pattern(), 'match');
if numel(numbers) >= 4
    nspoolge = str2double(numbers{4});
else
    nspoolge = 720;
end
if ~isempty(P.force_nspoolge)
    nspoolge = P.force_nspoolge;
end

startDay = double(P.fort63_start_day);
endDay = double(P.fort63_end_day);
if ~isfinite(startDay) || ~isfinite(endDay) || startDay < 0 || endDay <= startDay
    warning('Invalid fort.63 output window %.6g to %.6g days. fort.15 was not patched.', startDay, endDay);
    return;
end
if endDay > rnday + 1e-9
    warning('Requested fort.63 output end day %.6g exceeds RNDAY %.6g in %s. fort.15 was not patched.', ...
        endDay, rnday, fort15Path);
    return;
end

newLine = sprintf(' 1 %.6g %.6g %d ! NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE : GLOBAL ELEVATION OUTPUT INFO (UNIT  63)', ...
    startDay, endDay, round(nspoolge));

if string(oldLine) ~= string(newLine)
    lines(iOut) = string(newLine);
    writelines(lines, fort15Path);
    patched = true;
end
end

function rnday = parse_rnday(lines)
rnday = NaN;
iRnday = find(contains(lines, 'RNDAY'), 1, 'first');
if isempty(iRnday)
    return;
end
numbers = regexp(extract_before_comment(char(lines(iRnday))), number_pattern(), 'match');
if ~isempty(numbers)
    rnday = str2double(numbers{1});
end
end

function [sourcePath, sourceType, zipCasePrefix, status] = locate_case_source(scriptDir, runfileRoot, blockId, trackCase, selectedMaxeleFile)
candidates = strings(0, 1);
candidates(end + 1, 1) = string(fullfile(runfileRoot, blockId, trackCase));
candidates(end + 1, 1) = string(fullfile(scriptDir, blockId, trackCase));
candidates(end + 1, 1) = string(fullfile(scriptDir, 'ERA5_prerun', blockId, trackCase));

if strlength(string(selectedMaxeleFile)) > 0
    maxeleDir = string(fileparts(selectedMaxeleFile));
    candidates(end + 1, 1) = maxeleDir;

    marker = [filesep 'out_maxele' filesep];
    maxelePath = char(selectedMaxeleFile);
    markerIdx = strfind(maxelePath, marker);
    if ~isempty(markerIdx)
        prefix = maxelePath(1:markerIdx(1) - 1);
        suffix = maxelePath(markerIdx(1) + numel(marker):end);
        suffixDir = fileparts(suffix);
        candidates(end + 1, 1) = string(fullfile(prefix, 'runfile', suffixDir));
        candidates(end + 1, 1) = string(fullfile(prefix, suffixDir));
    end
end

candidates = unique(candidates, 'stable');
for i = 1:numel(candidates)
    candidate = candidates(i);
    if isfolder(candidate) && isfile(fullfile(candidate, 'fort.15')) && isfile(fullfile(candidate, 'fort.22'))
        sourcePath = candidate;
        sourceType = "directory";
        zipCasePrefix = "";
        status = "found";
        return;
    end
end

zipCandidates = unique([ ...
    string(fullfile(runfileRoot, [blockId '.zip'])); ...
    string(fullfile(scriptDir, [blockId '.zip'])); ...
    string(fullfile(scriptDir, 'ERA5_prerun', [blockId '.zip']))], 'stable');
prefixCandidates = [string(blockId) + "/" + string(trackCase), string(trackCase)];
for i = 1:numel(zipCandidates)
    zipCandidate = zipCandidates(i);
    if ~isfile(zipCandidate)
        continue;
    end
    try
        zipArchive = java.util.zip.ZipFile(java.io.File(char(zipCandidate)));
    catch javaError
        if ispc
            % Some valid ZIP archives are accepted by Windows tar/.NET but
            % rejected by MATLAB's stricter Java ZIP parser.  The copy step
            % validates and selectively extracts this exact case with tar.
            sourcePath = zipCandidate;
            sourceType = "zip";
            zipCasePrefix = prefixCandidates(1);
            status = "found";
            return;
        end
        rethrow(javaError);
    end
    zipCleanup = onCleanup(@() zipArchive.close());
    for j = 1:numel(prefixCandidates)
        prefix = prefixCandidates(j);
        fort15Entry = zipArchive.getEntry(char(prefix + "/fort.15"));
        fort22Entry = zipArchive.getEntry(char(prefix + "/fort.22"));
        if ~isempty(fort15Entry) && ~isempty(fort22Entry)
            sourcePath = zipCandidate;
            sourceType = "zip";
            zipCasePrefix = prefix;
            status = "found";
            clear zipCleanup;
            return;
        end
    end
    clear zipCleanup;
end

sourcePath = candidates(1);
sourceType = "missing";
zipCasePrefix = "";
status = "missing";
end

function values = cell_return_level_values(T, idx)
if ismember("water_level_m", string(T.Properties.VariableNames))
    values = T.water_level_m(idx);
elseif ismember("water_level_100yr_m", string(T.Properties.VariableNames))
    values = T.water_level_100yr_m(idx);
else
    values = NaN(size(idx));
end
end

function value = get_numeric_value(T, varName, idx, defaultValue)
if ismember(varName, string(T.Properties.VariableNames))
    col = T.(varName);
    if isnumeric(col)
        value = col(idx);
    else
        value = str2double(string(col(idx)));
    end
else
    value = defaultValue;
end
end

function values = get_numeric_vector(T, varName, idx)
if ismember(varName, string(T.Properties.VariableNames))
    col = T.(varName);
    if isnumeric(col)
        values = col(idx);
    else
        values = str2double(string(col(idx)));
    end
else
    values = NaN(size(idx));
end
end

function value = get_string_value(T, varName, idx, defaultValue)
if ismember(varName, string(T.Properties.VariableNames))
    value = string(T.(varName)(idx));
else
    value = string(defaultValue);
end
end

function write_block_manifests(allTable, outputRoot)
blockDir = fullfile(outputRoot, 'block_manifests');
ensure_dir(blockDir);

keys = unique(string(allTable.rp_tag) + "|" + string(allTable.block_id), 'stable');
for i = 1:numel(keys)
    parts = split(keys(i), "|");
    rpTag = parts(1);
    blockId = parts(2);
    keep = string(allTable.rp_tag) == rpTag & string(allTable.block_id) == blockId;
    blockTable = allTable(keep, :);
    outFile = fullfile(blockDir, sprintf('%s_%s_selected_tc_events_manifest.csv', rpTag, blockId));
    writetable(blockTable, outFile);
end
end

function print_copy_status_summary(allTable)
statuses = string(allTable.copy_status);
[G, names] = findgroups(statuses);
counts = splitapply(@sum, ones(size(statuses)), G);
fprintf('Copy status summary:\n');
for i = 1:numel(names)
    fprintf('  %s: %d\n', names(i), counts(i));
end
end

function write_submit_helper(P, rpTags)
helperFile = fullfile(P.output_root, 'P7_submit_selected_return_period_tc_events.sh');
lines = strings(0, 1);
lines(end + 1) = "#!/usr/bin/env bash";
lines(end + 1) = "set -euo pipefail";
lines(end + 1) = "";
lines(end + 1) = "BASE_DIR=""$(cd ""$(dirname ""${BASH_SOURCE[0]}"")"" && pwd)""";
lines(end + 1) = "if [[ -f ""$BASE_DIR/../submit_sub_intel_jobs.sh"" ]]; then";
lines(end + 1) = "  WORK_DIR=""$(cd ""$BASE_DIR/.."" && pwd)""";
lines(end + 1) = "elif [[ -f ""$BASE_DIR/../../submit_sub_intel_jobs.sh"" ]]; then";
lines(end + 1) = "  WORK_DIR=""$(cd ""$BASE_DIR/../.."" && pwd)""";
lines(end + 1) = "else";
lines(end + 1) = "  echo ""Could not find submit_sub_intel_jobs.sh from $BASE_DIR"" >&2";
lines(end + 1) = "  exit 2";
lines(end + 1) = "fi";
lines(end + 1) = "BLOCK_PATTERN=""${1:-ADC_*}""";
lines(end + 1) = "MAX_JOBS=""${MAX_JOBS:-490}""";
lines(end + 1) = "MAX_SUBMIT=""${MAX_SUBMIT:-0}""";
lines(end + 1) = "DRY_RUN_ARG=""${DRY_RUN_ARG:-}""";
lines(end + 1) = "";
for i = 1:numel(rpTags)
    rpTag = rpTags(i);
    lines(end + 1) = sprintf('if [[ -d "$BASE_DIR/%s" ]]; then', rpTag);
    lines(end + 1) = sprintf('  bash "$WORK_DIR/submit_sub_intel_jobs.sh" --root "$BASE_DIR/%s" --block-pattern "$BLOCK_PATTERN" --max-jobs "$MAX_JOBS" --max-submit "$MAX_SUBMIT" $DRY_RUN_ARG', rpTag);
    lines(end + 1) = "fi";
end

write_text_lines(helperFile, lines);
end

function write_readme(P, allTable)
readmeFile = fullfile(P.output_root, 'README_selected_event_reruns.txt');
rpTags = unique(string(allTable.rp_tag), 'stable');

lines = strings(0, 1);
lines(end + 1) = "Selected return-period ADCIRC rerun cases";
lines(end + 1) = "";
lines(end + 1) = "Directory layout:";
lines(end + 1) = "  <return_period>/<block_id>/<track_case>/fort.13 fort.14 fort.15 fort.19 fort.22 sub_intel.sh";
lines(end + 1) = "";
lines(end + 1) = "Main manifest:";
lines(end + 1) = "  all_selected_tc_events_manifest.csv";
lines(end + 1) = "";
lines(end + 1) = "Submit all selected events from the cluster login node:";
lines(end + 1) = "  cd <the directory containing this README>";
lines(end + 1) = "  bash P7_submit_selected_return_period_tc_events.sh";
lines(end + 1) = "";
lines(end + 1) = "Submit only WNP blocks:";
lines(end + 1) = "  bash P7_submit_selected_return_period_tc_events.sh 'ADC_WNP_*'";
lines(end + 1) = "";
lines(end + 1) = "Dry run:";
lines(end + 1) = "  DRY_RUN_ARG=--dry-run bash P7_submit_selected_return_period_tc_events.sh 'ADC_WNP_*'";
lines(end + 1) = "";
lines(end + 1) = "Return-period folders written:";
for i = 1:numel(rpTags)
    thisTag = rpTags(i);
    n = nnz(string(allTable.rp_tag) == thisTag);
    lines(end + 1) = sprintf("  %s: %d unique block/event rows", thisTag, n);
end
lines(end + 1) = "";
lines(end + 1) = sprintf("fort.15 fort.63 output-window patch enabled: %d", P.patch_fort15_for_fort63_window);
lines(end + 1) = sprintf("The NOUTGE line is set to write fort.63 from day %.3g to %.3g, preserving NSPOOLGE unless configured otherwise.", ...
    P.fort63_start_day, P.fort63_end_day);

write_text_lines(readmeFile, lines);
end

function write_text_lines(fileName, lines)
fid = fopen(fileName, 'w');
if fid < 0
    error('Could not open file for writing: %s', fileName);
end
cleanupObj = onCleanup(@() fclose(fid));
for i = 1:numel(lines)
    fprintf(fid, '%s\n', lines(i));
end
clear cleanupObj;
end

function s = extract_before_comment(line)
parts = split(string(line), "!");
s = char(parts(1));
end

function p = number_pattern()
p = '[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?';
end

function tag = return_period_tag(rp)
if abs(rp - round(rp)) < 1e-9
    tag = sprintf('%03dyr', round(rp));
else
    tag = sprintf('%06.1fyr', rp);
    tag = strrep(tag, '.', 'p');
end
end

function ensure_dir(pathName)
if ~isfolder(pathName)
    mkdir(pathName);
end
end

function value = env_flag(name, defaultValue)
raw = getenv(name);
if isempty(raw)
    value = defaultValue;
    return;
end
value = any(strcmpi(strtrim(raw), {'1', 'true', 'yes', 'on'}));
end

function value = env_path(name, defaultValue)
raw = getenv(name);
if isempty(raw)
    value = defaultValue;
else
    value = char(string(raw));
end
end

function tag = env_grid_tag(name, defaultValue)
raw = strip(string(getenv(name)));
if raw == ""
    tag = normalize_grid_tag(defaultValue);
    return;
end
tag = normalize_grid_tag(raw);
end

function tag = normalize_grid_tag(raw)
raw = lower(strip(string(raw)));
raw = replace(raw, " ", "");
if ~isempty(regexp(char(raw), '^\d{3}deg$', 'once')) || ...
        ~isempty(regexp(char(raw), '^\d+p\d+deg$', 'once'))
    tag = raw;
    return;
end
if ~isempty(regexp(char(raw), '^\d{3}$', 'once'))
    tag = raw + "deg";
    return;
end

valueText = erase(raw, "deg");
valueText = replace(valueText, "p", ".");
gridDeg = str2double(valueText);
if ~isfinite(gridDeg) || gridDeg <= 0
    error('Invalid %s value "%s". Use tags like 025deg, 050deg, 100deg, or numeric values like 0.25, 0.5, 1.0.', ...
        'PREPARE_RP_RERUN_GRID_TAG', raw);
end
tag = sprintf('%03ddeg', round(gridDeg .* 100));
end

function fileName = first_existing_file(candidates, defaultValue)
fileName = defaultValue;
for i = 1:numel(candidates)
    if isfile(candidates{i})
        fileName = candidates{i};
        return;
    end
end
end

function dirName = first_existing_dir(candidates, defaultValue)
dirName = defaultValue;
for i = 1:numel(candidates)
    if isfolder(candidates{i})
        dirName = candidates{i};
        return;
    end
end
end

function values = env_number_list(name, defaultValue)
raw = getenv(name);
if isempty(raw)
    values = defaultValue;
    return;
end
parts = split(string(raw), {',', ';', ' '});
parts(parts == "") = [];
values = str2double(parts).';
if any(~isfinite(values))
    error('Environment variable %s must be a comma/semicolon/space separated numeric list.', name);
end
end

function values = env_string_list(name, defaultValue)
raw = getenv(name);
if isempty(raw)
    values = string(defaultValue);
    return;
end
values = split(string(raw), {';', ','});
values = strip(values);
values(values == "") = [];
if isempty(values)
    values = string(defaultValue);
end
values = values(:).';
end

function rows = empty_manifest_rows()
rows = repmat(empty_manifest_row(), 0, 1);
end

function row = empty_manifest_row()
row = struct( ...
    'return_period_yr', NaN, ...
    'rp_tag', '', ...
    'block_id', '', ...
    'selected_track_case', '', ...
    'selected_track_index', NaN, ...
    'selected_event_index', NaN, ...
    'selected_event_rank', NaN, ...
    'selected_event_return_period_yr', NaN, ...
    'selected_event_peak_water_level_m', NaN, ...
    'controlled_cell_count', NaN, ...
    'max_cell_return_level_m', NaN, ...
    'source_case_dir', '', ...
    'dest_case_dir', '', ...
    'selected_maxele_file', '', ...
    'selection_csv', '', ...
    'copy_status', '', ...
    'copied_file_count', NaN, ...
    'skipped_file_count', NaN, ...
    'missing_required_files', '', ...
    'fort15_patched_fort63_window', false);
end

function row = empty_event_row()
row = struct( ...
    'block_id', '', ...
    'selected_track_case', '', ...
    'selected_track_index', NaN, ...
    'selected_event_index', NaN, ...
    'selected_event_rank', NaN, ...
    'selected_event_return_period_yr', NaN, ...
    'selected_event_peak_water_level_m', NaN, ...
    'controlled_cell_count', NaN, ...
    'max_cell_return_level_m', NaN, ...
    'selection_csv', '', ...
    'selected_maxele_file', '');
end

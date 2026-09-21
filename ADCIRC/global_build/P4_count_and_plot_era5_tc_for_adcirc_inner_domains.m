%% P4_count_and_plot_era5_tc_for_adcirc_inner_domains.m
% Count ERA5 downscaled TCs that affect each ADCIRC inner-refinement domain,
% then plot the ADCIRC outer/shelf domain, inner refined region, and the
% matching TC tracks for every block.
%
% Event definition follows P3_compare_tc_impact_return_period_era5_future.m
% except that the impact mask is the ADCIRC inner refined polygon:
%   - densify track segments with 6 substeps;
%   - a TC counts once per ADCIRC block when a densified point lies inside
%     the block's inner refined domain and vmax > 33 m/s.

clearvars;
clc;
close all;

SCRIPT_DIR = fileparts(mfilename('fullpath'));
if isempty(SCRIPT_DIR)
    SCRIPT_DIR = pwd;
end

P = struct();
P.block_dir = fullfile(SCRIPT_DIR, 'catalog');
P.blocks_csv = fullfile(P.block_dir, 'global_tc_adcirc_blocks.csv');
P.mesh_root = fullfile(SCRIPT_DIR, 'output', 'adcirc_fort14_meshes');
P.boundary_dir = fullfile(P.block_dir, 'block_boundaries');
P.track_nc = fullfile(SCRIPT_DIR, 'external', 'tracks_GL_era5_197501_201412.nc');

v = strtrim(string(getenv('ADCIRC_BLOCK_CATALOG_ROOT')));
if strlength(v) > 0
    P.block_dir = char(v);
    P.blocks_csv = fullfile(P.block_dir, 'global_tc_adcirc_blocks.csv');
    P.boundary_dir = fullfile(P.block_dir, 'block_boundaries');
end
v = strtrim(string(getenv('ADCIRC_MESH_OUTPUT_ROOT')));
if strlength(v) > 0, P.mesh_root = char(v); end
v = strtrim(string(getenv('ADCIRC_TC_TRACK_NC')));
if strlength(v) > 0, P.track_nc = char(v); end

P.start_year = 1975;
P.end_year = 2014;
P.min_inner_vmax_ms = 33.0;
P.min_record_lmi_ms = 18.0;
P.segment_substeps = 6;
P.nc_track_block_size = 4000;
P.impact_buffer_km = 0.0; % 0 = exact inner refined domain; >0 expands the impact mask in local km coordinates.

P.force_recount = false;
P.make_figures = true;
P.figure_visible = 'off';
P.figure_dpi = 350;
P.save_pdf = false;
P.output_summary_csv = fullfile(P.block_dir, 'era5_tc_count_by_adcirc_inner_domain_vmax33_matlab.csv');
P.output_events_csv = fullfile(P.block_dir, 'era5_tc_events_by_adcirc_inner_domain_vmax33_matlab.csv');
P.output_mat = fullfile(P.block_dir, 'era5_tc_events_by_adcirc_inner_domain_vmax33_matlab.mat');
P.figure_dir = fullfile(P.block_dir, 'era5_tc_track_maps_by_adcirc_inner_domain');

P.font_name = 'Arial';
P.axis_font_size = 8;
P.title_font_size = 12;
P.track_color = [0.08 0.32 0.86];
P.track_line_width = 0.35;
P.track_alpha = 0.16;
P.outer_face_color = [0.80 0.91 0.98];
P.outer_edge_color = [0.05 0.32 0.52];
P.inner_face_color = [1.00 0.73 0.28];
P.inner_edge_color = [0.86 0.25 0.05];
P.land_color = [0.86 0.84 0.77];
P.coast_edge_color = [0.40 0.38 0.32];
P.map_pad_fraction = 0.10;
P.map_min_pad_deg = 0.50;
P.figure_size_cm = [16.0 12.0];

bufferEnv = str2double(strtrim(string(getenv('ADCIRC_TC_IMPACT_BUFFER_KM'))));
if isfinite(bufferEnv) && bufferEnv >= 0
    P.impact_buffer_km = bufferEnv;
end
forceEnv = string(strtrim(getenv('ADCIRC_FORCE_RECOUNT')));
if any(strcmpi(forceEnv, ["1", "true", "yes", "y"]))
    P.force_recount = true;
elseif any(strcmpi(forceEnv, ["0", "false", "no", "n"]))
    P.force_recount = false;
end
figEnv = string(strtrim(getenv('ADCIRC_MAKE_FIGURES')));
if any(strcmpi(figEnv, ["1", "true", "yes", "y"]))
    P.make_figures = true;
elseif any(strcmpi(figEnv, ["0", "false", "no", "n"]))
    P.make_figures = false;
end

if P.impact_buffer_km > 0
    bufferSuffix = sprintf('_buffer%gkm', P.impact_buffer_km);
    bufferSuffix = regexprep(bufferSuffix, '\.', 'p');
    P.output_summary_csv = replace(P.output_summary_csv, '.csv', bufferSuffix + ".csv");
    P.output_events_csv = replace(P.output_events_csv, '.csv', bufferSuffix + ".csv");
    P.output_mat = replace(P.output_mat, '.mat', bufferSuffix + ".mat");
    P.figure_dir = fullfile(P.figure_dir, regexprep(bufferSuffix, '^_', ''));
end

% Optional debug limiter:
%   $env:ADCIRC_TC_MAX_BLOCKS='1'
envMaxBlocks = str2double(strtrim(string(getenv('ADCIRC_TC_MAX_BLOCKS'))));
if isfinite(envMaxBlocks) && envMaxBlocks > 0
    P.max_blocks = floor(envMaxBlocks);
    suffix = sprintf('_first%d', P.max_blocks);
    P.output_summary_csv = replace(P.output_summary_csv, '.csv', suffix + ".csv");
    P.output_events_csv = replace(P.output_events_csv, '.csv', suffix + ".csv");
    P.output_mat = replace(P.output_mat, '.mat', suffix + ".mat");
    P.figure_dir = fullfile(P.figure_dir, sprintf('first_%d_blocks', P.max_blocks));
else
    P.max_blocks = inf;
end

add_project_paths(SCRIPT_DIR);
ensure_dir(P.figure_dir);

assert(exist(P.blocks_csv, 'file') == 2, 'Missing block CSV: %s', P.blocks_csv);
assert(exist(P.track_nc, 'file') == 2, 'Missing ERA5 track NetCDF: %s', P.track_nc);

T_blocks = readtable(P.blocks_csv, 'TextType', 'string');
if isfinite(P.max_blocks) && height(T_blocks) > P.max_blocks
    T_blocks = T_blocks(1:P.max_blocks, :);
end

fprintf('\n============================================================\n');
fprintf('ERA5 TC count for ADCIRC inner domains\n');
fprintf('Blocks    : %d\n', height(T_blocks));
fprintf('Track file: %s\n', P.track_nc);
fprintf('Period    : %d-%d\n', P.start_year, P.end_year);
fprintf('Criterion : densified track inside inner domain + %.1f km buffer and vmax > %.1f m/s\n', ...
    P.impact_buffer_km, P.min_inner_vmax_ms);
fprintf('Outputs   :\n  %s\n  %s\n', P.output_summary_csv, P.output_events_csv);
fprintf('Figures   : %s\n', P.figure_dir);
fprintf('============================================================\n');

Domains = load_adcirc_block_domains(T_blocks, P);

if P.force_recount || exist(P.output_summary_csv, 'file') ~= 2 || exist(P.output_events_csv, 'file') ~= 2
    [Summary, Events] = count_era5_events_for_domains(P, T_blocks, Domains);
    writetable(Summary, P.output_summary_csv);
    writetable(Events, P.output_events_csv);
    save(P.output_mat, 'P', 'T_blocks', 'Domains', 'Summary', 'Events', '-v7.3');
else
    Summary = readtable(P.output_summary_csv, 'TextType', 'string');
    Events = readtable(P.output_events_csv, 'TextType', 'string');
end

if P.make_figures
    FigureIndex = plot_tc_track_maps_for_domains(P, T_blocks, Domains, Summary, Events);
    writetable(FigureIndex, fullfile(P.figure_dir, 'era5_tc_track_map_index.csv'));
end

fprintf('\nDone.\n');

%% ============================================================
function Domains = load_adcirc_block_domains(T_blocks, P)

n = height(T_blocks);
Domains = repmat(empty_domain_struct(), n, 1);
for i = 1:n
    blockId = string(T_blocks.block_id(i));
    innerRing = read_domain_ring(blockId, "inner", P);
    outerRing = read_domain_ring(blockId, "outer", P);
    assert(~isempty(innerRing), 'Missing inner ring for %s.', blockId);
    if isempty(outerRing)
        outerRing = innerRing;
        warning('Missing outer ring for %s; using inner ring for map extent.', blockId);
    end
    Domains(i).block_id = blockId;
    Domains(i).inner_ring = close_ring_if_needed(clean_lonlat_ring(innerRing));
    Domains(i).outer_ring = close_ring_if_needed(clean_lonlat_ring(outerRing));
    Domains(i).inner_bbox = ring_bbox(Domains(i).inner_ring);
    Domains(i).outer_bbox = ring_bbox(Domains(i).outer_ring);
    Domains(i).impact_ring = make_impact_ring(Domains(i).inner_ring, P.impact_buffer_km);
    Domains(i).impact_bbox = ring_bbox(Domains(i).impact_ring);
end
end

%% ============================================================
function D = empty_domain_struct()

D = struct();
D.block_id = "";
D.inner_ring = zeros(0, 2);
D.outer_ring = zeros(0, 2);
D.impact_ring = zeros(0, 2);
D.inner_bbox = [NaN NaN NaN NaN];
D.outer_bbox = [NaN NaN NaN NaN];
D.impact_bbox = [NaN NaN NaN NaN];
end

%% ============================================================
function ring = read_domain_ring(blockId, kind, P)

blockId = string(blockId);
kind = lower(string(kind));
paths = strings(0, 1);
if kind == "inner"
    paths(end + 1) = fullfile(P.mesh_root, char(blockId), 'inner_ring_ocean_split_used.csv');
    paths(end + 1) = fullfile(P.mesh_root, char(blockId), 'inner_ring.csv');
    paths(end + 1) = fullfile(P.boundary_dir, char(blockId + "_inner_ring.csv"));
else
    paths(end + 1) = fullfile(P.mesh_root, char(blockId), 'outer_ring.csv');
    paths(end + 1) = fullfile(P.boundary_dir, char(blockId + "_outer_ring.csv"));
end

ring = [];
for ip = 1:numel(paths)
    f = char(paths(ip));
    if exist(f, 'file') ~= 2
        continue;
    end
    T = readtable(f, 'TextType', 'string');
    vars = lower(string(T.Properties.VariableNames));
    ilon = find(vars == "lon" | vars == "longitude" | vars == "x", 1, 'first');
    ilat = find(vars == "lat" | vars == "latitude" | vars == "y", 1, 'first');
    if isempty(ilon) || isempty(ilat)
        continue;
    end
    ring = [double(T{:, ilon}), double(T{:, ilat})];
    return;
end
end

%% ============================================================
function ring = clean_lonlat_ring(ring)

ring = double(ring);
if size(ring, 2) < 2 && size(ring, 1) >= 2
    ring = ring.';
end
ring = ring(:, 1:2);
ok = all(isfinite(ring), 2);
ring = ring(ok, :);
ring(:, 1) = wrap_to_180_local(ring(:, 1));
end

%% ============================================================
function ring = close_ring_if_needed(ring)

if size(ring, 1) < 2
    return;
end
if any(abs(ring(1, :) - ring(end, :)) > 1e-10)
    ring(end + 1, :) = ring(1, :);
end
end

%% ============================================================
function bbox = ring_bbox(ring)

ring = double(ring);
if isempty(ring) || size(ring, 2) < 2
    bbox = [NaN NaN NaN NaN];
    return;
end
ok = isfinite(ring(:, 1)) & isfinite(ring(:, 2));
if ~any(ok)
    bbox = [NaN NaN NaN NaN];
else
    bbox = [min(ring(ok, 1)), max(ring(ok, 1)), min(ring(ok, 2)), max(ring(ok, 2))];
end
end

%% ============================================================
function impactRing = make_impact_ring(innerRing, bufferKm)

innerRing = close_ring_if_needed(clean_lonlat_ring(innerRing));
if ~isfinite(bufferKm) || bufferKm <= 0
    impactRing = innerRing;
    return;
end

ok = isfinite(innerRing(:, 1)) & isfinite(innerRing(:, 2));
ring = innerRing(ok, :);
if size(ring, 1) < 3
    impactRing = innerRing;
    return;
end

lon0 = mean(ring(:, 1), 'omitnan');
lat0 = mean(ring(:, 2), 'omitnan');
cosLat = max(cosd(lat0), 0.15);
x = wrap_to_180_local(ring(:, 1) - lon0) .* 111.32 .* cosLat;
y = (ring(:, 2) - lat0) .* 111.32;

try
    p = polyshape(x, y, 'Simplify', true);
    if area(p) <= 0
        impactRing = innerRing;
        return;
    end
    p = polybuffer(p, double(bufferKm));
    p = rmholes(p);
    [xb, yb] = boundary(p);
    if iscell(xb)
        areas = cellfun(@(xx, yy) abs(polyarea(xx, yy)), xb, yb);
        [~, imax] = max(areas);
        xb = xb{imax};
        yb = yb{imax};
    end
    lon = wrap_to_180_local(lon0 + xb(:) ./ (111.32 .* cosLat));
    lat = lat0 + yb(:) ./ 111.32;
    impactRing = close_ring_if_needed([lon, lat]);
catch ME
    warning('Impact buffer failed; using unbuffered inner ring: %s', ME.message);
    impactRing = innerRing;
end
end

%% ============================================================
function [Summary, Events] = count_era5_events_for_domains(P, T_blocks, Domains)

yearsAll = clean_fill_to_nan(double(ncread(P.track_nc, 'tc_years')));
yearsAll = yearsAll(:);
nTrack = numel(yearsAll);
selected = find(yearsAll >= P.start_year & yearsAll <= P.end_year);
assert(~isempty(selected), 'No tracks in requested period.');

[trackDim, nTime] = infer_track_dim_from_nc(P.track_nc, 'lon_trks', nTrack);
blocks = make_contiguous_blocks(selected, P.nc_track_block_size);
[~, stem, ~] = fileparts(P.track_nc);

fprintf('\nCounting events from %d selected tracks in %d chunks...\n', numel(selected), size(blocks, 1));

events = empty_event_struct();
nEv = 0;
domainCount = zeros(numel(Domains), 1);

for ib = 1:size(blocks, 1)
    firstTrack = blocks(ib, 1);
    nBlock = blocks(ib, 2);
    if trackDim == 1
        start = [firstTrack 1];
        count = [nBlock nTime];
        lonBlock = clean_fill_to_nan(double(ncread(P.track_nc, 'lon_trks', start, count)));
        latBlock = clean_fill_to_nan(double(ncread(P.track_nc, 'lat_trks', start, count)));
        vmaxBlock = clean_fill_to_nan(double(ncread(P.track_nc, 'vmax_trks', start, count)));
    else
        start = [1 firstTrack];
        count = [nTime nBlock];
        lonBlock = clean_fill_to_nan(double(ncread(P.track_nc, 'lon_trks', start, count))).';
        latBlock = clean_fill_to_nan(double(ncread(P.track_nc, 'lat_trks', start, count))).';
        vmaxBlock = clean_fill_to_nan(double(ncread(P.track_nc, 'vmax_trks', start, count))).';
    end

    for it = 1:nBlock
        trackIndex = firstTrack + it - 1;
        yy = yearsAll(trackIndex);
        vmax = double(vmaxBlock(it, :));
        goodV = isfinite(vmax);
        if ~any(goodV)
            continue;
        end
        lmi = max(vmax(goodV));
        if ~isfinite(lmi) || lmi < P.min_record_lmi_ms || lmi <= P.min_inner_vmax_ms
            continue;
        end

        [lonD, latD, windD, stepD] = densify_track_if_requested( ...
            lonBlock(it, :), latBlock(it, :), vmax, P);
        valid = isfinite(lonD) & isfinite(latD) & isfinite(windD) & ...
            latD >= -90 & latD <= 90 & windD > P.min_inner_vmax_ms;
        if ~any(valid)
            continue;
        end

        lonH = lonD(valid);
        latH = latD(valid);
        windH = windD(valid);
        stepH = stepD(valid);
        cand = candidate_domains_by_bbox(Domains, lonH, latH);
        if isempty(cand)
            continue;
        end

        trackId = sprintf('%s_%06d', stem, trackIndex);
        for ic = 1:numel(cand)
            id = cand(ic);
            ring = Domains(id).impact_ring;
            [in, on] = inpolygon(lonH, latH, ring(:, 1), ring(:, 2));
            inside = in | on;
            if ~any(inside)
                continue;
            end
            localIdx = find(inside);
            [mx, loc] = max(windH(inside));
            if ~isfinite(mx) || mx <= P.min_inner_vmax_ms
                continue;
            end
            k = localIdx(loc);

            nEv = nEv + 1;
            domainCount(id) = domainCount(id) + 1;
            events(nEv).block_id = Domains(id).block_id; %#ok<AGROW>
            events(nEv).track_id = string(trackId);
            events(nEv).track_index = double(trackIndex);
            events(nEv).year = double(yy);
            events(nEv).max_inner_vmax_ms = double(mx);
            events(nEv).lmi_ms = double(lmi);
            events(nEv).max_lon = double(lonH(k));
            events(nEv).max_lat = double(latH(k));
            events(nEv).max_time_index = double(stepH(k));
        end
    end

    fprintf('  chunk %02d/%02d processed; events so far: %d\n', ib, size(blocks, 1), nEv);
end

Events = event_struct_to_table(events);
Summary = make_summary_table(T_blocks, Domains, Events, domainCount, P);

fprintf('Total per-domain TC simulations: %d\n', sum(Summary.tc_count_vmax_gt33_ms));
fprintf('Unique TC tracks across all domains: %d\n', numel(unique(Events.track_id)));
end

%% ============================================================
function events = empty_event_struct()

events = struct('block_id', {}, 'track_id', {}, 'track_index', {}, 'year', {}, ...
    'max_inner_vmax_ms', {}, 'lmi_ms', {}, 'max_lon', {}, 'max_lat', {}, ...
    'max_time_index', {});
end

%% ============================================================
function T = event_struct_to_table(events)

if isempty(events)
    T = table(strings(0, 1), strings(0, 1), zeros(0, 1), zeros(0, 1), ...
        zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
        'VariableNames', {'block_id', 'track_id', 'track_index', 'year', ...
        'max_inner_vmax_ms', 'lmi_ms', 'max_lon', 'max_lat', 'max_time_index'});
    return;
end
T = struct2table(events);
T.block_id = string(T.block_id);
T.track_id = string(T.track_id);
end

%% ============================================================
function Summary = make_summary_table(T_blocks, Domains, Events, domainCount, P)

n = height(T_blocks);
rows = repmat(struct( ...
    'block_id', "", 'basin_id', "", 'basin_label', "", 'member_domains', "", ...
    'inner_area_km2', NaN, 'tc_count_vmax_gt33_ms', 0, 'max_inner_vmax_ms', NaN, ...
    'year_min', NaN, 'year_max', NaN), n, 1);

for i = 1:n
    bid = string(T_blocks.block_id(i));
    rows(i).block_id = bid;
    rows(i).basin_id = string(get_table_value_or_blank(T_blocks, i, "basin_id"));
    rows(i).basin_label = string(get_table_value_or_blank(T_blocks, i, "basin_label"));
    rows(i).member_domains = string(get_table_value_or_blank(T_blocks, i, "member_domains"));
    rows(i).inner_area_km2 = double(get_table_value_or_nan(T_blocks, i, "inner_area_km2"));
    rows(i).impact_buffer_km = double(P.impact_buffer_km);
    rows(i).tc_count_vmax_gt33_ms = double(domainCount(i));

    if ~isempty(Events)
        I = Events.block_id == bid;
        if any(I)
            rows(i).max_inner_vmax_ms = max(Events.max_inner_vmax_ms(I), [], 'omitnan');
            rows(i).year_min = min(Events.year(I), [], 'omitnan');
            rows(i).year_max = max(Events.year(I), [], 'omitnan');
        end
    end

    if strlength(rows(i).basin_id) == 0 && isfield(Domains(i), 'basin_id')
        rows(i).basin_id = Domains(i).basin_id;
    end
end

Summary = struct2table(rows);
end

%% ============================================================
function v = get_table_value_or_blank(T, i, name)

name = string(name);
if any(string(T.Properties.VariableNames) == name)
    v = T{i, char(name)};
else
    v = "";
end
end

%% ============================================================
function v = get_table_value_or_nan(T, i, name)

name = string(name);
if any(string(T.Properties.VariableNames) == name)
    v = T{i, char(name)};
else
    v = NaN;
end
end

%% ============================================================
function cand = candidate_domains_by_bbox(Domains, lon, lat)

if isempty(lon) || isempty(lat)
    cand = [];
    return;
end
minLon = min(lon);
maxLon = max(lon);
minLat = min(lat);
maxLat = max(lat);
cand = [];
for i = 1:numel(Domains)
    b = Domains(i).impact_bbox;
    if b(2) < minLon || b(1) > maxLon || b(4) < minLat || b(3) > maxLat
        continue;
    end
    cand(end + 1) = i; %#ok<AGROW>
end
end

%% ============================================================
function [lonD, latD, windD, stepD] = densify_track_if_requested(lon, lat, wind, P)

n = min([numel(lon), numel(lat), numel(wind)]);
lon = wrap_to_180_local(double(lon(1:n)));
lat = double(lat(1:n));
wind = double(wind(1:n));

valid = isfinite(lon) & isfinite(lat) & isfinite(wind);
orig = find(valid);
lon = lon(valid);
lat = lat(valid);
wind = wind(valid);

lonD = [];
latD = [];
windD = [];
stepD = [];

if isempty(lon)
    return;
end
if numel(lon) < 2
    lonD = lon(:);
    latD = lat(:);
    windD = wind(:);
    stepD = orig(:);
    return;
end

substeps = max(1, round(P.segment_substeps));
for i = 1:(numel(lon) - 1)
    lon0 = lon(i);
    lon1 = lon0 + wrap_to_180_local(lon(i + 1) - lon0);
    lat0 = lat(i);
    lat1 = lat(i + 1);
    w0 = wind(i);
    w1 = wind(i + 1);

    t = (0:substeps-1).' / substeps;
    lonSeg = wrap_to_180_local(lon0 + t .* (lon1 - lon0));
    latSeg = lat0 + t .* (lat1 - lat0);
    wSeg = w0 + t .* (w1 - w0);
    sSeg = round(orig(i) + t .* (orig(i + 1) - orig(i)));

    lonD = [lonD; lonSeg]; %#ok<AGROW>
    latD = [latD; latSeg]; %#ok<AGROW>
    windD = [windD; wSeg]; %#ok<AGROW>
    stepD = [stepD; sSeg]; %#ok<AGROW>
end

lonD = [lonD; lon(end)];
latD = [latD; lat(end)];
windD = [windD; wind(end)];
stepD = [stepD; orig(end)];
end

%% ============================================================
function FigureIndex = plot_tc_track_maps_for_domains(P, T_blocks, Domains, Summary, Events)

ensure_dir(P.figure_dir);
[trackDim, nTime] = infer_track_dim_from_nc(P.track_nc, 'lon_trks', numel(ncread(P.track_nc, 'tc_years')));

n = numel(Domains);
rows = repmat(struct('block_id', "", 'tc_count', 0, 'png', "", 'pdf', ""), n, 1);
fprintf('\nPlotting TC track maps for %d ADCIRC blocks...\n', n);

for i = 1:n
    blockId = Domains(i).block_id;
    I = Events.block_id == blockId;
    E = Events(I, :);
    trackIdx = unique(E.track_index);
    [trackLon, trackLat] = read_track_lines_for_indices(P.track_nc, trackDim, nTime, trackIdx);

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 P.figure_size_cm], 'Visible', P.figure_visible);
    ax = axes('Parent', fig);
    hold(ax, 'on');
    set(ax, 'FontName', P.font_name, 'FontSize', P.axis_font_size, ...
        'Layer', 'top', 'Box', 'on');

    [lonLim, latLim] = map_limits_from_domain(Domains(i), P);
    useMMap = setup_map_axes(lonLim, latLim);

    draw_land_background(useMMap, P);
    hTracks = draw_track_lines(trackLon, trackLat, P, useMMap);
    hInner = draw_polygon(Domains(i).inner_ring, P.inner_face_color, P.inner_edge_color, 1.70, 0.38, useMMap);

    title(ax, sprintf('%s', blockId), ...
        'FontName', P.font_name, 'FontSize', P.title_font_size, ...
        'FontWeight', 'normal', 'Interpreter', 'none');

    tcCount = height(E);
    text(ax, 0.985, 0.975, sprintf('TC count: %d', tcCount), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'FontName', P.font_name, ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.05 0.05 0.05], ...
        'BackgroundColor', [1 1 1], 'Margin', 4, 'Interpreter', 'none');

    legendHandles = [hInner, hTracks];
    legendLabels = {'Inner refined domain', 'ERA5 TC tracks'};
    ok = isgraphics(legendHandles);
    if any(ok)
        legend(ax, legendHandles(ok), legendLabels(ok), 'Location', 'southwest', ...
            'FontName', P.font_name, 'FontSize', 7, 'Box', 'on');
    end

    base = sprintf('%s_era5_tc_tracks_inner_vmax_gt33', sanitize_filename(blockId));
    outPng = fullfile(P.figure_dir, [base '.png']);
    exportgraphics(fig, outPng, 'Resolution', P.figure_dpi, 'BackgroundColor', 'white');
    outPdf = "";
    if P.save_pdf
        outPdf = fullfile(P.figure_dir, [base '.pdf']);
        exportgraphics(fig, outPdf, 'ContentType', 'vector', 'BackgroundColor', 'white');
    end
    close(fig);

    rows(i).block_id = blockId;
    rows(i).tc_count = tcCount;
    rows(i).png = string(outPng);
    rows(i).pdf = string(outPdf);
    fprintf('  [%02d/%02d] %s: %d tracks\n', i, n, blockId, tcCount);
end

FigureIndex = struct2table(rows);
end

%% ============================================================
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

%% ============================================================
function [trackLonLine, trackLatLine] = read_track_lines_for_indices(ncFile, trackDim, nTime, trackIdx)

trackLonLine = [];
trackLatLine = [];
trackIdx = unique(double(trackIdx(:)));
trackIdx = trackIdx(isfinite(trackIdx) & trackIdx >= 1);
if isempty(trackIdx)
    return;
end

blocks = make_contiguous_blocks(trackIdx, inf);
for ib = 1:size(blocks, 1)
    firstTrack = blocks(ib, 1);
    nBlock = blocks(ib, 2);
    if trackDim == 1
        start = [firstTrack 1];
        count = [nBlock nTime];
        lonBlock = clean_fill_to_nan(double(ncread(ncFile, 'lon_trks', start, count)));
        latBlock = clean_fill_to_nan(double(ncread(ncFile, 'lat_trks', start, count)));
    else
        start = [1 firstTrack];
        count = [nTime nBlock];
        lonBlock = clean_fill_to_nan(double(ncread(ncFile, 'lon_trks', start, count))).';
        latBlock = clean_fill_to_nan(double(ncread(ncFile, 'lat_trks', start, count))).';
    end
    for it = 1:nBlock
        lon = wrap_to_180_local(lonBlock(it, :));
        lat = latBlock(it, :);
        ok = isfinite(lon) & isfinite(lat) & lat >= -90 & lat <= 90;
        lon = lon(ok);
        lat = lat(ok);
        if numel(lon) < 2
            continue;
        end
        [lon, lat] = break_dateline_segments(lon(:), lat(:));
        trackLonLine = [trackLonLine; lon(:); NaN]; %#ok<AGROW>
        trackLatLine = [trackLatLine; lat(:); NaN]; %#ok<AGROW>
    end
end
end

%% ============================================================
function [lon, lat] = break_dateline_segments(lon, lat)

if numel(lon) < 2
    return;
end
d = abs(diff(lon));
breaks = find(d > 180);
for k = numel(breaks):-1:1
    j = breaks(k);
    lon = [lon(1:j); NaN; lon(j+1:end)]; %#ok<AGROW>
    lat = [lat(1:j); NaN; lat(j+1:end)]; %#ok<AGROW>
end
end

%% ============================================================
function [lonLim, latLim] = map_limits_from_domain(Domain, P)

ring = Domain.outer_ring;
if isempty(ring)
    ring = Domain.inner_ring;
end
b = ring_bbox(ring);
dx = max(b(2) - b(1), 1);
dy = max(b(4) - b(3), 1);
padX = max(P.map_min_pad_deg, P.map_pad_fraction * dx);
padY = max(P.map_min_pad_deg, P.map_pad_fraction * dy);
lonLim = [b(1) - padX, b(2) + padX];
latLim = [b(3) - padY, b(4) + padY];
lonLim = max([-180 -180], min([180 180], lonLim));
latLim = max([-89 -89], min([89 89], latLim));
if diff(lonLim) < 2
    lonLim = mean(lonLim) + [-1 1];
end
if diff(latLim) < 2
    latLim = mean(latLim) + [-1 1];
end
end

%% ============================================================
function useMMap = setup_map_axes(lonLim, latLim)

useMMap = exist('m_proj', 'file') == 2 && exist('m_grid', 'file') == 2;
if useMMap
    try
        m_proj('miller', 'lon', lonLim, 'lat', latLim);
        m_grid('box', 'fancy', 'tickdir', 'out', 'linestyle', ':', ...
            'fontsize', 8, 'fontname', 'Arial');
        return;
    catch ME
        warning('m_map setup failed; falling back to lon/lat axes: %s', ME.message);
        useMMap = false;
    end
end
axis equal;
xlim(lonLim);
ylim(latLim);
grid on;
xlabel('Longitude');
ylabel('Latitude');
end

%% ============================================================
function draw_land_background(useMMap, P)

if useMMap && exist('m_gshhs_i', 'file') == 2
    try
        m_gshhs_i('patch', P.land_color, 'edgecolor', P.coast_edge_color, 'linewidth', 0.35);
    catch
    end
end
end

%% ============================================================
function h = draw_polygon(ring, faceColor, edgeColor, lineWidth, faceAlpha, useMMap)

h = gobjects(1, 1);
if isempty(ring) || size(ring, 1) < 3
    return;
end
lon = ring(:, 1);
lat = ring(:, 2);
hold on;
try
    if useMMap
        h = m_patch(lon, lat, faceColor, 'EdgeColor', edgeColor, 'LineWidth', lineWidth);
    else
        h = patch(lon, lat, faceColor, 'EdgeColor', edgeColor, 'LineWidth', lineWidth);
    end
    set(h, 'FaceAlpha', faceAlpha);
catch
    if useMMap
        h = m_plot(lon, lat, '-', 'Color', edgeColor, 'LineWidth', lineWidth);
    else
        h = plot(lon, lat, '-', 'Color', edgeColor, 'LineWidth', lineWidth);
    end
end
end

%% ============================================================
function h = draw_track_lines(lon, lat, P, useMMap)

h = gobjects(1, 1);
if isempty(lon) || numel(lon) < 2
    return;
end
hold on;
try
    if useMMap
        [x, y] = m_ll2xy(lon, lat, 'clip', 'off');
    else
        x = lon;
        y = lat;
    end
    h = patch('XData', x(:), 'YData', y(:), ...
        'FaceColor', 'none', ...
        'EdgeColor', P.track_color, ...
        'EdgeAlpha', P.track_alpha, ...
        'LineWidth', P.track_line_width, ...
        'Tag', 'era5_tc_track_overlay');
catch
    if useMMap
        h = m_plot(lon, lat, '-', 'Color', P.track_color, 'LineWidth', P.track_line_width);
    else
        h = plot(lon, lat, '-', 'Color', P.track_color, 'LineWidth', P.track_line_width);
    end
end
end

%% ============================================================
function blocks = make_contiguous_blocks(indices, maxBlock)

indices = double(indices(:));
if isempty(indices)
    blocks = zeros(0, 2);
    return;
end
indices = sort(indices);
gapAfter = find(diff(indices) > 1);
starts = [1; gapAfter + 1];
ends = [gapAfter; numel(indices)];
blocks = zeros(0, 2);
if ~isfinite(maxBlock)
    maxBlock = inf;
end
for ir = 1:numel(starts)
    a = indices(starts(ir));
    b = indices(ends(ir));
    s = a;
    while s <= b
        c = min(maxBlock, b - s + 1);
        blocks(end + 1, :) = [s c]; %#ok<AGROW>
        s = s + c;
    end
end
end

%% ============================================================
function x = clean_fill_to_nan(x)

x = double(x);
x(abs(x) > 1.0e20) = NaN;
x(x <= -9998) = NaN;
end

%% ============================================================
function lon = wrap_to_180_local(lon)

lon = mod(double(lon) + 180, 360) - 180;
lon(lon == -180) = 180;
end

%% ============================================================
function name = sanitize_filename(x)

name = regexprep(char(string(x)), '[^A-Za-z0-9_\-]+', '_');
name = regexprep(name, '_+', '_');
end

%% ============================================================
function ensure_dir(d)

if exist(d, 'dir') ~= 7
    mkdir(d);
end
end

%% ============================================================
function add_project_paths(scriptDir)

repoDir = find_parent_dir_named(scriptDir, 'OceanMesh2D-Projection');
if ~isempty(repoDir)
    addpath(genpath(fullfile(repoDir, 'm_map')));
    addpath(genpath(fullfile(repoDir, 'utilities')));
end
end

%% ============================================================
function root = find_parent_dir_named(startDir, targetName)

root = '';
d = char(startDir);
targetName = char(targetName);
while true
    [parent, name] = fileparts(d);
    if strcmpi(name, targetName)
        root = d;
        return;
    end
    if isempty(parent) || strcmp(parent, d)
        return;
    end
    d = parent;
end
end

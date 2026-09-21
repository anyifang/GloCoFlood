%% P1_build_adcirc_model_blocks_from_partitions.m
% Build ADCIRC/OceanMesh2D modelling blocks from the TC-exposed coastal
% partitions produced by partition_global_tc_coastal_model_domains.m.
%
% This script prepares block-scale ADCIRC domains without requiring Mapping
% Toolbox.  It writes GeoJSON/CSV block boundaries and a separate optional
% OceanMesh2D runner.  The runner follows the Global_auto.m workflow when
% OceanMesh2D classes are available.

clearvars;
clc;
close all;
warning('off', 'MATLAB:polyshape:repairedBySimplify');

%% ------------------------- user config -------------------------
SCRIPT_DIR = fileparts(mfilename('fullpath'));
if isempty(SCRIPT_DIR)
    SCRIPT_DIR = pwd;
end

P = struct();
P.partition_dir = fullfile(SCRIPT_DIR, 'external', 'global_tc_coastal_model_partitions');
P.domain_csv = fullfile(P.partition_dir, 'global_tc_coastal_model_domains.csv');
P.cell_csv = fullfile(P.partition_dir, 'global_tc_coastal_model_domain_cells.csv');
P.domain_geojson = fullfile(P.partition_dir, 'global_tc_coastal_model_domain_boundaries.geojson');

P.output_dir = fullfile(SCRIPT_DIR, 'output', 'global_tc_adcirc_model_blocks');
P.block_boundary_dir = fullfile(P.output_dir, 'block_boundaries');
P.quicklook_dir = fullfile(P.output_dir, 'quicklooks');
P.selected_block_ids = strings(0, 1);

v = strtrim(string(getenv('ADCIRC_PARTITION_INPUT_ROOT')));
if strlength(v) > 0
    P.partition_dir = char(v);
    P.domain_csv = fullfile(P.partition_dir, 'global_tc_coastal_model_domains.csv');
    P.cell_csv = fullfile(P.partition_dir, 'global_tc_coastal_model_domain_cells.csv');
    P.domain_geojson = fullfile(P.partition_dir, 'global_tc_coastal_model_domain_boundaries.geojson');
end
v = strtrim(string(getenv('ADCIRC_BLOCK_OUTPUT_ROOT')));
if strlength(v) > 0
    P.output_dir = char(v);
    P.block_boundary_dir = fullfile(P.output_dir, 'block_boundaries');
    P.quicklook_dir = fullfile(P.output_dir, 'quicklooks');
end

% Keep ADCIRC domains regional, not estuary-by-estuary.  Nearby coastal
% partitions are grouped into one block until the alongshore/cell span would
% make the open-ocean domain too broad.
P.max_block_alongshore_km = 1500;
P.max_block_span_km = 2600;
P.max_domains_per_block = 10;
% OceanMesh2D concentrates most elements inside the inner refinement domain.
% Keep merged blocks from becoming too expensive by limiting this area.
P.max_inner_refinement_area_km2 = 5.0e5;
P.centroid_gap_split_km = 500;
P.tiny_block_min_domains = 1;
P.tiny_block_min_cells = 0;

% Second-pass block merging.  The first pass keeps the coastal segmentation
% conservative; this pass merges neighbouring blocks whose final outer
% ADCIRC domains overlap substantially, while still preventing basin-scale
% open-ocean domains.
P.merge_overlapping_blocks = true;
P.merge_outer_overlap_min_fraction = 0.22;  % intersection / smaller outer domain
P.merge_outer_union_min_fraction = 0.065;   % intersection / union of two outer domains
P.merge_outer_gap_km = 100;                 % also merge nearly touching outer domains
P.merge_candidate_max_outer_area_km2 = 8.0e6;
P.merge_candidate_max_span_km = 2000;
P.merge_candidate_max_domains = 10;
P.coastal_connect_gap_km = 180;
P.merge_coastal_gap_km = 180;

% Final post-check for very small blocks.  Tiny island/coastal remnants are
% absorbed into a nearby major block when their refined coastal domains are
% close, or merged with neighbouring tiny remnants that form the same island
% chain.  This runs after the conservative connectivity split above.
P.absorb_small_blocks = true;
P.small_block_max_cells = 80;
P.small_block_cluster_max_cells = 140;
P.small_block_cluster_max_domains = 6;
P.small_block_cluster_max_inner_area_km2 = 450000;
P.small_block_max_domains = 2;
P.small_block_max_inner_area_km2 = 220000;
P.small_major_min_cells = 150;
P.small_merge_to_major_gap_km = 240;
P.small_merge_small_gap_km = 650;
P.small_merge_outer_overlap_min_fraction = 0.35;
P.small_merge_outer_union_min_fraction = 0.10;
P.small_merge_candidate_max_outer_area_km2 = 9.0e6;
P.small_merge_candidate_max_span_km = 2800;
P.small_merge_candidate_max_domains = 55;

% Geometry generation.  The inner refinement domain is the original
% partition boundary expanded by this distance.  The outer ADCIRC domain is
% generated from the coastal research segment using morphological bay
% closing/opening ideas from Examples/Global_autofunction.
P.inner_refinement_buffer_km = 50;
P.outer_base_buffer_km = 300;
P.outer_shelf_margin_km = 180;
P.outer_min_buffer_km = 300;
P.outer_max_buffer_km = 500;
P.outer_bay_close_km = 160;
P.outer_smooth_km = 35;
P.use_gebco_shelf_width = true;
P.shelf_isobath_m = -200;
P.shelf_min_width_km = 250;
P.shelf_max_width_km = 760;
P.shelf_sample_step = 8;
P.outer_open_ocean_max_area_km2 = 8.0e6;
P.fast_blocking = true;

% MATLAB-only sea-side geometry guard.  The ocean-side domain remains the
% full shelf/bay buffer; only the landward/backshore side is clipped so a
% domain does not spill through a narrow land barrier into the opposite sea.
P.outer_backshore_limit_km = 20;
P.keep_largest_outer_region = true;
P.force_single_outer_region = true;
P.final_keep_only_inner_touching_regions = true;
P.final_inner_touch_buffer_km = 8;
P.final_inner_touch_min_area_km2 = 1;
P.final_enforce_single_ocean_component = true;
P.final_ocean_seed_buffer_km = 160;
P.final_ocean_land_keep_margin_km = 80;
P.final_ocean_min_component_area_km2 = 5000;
P.final_ocean_keep_mask_buffer_km = 10;
P.final_ocean_connectivity_report = true;
P.final_ocean_prefer_inner_touching = true;
P.final_strict_single_ocean_component = true;
P.final_strict_ocean_min_component_area_km2 = 10;
P.final_strict_ocean_land_margin_km = 80;
P.final_strict_ocean_boundary_smooth_km = 55;
P.final_strict_ocean_report = true;
P.final_smooth_open_ocean_boundary = true;
P.final_open_boundary_fill_only = true;
P.final_open_boundary_protect_km = 150;
P.final_open_boundary_close_km = 460;
P.final_open_boundary_open_km = 0;
P.final_open_boundary_round_km = 180;
P.final_open_boundary_round_inset_fraction = 0.84;
P.final_open_boundary_edge_round_km = 220;
P.final_open_boundary_edge_round_inset_fraction = 0.86;
P.final_open_boundary_edge_round_cover_min = 0.965;
P.final_open_boundary_curve_smooth_enable = true;
P.final_open_boundary_curve_step_km = 35;
P.final_open_boundary_curve_window_km = 520;
P.final_open_boundary_curve_cover_pad_km = 90;
P.final_open_boundary_curve_cover_min = 0.965;
P.final_open_boundary_curve_max_area_increase_fraction = 0.26;
P.final_open_boundary_envelope_enable = true;
P.final_open_boundary_envelope_shrink = 0.52;
P.final_open_boundary_envelope_pad_km = 80;
P.final_open_boundary_envelope_round_km = 180;
P.final_open_boundary_envelope_round_inset_fraction = 0.90;
P.final_open_boundary_envelope_cover_min = 0.965;
P.final_open_boundary_envelope_max_area_increase_fraction = 0.26;
P.final_open_boundary_max_extra_km = 380;
P.final_open_boundary_min_area_fraction = 0.995;
P.final_open_boundary_max_area_increase_fraction = 0.38;
P.final_open_boundary_report = true;
P.final_output_outer_smooth_enable = true;
P.final_output_outer_protect_km = 130;
P.final_output_outer_round_km = 160;
P.final_output_outer_round_inset_fraction = 0.90;
P.final_output_outer_curve_step_km = 35;
P.final_output_outer_curve_window_km = 720;
P.final_output_outer_curve_cover_pad_km = 40;
P.final_output_outer_curve_cover_min = 0.88;
P.final_output_outer_min_area_fraction = 0.80;
P.final_output_outer_max_area_increase_fraction = 0.35;
P.limit_outer_seed_alongshore = true;
P.outer_seed_alongshore_pad_km = 120;
P.use_gebco_topology = true;
P.topology_grid_km = 18;
P.topology_max_cells = 2.2e5;
P.topology_seed_buffer_km = 140;
P.topology_bay_extra_km = 1000;
P.topology_major_land_min_area_km2 = 25000;
P.topology_capsule_style = true;
P.topology_capsule_inland_km = 1.5 * 111.32;
P.topology_capsule_bay_close_km = 3.0 * 111.32;
P.topology_capsule_razor_km = 90;
P.topology_capsule_fill_km = 70;
P.topology_capsule_min_component_area_km2 = 9000;
P.topology_capsule_other_ocean_margin_cells = 1;
P.topology_strait_barrier_km = 150;
P.topology_strait_barrier_smooth_km = 45;
P.topology_strait_barrier_max_land_components = 25;
% GEBCO topology adjustment is split into two physical steps:
%   1) shelf expansion: follow important shallow continental shelves;
%   2) bay completion: only add water bodies that become closed by nearby
%      major land, e.g. the Gulf of Mexico, not ordinary open coast arcs.
P.topology_shelf_expand = false;     % raw shelf crawling is too branchy; keep shelf only for gap-fill below
P.topology_shelf_isobath_m = P.shelf_isobath_m;
P.topology_shelf_extra_km = 720;
P.topology_shelf_touch_km = 180;
P.topology_shelf_reach_km = 620;
P.topology_shelf_max_area_km2 = 3.0e6;
P.topology_shelf_total_max_area_km2 = 4.2e6;
P.topology_shelf_search_edge_touch_max_fraction = 0.08;
P.topology_shelf_razor_km = 95;
P.topology_shelf_fill_km = 80;
P.topology_shelf_smooth_km = 70;
P.topology_shelf_min_minor_width_km = 140;
P.topology_shelf_min_fill_fraction = 0.12;
P.topology_shelf_gap_fill_km = 460;
P.topology_shelf_gap_search_km = 700;
P.topology_shelf_gap_shallow_only = false;
P.topology_shelf_gap_max_area_km2 = 2.4e6;
P.topology_shelf_gap_total_max_area_km2 = 3.2e6;
P.topology_shelf_gap_edge_touch_max_fraction = 0.040;
P.topology_land_bay_bridge_km = 320;
P.topology_land_bay_search_km = 700;
P.topology_land_bay_touch_km = 220;
P.topology_land_bay_min_area_km2 = 30000;
P.topology_land_bay_max_area_km2 = 2.2e6;
P.topology_land_bay_total_max_area_km2 = 3.2e6;
P.topology_land_bay_edge_touch_max_fraction = 0.050;
P.topology_land_bay_min_minor_width_km = 90;
P.topology_land_bay_min_fill_fraction = 0.080;
P.topology_land_bay_min_touch_fraction = 0.075;
P.topology_land_bay_enclosure_km = 650;
P.topology_land_bay_min_land_sectors = 3;
P.topology_land_bay_min_opposite_pairs = 1;
P.topology_land_bay_max_major_span_km = 2600;
P.topology_major_bay_enable = true;
P.topology_major_bay_close_km = 520;
P.topology_major_bay_recover_km = 820;
P.topology_major_bay_search_km = 1200;
P.topology_major_bay_touch_km = 420;
P.topology_major_bay_min_area_km2 = 250000;
P.topology_major_bay_max_area_km2 = 4.2e6;
P.topology_major_bay_total_max_area_km2 = 5.2e6;
P.topology_major_bay_edge_touch_max_fraction = 0.035;
P.topology_major_bay_min_land_sector_fraction = 0.35;
P.topology_major_bay_min_touch_fraction = 0.025;
P.topology_major_bay_min_touch_sectors = 2;
P.topology_named_embayment_enable = true;
P.topology_named_embayment_extra_km = 2200;
P.topology_named_embayment_window_pad_km = 120;
P.topology_named_embayment_window_use_superellipse = true;
P.topology_named_embayment_window_scale = 1.10;
P.topology_named_embayment_window_exponent = 2.0;
P.topology_named_embayment_touch_km = 650;
P.topology_named_embayment_smooth_km = 260;
P.topology_named_embayment_trigger_km = 180;
P.topology_named_caribbean_embayments_enable = true;
P.topology_named_wpac_embayments_enable = true;
P.topology_bay_ring_enable = false;
P.topology_bay_ring_close_km = 380;
P.topology_bay_ring_recover_km = 600;
P.topology_bay_ring_touch_km = 280;
P.topology_bay_ring_min_area_km2 = 50000;
P.topology_bay_ring_max_area_km2 = 3.8e6;
P.topology_bay_ring_total_max_area_km2 = 4.8e6;
P.topology_bay_ring_search_edge_touch_max_fraction = 0.045;
P.topology_bay_ring_current_touch_min_cells = 6;
P.topology_bay_domain_smooth_km = 220;
P.topology_bay_domain_smooth_ocean_only = true;
P.topology_bay_domain_land_protect_km = 90;
P.topology_bay_domain_smooth_inset_fraction = 0.70;
P.topology_bay_domain_smooth_cover_min = 0.88;
P.topology_bay_domain_smooth_max_area_increase_fraction = 0.50;
P.topology_bay_domain_curve_smooth_enable = true;
P.topology_bay_domain_curve_step_km = 30;
P.topology_bay_domain_curve_window_km = 560;
P.topology_bay_domain_curve_cover_min = 0.86;
P.topology_bay_domain_curve_max_area_increase_fraction = 0.55;
P.topology_bay_domain_final_smooth_enable = true;
P.topology_bay_domain_final_smooth_km = 360;
P.topology_bay_domain_final_land_protect_km = 35;
P.topology_bay_domain_final_inset_fraction = 0.88;
P.topology_bay_domain_final_cover_min = 0.72;
P.topology_bay_domain_final_min_area_fraction = 0.60;
P.topology_bay_domain_final_max_area_increase_fraction = 0.85;
P.topology_bay_domain_final_curve_enable = true;
P.topology_bay_domain_final_curve_step_km = 35;
P.topology_bay_domain_final_curve_window_km = 900;
P.topology_bay_domain_final_curve_cover_min = 0.70;
P.topology_bay_domain_final_curve_max_area_increase_fraction = 0.90;
P.topology_ocean_zmax_m = -1.0;
P.topology_boundary_smooth_km = 260;
P.topology_boundary_inset_fraction = 0.95;
P.topology_land_other_ocean_margin_cells = 1;

% Split final block groupings when the combined inner-refinement buffer spans
% disconnected sea basins. This keeps each ADCIRC block focused on one sea
% without clipping the final black inner-refinement polygon to the coastline.
P.split_blocks_by_inner_ocean_components = true;
P.inner_ocean_split_grid_km = 10;
P.inner_ocean_split_component_pad_km = 12;
P.inner_ocean_split_smooth_km = 8;
P.inner_ocean_split_min_area_km2 = 500;
P.inner_ocean_split_max_grid_cells = 2.5e5;
P.inner_ocean_split_max_components = 12;

% Global_auto/OceanMesh2D-style parameters exported for each block.
P.auto_inland_deg = 1.5;
P.auto_bay_close_deg = 3.0;
P.auto_inner_radius_min_deg = 1.5;
P.auto_inner_radius_max_deg = 7.0;
P.auto_alongshore_min_deg = 2.5;
P.auto_alongshore_max_deg = 12.0;
P.auto_offshore_min_deg = 3.0;
P.auto_offshore_max_deg = 7.0;

P.figure_dpi = 240;

P.repo_dir = find_parent_dir_named(SCRIPT_DIR, 'OceanMesh2D-Projection');
if isempty(P.repo_dir)
    P.repo_dir = fileparts(SCRIPT_DIR);
end
P.gebco_path = fullfile(P.repo_dir, 'datasets', 'GEBCO', 'GEBCO_2025_sub_ice.nc');
add_existing_genpath({fullfile(P.repo_dir, 'm_map')});
add_existing_genpath({fullfile(P.repo_dir, 'Examples', 'Global_autofunction')});

% Runtime overrides:
%   setenv('ADCIRC_FAST_BLOCKING','1')       % default, fast preliminary blocks
%   setenv('ADCIRC_FAST_BLOCKING','0')       % enable GEBCO shelf-width estimate
%   setenv('ADCIRC_USE_GEBCO_SHELF','1')    % force GEBCO shelf-width estimate
%   setenv('ADCIRC_MAX_INNER_REFINEMENT_AREA_KM2','1000000')
%   setenv('ADCIRC_OUTPUT_DIR','...\debug_blocks')
%   setenv('ADCIRC_SELECTED_BLOCKS','ADC_WNP_06,ADC_NATL_05,ADC_NATL_01')
fastEnv = string(strtrim(getenv('ADCIRC_FAST_BLOCKING')));
if any(strcmpi(fastEnv, ["0", "false", "no", "n"]))
    P.fast_blocking = false;
elseif any(strcmpi(fastEnv, ["1", "true", "yes", "y"]))
    P.fast_blocking = true;
end
if P.fast_blocking
    P.use_gebco_shelf_width = false;
    P.figure_dpi = min(P.figure_dpi, 180);
end
gebcoEnv = string(strtrim(getenv('ADCIRC_USE_GEBCO_SHELF')));
if any(strcmpi(gebcoEnv, ["1", "true", "yes", "y"]))
    P.use_gebco_shelf_width = true;
elseif any(strcmpi(gebcoEnv, ["0", "false", "no", "n"]))
    P.use_gebco_shelf_width = false;
end
innerAreaCapEnv = str2double(getenv('ADCIRC_MAX_INNER_REFINEMENT_AREA_KM2'));
if isfinite(innerAreaCapEnv) && innerAreaCapEnv > 0
    P.max_inner_refinement_area_km2 = innerAreaCapEnv;
end
outputDirEnv = string(strtrim(getenv('ADCIRC_OUTPUT_DIR')));
if strlength(outputDirEnv) > 0
    P.output_dir = char(outputDirEnv);
    P.block_boundary_dir = fullfile(P.output_dir, 'block_boundaries');
    P.quicklook_dir = fullfile(P.output_dir, 'quicklooks');
end
selectedBlocksEnv = string(strtrim(getenv('ADCIRC_SELECTED_BLOCKS')));
if strlength(selectedBlocksEnv) > 0 && ~any(strcmpi(selectedBlocksEnv, ["all", "none", "0", "false"]))
    parts = regexp(char(selectedBlocksEnv), '[,;\s]+', 'split');
    parts = string(parts(:));
    parts = upper(strtrim(parts(parts ~= "")));
    P.selected_block_ids = unique(parts, 'stable');
end

%% ------------------------- init -------------------------
ensure_dir(P.output_dir);
ensure_dir(P.block_boundary_dir);
ensure_dir(P.quicklook_dir);
cleanup_previous_block_outputs(P);

fprintf('\n============================================================\n');
fprintf('Build ADCIRC model blocks from TC coastal partitions\n');
fprintf('Domain CSV : %s\n', P.domain_csv);
fprintf('Cell CSV   : %s\n', P.cell_csv);
fprintf('Output dir : %s\n', P.output_dir);
fprintf('Max block  : %.0f km alongshore, %.0f km total span, %d domains\n', ...
    P.max_block_alongshore_km, P.max_block_span_km, P.max_domains_per_block);
fprintf('Max inner refinement area after merge: %.0f km^2\n', ...
    P.max_inner_refinement_area_km2);
fprintf('Fast mode  : %d | GEBCO shelf: %d\n', ...
    P.fast_blocking, P.use_gebco_shelf_width);
fprintf('============================================================\n');

assert(exist(P.domain_csv, 'file') == 2, 'Missing domain CSV: %s', P.domain_csv);
assert(exist(P.cell_csv, 'file') == 2, 'Missing cell CSV: %s', P.cell_csv);

T_domain = read_table_as_strings(P.domain_csv);
T_cell = read_table_as_strings(P.cell_csv);
T_domain = ensure_domain_basins(T_domain);
DomainPolys = read_domain_boundary_geojson(P.domain_geojson);

fprintf('\n[1/5] Grouping partition domains into ADCIRC blocks...\n');
[T_block, T_member, BlockGeom] = build_adcirc_blocks(T_domain, T_cell, P, DomainPolys);
writetable(T_block, fullfile(P.output_dir, 'global_tc_adcirc_blocks.csv'));
writetable(T_member, fullfile(P.output_dir, 'global_tc_adcirc_block_members.csv'));
write_block_range_summary(fullfile(P.output_dir, 'global_tc_adcirc_block_ranges.csv'), T_block, P);
fprintf('  -> blocks: %d\n', height(T_block));

fprintf('\n[2/5] Writing block boundary GeoJSON and per-block ring CSV files...\n');
write_block_geojson(fullfile(P.output_dir, 'global_tc_adcirc_block_outer_domains.geojson'), ...
    T_block, BlockGeom, 'outer');
write_block_geojson(fullfile(P.output_dir, 'global_tc_adcirc_block_inner_domains.geojson'), ...
    T_block, BlockGeom, 'inner');
write_block_ring_csvs(T_block, BlockGeom, P);

fprintf('\n[3/5] Writing optional OceanMesh2D runner...\n');
write_oceanmesh2d_runner(fullfile(P.output_dir, 'P2_run_global_tc_adcirc_block_meshes.m'), P);

fprintf('\n[4/5] Making block quicklook maps...\n');
make_global_block_quicklook(T_block, BlockGeom, P);
make_basin_block_quicklooks(T_block, BlockGeom, P);
make_block_panel_quicklook(T_block, BlockGeom, P);

fprintf('\n[5/5] Saving MAT summary...\n');
Results = struct();
Results.config = P;
Results.domainTable = T_domain;
Results.cellTable = T_cell;
Results.domainPolygons = DomainPolys;
Results.blockTable = T_block;
Results.memberTable = T_member;
Results.blockGeometry = BlockGeom;
save(fullfile(P.output_dir, 'global_tc_adcirc_model_blocks.mat'), 'Results', '-v7.3');
write_run_metadata(fullfile(P.output_dir, '00_run_metadata.txt'), P, T_block);

fprintf('\nDone. ADCIRC block products saved to:\n  %s\n', P.output_dir);


%% ============================================================
function [T_block, T_member, BlockGeom] = build_adcirc_blocks(T_domain, T_cell, P, DomainPolys)

basins = unique(T_domain.basin_id, 'stable');
rawBlocks = struct('basin_id', {}, 'basin_label', {}, 'domain_ids', {});

for ib = 1:numel(basins)
    bid = basins(ib);
    D = T_domain(T_domain.basin_id == bid, :);
    if isempty(D)
        continue;
    end

    [x, y, lon0, lat0] = lonlat_to_local_km(D.lon, D.lat);
    w = table_weights(D);
    axisU = principal_axis_2d(x, y, w);
    s = x .* axisU(1) + y .* axisU(2);
    [~, order] = sort(s);
    D = D(order, :);
    x = x(order);
    y = y(order);
    s = s(order);

    current = strings(0, 1);
    currentRows = [];
    for i = 1:height(D)
        proposed = [current; D.model_domain_id(i)];
        proposedRows = [currentRows; i]; %#ok<AGROW>
        splitNow = false;

        if ~isempty(current)
            gapKm = hypot(x(i) - x(i - 1), y(i) - y(i - 1));
            splitNow = gapKm > P.centroid_gap_split_km || ...
                numel(proposed) > P.max_domains_per_block || ...
                along_range_from_indices(s, proposedRows) > P.max_block_alongshore_km || ...
                span_from_indices(x, y, proposedRows) > P.max_block_span_km || ...
                cell_span_for_domain_ids(T_cell, proposed) > P.max_block_span_km;
            if ~splitNow && inner_refinement_area_limit_enabled(P)
                proposedInnerArea = candidate_inner_refinement_area_km2(proposed, T_domain, T_cell, P, DomainPolys);
                splitNow = proposedInnerArea > P.max_inner_refinement_area_km2;
            end
        end

        if splitNow
            rawBlocks(end + 1) = make_raw_block(D.basin_id(1), D.basin_label(1), current); %#ok<AGROW>
            current = D.model_domain_id(i);
            currentRows = i;
        else
            current = proposed;
            currentRows = proposedRows;
        end
    end

    if ~isempty(current)
        rawBlocks(end + 1) = make_raw_block(D.basin_id(1), D.basin_label(1), current); %#ok<AGROW>
    end
end

rawBlocks = merge_tiny_blocks(rawBlocks, T_domain, T_cell, P, DomainPolys);
rawBlocks = split_blocks_by_coastal_connectivity(rawBlocks, T_domain, DomainPolys, P);
rawBlocks = merge_overlapping_adc_blocks(rawBlocks, T_domain, T_cell, P, DomainPolys);
rawBlocks = absorb_small_nearby_adc_blocks(rawBlocks, T_domain, T_cell, P, DomainPolys);
rawBlocks = split_blocks_by_inner_area_limit(rawBlocks, T_domain, T_cell, P, DomainPolys);
rawBlocks = split_blocks_by_inner_ocean_components(rawBlocks, T_domain, P, DomainPolys);

preassignedBlockIds = assigned_adc_block_ids(rawBlocks);
if isfield(P, 'selected_block_ids') && ~isempty(P.selected_block_ids)
    keep = ismember(preassignedBlockIds, upper(string(P.selected_block_ids(:))));
    missing = setdiff(upper(string(P.selected_block_ids(:))), preassignedBlockIds, 'stable');
    if ~isempty(missing)
        fprintf('  !! selected block(s) not found after grouping: %s\n', strjoin(missing, ', '));
    end
    rawBlocks = rawBlocks(keep);
    preassignedBlockIds = preassignedBlockIds(keep);
    fprintf('  -> selected block debug mode: keeping %d block(s): %s\n', ...
        numel(rawBlocks), strjoin(preassignedBlockIds, ', '));
end

n = numel(rawBlocks);
BlockGeom = struct('block_id', {}, 'outer_lon', {}, 'outer_lat', {}, ...
    'inner_lon', {}, 'inner_lat', {}, 'shelf_lon', {}, 'shelf_lat', {}, ...
    'bay_lon', {}, 'bay_lat', {});

block_id = strings(n, 1);
basin_id = strings(n, 1);
basin_label = strings(n, 1);
rank_block = (1:n).';
domain_count = zeros(n, 1);
grid_cell_count = zeros(n, 1);
lowlying_area_km2 = zeros(n, 1);
priority_max = nan(n, 1);
tc_landfall_count_max = nan(n, 1);
population_sum_2020 = nan(n, 1);
lon = nan(n, 1);
lat = nan(n, 1);
inner_area_km2 = nan(n, 1);
outer_area_km2 = nan(n, 1);
inner_along_half_km = nan(n, 1);
inner_cross_half_km = nan(n, 1);
outer_along_half_km = nan(n, 1);
outer_cross_half_km = nan(n, 1);
shelf_width_km = nan(n, 1);
outer_buffer_km = nan(n, 1);
shelf_add_area_km2 = nan(n, 1);
shelf_add_cell_count = zeros(n, 1);
shelf_add_component_count = zeros(n, 1);
bay_add_area_km2 = nan(n, 1);
bay_add_cell_count = zeros(n, 1);
bay_add_core_cell_count = zeros(n, 1);
bay_add_component_count = zeros(n, 1);
block_span_km = nan(n, 1);
open_ocean_area_risk = false(n, 1);
auto_alongshore_radius_deg = nan(n, 1);
auto_offshore_deg = nan(n, 1);
auto_inland_deg = nan(n, 1);
auto_bay_close_deg = nan(n, 1);
auto_ellipse_ratio = nan(n, 1);
auto_inner_radius = nan(n, 1);
member_domains = strings(n, 1);

memberRows = struct('block_id', {}, 'basin_id', {}, 'basin_label', {}, ...
    'model_domain_id', {}, 'rank_domain', {}, 'domain_lon', {}, 'domain_lat', {}, ...
    'domain_priority', {}, 'domain_grid_cell_count', {});

for iblk = 1:n
    did = rawBlocks(iblk).domain_ids(:);
    D = T_domain(ismember(T_domain.model_domain_id, did), :);
    C = T_cell(ismember(T_cell.model_domain_id, did), :);

    bid = preassignedBlockIds(iblk);

    fprintf('  -> final geometry %d / %d: %s\n', iblk, n, char(bid));
    [geom, metrics] = make_block_geometry_cached(D, C, P, DomainPolys);
    if inner_refinement_area_limit_enabled(P) && metrics.inner_area_km2 > P.max_inner_refinement_area_km2
        fprintf(['  !! %s inner refinement area %.0f km^2 exceeds cap %.0f km^2; ', ...
            'this usually means at least one source partition is already too large.\n'], ...
            char(bid), metrics.inner_area_km2, P.max_inner_refinement_area_km2);
    end

    block_id(iblk) = bid;
    basin_id(iblk) = rawBlocks(iblk).basin_id;
    basin_label(iblk) = rawBlocks(iblk).basin_label;
    domain_count(iblk) = height(D);
    grid_cell_count(iblk) = height(C);
    lowlying_area_km2(iblk) = table_sum_if_present(D, 'lowlying_area_km2');
    priority_max(iblk) = table_max_if_present(D, 'priority');
    tc_landfall_count_max(iblk) = table_max_if_present(D, 'tc_landfall_count_max');
    population_sum_2020(iblk) = table_sum_if_present(D, 'population_sum_2020');
    lon(iblk) = metrics.lon_center;
    lat(iblk) = metrics.lat_center;
    inner_area_km2(iblk) = metrics.inner_area_km2;
    outer_area_km2(iblk) = metrics.outer_area_km2;
    inner_along_half_km(iblk) = metrics.inner_along_half_km;
    inner_cross_half_km(iblk) = metrics.inner_cross_half_km;
    outer_along_half_km(iblk) = metrics.outer_along_half_km;
    outer_cross_half_km(iblk) = metrics.outer_cross_half_km;
    shelf_width_km(iblk) = metrics.shelf_width_km;
    outer_buffer_km(iblk) = metrics.outer_buffer_km;
    shelf_add_area_km2(iblk) = metrics.shelf_add_area_km2;
    shelf_add_cell_count(iblk) = metrics.shelf_add_cell_count;
    shelf_add_component_count(iblk) = metrics.shelf_add_component_count;
    bay_add_area_km2(iblk) = metrics.bay_add_area_km2;
    bay_add_cell_count(iblk) = metrics.bay_add_cell_count;
    bay_add_core_cell_count(iblk) = metrics.bay_add_core_cell_count;
    bay_add_component_count(iblk) = metrics.bay_add_component_count;
    block_span_km(iblk) = metrics.block_span_km;
    open_ocean_area_risk(iblk) = metrics.outer_area_km2 > P.outer_open_ocean_max_area_km2;
    auto_alongshore_radius_deg(iblk) = metrics.auto_alongshore_radius_deg;
    auto_offshore_deg(iblk) = metrics.auto_offshore_deg;
    auto_inland_deg(iblk) = P.auto_inland_deg;
    auto_bay_close_deg(iblk) = P.auto_bay_close_deg;
    auto_ellipse_ratio(iblk) = metrics.auto_ellipse_ratio;
    auto_inner_radius(iblk) = metrics.auto_inner_radius;
    member_domains(iblk) = strjoin(did, ';');

    BlockGeom(iblk).block_id = bid; %#ok<AGROW>
    BlockGeom(iblk).outer_lon = geom.outer_lon(:);
    BlockGeom(iblk).outer_lat = geom.outer_lat(:);
    BlockGeom(iblk).inner_lon = geom.inner_lon(:);
    BlockGeom(iblk).inner_lat = geom.inner_lat(:);
    BlockGeom(iblk).shelf_lon = geom.shelf_lon(:);
    BlockGeom(iblk).shelf_lat = geom.shelf_lat(:);
    BlockGeom(iblk).bay_lon = geom.bay_lon(:);
    BlockGeom(iblk).bay_lat = geom.bay_lat(:);

    for i = 1:height(D)
        r = struct();
        r.block_id = bid;
        r.basin_id = rawBlocks(iblk).basin_id;
        r.basin_label = rawBlocks(iblk).basin_label;
        r.model_domain_id = D.model_domain_id(i);
        r.rank_domain = D.rank_domain(i);
        r.domain_lon = D.lon(i);
        r.domain_lat = D.lat(i);
        if ismember('priority', D.Properties.VariableNames)
            r.domain_priority = D.priority(i);
        else
            r.domain_priority = NaN;
        end
        if ismember('grid_cell_count', D.Properties.VariableNames)
            r.domain_grid_cell_count = D.grid_cell_count(i);
        else
            r.domain_grid_cell_count = NaN;
        end
        memberRows(end + 1) = r; %#ok<AGROW>
    end
end

T_block = table(block_id, rank_block, basin_id, basin_label, lon, lat, ...
    domain_count, grid_cell_count, lowlying_area_km2, priority_max, ...
    tc_landfall_count_max, population_sum_2020, ...
    inner_area_km2, outer_area_km2, block_span_km, ...
    inner_along_half_km, inner_cross_half_km, ...
    outer_along_half_km, outer_cross_half_km, shelf_width_km, outer_buffer_km, ...
    shelf_add_area_km2, shelf_add_cell_count, shelf_add_component_count, ...
    bay_add_area_km2, bay_add_cell_count, bay_add_core_cell_count, bay_add_component_count, ...
    open_ocean_area_risk, ...
    auto_alongshore_radius_deg, auto_offshore_deg, auto_inland_deg, ...
    auto_bay_close_deg, auto_ellipse_ratio, auto_inner_radius, member_domains);

[T_block, sortIdx] = sortrows(T_block, {'basin_id', 'rank_block'});
BlockGeom = BlockGeom(sortIdx);
T_block.rank_block = (1:height(T_block)).';

if isempty(memberRows)
    T_member = table();
else
    T_member = struct2table(memberRows);
    T_member.block_id = string(T_member.block_id);
    T_member.basin_id = string(T_member.basin_id);
    T_member.basin_label = string(T_member.basin_label);
    T_member.model_domain_id = string(T_member.model_domain_id);
end
end

%% ============================================================
function raw = make_raw_block(basin_id, basin_label, domain_ids)

raw = struct();
raw.basin_id = string(basin_id);
raw.basin_label = string(basin_label);
raw.domain_ids = string(domain_ids(:));
end

%% ============================================================
function blockIds = assigned_adc_block_ids(rawBlocks)

n = numel(rawBlocks);
blockIds = strings(n, 1);
for i = 1:n
    short = basin_short_name(rawBlocks(i).basin_id);
    seq = 1 + nnz(string({rawBlocks(1:i-1).basin_id}) == rawBlocks(i).basin_id);
    blockIds(i) = sprintf("ADC_%s_%02d", short, seq);
end
end

%% ============================================================
function rawBlocks = merge_tiny_blocks(rawBlocks, T_domain, T_cell, P, DomainPolys)

if numel(rawBlocks) < 2
    return;
end

changed = true;
while changed
    changed = false;
    for i = 1:numel(rawBlocks)
        did = rawBlocks(i).domain_ids;
        cells = nnz(ismember(T_cell.model_domain_id, did));
        if numel(did) >= P.tiny_block_min_domains && cells >= P.tiny_block_min_cells
            continue;
        end

        same = find(string({rawBlocks.basin_id}) == rawBlocks(i).basin_id);
        same = same(same ~= i);
        if isempty(same)
            continue;
        end

        ci = block_centroid(rawBlocks(i), T_domain);
        best = 0;
        bestD = inf;
        for j = same(:).'
            cj = block_centroid(rawBlocks(j), T_domain);
            d = local_distance_km(ci(1), ci(2), cj(1), cj(2));
            candidateIds = unique([rawBlocks(j).domain_ids; rawBlocks(i).domain_ids], 'stable');
            if inner_refinement_area_limit_enabled(P)
                candidateArea = candidate_inner_refinement_area_km2(candidateIds, T_domain, T_cell, P, DomainPolys);
                if candidateArea > P.max_inner_refinement_area_km2
                    continue;
                end
            end
            if d < bestD
                bestD = d;
                best = j;
            end
        end
        if best == 0
            continue;
        end

        rawBlocks(best).domain_ids = unique([rawBlocks(best).domain_ids; rawBlocks(i).domain_ids], 'stable');
        rawBlocks(i) = [];
        changed = true;
        break;
    end
end
end

%% ============================================================
function rawBlocks = merge_overlapping_adc_blocks(rawBlocks, T_domain, T_cell, P, DomainPolys)

if ~isfield(P, 'merge_overlapping_blocks') || ~P.merge_overlapping_blocks || numel(rawBlocks) < 2
    return;
end

mergeCount = 0;
while true
    n = numel(rawBlocks);
    geomCache = cell(n, 1);
    metricsCache = cell(n, 1);

    for i = 1:n
        did = rawBlocks(i).domain_ids(:);
        D = T_domain(ismember(T_domain.model_domain_id, did), :);
        C = T_cell(ismember(T_cell.model_domain_id, did), :);
        [geomCache{i}, metricsCache{i}] = make_block_geometry_cached_for_merge(D, C, P, DomainPolys);
    end

    bestPair = [];
    bestScore = -inf;
    bestStats = struct();
    bestMetrics = struct();

    for i = 1:n-1
        for j = i+1:n
            if rawBlocks(i).basin_id ~= rawBlocks(j).basin_id
                continue;
            end

            coastalGapKm = raw_blocks_coastal_gap_km(rawBlocks(i), rawBlocks(j), T_domain, DomainPolys);
            if coastalGapKm > P.merge_coastal_gap_km
                continue;
            end

            stats = outer_domain_relation(geomCache{i}, geomCache{j}, metricsCache{i}, metricsCache{j});
            stats.coastal_gap_km = coastalGapKm;
            mergeByOverlap = stats.overlap_min_fraction >= P.merge_outer_overlap_min_fraction || ...
                stats.overlap_union_fraction >= P.merge_outer_union_min_fraction;
            mergeByTouch = stats.gap_km <= P.merge_outer_gap_km;
            if ~(mergeByOverlap || mergeByTouch)
                continue;
            end

            candidateIds = unique([rawBlocks(i).domain_ids; rawBlocks(j).domain_ids], 'stable');
            if numel(candidateIds) > P.merge_candidate_max_domains
                continue;
            end
            if cell_span_for_domain_ids(T_cell, candidateIds) > P.merge_candidate_max_span_km
                continue;
            end

            Dm = T_domain(ismember(T_domain.model_domain_id, candidateIds), :);
            Cm = T_cell(ismember(T_cell.model_domain_id, candidateIds), :);
            [~, candidateMetrics] = make_block_geometry_cached_for_merge(Dm, Cm, P, DomainPolys);
            if candidateMetrics.inner_area_km2 > P.max_inner_refinement_area_km2 || ...
                    candidateMetrics.outer_area_km2 > P.merge_candidate_max_outer_area_km2 || ...
                    candidateMetrics.block_span_km > P.merge_candidate_max_span_km
                continue;
            end

            touchScore = max(0, 1 - stats.gap_km / max(P.merge_outer_gap_km, eps));
            score = 100 * stats.overlap_min_fraction + ...
                70 * stats.overlap_union_fraction + ...
                20 * touchScore - ...
                0.08 * coastalGapKm - ...
                2 * candidateMetrics.outer_area_km2 / P.merge_candidate_max_outer_area_km2 - ...
                1.2 * candidateMetrics.inner_area_km2 / P.max_inner_refinement_area_km2 - ...
                0.01 * stats.center_distance_km;

            if score > bestScore
                bestScore = score;
                bestPair = [i, j]; %#ok<AGROW>
                bestStats = stats;
                bestMetrics = candidateMetrics;
            end
        end
    end

    if isempty(bestPair)
        break;
    end

    i = bestPair(1);
    j = bestPair(2);
    fprintf(['  -> merging %s-side adjacent blocks in %s: %d + %d ', ...
        '(outer overlap %.2f, union overlap %.2f, outer gap %.0f km, coastal gap %.0f km, merged inner %.2f Mkm2, outer %.2f Mkm2)\n'], ...
        char(basin_short_name(rawBlocks(i).basin_id)), char(rawBlocks(i).basin_id), i, j, ...
        bestStats.overlap_min_fraction, bestStats.overlap_union_fraction, ...
        bestStats.gap_km, bestStats.coastal_gap_km, ...
        bestMetrics.inner_area_km2 / 1e6, bestMetrics.outer_area_km2 / 1e6);

    rawBlocks(i).domain_ids = unique([rawBlocks(i).domain_ids; rawBlocks(j).domain_ids], 'stable');
    rawBlocks(j) = [];
    mergeCount = mergeCount + 1;
end

if mergeCount > 0
    fprintf('  -> second-pass outer-domain merges: %d, remaining blocks before final numbering: %d\n', ...
        mergeCount, numel(rawBlocks));
end
end

%% ============================================================
function rawBlocks = split_blocks_by_coastal_connectivity(rawBlocks, T_domain, DomainPolys, P)

if numel(rawBlocks) < 1
    return;
end

newBlocks = struct('basin_id', {}, 'basin_label', {}, 'domain_ids', {});
splitCount = 0;

for ib = 1:numel(rawBlocks)
    ids = string(rawBlocks(ib).domain_ids(:));
    if numel(ids) <= 1
        newBlocks(end + 1) = rawBlocks(ib); %#ok<AGROW>
        continue;
    end

    comps = connected_domain_components(ids, T_domain, DomainPolys, P.coastal_connect_gap_km);
    if numel(comps) <= 1
        newBlocks(end + 1) = rawBlocks(ib); %#ok<AGROW>
        continue;
    end

    splitCount = splitCount + numel(comps) - 1;
    for ic = 1:numel(comps)
        newBlocks(end + 1) = make_raw_block(rawBlocks(ib).basin_id, rawBlocks(ib).basin_label, comps{ic}); %#ok<AGROW>
    end
end

rawBlocks = newBlocks;
if splitCount > 0
    fprintf('  -> split non-contiguous preliminary blocks by coastal gap: +%d blocks\n', splitCount);
end
end

%% ============================================================
function comps = connected_domain_components(ids, T_domain, DomainPolys, maxGapKm)

ids = string(ids(:));
n = numel(ids);
if n <= 1
    comps = {ids};
    return;
end

D = T_domain(ismember(T_domain.model_domain_id, ids), :);
w = table_weights(D);
lon0 = weighted_circular_mean_lon(D.lon, w);
lat0 = weighted_mean(D.lat, w);

Pdom = cell(n, 1);
centers = nan(n, 2);
for i = 1:n
    Di = T_domain(T_domain.model_domain_id == ids(i), :);
    centers(i, :) = [weighted_circular_mean_lon(Di.lon, table_weights(Di)), weighted_mean(Di.lat, table_weights(Di))];
    Pdom{i} = domain_polygons_to_local_polyshape(ids(i), DomainPolys, lon0, lat0);
end

adj = false(n, n);
for i = 1:n-1
    for j = i+1:n
        gapKm = inf;
        if area(Pdom{i}) > 0 && area(Pdom{j}) > 0
            gapKm = polyshape_gap_km(Pdom{i}, Pdom{j});
        end
        if ~isfinite(gapKm)
            gapKm = local_distance_km(centers(i, 1), centers(i, 2), centers(j, 1), centers(j, 2));
        end
        if gapKm <= maxGapKm
            adj(i, j) = true;
            adj(j, i) = true;
        end
    end
end

seen = false(n, 1);
comps = {};
for i = 1:n
    if seen(i)
        continue;
    end
    queue = i;
    seen(i) = true;
    members = [];
    while ~isempty(queue)
        q = queue(1);
        queue(1) = [];
        members(end + 1, 1) = q; %#ok<AGROW>
        nb = find(adj(q, :) & ~seen.');
        if ~isempty(nb)
            seen(nb) = true;
            queue = [queue; nb(:)]; %#ok<AGROW>
        end
    end
    comps{end + 1} = ids(members); %#ok<AGROW>
end
end

function gapKm = raw_blocks_coastal_gap_km(raw1, raw2, T_domain, DomainPolys)

ids1 = string(raw1.domain_ids(:));
ids2 = string(raw2.domain_ids(:));
D = T_domain(ismember(T_domain.model_domain_id, [ids1; ids2]), :);
w = table_weights(D);
lon0 = weighted_circular_mean_lon(D.lon, w);
lat0 = weighted_mean(D.lat, w);

p1 = domain_polygons_to_local_polyshape(ids1, DomainPolys, lon0, lat0);
p2 = domain_polygons_to_local_polyshape(ids2, DomainPolys, lon0, lat0);
if area(p1) > 0 && area(p2) > 0
    gapKm = polyshape_gap_km(p1, p2);
else
    c1 = block_centroid(raw1, T_domain);
    c2 = block_centroid(raw2, T_domain);
    gapKm = local_distance_km(c1(1), c1(2), c2(1), c2(2));
end
end

%% ============================================================
function stats = outer_domain_relation(geom1, geom2, metrics1, metrics2)

lon0 = weighted_circular_mean_lon([metrics1.lon_center; metrics2.lon_center], [1; 1]);
lat0 = mean([metrics1.lat_center; metrics2.lat_center], 'omitnan');

p1 = lonlat_vectors_to_local_polyshape(geom1.outer_lon, geom1.outer_lat, lon0, lat0);
p2 = lonlat_vectors_to_local_polyshape(geom2.outer_lon, geom2.outer_lat, lon0, lat0);
a1 = area(p1);
a2 = area(p2);

stats = struct();
stats.overlap_min_fraction = 0;
stats.overlap_union_fraction = 0;
stats.gap_km = inf;
stats.center_distance_km = local_distance_km(metrics1.lon_center, metrics1.lat_center, ...
    metrics2.lon_center, metrics2.lat_center);

if a1 <= 0 || a2 <= 0
    return;
end

try
    pInt = intersect(p1, p2);
    aInt = area(pInt);
    pUnion = union(p1, p2);
    aUnion = area(pUnion);
catch
    aInt = 0;
    aUnion = a1 + a2;
end

if aInt > 0
    stats.overlap_min_fraction = aInt / max(eps, min(a1, a2));
    stats.overlap_union_fraction = aInt / max(eps, aUnion);
    stats.gap_km = 0;
else
    stats.gap_km = polyshape_gap_km(p1, p2);
end
end

%% ============================================================
function d = polyshape_gap_km(p1, p2)

try
    if area(intersect(p1, p2)) > 0
        d = 0;
        return;
    end
catch
end

[x1, y1] = boundary(p1);
[x2, y2] = boundary(p2);
v1 = isfinite(x1) & isfinite(y1);
v2 = isfinite(x2) & isfinite(y2);
if nnz(v1) < 1 || nnz(v2) < 1
    d = inf;
    return;
end
d12 = min_distance_to_vertices_km(x1(v1), y1(v1), x2(v2), y2(v2));
d21 = min_distance_to_vertices_km(x2(v2), y2(v2), x1(v1), y1(v1));
d = min([d12(:); d21(:)], [], 'omitnan');
if ~isfinite(d)
    d = inf;
end
end

%% ============================================================
function rawBlocks = absorb_small_nearby_adc_blocks(rawBlocks, T_domain, T_cell, P, DomainPolys)

if ~isfield(P, 'absorb_small_blocks') || ~P.absorb_small_blocks || numel(rawBlocks) < 2
    return;
end

absorbCount = 0;
while true
    n = numel(rawBlocks);
    geomCache = cell(n, 1);
    metricsCache = cell(n, 1);
    cells = zeros(n, 1);
    domains = zeros(n, 1);
    isSmall = false(n, 1);
    isMajor = false(n, 1);

    for i = 1:n
        did = rawBlocks(i).domain_ids(:);
        D = T_domain(ismember(T_domain.model_domain_id, did), :);
        C = T_cell(ismember(T_cell.model_domain_id, did), :);
        [geomCache{i}, metricsCache{i}] = make_block_geometry_cached_for_merge(D, C, P, DomainPolys);
        cells(i) = height(C);
        domains(i) = numel(did);
        isSmall(i) = cells(i) <= P.small_block_max_cells || ...
            (domains(i) <= P.small_block_max_domains && ...
            metricsCache{i}.inner_area_km2 <= P.small_block_max_inner_area_km2) || ...
            (cells(i) <= P.small_block_cluster_max_cells && ...
            domains(i) <= P.small_block_cluster_max_domains && ...
            metricsCache{i}.inner_area_km2 <= P.small_block_cluster_max_inner_area_km2);
        isMajor(i) = cells(i) >= P.small_major_min_cells;
    end

    bestPair = [];
    bestScore = -inf;
    bestStats = struct();
    bestMetrics = struct();
    bestMode = "";

    for i = 1:n-1
        for j = i+1:n
            if rawBlocks(i).basin_id ~= rawBlocks(j).basin_id
                continue;
            end
            if ~(isSmall(i) || isSmall(j))
                continue;
            end

            mode = "";
            if (isSmall(i) && isMajor(j)) || (isSmall(j) && isMajor(i))
                mode = "small_to_major";
            elseif isSmall(i) && isSmall(j)
                mode = "small_to_small";
            else
                continue;
            end

            outerStats = outer_domain_relation(geomCache{i}, geomCache{j}, metricsCache{i}, metricsCache{j});
            innerGapKm = inner_domain_gap_km(geomCache{i}, geomCache{j}, metricsCache{i}, metricsCache{j});
            coastalGapKm = raw_blocks_coastal_gap_km(rawBlocks(i), rawBlocks(j), T_domain, DomainPolys);

            if mode == "small_to_major"
                ok = innerGapKm <= P.small_merge_to_major_gap_km || ...
                    coastalGapKm <= P.small_merge_to_major_gap_km || ...
                    outerStats.overlap_min_fraction >= P.small_merge_outer_overlap_min_fraction || ...
                    outerStats.overlap_union_fraction >= P.small_merge_outer_union_min_fraction;
            else
                ok = innerGapKm <= P.small_merge_small_gap_km && ...
                    (outerStats.overlap_min_fraction >= 0.08 || ...
                    outerStats.overlap_union_fraction >= 0.035 || ...
                    outerStats.gap_km <= P.merge_outer_gap_km);
            end
            if ~ok
                continue;
            end

            candidateIds = unique([rawBlocks(i).domain_ids; rawBlocks(j).domain_ids], 'stable');
            if numel(candidateIds) > P.small_merge_candidate_max_domains
                continue;
            end
            if cell_span_for_domain_ids(T_cell, candidateIds) > P.small_merge_candidate_max_span_km
                continue;
            end

            Dm = T_domain(ismember(T_domain.model_domain_id, candidateIds), :);
            Cm = T_cell(ismember(T_cell.model_domain_id, candidateIds), :);
            [~, candidateMetrics] = make_block_geometry_cached_for_merge(Dm, Cm, P, DomainPolys);
            if candidateMetrics.inner_area_km2 > P.max_inner_refinement_area_km2 || ...
                    candidateMetrics.outer_area_km2 > P.small_merge_candidate_max_outer_area_km2 || ...
                    candidateMetrics.block_span_km > P.small_merge_candidate_max_span_km
                continue;
            end

            if mode == "small_to_major"
                score = 80 + 80 * outerStats.overlap_min_fraction + ...
                    60 * outerStats.overlap_union_fraction - ...
                    0.05 * min(innerGapKm, coastalGapKm) - ...
                    1.5 * candidateMetrics.outer_area_km2 / P.small_merge_candidate_max_outer_area_km2 - ...
                    1.0 * candidateMetrics.inner_area_km2 / P.max_inner_refinement_area_km2;
            else
                score = 45 + 70 * outerStats.overlap_min_fraction + ...
                    55 * outerStats.overlap_union_fraction - ...
                    0.035 * innerGapKm - ...
                    1.5 * candidateMetrics.outer_area_km2 / P.small_merge_candidate_max_outer_area_km2 - ...
                    1.0 * candidateMetrics.inner_area_km2 / P.max_inner_refinement_area_km2;
            end

            if score > bestScore
                bestPair = [i, j]; %#ok<AGROW>
                bestScore = score;
                bestStats = outerStats;
                bestStats.inner_gap_km = innerGapKm;
                bestStats.coastal_gap_km = coastalGapKm;
                bestMetrics = candidateMetrics;
                bestMode = mode;
            end
        end
    end

    if isempty(bestPair)
        break;
    end

    i = bestPair(1);
    j = bestPair(2);
    fprintf(['  -> absorbing small nearby block in %s: %d + %d (%s, ', ...
        'outer overlap %.2f, inner gap %.0f km, coastal gap %.0f km, merged inner %.2f Mkm2, outer %.2f Mkm2)\n'], ...
        char(rawBlocks(i).basin_id), i, j, char(bestMode), ...
        bestStats.overlap_min_fraction, bestStats.inner_gap_km, bestStats.coastal_gap_km, ...
        bestMetrics.inner_area_km2 / 1e6, bestMetrics.outer_area_km2 / 1e6);

    rawBlocks(i).domain_ids = unique([rawBlocks(i).domain_ids; rawBlocks(j).domain_ids], 'stable');
    rawBlocks(j) = [];
    absorbCount = absorbCount + 1;
end

if absorbCount > 0
    fprintf('  -> final small-block absorptions: %d, remaining blocks before final numbering: %d\n', ...
        absorbCount, numel(rawBlocks));
end
end

%% ============================================================
function rawBlocks = split_blocks_by_inner_area_limit(rawBlocks, T_domain, T_cell, P, DomainPolys)

if ~inner_refinement_area_limit_enabled(P) || numel(rawBlocks) < 1
    return;
end

newBlocks = struct('basin_id', {}, 'basin_label', {}, 'domain_ids', {});
splitCount = 0;

for ib = 1:numel(rawBlocks)
    ids = string(rawBlocks(ib).domain_ids(:));
    if numel(ids) <= 1
        newBlocks(end + 1) = rawBlocks(ib); %#ok<AGROW>
        continue;
    end

    innerArea = candidate_inner_refinement_area_km2(ids, T_domain, T_cell, P, DomainPolys);
    if innerArea <= P.max_inner_refinement_area_km2
        newBlocks(end + 1) = rawBlocks(ib); %#ok<AGROW>
        continue;
    end

    D = T_domain(ismember(T_domain.model_domain_id, ids), :);
    [x, y] = lonlat_to_local_km(D.lon, D.lat);
    w = table_weights(D);
    axisU = principal_axis_2d(x, y, w);
    s = x .* axisU(1) + y .* axisU(2);
    [~, order] = sort(s);
    D = D(order, :);

    parts = {};
    current = strings(0, 1);
    for k = 1:height(D)
        proposed = [current; D.model_domain_id(k)];
        if ~isempty(current)
            proposedArea = candidate_inner_refinement_area_km2(proposed, T_domain, T_cell, P, DomainPolys);
            if proposedArea > P.max_inner_refinement_area_km2
                parts{end + 1, 1} = current; %#ok<AGROW>
                current = D.model_domain_id(k);
            else
                current = proposed;
            end
        else
            current = proposed;
        end
    end
    if ~isempty(current)
        parts{end + 1, 1} = current; %#ok<AGROW>
    end

    if numel(parts) <= 1
        newBlocks(end + 1) = rawBlocks(ib); %#ok<AGROW>
        fprintf(['  !! inner refinement area %.0f km^2 exceeds cap %.0f km^2 ', ...
            'but block has no splittable multi-domain grouping in %s\n'], ...
            innerArea, P.max_inner_refinement_area_km2, char(rawBlocks(ib).basin_id));
        continue;
    end

    splitCount = splitCount + numel(parts) - 1;
    fprintf('  -> split oversized inner-refinement block in %s: %.2f Mkm2 -> %d blocks\n', ...
        char(rawBlocks(ib).basin_id), innerArea / 1e6, numel(parts));
    for ip = 1:numel(parts)
        newBlocks(end + 1) = make_raw_block(rawBlocks(ib).basin_id, rawBlocks(ib).basin_label, parts{ip}); %#ok<AGROW>
    end
end

if splitCount > 0
    fprintf('  -> final inner-refinement area splits: +%d blocks, final preliminary blocks: %d\n', ...
        splitCount, numel(newBlocks));
end
rawBlocks = newBlocks;
end

%% ============================================================
function rawBlocks = split_blocks_by_inner_ocean_components(rawBlocks, T_domain, P, DomainPolys)

if ~isfield(P, 'split_blocks_by_inner_ocean_components') || ~P.split_blocks_by_inner_ocean_components || ...
        ~isfield(P, 'gebco_path') || exist(P.gebco_path, 'file') ~= 2 || isempty(rawBlocks)
    return;
end

newBlocks = struct('basin_id', {}, 'basin_label', {}, 'domain_ids', {});
splitCount = 0;

for ib = 1:numel(rawBlocks)
    ids = string(rawBlocks(ib).domain_ids(:));
    if numel(ids) <= 1
        newBlocks(end + 1) = rawBlocks(ib); %#ok<AGROW>
        continue;
    end

    D = T_domain(ismember(T_domain.model_domain_id, ids), :);
    if isempty(D)
        newBlocks(end + 1) = rawBlocks(ib); %#ok<AGROW>
        continue;
    end

    wD = table_weights(D);
    lon0 = weighted_circular_mean_lon(D.lon, wD);
    lat0 = weighted_mean(D.lat, wD);
    studyPoly = domain_polygons_to_local_polyshape(D.model_domain_id, DomainPolys, lon0, lat0);
    if area(studyPoly) == 0
        newBlocks(end + 1) = rawBlocks(ib); %#ok<AGROW>
        continue;
    end

    innerPoly = rmholes(union(safe_polybuffer(studyPoly, P.inner_refinement_buffer_km)));
    components = inner_ocean_component_polyshapes(innerPoly, lon0, lat0, P);
    if numel(components) <= 1
        newBlocks(end + 1) = rawBlocks(ib); %#ok<AGROW>
        continue;
    end

    compIdx = assign_domains_to_ocean_components(D, DomainPolys, components, lon0, lat0);
    usedComponents = unique(compIdx(compIdx > 0), 'stable');
    if numel(usedComponents) <= 1
        newBlocks(end + 1) = rawBlocks(ib); %#ok<AGROW>
        continue;
    end

    splitCount = splitCount + numel(usedComponents) - 1;
    fprintf('  -> split cross-sea inner block in %s: %d domains -> %d one-sea blocks\n', ...
        char(rawBlocks(ib).basin_id), numel(ids), numel(usedComponents));
    for ic = usedComponents(:).'
        groupIds = D.model_domain_id(compIdx == ic);
        if isempty(groupIds)
            continue;
        end
        newBlocks(end + 1) = make_raw_block(rawBlocks(ib).basin_id, rawBlocks(ib).basin_label, groupIds); %#ok<AGROW>
    end
end

if splitCount > 0
    fprintf('  -> final cross-sea inner splits: +%d blocks, final preliminary blocks: %d\n', ...
        splitCount, numel(newBlocks));
end
rawBlocks = newBlocks;
end

%% ============================================================
function compIdx = assign_domains_to_ocean_components(D, DomainPolys, components, lon0, lat0)

compIdx = zeros(height(D), 1);
for i = 1:height(D)
    did = string(D.model_domain_id(i));
    p = domain_polygons_to_local_polyshape(did, DomainPolys, lon0, lat0);
    scores = zeros(numel(components), 1);
    if area(p) > 0
        for ic = 1:numel(components)
            try
                scores(ic) = area(intersect(p, components{ic}));
            catch
                scores(ic) = 0;
            end
        end
    end
    [bestScore, best] = max(scores);
    if bestScore > 0
        compIdx(i) = best;
        continue;
    end

    [x, y] = lonlat_to_local_km(D.lon(i), D.lat(i), lon0, lat0);
    inside = false(numel(components), 1);
    for ic = 1:numel(components)
        try
            inside(ic) = isinterior(components{ic}, x, y);
        catch
            inside(ic) = false;
        end
    end
    if any(inside)
        compIdx(i) = find(inside, 1, 'first');
    else
        compIdx(i) = nearest_polyshape_component(x, y, components);
    end
end
end

%% ============================================================
function idx = nearest_polyshape_component(x, y, components)

idx = 0;
bestD = inf;
for ic = 1:numel(components)
    [xb, yb] = boundary(components{ic});
    v = isfinite(xb) & isfinite(yb);
    if ~any(v)
        continue;
    end
    d = min(hypot(xb(v) - x, yb(v) - y), [], 'omitnan');
    if d < bestD
        bestD = d;
        idx = ic;
    end
end
if idx == 0 && ~isempty(components)
    idx = 1;
end
end

%% ============================================================
function tf = inner_refinement_area_limit_enabled(P)

tf = isfield(P, 'max_inner_refinement_area_km2') && ...
    isfinite(P.max_inner_refinement_area_km2) && P.max_inner_refinement_area_km2 > 0;
end

%% ============================================================
function areaKm2 = candidate_inner_refinement_area_km2(domainIds, T_domain, T_cell, P, DomainPolys)

domainIds = string(domainIds(:));
D = T_domain(ismember(T_domain.model_domain_id, domainIds), :);
C = T_cell(ismember(T_cell.model_domain_id, domainIds), :);
if isempty(D)
    areaKm2 = 0;
    return;
end
[~, metrics] = make_block_geometry_cached_for_merge(D, C, P, DomainPolys);
areaKm2 = metrics.inner_area_km2;
end

%% ============================================================
function gapKm = inner_domain_gap_km(geom1, geom2, metrics1, metrics2)

lon0 = weighted_circular_mean_lon([metrics1.lon_center; metrics2.lon_center], [1; 1]);
lat0 = mean([metrics1.lat_center; metrics2.lat_center], 'omitnan');
p1 = lonlat_vectors_to_local_polyshape(geom1.inner_lon, geom1.inner_lat, lon0, lat0);
p2 = lonlat_vectors_to_local_polyshape(geom2.inner_lon, geom2.inner_lat, lon0, lat0);
if area(p1) > 0 && area(p2) > 0
    gapKm = polyshape_gap_km(p1, p2);
else
    gapKm = local_distance_km(metrics1.lon_center, metrics1.lat_center, ...
        metrics2.lon_center, metrics2.lat_center);
end
end

%% ============================================================
function [geom, metrics] = make_block_geometry_cached(D, C, P, DomainPolys)

persistent GeometryCache
if isempty(GeometryCache)
    GeometryCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end

key = geometry_cache_key(D, P);

if isKey(GeometryCache, key)
    rec = GeometryCache(key);
    geom = rec.geom;
    metrics = rec.metrics;
    return;
end

[geom, metrics] = make_block_geometry(D, C, P, DomainPolys);
GeometryCache(key) = struct('geom', geom, 'metrics', metrics);
end

%% ============================================================
function [geom, metrics] = make_block_geometry_cached_for_merge(D, C, P, DomainPolys)

P.use_gebco_topology = false;
[geom, metrics] = make_block_geometry_cached(D, C, P, DomainPolys);
end

%% ============================================================
function key = geometry_cache_key(D, P)

ids = sort(string(D.model_domain_id(:)));
key = char(strjoin(ids, ',') + "|cfg=" + string(config_signature(P)));
end

%% ============================================================
function sig = config_signature(P)

try
    txt = jsonencode(orderfields(P));
catch
    txt = evalc('disp(P)');
end
sig = simple_string_hash(txt);
end

%% ============================================================
function h = simple_string_hash(txt)

bytes = uint8(char(txt));
hash = uint32(2166136261);
prime = uint64(16777619);
modulus = uint64(4294967296);
for i = 1:numel(bytes)
    hash = bitxor(hash, uint32(bytes(i)));
    hash = uint32(mod(uint64(hash) * prime, modulus));
end
h = dec2hex(double(hash), 8);
end

%% ============================================================
function c = block_centroid(raw, T_domain)

D = T_domain(ismember(T_domain.model_domain_id, raw.domain_ids), :);
w = table_weights(D);
c = [weighted_circular_mean_lon(D.lon, w), weighted_mean(D.lat, w)];
end

%% ============================================================
function [geom, metrics] = make_block_geometry(D, C, P, DomainPolys)

wD = table_weights(D);
lon0 = weighted_circular_mean_lon(D.lon, wD);
lat0 = weighted_mean(D.lat, wD);

studyPoly = domain_polygons_to_local_polyshape(D.model_domain_id, DomainPolys, lon0, lat0);
if area(studyPoly) == 0
    studyPoly = cells_to_local_polyshape(C, lon0, lat0);
end
studyPoly = rmholes(union(studyPoly));

[vx, vy] = boundary(studyPoly);
valid = isfinite(vx) & isfinite(vy);
if nnz(valid) < 3
    [xFallback, yFallback] = lonlat_to_local_km(D.lon, D.lat, lon0, lat0);
    vx = xFallback;
    vy = yFallback;
    valid = isfinite(vx) & isfinite(vy);
end
axisU = principal_axis_2d(vx(valid), vy(valid), ones(nnz(valid), 1));
axisV = [-axisU(2); axisU(1)];

% Inner refinement domain: exact research partition geometry expanded 50 km
% in local kilometre coordinates.
innerPoly = safe_polybuffer(studyPoly, P.inner_refinement_buffer_km);
innerPoly = rmholes(union(innerPoly));

% Outer ADCIRC domain: morphological bay closure + shelf-aware offshore
% expansion.  This follows the idea of Global_autofunction: close small bays
% and gaps before expanding the ocean-side/open-boundary domain, but it uses
% the already-screened research coastal segment as the seed geometry.
macroPoly = safe_polybuffer(studyPoly, P.outer_bay_close_km);
macroPoly = safe_polybuffer(macroPoly, -P.outer_bay_close_km);
macroPoly = union(macroPoly, studyPoly);
macroPoly = rmholes(union(macroPoly));

shelfWidthKm = estimate_shelf_width_km(P, studyPoly, lon0, lat0);
outerBufferKm = max(P.outer_min_buffer_km, max(P.outer_base_buffer_km, shelfWidthKm + P.outer_shelf_margin_km));
outerBufferKm = min(P.outer_max_buffer_km, outerBufferKm);

outerPoly = safe_polybuffer(macroPoly, outerBufferKm);
outerPoly = limit_initial_outer_seed_alongshore(outerPoly, studyPoly, innerPoly, outerBufferKm, P);
outerPoly = bias_outer_seed_to_tc_approach(outerPoly, studyPoly, innerPoly, ...
    lon0, lat0, D.basin_id(1), outerBufferKm);
if P.outer_smooth_km > 0
    smoothed = safe_polybuffer(outerPoly, P.outer_smooth_km);
    smoothed = safe_polybuffer(smoothed, -P.outer_smooth_km);
    outerPoly = union(outerPoly, smoothed); % never shrink the computed domain
end
% Before land/topology clipping, force the preliminary outer domain to
% fully cover the inner refinement domain so the black region is never
% dropped by later sea-side partitioning.
outerPoly = ensure_outer_covers_inner_refinement(outerPoly, innerPoly);

topoStats = empty_topology_stats();

% Keep the offshore expansion broad, but use GEBCO topology to avoid
% crossing through land into the wrong ocean.  Fallback to the geometric
% backshore guard when local bathymetry cannot be sampled.
[topoOuter, topoOk, topoStats] = apply_gebco_topology_outer_adjustment(outerPoly, studyPoly, innerPoly, ...
    lon0, lat0, outerBufferKm, D.basin_id(1), P);
if topoOk
    outerPoly = topoOuter;
else
    corridorPoly = make_backshore_clip_corridor(outerPoly, studyPoly, innerPoly, D.basin_id(1), ...
        lon0, lat0, axisV, outerBufferKm, P);
    if area(corridorPoly) > 0
        clippedOuter = intersect(outerPoly, corridorPoly);
        if area(clippedOuter) > max(area(innerPoly) * 1.08, area(studyPoly) + 1)
            outerPoly = clippedOuter;
        end
    end
end

outerPoly = ensure_outer_covers_inner_refinement(outerPoly, innerPoly);
outerPoly = rmholes(union(outerPoly));
outerPoly = keep_largest_relevant_region(outerPoly, innerPoly, studyPoly, P);
outerPoly = smooth_final_output_outer_boundary(outerPoly, innerPoly, studyPoly, P);
outerPoly = ensure_outer_covers_inner_refinement(outerPoly, innerPoly);
outerPoly = rmholes(union(outerPoly));

[innerLon, innerLat] = local_polyshape_to_lonlat_vectors(innerPoly, lon0, lat0);
[outerLon, outerLat] = local_polyshape_to_lonlat_vectors(outerPoly, lon0, lat0);
[shelfLon, shelfLat] = local_polyshape_to_lonlat_vectors(topoStats.shelf_domain, lon0, lat0);
[bayLon, bayLat] = local_polyshape_to_lonlat_vectors(topoStats.bay_domain, lon0, lat0);

[xStudy, yStudy] = boundary(studyPoly);
[xInner, yInner] = boundary(innerPoly);
[xOuter, yOuter] = boundary(outerPoly);
spanAll = pairwise_span_km(xStudy(isfinite(xStudy)), yStudy(isfinite(yStudy)));
[innerAlong, innerCross] = half_extents_along_axis(xInner, yInner, axisU, axisV);
[outerAlong, outerCross] = half_extents_along_axis(xOuter, yOuter, axisU, axisV);

geom = struct();
geom.inner_lon = innerLon;
geom.inner_lat = innerLat;
geom.outer_lon = outerLon;
geom.outer_lat = outerLat;
geom.shelf_lon = shelfLon;
geom.shelf_lat = shelfLat;
geom.bay_lon = bayLon;
geom.bay_lat = bayLat;

metrics = struct();
metrics.lon_center = lon0;
metrics.lat_center = lat0;
metrics.inner_along_half_km = innerAlong;
metrics.inner_cross_half_km = innerCross;
metrics.outer_along_half_km = outerAlong;
metrics.outer_cross_half_km = outerCross;
metrics.block_span_km = spanAll;
metrics.inner_area_km2 = area(innerPoly);
metrics.outer_area_km2 = area(outerPoly);
metrics.shelf_width_km = shelfWidthKm;
metrics.outer_buffer_km = outerBufferKm;
metrics.shelf_add_area_km2 = topoStats.shelf_add_area_km2;
metrics.shelf_add_cell_count = topoStats.shelf_add_cell_count;
metrics.shelf_add_component_count = topoStats.shelf_add_component_count;
metrics.bay_add_area_km2 = topoStats.bay_add_area_km2;
metrics.bay_add_cell_count = topoStats.bay_add_cell_count;
metrics.bay_add_core_cell_count = topoStats.bay_add_core_cell_count;
metrics.bay_add_component_count = topoStats.bay_add_component_count;

metrics.auto_alongshore_radius_deg = max(P.auto_alongshore_min_deg, ...
    min(P.auto_alongshore_max_deg, outerAlong / 111.32));
metrics.auto_offshore_deg = max(P.auto_offshore_min_deg, ...
    min(P.auto_offshore_max_deg, outerCross / 111.32));
metrics.auto_inner_radius = max(P.auto_inner_radius_min_deg, ...
    min(P.auto_inner_radius_max_deg, innerAlong / 111.32));
metrics.auto_ellipse_ratio = max(0.65, min(1.55, outerCross / max(outerAlong, eps)));
end

%% ============================================================
function componentsPoly = inner_ocean_component_polyshapes(innerPoly, lon0, lat0, P)

componentsPoly = {};
if area(innerPoly) <= 0 || ~isfield(P, 'gebco_path') || exist(P.gebco_path, 'file') ~= 2
    return;
end

[xb, yb] = boundary(innerPoly);
v = isfinite(xb) & isfinite(yb);
if nnz(v) < 3
    return;
end

dx = 10;
if isfield(P, 'inner_ocean_split_grid_km')
    dx = max(2, double(P.inner_ocean_split_grid_km));
end
xmin = floor((min(xb(v)) - 2 * dx) / dx) * dx;
xmax = ceil((max(xb(v)) + 2 * dx) / dx) * dx;
ymin = floor((min(yb(v)) - 2 * dx) / dx) * dx;
ymax = ceil((max(yb(v)) + 2 * dx) / dx) * dx;
if ~(isfinite(xmin) && isfinite(xmax) && isfinite(ymin) && isfinite(ymax) && xmax > xmin && ymax > ymin)
    return;
end

maxCells = 2.5e5;
if isfield(P, 'inner_ocean_split_max_grid_cells')
    maxCells = max(1e4, double(P.inner_ocean_split_max_grid_cells));
end
nx0 = max(2, ceil((xmax - xmin) / dx) + 1);
ny0 = max(2, ceil((ymax - ymin) / dx) + 1);
if nx0 * ny0 > maxCells
    dx = dx * sqrt((nx0 * ny0) / maxCells);
end
xv = xmin:dx:xmax;
yv = ymin:dx:ymax;
if numel(xv) < 3 || numel(yv) < 3
    return;
end

[X, Y] = meshgrid(xv, yv);
innerMask = polyshape_grid_mask_local(innerPoly, X, Y);
if nnz(innerMask) < 4
    return;
end

[lonQ, latQ] = local_km_to_lonlat(X(:), Y(:), lon0, lat0);
z = sample_gebco_nearest(P.gebco_path, lonQ, latQ, dx);
if isempty(z)
    return;
end
z = reshape(z, size(X));
zmax = P.topology_ocean_zmax_m;
if isfield(P, 'inner_ocean_split_ocean_zmax_m')
    zmax = P.inner_ocean_split_ocean_zmax_m;
end
waterMask = innerMask & isfinite(z) & z <= zmax;

minArea = 500;
if isfield(P, 'inner_ocean_split_min_area_km2')
    minArea = max(0, double(P.inner_ocean_split_min_area_km2));
end
minCells = max(2, ceil(minArea / max(dx * dx, eps)));
waterMask = areaopen_mask_local(waterMask, minCells);
components = connected_components_pixelidx(waterMask);
if isempty(components)
    return;
end

componentSizes = cellfun(@numel, components(:));
[~, order] = sort(componentSizes, 'descend');
maxComponents = 12;
if isfield(P, 'inner_ocean_split_max_components')
    maxComponents = max(1, round(double(P.inner_ocean_split_max_components)));
end
order = order(1:min(numel(order), maxComponents));

padKm = 12;
if isfield(P, 'inner_ocean_split_component_pad_km')
    padKm = max(0, double(P.inner_ocean_split_component_pad_km));
end
padCells = max(0, round(padKm / dx));
smoothKm = 8;
if isfield(P, 'inner_ocean_split_smooth_km')
    smoothKm = max(0, double(P.inner_ocean_split_smooth_km));
end

for ii = 1:numel(order)
    compMask = false(size(waterMask));
    compMask(components{order(ii)}) = true;
    if padCells > 0
        compMask = dilate_mask_local(compMask, padCells) & innerMask;
    end
    compPoly = grid_mask_to_polyshape_local(compMask, xv, yv);
    if area(compPoly) <= 0
        continue;
    end
    compPoly = intersect(compPoly, innerPoly);
    if smoothKm > 0 && area(compPoly) > 0
        smoothed = safe_polybuffer(compPoly, smoothKm);
        smoothed = safe_polybuffer(smoothed, -0.75 * smoothKm);
        smoothed = intersect(smoothed, innerPoly);
        if area(smoothed) >= 0.55 * area(compPoly)
            compPoly = smoothed;
        end
    end
    regs = regions(compPoly);
    for ir = 1:numel(regs)
        if area(regs(ir)) >= minArea
            componentsPoly{end + 1, 1} = regs(ir); %#ok<AGROW>
        end
    end
end
end

%% ============================================================
function Pout = limit_initial_outer_seed_alongshore(outerPoly, studyPoly, innerPoly, outerBufferKm, P)

Pout = outerPoly;
if ~isfield(P, 'limit_outer_seed_alongshore') || ~P.limit_outer_seed_alongshore || ...
        area(outerPoly) <= 0 || area(studyPoly) <= 0
    return;
end

basePoly = rmholes(union(studyPoly, innerPoly));
[xb, yb] = boundary(basePoly);
v = isfinite(xb) & isfinite(yb);
if nnz(v) < 3
    return;
end

axisU = principal_axis_2d(xb(v), yb(v), ones(nnz(v), 1));
axisU = axisU(:) / max(norm(axisU), eps);
axisV = [-axisU(2); axisU(1)];

sBase = xb(v) .* axisU(1) + yb(v) .* axisU(2);
tBase = xb(v) .* axisV(1) + yb(v) .* axisV(2);
padS = 120;
if isfield(P, 'outer_seed_alongshore_pad_km')
    padS = max(0, P.outer_seed_alongshore_pad_km);
end

% This changes the seed used for ordinary coastline expansion.  A simple
% isotropic buffer around a long coastal strip produces large rounded end
% caps.  Build an anisotropic buffer instead: keep the full cross-shore
% shelf/bay width, but shorten only the alongshore end-cap radius.  The
% large bay-search step runs later and can still add a bay where the GEBCO
% topology says it is a real connected embayment.
s0 = min(sBase, [], 'omitnan');
s1 = max(sBase, [], 'omitnan');
t0 = min(tBase, [], 'omitnan');
t1 = max(tBase, [], 'omitnan');
if ~(isfinite(s0) && isfinite(s1) && isfinite(t0) && isfinite(t1) && s1 > s0 && t1 > t0)
    return;
end

try
    sCtr = 0.5 * (s0 + s1);
    alongScale = max(2.2, min(7.5, max(outerBufferKm, P.outer_min_buffer_km) / max(padS, 50)));
    baseScaled = polyshape();
    rr = regions(basePoly);
    for i = 1:numel(rr)
        [xr, yr] = boundary(rr(i));
        vr = isfinite(xr) & isfinite(yr);
        if nnz(vr) < 3
            continue;
        end
        sr = xr(vr) .* axisU(1) + yr(vr) .* axisU(2);
        tr = xr(vr) .* axisV(1) + yr(vr) .* axisV(2);
        srScaled = sCtr + (sr - sCtr) .* alongScale;
        baseScaled = union(baseScaled, polyshape(srScaled, tr, 'Simplify', true));
    end
    if area(baseScaled) <= 0
        return;
    end
    seedScaled = safe_polybuffer(baseScaled, max(outerBufferKm, P.outer_min_buffer_km));
    seedLimit = polyshape();
    rrSeed = regions(seedScaled);
    for i = 1:numel(rrSeed)
        [ss, tt] = boundary(rrSeed(i));
        vs = isfinite(ss) & isfinite(tt);
        if nnz(vs) < 3
            continue;
        end
        sBack = sCtr + (ss(vs) - sCtr) ./ alongScale;
        xBack = sBack .* axisU(1) + tt(vs) .* axisV(1);
        yBack = sBack .* axisU(2) + tt(vs) .* axisV(2);
        seedLimit = union(seedLimit, polyshape(xBack, yBack, 'Simplify', true));
    end
    seedLimit = rmholes(union(seedLimit));
    if area(seedLimit) <= 0
        return;
    end
    limited = intersect(outerPoly, seedLimit);
    limited = union(limited, innerPoly);
    limited = rmholes(union(limited));
catch
    return;
end

if area(limited) <= max(area(innerPoly) * 1.05, area(studyPoly) + 1)
    return;
end
try
    innerCoverage = area(intersect(limited, innerPoly)) / max(area(innerPoly), eps);
catch
    innerCoverage = 0;
end
if innerCoverage < 0.98
    return;
end

Pout = limited;
end

%% ============================================================
function Pout = bias_outer_seed_to_tc_approach(outerPoly, studyPoly, innerPoly, lon0, lat0, basinId, outerBufferKm)

Pout = outerPoly;
if area(outerPoly) <= 0 || area(studyPoly) <= 0
    return;
end

[refLon, refLat, known] = basin_tc_approach_reference_lonlat(basinId, lon0, lat0);
if ~known
    return;
end

[rx, ry] = lonlat_to_local_km(refLon, refLat, lon0, lat0);
approachVec = [rx; ry];
if ~all(isfinite(approachVec)) || norm(approachVec) < 1
    return;
end
approachVec = approachVec / norm(approachVec);
sideVec = [-approachVec(2); approachVec(1)];

basePoly = rmholes(union(studyPoly, innerPoly));
if area(basePoly) <= 0
    return;
end

shiftFractions = [0.35, 0.80, 1.25];
bufferFractions = [0.92, 0.82, 0.72];
lateralFractions = [0.00, 0.22, -0.22];
approachPoly = polyshape();
for i = 1:numel(shiftFractions)
    d = shiftFractions(i) * outerBufferKm;
    for j = 1:numel(lateralFractions)
        latShift = lateralFractions(j) * d;
        shifted = translate_polyshape_local(basePoly, ...
            approachVec(1) * d + sideVec(1) * latShift, ...
            approachVec(2) * d + sideVec(2) * latShift);
        if area(shifted) <= 0
            continue;
        end
        approachPoly = union(approachPoly, safe_polybuffer(shifted, bufferFractions(i) * outerBufferKm));
    end
end
approachPoly = rmholes(union(approachPoly));
if area(approachPoly) <= 0
    return;
end

Pout = rmholes(union(Pout, approachPoly));
end

%% ============================================================
function Pout = ensure_outer_covers_inner_refinement(Pin, innerPoly)

Pout = Pin;
if area(Pin) <= 0 || area(innerPoly) <= 0
    return;
end

try
    missingInner = subtract(innerPoly, Pin);
catch
    missingInner = polyshape();
end
if area(missingInner) <= 0
    return;
end

Pout = rmholes(union(Pin, missingInner));
end

%% ============================================================
function stats = empty_topology_stats()

stats = struct();
stats.shelf_domain = polyshape();
stats.shelf_add_area_km2 = 0;
stats.shelf_add_cell_count = 0;
stats.shelf_add_component_count = 0;
stats.bay_domain = polyshape();
stats.bay_add_area_km2 = 0;
stats.bay_add_cell_count = 0;
stats.bay_add_core_cell_count = 0;
stats.bay_add_component_count = 0;
end

%% ============================================================
function searchKm = topology_named_embayment_search_km(P, basinId)

searchKm = max(P.topology_major_bay_search_km, P.topology_bay_extra_km);
if isfield(P, 'topology_named_embayment_enable') && P.topology_named_embayment_enable && ...
        isfield(P, 'topology_named_embayment_extra_km') && ...
        P.topology_named_embayment_extra_km > 0 && ~isempty(named_embayment_specs(basinId, P))
    searchKm = max(searchKm, P.topology_named_embayment_extra_km);
end
end

%% ============================================================
function [Pout, ok] = build_named_embayment_search_poly(outerPoly, innerPoly, basinId, lon0, lat0, P)

Pout = safe_polybuffer(outerPoly, P.topology_bay_extra_km);
ok = false;
Specs = named_embayment_specs(basinId, P);
if isempty(Specs)
    return;
end

maxGapKm = topology_named_embayment_search_km(P, basinId);
padKm = 120;
if isfield(P, 'topology_named_embayment_window_pad_km')
    padKm = max(0, P.topology_named_embayment_window_pad_km);
end
triggerKm = 180;
if isfield(P, 'topology_named_embayment_trigger_km')
    triggerKm = max(0, P.topology_named_embayment_trigger_km);
end
innerTrigger = safe_polybuffer(innerPoly, triggerKm);

baseSearch = Pout;
for ispec = 1:numel(Specs)
    specPoly = named_embayment_spec_poly(Specs(ispec), lon0, lat0, P);
    triggerPoly = named_embayment_trigger_poly(Specs(ispec), lon0, lat0);
    triggerMode = named_embayment_trigger_mode(Specs(ispec));
    if area(specPoly) <= 0
        continue;
    end
    nearNamedBay = true;
    if triggerMode == "inner_only" || triggerMode == "inner_or_current"
        if area(innerTrigger) <= 0 || area(triggerPoly) <= 0
            nearNamedBay = false;
        else
            try
                nearNamedBay = area(intersect(innerTrigger, safe_polybuffer(triggerPoly, padKm))) > 0;
            catch
                nearNamedBay = false;
            end
        end
    end
    if ~nearNamedBay
        continue;
    end
    gapKm = polyshape_gap_km(baseSearch, specPoly);
    if isfinite(gapKm) && gapKm <= maxGapKm
        Pout = union(Pout, safe_polybuffer(specPoly, padKm));
        ok = true;
    end
end
if ok
    Pout = rmholes(union(Pout));
end
end

%% ============================================================
function Pout = named_embayment_spec_poly(S, lon0, lat0, P)

if nargin < 4
    P = struct();
end

if isfield(S, 'poly_lon') && isfield(S, 'poly_lat') && ...
        numel(S.poly_lon) >= 4 && numel(S.poly_lat) == numel(S.poly_lon)
    lon = S.poly_lon(:);
    lat = S.poly_lat(:);
else
    useSuperellipse = true;
    if isfield(P, 'topology_named_embayment_window_use_superellipse')
        useSuperellipse = logical(P.topology_named_embayment_window_use_superellipse);
    end
    if useSuperellipse
        [lon, lat] = named_embayment_superellipse_lonlat(S, P);
    else
        lon = [S.lon_min; S.lon_max; S.lon_max; S.lon_min; S.lon_min];
        lat = [S.lat_min; S.lat_min; S.lat_max; S.lat_max; S.lat_min];
    end
end
[x, y] = lonlat_to_local_km(lon, lat, lon0, lat0);
try
    Pout = polyshape(x, y, 'Simplify', true);
catch
    Pout = polyshape();
end
end

%% ============================================================
function [lon, lat] = named_embayment_superellipse_lonlat(S, P)

n = 144;
theta = linspace(0, 2 * pi, n + 1).';
theta(end) = [];
cx = S.lon_center;
cy = S.lat_center;
rx = max(0.05, 0.5 * abs(S.lon_max - S.lon_min));
ry = max(0.05, 0.5 * abs(S.lat_max - S.lat_min));

scale = 1.0;
expo = 2.0;
if nargin >= 2 && isfield(P, 'topology_named_embayment_window_scale')
    scale = P.topology_named_embayment_window_scale;
end
if nargin >= 2 && isfield(P, 'topology_named_embayment_window_exponent')
    expo = P.topology_named_embayment_window_exponent;
end
scale = max(0.85, min(1.25, scale));
expo = max(2.0, min(6.0, expo));

c = cos(theta);
s = sin(theta);
lon = cx + scale * rx * sign(c) .* abs(c).^(2 / expo);
lat = cy + scale * ry * sign(s) .* abs(s).^(2 / expo);
lon = lon(:);
lat = lat(:);
end

%% ============================================================
function Pout = named_embayment_trigger_poly(S, lon0, lat0)

if isfield(S, 'trigger_poly_lon') && isfield(S, 'trigger_poly_lat') && ...
        numel(S.trigger_poly_lon) >= 4 && numel(S.trigger_poly_lat) == numel(S.trigger_poly_lon)
    lon = S.trigger_poly_lon(:);
    lat = S.trigger_poly_lat(:);
elseif isfield(S, 'trigger_lon_min') && isfield(S, 'trigger_lon_max') && ...
        isfield(S, 'trigger_lat_min') && isfield(S, 'trigger_lat_max')
    lon = [S.trigger_lon_min; S.trigger_lon_max; S.trigger_lon_max; S.trigger_lon_min; S.trigger_lon_min];
    lat = [S.trigger_lat_min; S.trigger_lat_min; S.trigger_lat_max; S.trigger_lat_max; S.trigger_lat_min];
else
    Pout = named_embayment_spec_poly(S, lon0, lat0);
    return;
end

[x, y] = lonlat_to_local_km(lon, lat, lon0, lat0);
try
    Pout = polyshape(x, y, 'Simplify', true);
catch
    Pout = polyshape();
end
end

%% ============================================================
function mode = named_embayment_trigger_mode(S)

mode = "inner_or_current";
if isfield(S, 'trigger_mode') && strlength(string(S.trigger_mode)) > 0
    mode = string(S.trigger_mode);
end
mode = lower(strtrim(mode));
end

%% ============================================================
function mask = named_embayment_window_mask(S, X, Y, ctx, P)

mask = false(size(X));
try
    specPoly = named_embayment_spec_poly(S, ctx.lon0, ctx.lat0, P);
    if area(specPoly) <= 0
        return;
    end
    padKm = 0;
    if isfield(P, 'topology_named_embayment_window_pad_km')
        padKm = max(0, 0.25 * P.topology_named_embayment_window_pad_km);
    end
    if padKm > 0
        specPoly = safe_polybuffer(specPoly, padKm);
    end
    mask = polyshape_grid_mask_local(specPoly, X, Y) & ctx.domainMask;
catch
    [lonG, latG] = local_km_to_lonlat(X, Y, ctx.lon0, ctx.lat0);
    lonG = reshape(wrapTo180_local(lonG), size(X));
    latG = reshape(latG, size(X));
    mask = lonG >= S.lon_min & lonG <= S.lon_max & ...
        latG >= S.lat_min & latG <= S.lat_max & ctx.domainMask;
end
end

%% ============================================================
function [Pout, ok, stats] = apply_gebco_topology_outer_adjustment(outerPoly, studyPoly, innerPoly, lon0, lat0, outerBufferKm, basinId, P)

Pout = outerPoly;
ok = false;
stats = empty_topology_stats();
if ~isfield(P, 'use_gebco_topology') || ~P.use_gebco_topology || ...
        exist(P.gebco_path, 'file') ~= 2 || area(outerPoly) <= 0
    return;
end

searchPoly = safe_polybuffer(outerPoly, max(P.topology_bay_extra_km, P.outer_backshore_limit_km));
ctx = build_gebco_topology_context(searchPoly, studyPoly, innerPoly, lon0, lat0, outerBufferKm, basinId, P);
if ~ctx.ok
    return;
end

[clipped, capsuleOk] = build_capsule_style_topology_domain(ctx, outerPoly, innerPoly, P);
if ~capsuleOk
    backCells = max(1, round(P.outer_backshore_limit_km / ctx.dx));
    otherOceanMask = ctx.oceanMask & ~ctx.targetOceanMask;
    targetDist = grid_distance_steps_from_sources(ctx.targetOceanMask, backCells);
    otherDist = grid_distance_steps_from_sources(otherOceanMask, backCells);
    otherMargin = 1;
    if isfield(P, 'topology_land_other_ocean_margin_cells')
        otherMargin = max(0, round(P.topology_land_other_ocean_margin_cells));
    end
    coastalLandMask = ctx.majorLandMask | (ctx.landMask & ctx.innerMask);
    landNearTarget = coastalLandMask & targetDist <= backCells & ...
        (targetDist + otherMargin <= otherDist | ~isfinite(otherDist));
    landNearTarget = landNearTarget | (ctx.landMask & ctx.innerMask);
    innerTopoMask = ctx.innerMask & (ctx.landMask | ctx.targetOceanMask | dilate_mask_local(ctx.targetOceanMask, 1));
    topoDomainMask = ctx.targetOceanMask | landNearTarget | innerTopoMask;
    topoDomainMask = topoDomainMask & ctx.domainMask;

    topoDomainPoly = grid_mask_to_polyshape_local(topoDomainMask, ctx.xv, ctx.yv);
    if area(topoDomainPoly) <= 0
        return;
    end

    clipped = intersect(outerPoly, topoDomainPoly);
    clipped = union(clipped, innerPoly);
    clipped = rmholes(union(clipped));
    if area(clipped) <= max(area(innerPoly) * 1.08, area(studyPoly) + 1)
        return;
    end
end

shelfDomain = polyshape();
[shelfAddMask, shelfInfo] = select_topologic_shelf_additions(ctx, clipped, outerPoly, P);
[shelfGapMask, shelfGapInfo] = select_topologic_shelf_corner_gap_additions(ctx, clipped, outerPoly, shelfAddMask, P);
shelfAddMask = shelfAddMask | shelfGapMask;
stats.shelf_add_component_count = shelfInfo.accepted_component_count + shelfGapInfo.accepted_component_count;
if any(shelfAddMask(:))
    shelfPoly = grid_mask_to_polyshape_local(shelfAddMask, ctx.xv, ctx.yv);
    if area(shelfPoly) > 0
        shelfSmooth = max(0, P.topology_shelf_smooth_km);
        if shelfSmooth > 0
            shelfSmoothed = safe_polybuffer(shelfPoly, shelfSmooth);
            shelfSmoothed = safe_polybuffer(shelfSmoothed, -0.78 * shelfSmooth);
            if area(shelfSmoothed) > 0
                shelfPoly = rmholes(union(shelfPoly, shelfSmoothed));
            end
        end
        shelfDomain = safe_polybuffer(shelfPoly, 0.35 * ctx.dx);
        shelfDomain = intersect(shelfDomain, searchPoly);
        stats.shelf_domain = shelfDomain;
        stats.shelf_add_area_km2 = area(shelfDomain);
        stats.shelf_add_cell_count = nnz(shelfAddMask);
        clipped = union(clipped, shelfDomain);
        clipped = rmholes(union(clipped));
    end
end

[landBayAddMask, landBayInfo] = select_topologic_land_closed_bay_additions(ctx, clipped, outerPoly, P);
clippedForRingBay = clipped;
if any(landBayAddMask(:))
    landBayDomain0 = bay_mask_to_smoothed_domain(landBayAddMask, searchPoly, ctx, P);
    if area(landBayDomain0) > 0
        clippedForRingBay = rmholes(union(clippedForRingBay, landBayDomain0));
    end
end

[majorBayAddMask, majorBayInfo] = select_topologic_major_embayment_additions(ctx, clippedForRingBay, outerPoly, P);
if any(majorBayAddMask(:))
    majorBayDomain0 = bay_mask_to_smoothed_domain(majorBayAddMask, searchPoly, ctx, P);
    if area(majorBayDomain0) > 0
        clippedForRingBay = rmholes(union(clippedForRingBay, majorBayDomain0));
    end
else
    majorBayInfo = struct('core_cell_count', 0, 'accepted_component_count', 0);
end

namedCtx = ctx;
namedSearchPoly = searchPoly;
namedSearchKm = topology_named_embayment_search_km(P, basinId);
[namedSearchPoly0, namedSearchOk] = build_named_embayment_search_poly(outerPoly, innerPoly, basinId, lon0, lat0, P);
if namedSearchOk && namedSearchKm > P.topology_bay_extra_km
    namedCtx0 = build_gebco_topology_context(namedSearchPoly0, studyPoly, innerPoly, ...
        lon0, lat0, outerBufferKm, basinId, P);
    if namedCtx0.ok
        namedCtx = namedCtx0;
        namedSearchPoly = namedSearchPoly0;
    end
end

namedBayDomain = polyshape();
namedBayCellCount = 0;
[namedBayAddMask, namedBayInfo] = select_named_embayment_additions(namedCtx, clippedForRingBay, innerPoly, outerPoly, basinId, P);
if any(namedBayAddMask(:))
    namedBayDomain0 = bay_mask_to_smoothed_domain(namedBayAddMask, namedSearchPoly, namedCtx, P);
    if area(namedBayDomain0) > 0
        namedBayDomain = namedBayDomain0;
        namedBayCellCount = nnz(namedBayAddMask);
        clippedForRingBay = rmholes(union(clippedForRingBay, namedBayDomain0));
    end
else
    namedBayInfo = struct('core_cell_count', 0, 'accepted_component_count', 0);
end

ringBayAddMask = false(size(landBayAddMask));
bayInfo = struct('core_cell_count', 0, 'accepted_component_count', 0);
if ~isfield(P, 'topology_bay_ring_enable') || P.topology_bay_ring_enable
    [ringBayAddMask, bayInfo] = select_topologic_bay_additions(ctx, clippedForRingBay, outerPoly, P);
end
bayAddMask = landBayAddMask | majorBayAddMask | ringBayAddMask;
stats.bay_add_core_cell_count = bayInfo.core_cell_count + landBayInfo.core_cell_count + ...
    majorBayInfo.core_cell_count + namedBayInfo.core_cell_count;
stats.bay_add_component_count = bayInfo.accepted_component_count + landBayInfo.accepted_component_count + ...
    majorBayInfo.accepted_component_count + namedBayInfo.accepted_component_count;
bayDomain = namedBayDomain;
bayCellCount = namedBayCellCount;
if any(bayAddMask(:))
    bayDomain0 = bay_mask_to_smoothed_domain(bayAddMask, searchPoly, ctx, P);
    if area(bayDomain0) > 0
        bayDomain = rmholes(union(bayDomain, bayDomain0));
        bayCellCount = bayCellCount + nnz(bayAddMask);
    end
end
if area(bayDomain) > 0
    bayDomain = smooth_bay_domain_final_outline(bayDomain, searchPoly, ctx, P);
    stats.bay_domain = bayDomain;
    stats.bay_add_area_km2 = area(bayDomain);
    stats.bay_add_cell_count = bayCellCount;
    clipped = union(clipped, bayDomain);
    clipped = rmholes(union(clipped));
end

clipped = smooth_topology_polyshape(clipped, P);
if area(clipped) <= 0
    return;
end

[postBayAddMask, postBayInfo] = select_topologic_land_closed_bay_additions(ctx, clipped, outerPoly, P);
if any(postBayAddMask(:))
    postBayDomain = bay_mask_to_smoothed_domain(postBayAddMask, searchPoly, ctx, P);
    if area(postBayDomain) > 0
        stats.bay_domain = rmholes(union(stats.bay_domain, postBayDomain));
        stats.bay_domain = smooth_bay_domain_final_outline(stats.bay_domain, searchPoly, ctx, P);
        stats.bay_add_area_km2 = area(stats.bay_domain);
        stats.bay_add_cell_count = stats.bay_add_cell_count + nnz(postBayAddMask);
        stats.bay_add_core_cell_count = stats.bay_add_core_cell_count + postBayInfo.core_cell_count;
        stats.bay_add_component_count = stats.bay_add_component_count + postBayInfo.accepted_component_count;
        clipped = rmholes(union(clipped, stats.bay_domain));
        clipped = smooth_topology_polyshape(clipped, P);
    end
end

clipped = enforce_final_single_ocean_component(clipped, innerPoly, studyPoly, ctx, P);
clipped = smooth_final_open_ocean_boundary(clipped, innerPoly, studyPoly, ctx, P);
clipped = enforce_final_single_ocean_component(clipped, innerPoly, studyPoly, ctx, P);
clipped = enforce_strict_final_single_ocean_component(clipped, innerPoly, studyPoly, ctx, P);
clipped = keep_largest_relevant_region(clipped, innerPoly, studyPoly, P);
clipped = enforce_strict_final_single_ocean_component(clipped, innerPoly, studyPoly, ctx, P);
if area(namedBayDomain) > 0
    clipped = rmholes(union(clipped, namedBayDomain));
end
Pout = clipped;
ok = true;
end

%% ============================================================
function [Pout, ok] = build_capsule_style_topology_domain(ctx, outerPoly, innerPoly, P)

Pout = polyshape();
ok = false;
if ~isfield(P, 'topology_capsule_style') || ~P.topology_capsule_style || area(outerPoly) <= 0
    return;
end

[X, Y] = meshgrid(ctx.xv, ctx.yv);
baseMask = polyshape_grid_mask_local(outerPoly, X, Y) & ctx.domainMask;
targetBase = ctx.targetOceanMask & baseMask;
if nnz(targetBase) < 20
    return;
end

otherOcean = ctx.oceanMask & baseMask & ~ctx.targetOceanMask;
otherOcean = areaopen_mask_local(otherOcean, max(20, topology_min_component_cells(P)));

bayCloseCells = max(1, round(P.topology_capsule_bay_close_km / ctx.dx));
oceanDomain = close_mask_local(targetBase, bayCloseCells);
oceanDomain = oceanDomain & baseMask & ~otherOcean;
oceanDomain = flood_fill_from_seed_mask(oceanDomain, targetBase);
if nnz(oceanDomain) < 20
    return;
end

inlandCells = max(1, round(P.topology_capsule_inland_km / ctx.dx));
landReachCells = inlandCells + max(1, round(35 / ctx.dx));
maxDistSteps = ceil(hypot(size(baseMask, 1), size(baseMask, 2))) + 5;
targetDist = grid_distance_steps_from_sources(targetBase, maxDistSteps);
if any(otherOcean(:))
    otherDist = grid_distance_steps_from_sources(otherOcean, maxDistSteps);
else
    otherDist = inf(size(baseMask));
end
otherMargin = max(0, round(P.topology_capsule_other_ocean_margin_cells));

nearTargetOcean = dilate_mask_local(oceanDomain, landReachCells);
landDomain = ctx.landMask & ctx.domainMask & nearTargetOcean & ...
    targetDist <= inlandCells & ...
    (targetDist + otherMargin < otherDist | ~isfinite(otherDist));

protectedMask = ctx.innerMask & ctx.domainMask;
protectedTargetMask = protectedMask & ...
    (ctx.landMask | ctx.targetOceanMask | dilate_mask_local(ctx.targetOceanMask, 1));
blockedOcean = ctx.oceanMask & ~ctx.targetOceanMask;
finalMask = (oceanDomain | landDomain | protectedTargetMask) & ctx.domainMask;
finalMask = fill_holes_mask_local(finalMask);
finalMask = (finalMask & ~blockedOcean) | protectedTargetMask;

razorCells = max(0, round(P.topology_capsule_razor_km / ctx.dx));
if razorCells > 0
    opened = open_mask_local(finalMask, razorCells);
    opened = opened | protectedTargetMask;
    if nnz(opened) > 0.60 * nnz(finalMask)
        finalMask = opened;
    end
end

minCells = max(10, round(P.topology_capsule_min_component_area_km2 / max(ctx.dx * ctx.dx, eps)));
finalMask = areaopen_mask_local(finalMask, minCells);
finalMask = finalMask | protectedTargetMask;
finalMask = fill_holes_mask_local(finalMask);
finalMask = (finalMask & ~blockedOcean) | protectedTargetMask;

fillCells = max(0, round(P.topology_capsule_fill_km / ctx.dx));
if fillCells > 0
    finalMask = close_mask_local(finalMask, fillCells);
end
finalMask = ((finalMask & ~blockedOcean) | protectedTargetMask) & ctx.domainMask;
finalMask = flood_fill_from_seed_mask(finalMask, targetBase | protectedTargetMask);
finalMask = fill_holes_mask_local(finalMask);
finalMask = (finalMask & ~blockedOcean) | protectedTargetMask;

if nnz(finalMask) < 20
    return;
end

Pmask = grid_mask_to_polyshape_local(finalMask, ctx.xv, ctx.yv);
if area(Pmask) <= 0
    return;
end
Pmask = union(Pmask, innerPoly);
Pmask = rmholes(union(Pmask));

smoothR = max(0, 0.70 * P.topology_boundary_smooth_km);
if smoothR > 0
    try
        sm = safe_polybuffer(Pmask, smoothR);
        sm = safe_polybuffer(sm, -0.74 * smoothR);
        if area(sm) > 0
            cover = area(intersect(sm, Pmask)) / max(area(Pmask), eps);
            if cover >= 0.97
                Pmask = rmholes(union(sm, innerPoly));
            else
                Pmask = rmholes(union(Pmask, sm));
            end
        end
    catch
    end
end

Pout = Pmask;
ok = area(Pout) > 0;
end

%% ============================================================
function Pout = smooth_topology_polyshape(Pin, P)

Pout = Pin;
if ~isfield(P, 'topology_boundary_smooth_km') || P.topology_boundary_smooth_km <= 0 || area(Pin) <= 0
    return;
end
r = P.topology_boundary_smooth_km;
try
    insetFraction = 0.28;
    if isfield(P, 'topology_boundary_inset_fraction')
        insetFraction = max(0.05, min(0.85, P.topology_boundary_inset_fraction));
    end
    sm = safe_polybuffer(Pin, r);
    sm = safe_polybuffer(sm, -insetFraction * r);
    if area(sm) > 0.80 * area(Pin)
        try
            cover = area(intersect(sm, Pin)) / max(area(Pin), eps);
        catch
            cover = 0;
        end
        if cover >= 0.985
            % Use the rounded topology envelope directly when it still
            % covers the computed footprint; this removes GEBCO-grid and
            % seed-envelope straight edges from the open-sea boundary.
            Pout = rmholes(sm);
        else
            Pout = rmholes(union(sm, Pin));
        end
    end
catch
end
end

%% ============================================================
function Pout = smooth_final_output_outer_boundary(Pin, innerPoly, studyPoly, P)

Pout = Pin;
if area(Pin) <= 0 || ~isfield(P, 'final_output_outer_smooth_enable') || ...
        ~P.final_output_outer_smooth_enable
    return;
end

oldArea = area(Pin);
protectKm = 130;
if isfield(P, 'final_output_outer_protect_km')
    protectKm = max(0, P.final_output_outer_protect_km);
end
protectedPoly = rmholes(union(innerPoly, studyPoly));
if protectKm > 0
    protectedPoly = safe_polybuffer(protectedPoly, protectKm);
end

candidate = Pin;
roundKm = 0;
if isfield(P, 'final_output_outer_round_km')
    roundKm = max(0, P.final_output_outer_round_km);
end
if roundKm > 0
    insetFraction = 0.90;
    if isfield(P, 'final_output_outer_round_inset_fraction')
        insetFraction = max(0.55, min(1.05, P.final_output_outer_round_inset_fraction));
    end
    try
        rounded = safe_polybuffer(candidate, roundKm);
        rounded = safe_polybuffer(rounded, -insetFraction * roundKm);
        rounded = rmholes(union(rounded, protectedPoly));
        rounded = union(rounded, innerPoly);
        if output_outer_candidate_ok(rounded, Pin, innerPoly, oldArea, P)
            candidate = rounded;
        end
    catch
    end
end

curveP = P;
curveP.final_open_boundary_curve_smooth_enable = true;
curveP.final_open_boundary_curve_step_km = 35;
curveP.final_open_boundary_curve_window_km = 720;
curveP.final_open_boundary_curve_cover_pad_km = 40;
curveP.final_open_boundary_curve_cover_min = 0.88;
curveP.final_open_boundary_curve_max_area_increase_fraction = 0.35;
if isfield(P, 'final_output_outer_curve_step_km')
    curveP.final_open_boundary_curve_step_km = max(12, P.final_output_outer_curve_step_km);
end
if isfield(P, 'final_output_outer_curve_window_km')
    curveP.final_open_boundary_curve_window_km = max(curveP.final_open_boundary_curve_step_km, ...
        P.final_output_outer_curve_window_km);
end
if isfield(P, 'final_output_outer_curve_cover_pad_km')
    curveP.final_open_boundary_curve_cover_pad_km = max(0, P.final_output_outer_curve_cover_pad_km);
end
if isfield(P, 'final_output_outer_curve_cover_min')
    curveP.final_open_boundary_curve_cover_min = max(0.75, min(0.98, ...
        P.final_output_outer_curve_cover_min));
end
if isfield(P, 'final_output_outer_max_area_increase_fraction')
    curveP.final_open_boundary_curve_max_area_increase_fraction = max(0, min(0.70, ...
        P.final_output_outer_max_area_increase_fraction));
end

try
    curved = smooth_outer_boundary_curve(candidate, innerPoly, protectedPoly, curveP);
    curved = rmholes(union(curved, protectedPoly));
    curved = union(curved, innerPoly);
    if output_outer_candidate_ok(curved, Pin, innerPoly, oldArea, P)
        candidate = curved;
    end
catch
end

if output_outer_candidate_ok(candidate, Pin, innerPoly, oldArea, P)
    Pout = rmholes(union(candidate));
end
end

%% ============================================================
function ok = output_outer_candidate_ok(candidate, original, innerPoly, oldArea, P)

ok = false;
try
    newArea = area(candidate);
    innerCover = area(intersect(candidate, innerPoly)) / max(area(innerPoly), eps);
    coverFrac = area(intersect(candidate, original)) / max(oldArea, eps);
catch
    return;
end
if newArea <= 0 || innerCover < 0.995
    return;
end
minFrac = 0.80;
if isfield(P, 'final_output_outer_min_area_fraction')
    minFrac = max(0.50, min(0.98, P.final_output_outer_min_area_fraction));
end
maxIncrease = 0.35;
if isfield(P, 'final_output_outer_max_area_increase_fraction')
    maxIncrease = max(0, min(0.80, P.final_output_outer_max_area_increase_fraction));
end
if newArea < minFrac * oldArea || newArea > (1 + maxIncrease) * oldArea
    return;
end
if coverFrac < minFrac
    return;
end
ok = true;
end

%% ============================================================
function Pout = smooth_final_open_ocean_boundary(Pin, innerPoly, studyPoly, ctx, P)

Pout = Pin;
if ~isfield(P, 'final_smooth_open_ocean_boundary') || ...
        ~P.final_smooth_open_ocean_boundary || area(Pin) <= 0 || ...
        ~isfield(ctx, 'ok') || ~ctx.ok || isempty(ctx.oceanMask)
    return;
end

try
    [X, Y] = meshgrid(ctx.xv, ctx.yv);
    outerMask = polyshape_grid_mask_local(Pin, X, Y) & ctx.domainMask;
catch
    return;
end
if nnz(outerMask) < 20
    return;
end

protectKm = 180;
if isfield(P, 'final_open_boundary_protect_km')
    protectKm = max(0, P.final_open_boundary_protect_km);
end
protectCells = max(0, round(protectKm / max(ctx.dx, eps)));
innerStudyMask = outerMask & (ctx.innerMask | polyshape_grid_mask_local(studyPoly, X, Y));
coastalSource = ctx.majorLandMask | ctx.innerMask | innerStudyMask;
if protectCells > 0
    protectedMask = outerMask & dilate_mask_local(coastalSource, protectCells);
else
    protectedMask = innerStudyMask;
end
protectedMask = protectedMask | innerStudyMask;

fillOnly = true;
if isfield(P, 'final_open_boundary_fill_only')
    fillOnly = logical(P.final_open_boundary_fill_only);
end

closeKm = 0;
if isfield(P, 'final_open_boundary_close_km')
    closeKm = max(0, P.final_open_boundary_close_km);
end
closeCells = max(0, round(closeKm / max(ctx.dx, eps)));
if closeCells <= 0
    return;
end

closedMask = close_mask_local(outerMask, closeCells) & ctx.domainMask;
if ~any(closedMask(:))
    return;
end

if fillOnly
    fillMask = closedMask & ~outerMask & ~protectedMask;
else
    fillMask = closedMask & ~outerMask;
end

maxExtraKm = 180;
if isfield(P, 'final_open_boundary_max_extra_km')
    maxExtraKm = max(0, P.final_open_boundary_max_extra_km);
end
if maxExtraKm > 0
    extraMask = polyshape_grid_mask_local(safe_polybuffer(Pin, maxExtraKm), X, Y);
    fillMask = fillMask & extraMask;
end
if nnz(fillMask) < 4
    return;
end

fillPoly = grid_mask_to_polyshape_local(fillMask, ctx.xv, ctx.yv);
if area(fillPoly) <= 0
    return;
end

roundKm = 0;
if isfield(P, 'final_open_boundary_round_km')
    roundKm = max(0, P.final_open_boundary_round_km);
end
if roundKm > 0
    insetFraction = 0.82;
    if isfield(P, 'final_open_boundary_round_inset_fraction')
        insetFraction = max(0.55, min(1.05, P.final_open_boundary_round_inset_fraction));
    end
    try
        rounded = safe_polybuffer(fillPoly, roundKm);
        rounded = safe_polybuffer(rounded, -insetFraction * roundKm);
        if area(rounded) > 0
            fillPoly = rmholes(union(fillPoly, rounded));
        end
    catch
    end
end

protectedPoly = grid_mask_to_polyshape_local(protectedMask, ctx.xv, ctx.yv);
candidate = rmholes(union(Pin, fillPoly));
candidate = rmholes(union(candidate, protectedPoly));
candidate = union(candidate, innerPoly);
if maxExtraKm > 0
    candidate = intersect(candidate, safe_polybuffer(Pin, maxExtraKm));
    candidate = union(candidate, innerPoly);
end

edgeRoundKm = 0;
if isfield(P, 'final_open_boundary_edge_round_km')
    edgeRoundKm = max(0, P.final_open_boundary_edge_round_km);
end
if edgeRoundKm > 0
    edgeInsetFraction = 1.00;
    if isfield(P, 'final_open_boundary_edge_round_inset_fraction')
        edgeInsetFraction = max(0.65, min(1.00, P.final_open_boundary_edge_round_inset_fraction));
    end
    edgeCoverMin = 0.985;
    if isfield(P, 'final_open_boundary_edge_round_cover_min')
        edgeCoverMin = max(0.90, min(0.999, P.final_open_boundary_edge_round_cover_min));
    end
    try
        roundedCandidate = safe_polybuffer(candidate, edgeRoundKm);
        roundedCandidate = safe_polybuffer(roundedCandidate, -edgeInsetFraction * edgeRoundKm);
        if area(roundedCandidate) > 0
            coverFrac = area(intersect(roundedCandidate, candidate)) / max(area(candidate), eps);
            if coverFrac >= edgeCoverMin
                candidate = roundedCandidate;
            else
                candidate = union(candidate, roundedCandidate);
            end
            candidate = rmholes(union(candidate, protectedPoly));
            candidate = union(candidate, innerPoly);
        end
    catch
    end
end

if isfield(P, 'final_open_boundary_curve_smooth_enable') && ...
        P.final_open_boundary_curve_smooth_enable
    candidate = smooth_outer_boundary_curve(candidate, innerPoly, protectedPoly, P);
end

if isfield(P, 'final_open_boundary_envelope_enable') && ...
        P.final_open_boundary_envelope_enable
    candidate = smooth_outer_boundary_with_envelope(candidate, innerPoly, protectedPoly, P);
end

try
    oldArea = area(Pin);
    newArea = area(candidate);
catch
    return;
end
if newArea <= 0
    return;
end

minFrac = 0.72;
if isfield(P, 'final_open_boundary_min_area_fraction')
    minFrac = max(0.35, min(0.98, P.final_open_boundary_min_area_fraction));
end
maxIncrease = 0.12;
if isfield(P, 'final_open_boundary_max_area_increase_fraction')
    maxIncrease = max(0, min(0.80, P.final_open_boundary_max_area_increase_fraction));
end
innerCover = area(intersect(candidate, innerPoly)) / max(area(innerPoly), eps);
if innerCover < 0.98 || newArea < minFrac * oldArea || newArea > (1 + maxIncrease) * oldArea
    return;
end

Pout = rmholes(union(candidate));

doReport = false;
if isfield(P, 'final_open_boundary_report')
    doReport = logical(P.final_open_boundary_report);
end
deltaArea = area(Pout) - oldArea;
if doReport && abs(deltaArea) > max(1000, 0.001 * oldArea)
    fprintf('  -> final open-boundary smoothing near %.2fE %.2fN: area %+.0f km2\n', ...
        ctx.lon0, ctx.lat0, deltaArea);
end
end

%% ============================================================
function Pout = smooth_outer_boundary_curve(Pin, innerPoly, protectedPoly, P)

Pout = Pin;
if area(Pin) <= 0
    return;
end

stepKm = 45;
if isfield(P, 'final_open_boundary_curve_step_km')
    stepKm = max(12, P.final_open_boundary_curve_step_km);
end
windowKm = 360;
if isfield(P, 'final_open_boundary_curve_window_km')
    windowKm = max(stepKm, P.final_open_boundary_curve_window_km);
end
coverPadKm = 65;
if isfield(P, 'final_open_boundary_curve_cover_pad_km')
    coverPadKm = max(0, P.final_open_boundary_curve_cover_pad_km);
end
coverMin = 0.982;
if isfield(P, 'final_open_boundary_curve_cover_min')
    coverMin = max(0.90, min(0.999, P.final_open_boundary_curve_cover_min));
end
maxIncrease = 0.16;
if isfield(P, 'final_open_boundary_curve_max_area_increase_fraction')
    maxIncrease = max(0, min(0.45, P.final_open_boundary_curve_max_area_increase_fraction));
end

try
    [xb, yb] = boundary(Pin);
    segments = nan_split_segments(xb, yb);
catch
    return;
end
if isempty(segments)
    return;
end

smoothPoly = polyshape();
for iseg = 1:numel(segments)
    seg = segments{iseg};
    if size(seg, 1) < 4
        continue;
    end
    ring = resample_closed_ring_local(seg, stepKm);
    if size(ring, 1) < 12
        continue;
    end
    ring = circular_smooth_ring_local(ring, max(3, round(windowKm / stepKm)));
    try
        part = polyshape(ring(:, 1), ring(:, 2), 'Simplify', true);
        if area(part) > 0
            smoothPoly = union(smoothPoly, part);
        end
    catch
    end
end
if area(smoothPoly) <= 0
    return;
end

try
    if coverPadKm > 0
        smoothPoly = safe_polybuffer(smoothPoly, coverPadKm);
        smoothPoly = safe_polybuffer(smoothPoly, -0.25 * coverPadKm);
    end
    smoothPoly = rmholes(union(smoothPoly));
    smoothPoly = rmholes(union(smoothPoly, protectedPoly));
    smoothPoly = union(smoothPoly, innerPoly);
catch
    return;
end

try
    coverFrac = area(intersect(smoothPoly, Pin)) / max(area(Pin), eps);
    smArea = area(smoothPoly);
    oldArea = area(Pin);
catch
    return;
end

if coverFrac < coverMin || smArea <= 0 || smArea > (1 + maxIncrease) * oldArea
    return;
end

Pout = rmholes(smoothPoly);
end

%% ============================================================
function ring = resample_closed_ring_local(seg, stepKm)

seg = double(seg(:, 1:2));
seg = seg(all(isfinite(seg), 2), :);
if size(seg, 1) < 4
    ring = seg;
    return;
end
if norm(seg(1, :) - seg(end, :)) == 0
    seg = seg(1:end-1, :);
end
if size(seg, 1) < 3
    ring = [seg; seg(1, :)];
    return;
end

x = seg(:, 1);
y = seg(:, 2);
xClosed = [x; x(1)];
yClosed = [y; y(1)];
ds = hypot(diff(xClosed), diff(yClosed));
s = [0; cumsum(ds)];
perim = s(end);
if perim <= 0
    ring = [seg; seg(1, :)];
    return;
end
n = max(18, ceil(perim / max(stepKm, eps)));
tq = linspace(0, perim, n + 1).';
tq(end) = [];
xq = interp1(s, xClosed, tq, 'linear');
yq = interp1(s, yClosed, tq, 'linear');
ring = [xq, yq; xq(1), yq(1)];
end

%% ============================================================
function ring = circular_smooth_ring_local(ring, windowN)

ring = double(ring(:, 1:2));
if norm(ring(1, :) - ring(end, :)) == 0
    pts = ring(1:end-1, :);
else
    pts = ring;
end
n = size(pts, 1);
if n < 8
    ring = [pts; pts(1, :)];
    return;
end
windowN = max(3, min(n - 1, round(windowN)));
if mod(windowN, 2) == 0
    windowN = windowN + 1;
end
halfN = floor(windowN / 2);
idx = mod((-halfN:n+halfN-1), n) + 1;
pad = pts(idx, :);
kernel = ones(windowN, 1) / windowN;
xs = conv(pad(:, 1), kernel, 'same');
ys = conv(pad(:, 2), kernel, 'same');
pts2 = [xs(halfN + 1:halfN + n), ys(halfN + 1:halfN + n)];
ring = [pts2; pts2(1, :)];
end

%% ============================================================
function Pout = smooth_outer_boundary_with_envelope(Pin, innerPoly, protectedPoly, P)

Pout = Pin;
if area(Pin) <= 0
    return;
end

try
    [xb, yb] = boundary(Pin);
    valid = isfinite(xb) & isfinite(yb);
    xb = xb(valid);
    yb = yb(valid);
catch
    return;
end
if numel(xb) < 4
    return;
end

shrink = 0.30;
if isfield(P, 'final_open_boundary_envelope_shrink')
    shrink = max(0.02, min(0.90, P.final_open_boundary_envelope_shrink));
end
padKm = 0;
if isfield(P, 'final_open_boundary_envelope_pad_km')
    padKm = max(0, P.final_open_boundary_envelope_pad_km);
end
roundKm = 0;
if isfield(P, 'final_open_boundary_envelope_round_km')
    roundKm = max(0, P.final_open_boundary_envelope_round_km);
end
roundInset = 0.96;
if isfield(P, 'final_open_boundary_envelope_round_inset_fraction')
    roundInset = max(0.70, min(1.00, P.final_open_boundary_envelope_round_inset_fraction));
end
coverMin = 0.98;
if isfield(P, 'final_open_boundary_envelope_cover_min')
    coverMin = max(0.90, min(0.999, P.final_open_boundary_envelope_cover_min));
end
maxIncrease = 0.18;
if isfield(P, 'final_open_boundary_envelope_max_area_increase_fraction')
    maxIncrease = max(0, min(0.60, P.final_open_boundary_envelope_max_area_increase_fraction));
end

try
    k = boundary(xb(:), yb(:), shrink);
    env = polyshape(xb(k), yb(k), 'Simplify', true);
    if area(env) <= 0
        return;
    end
    if padKm > 0
        env = safe_polybuffer(env, padKm);
    end
    if roundKm > 0
        env = safe_polybuffer(env, roundKm);
        env = safe_polybuffer(env, -roundInset * roundKm);
    end
    env = rmholes(union(env));
    env = rmholes(union(env, protectedPoly));
    env = union(env, innerPoly);
catch
    return;
end

try
    coverFrac = area(intersect(env, Pin)) / max(area(Pin), eps);
    envArea = area(env);
    oldArea = area(Pin);
catch
    return;
end

if coverFrac < coverMin || envArea <= 0 || envArea > (1 + maxIncrease) * oldArea
    return;
end

Pout = rmholes(env);
end

%% ============================================================
function Pout = enforce_final_single_ocean_component(Pin, innerPoly, studyPoly, ctx, P)

Pout = Pin;
if ~isfield(P, 'final_enforce_single_ocean_component') || ...
        ~P.final_enforce_single_ocean_component || area(Pin) <= 0 || ...
        ~isfield(ctx, 'ok') || ~ctx.ok || isempty(ctx.oceanMask)
    return;
end

try
    [X, Y] = meshgrid(ctx.xv, ctx.yv);
    outerMask = polyshape_grid_mask_local(Pin, X, Y) & ctx.domainMask;
catch
    return;
end
if nnz(outerMask) < 20
    return;
end

oceanInOuter = ctx.oceanMask & outerMask;
components = connected_components_pixelidx(oceanInOuter);
if numel(components) <= 1
    return;
end

minAreaKm2 = 5000;
if isfield(P, 'final_ocean_min_component_area_km2')
    minAreaKm2 = max(0, P.final_ocean_min_component_area_km2);
end
minCells = max(1, round(minAreaKm2 / max(ctx.dx * ctx.dx, eps)));

seedBufferKm = 160;
if isfield(P, 'final_ocean_seed_buffer_km')
    seedBufferKm = max(0, P.final_ocean_seed_buffer_km);
end
seedPoly = safe_polybuffer(innerPoly, seedBufferKm);
seedMask = oceanInOuter & polyshape_grid_mask_local(seedPoly, X, Y);
if ~any(seedMask(:))
    seedPoly = safe_polybuffer(studyPoly, max(seedBufferKm, 0.5 * seedBufferKm + 50));
    seedMask = oceanInOuter & polyshape_grid_mask_local(seedPoly, X, Y);
end
targetMask = oceanInOuter & ctx.targetOceanMask;
studyMask = polyshape_grid_mask_local(studyPoly, X, Y);
innerStudyMask = outerMask & (ctx.innerMask | studyMask);
innerTouchMask = oceanInOuter & dilate_mask_local(innerStudyMask, final_inner_touch_cells(P, ctx));

compCells = zeros(numel(components), 1);
seedCells = zeros(numel(components), 1);
targetCells = zeros(numel(components), 1);
innerCells = zeros(numel(components), 1);
for i = 1:numel(components)
    pix = components{i};
    compCells(i) = numel(pix);
    seedCells(i) = nnz(seedMask(pix));
    targetCells(i) = nnz(targetMask(pix));
    innerCells(i) = nnz(innerTouchMask(pix));
end

[best, significant, rejectCount, keepMode] = select_final_ocean_component( ...
    compCells, seedCells, targetCells, innerCells, minCells, P);
if best <= 0
    return;
end

keptOcean = false(size(oceanInOuter));
keptOcean(components{best}) = true;
rejectOcean = false(size(oceanInOuter));
for i = find(significant(:).')
    if i == best
        continue;
    end
    rejectOcean(components{i}) = true;
end
if ~any(rejectOcean(:))
    return;
end

maxSteps = ceil(hypot(size(outerMask, 1), size(outerMask, 2))) + 5;
distKeep = grid_distance_steps_from_sources(keptOcean, maxSteps);
distReject = grid_distance_steps_from_sources(rejectOcean, maxSteps);

landKeepMarginKm = 80;
if isfield(P, 'final_ocean_land_keep_margin_km')
    landKeepMarginKm = max(0, P.final_ocean_land_keep_margin_km);
end
landKeepCells = max(0, round(landKeepMarginKm / max(ctx.dx, eps)));
if landKeepCells > 0
    nearKeptOcean = dilate_mask_local(keptOcean, landKeepCells);
else
    nearKeptOcean = keptOcean;
end

nearKeptSide = distKeep <= distReject | ~isfinite(distReject);
protectedNearKept = innerStudyMask & (nearKeptSide | nearKeptOcean);

keepMask = outerMask & (nearKeptSide | nearKeptOcean | protectedNearKept);
keepMask = (keepMask & ~rejectOcean) | keptOcean | protectedNearKept;
if nnz(keepMask) < 20
    return;
end

keepPoly = grid_mask_to_polyshape_local(keepMask, ctx.xv, ctx.yv);
if area(keepPoly) <= 0
    return;
end

maskBufferKm = 0;
if isfield(P, 'final_ocean_keep_mask_buffer_km')
    maskBufferKm = max(0, P.final_ocean_keep_mask_buffer_km);
end
if maskBufferKm > 0
    keepPoly = safe_polybuffer(keepPoly, maskBufferKm);
end

try
    trimmed = intersect(Pin, keepPoly);
    if area(trimmed) <= 0
        return;
    end
    trimmed = rmholes(union(trimmed));
catch
    return;
end

oldArea = area(Pin);
newArea = area(trimmed);
if newArea <= 0 || newArea < 0.35 * max(area(innerPoly), area(studyPoly))
    return;
end
Pout = trimmed;

doReport = false;
if isfield(P, 'final_ocean_connectivity_report')
    doReport = logical(P.final_ocean_connectivity_report);
end
removedArea = max(0, oldArea - newArea);
if doReport && removedArea > max(1000, 0.002 * oldArea)
    fprintf(['  -> final ocean connectivity trim near %.2fE %.2fN: ', ...
        'kept %s ocean component, removed %d secondary components, outer -%.0f km2\n'], ...
        ctx.lon0, ctx.lat0, char(keepMode), rejectCount, removedArea);
end
end

%% ============================================================
function Pout = enforce_strict_final_single_ocean_component(Pin, innerPoly, studyPoly, ctx, P)

Pout = Pin;
if ~isfield(P, 'final_strict_single_ocean_component') || ...
        ~P.final_strict_single_ocean_component || area(Pin) <= 0 || ...
        ~isfield(ctx, 'ok') || ~ctx.ok || isempty(ctx.oceanMask)
    return;
end

try
    [X, Y] = meshgrid(ctx.xv, ctx.yv);
    outerMask = polyshape_grid_mask_local(Pin, X, Y) & ctx.domainMask;
catch
    return;
end
if nnz(outerMask) < 20
    return;
end

oceanInOuter = ctx.oceanMask & outerMask;
components = connected_components_pixelidx(oceanInOuter);
if numel(components) <= 1
    return;
end

minAreaKm2 = 0;
if isfield(P, 'final_strict_ocean_min_component_area_km2')
    minAreaKm2 = max(0, P.final_strict_ocean_min_component_area_km2);
end
minCells = max(1, round(minAreaKm2 / max(ctx.dx * ctx.dx, eps)));

seedBufferKm = 160;
if isfield(P, 'final_ocean_seed_buffer_km')
    seedBufferKm = max(0, P.final_ocean_seed_buffer_km);
end
seedPoly = safe_polybuffer(innerPoly, seedBufferKm);
seedMask = oceanInOuter & polyshape_grid_mask_local(seedPoly, X, Y);
if ~any(seedMask(:))
    seedPoly = safe_polybuffer(studyPoly, max(seedBufferKm, 0.5 * seedBufferKm + 50));
    seedMask = oceanInOuter & polyshape_grid_mask_local(seedPoly, X, Y);
end
studyMask = polyshape_grid_mask_local(studyPoly, X, Y);
innerStudyMask = outerMask & (ctx.innerMask | studyMask);
targetMask = oceanInOuter & ctx.targetOceanMask;
innerTouchMask = oceanInOuter & dilate_mask_local(innerStudyMask, final_inner_touch_cells(P, ctx));

compCells = zeros(numel(components), 1);
seedCells = zeros(numel(components), 1);
targetCells = zeros(numel(components), 1);
innerCells = zeros(numel(components), 1);
for i = 1:numel(components)
    pix = components{i};
    compCells(i) = numel(pix);
    seedCells(i) = nnz(seedMask(pix));
    targetCells(i) = nnz(targetMask(pix));
    innerCells(i) = nnz(innerTouchMask(pix));
end

[best, eligible, rejectCount, keepMode] = select_final_ocean_component( ...
    compCells, seedCells, targetCells, innerCells, minCells, P);
if best <= 0
    return;
end

keptOcean = false(size(oceanInOuter));
keptOcean(components{best}) = true;
rejectOcean = false(size(oceanInOuter));
for i = find(eligible(:).')
    if i == best
        continue;
    end
    rejectOcean(components{i}) = true;
end
if rejectCount == 0
    return;
end
allowedOcean = oceanInOuter & ~rejectOcean;

landMarginKm = 80;
if isfield(P, 'final_strict_ocean_land_margin_km')
    landMarginKm = max(0, P.final_strict_ocean_land_margin_km);
end
landCells = max(0, round(landMarginKm / max(ctx.dx, eps)));
if landCells > 0
    landAroundKept = dilate_mask_local(keptOcean, landCells) & outerMask & ~ctx.oceanMask;
else
    landAroundKept = outerMask & ~ctx.oceanMask & innerStudyMask;
end
keepMask = allowedOcean | landAroundKept;

keepPoly = grid_mask_to_polyshape_local(keepMask, ctx.xv, ctx.yv);
if area(keepPoly) <= 0
    return;
end

smoothKm = 0;
if isfield(P, 'final_strict_ocean_boundary_smooth_km')
    smoothKm = max(0, P.final_strict_ocean_boundary_smooth_km);
end
if smoothKm > 0
    try
        sm = safe_polybuffer(keepPoly, smoothKm);
        sm = safe_polybuffer(sm, -0.45 * smoothKm);
        if area(sm) > 0
            keepPoly = rmholes(union(keepPoly, sm));
        end
    catch
    end
end

try
    trimmed = intersect(Pin, keepPoly);
    if area(trimmed) <= 0
        return;
    end
    trimmed = rmholes(union(trimmed));
catch
    return;
end

try
    trimMask = polyshape_grid_mask_local(trimmed, X, Y) & ctx.domainMask;
    trimOcean = trimMask & ctx.oceanMask;
    finalComponents = connected_components_pixelidx(trimOcean);
    finalBadOcean = false(size(trimOcean));
    finalRejectCount = 0;
    for i = 1:numel(finalComponents)
        idx = finalComponents{i};
        if numel(idx) < minCells || any(keptOcean(idx))
            continue;
        end
        finalBadOcean(idx) = true;
        finalRejectCount = finalRejectCount + 1;
    end
    if finalRejectCount > 0
        finalKeepMask = trimMask & ~finalBadOcean;
        finalPoly = grid_mask_to_polyshape_local(finalKeepMask, ctx.xv, ctx.yv);
        if area(finalPoly) > 0
            trimmed = rmholes(union(intersect(trimmed, finalPoly)));
        end
    end
catch
end

oldArea = area(Pin);
newArea = area(trimmed);
if newArea <= 0 || newArea < 0.20 * max(area(innerPoly), area(studyPoly))
    return;
end

Pout = trimmed;

doReport = false;
if isfield(P, 'final_strict_ocean_report')
    doReport = logical(P.final_strict_ocean_report);
end
removedArea = max(0, oldArea - newArea);
if doReport && (rejectCount > 0 || removedArea > max(1000, 0.001 * oldArea))
    fprintf(['  -> strict final ocean connectivity near %.2fE %.2fN: ', ...
        'kept %s water body, removed %d disconnected water bodies, outer -%.0f km2\n'], ...
        ctx.lon0, ctx.lat0, char(keepMode), rejectCount, removedArea);
end
end

%% ============================================================
function n = final_inner_touch_cells(P, ctx)

bufferKm = 0;
if isfield(P, 'final_ocean_inner_touch_buffer_km')
    bufferKm = max(0, P.final_ocean_inner_touch_buffer_km);
elseif isfield(P, 'final_inner_touch_buffer_km')
    bufferKm = max(0, P.final_inner_touch_buffer_km);
end
n = max(1, round(bufferKm / max(ctx.dx, eps)));
end

%% ============================================================
function [best, eligible, rejectCount, modeName] = select_final_ocean_component( ...
    compCells, seedCells, targetCells, innerCells, minCells, P)

best = 0;
rejectCount = 0;
modeName = "selected";
eligible = compCells >= minCells;
if nnz(eligible) <= 1
    return;
end

preferInner = true;
if isfield(P, 'final_ocean_prefer_inner_touching')
    preferInner = logical(P.final_ocean_prefer_inner_touching);
end

scores = -inf(size(compCells));
innerRelevant = eligible & (seedCells > 0 | innerCells > 0);
if preferInner && any(innerRelevant)
    candidates = innerRelevant;
    modeName = "inner-touching";
    for i = find(candidates(:).')
        scores(i) = 1e7 * double(seedCells(i) > 0) + ...
            7e6 * double(innerCells(i) > 0) + ...
            5000 * double(seedCells(i)) + 2500 * double(innerCells(i)) + ...
            0.01 * double(targetCells(i)) + 1e-3 * double(compCells(i));
    end
else
    candidates = eligible & (targetCells > 0 | seedCells > 0 | innerCells > 0);
    if ~any(candidates)
        candidates = eligible;
        modeName = "largest";
    else
        modeName = "target";
    end
    for i = find(candidates(:).')
        scores(i) = 1e6 * double(targetCells(i) > 0) + ...
            5e5 * double(seedCells(i) > 0) + ...
            3e5 * double(innerCells(i) > 0) + ...
            50 * log(1 + double(seedCells(i) + innerCells(i))) + ...
            log(1 + double(compCells(i)));
    end
end

[bestScore, best] = max(scores);
if ~isfinite(bestScore) || best <= 0 || ~eligible(best)
    best = 0;
    return;
end
rejectCount = nnz(eligible) - 1;
end

%% ============================================================
function ctx = build_gebco_topology_context(searchPoly, studyPoly, innerPoly, lon0, lat0, outerBufferKm, basinId, P)

ctx = struct('ok', false, 'xv', [], 'yv', [], 'dx', NaN, ...
    'lon0', NaN, 'lat0', NaN, ...
    'domainMask', [], 'landMask', [], 'majorLandMask', [], 'oceanMask', [], ...
    'straitBarrierMask', [], 'targetOceanMask', [], ...
    'innerMask', [], 'z', []);

[xb, yb] = boundary(searchPoly);
v = isfinite(xb) & isfinite(yb);
if nnz(v) < 3
    return;
end

pad = max(P.topology_grid_km, 0.04 * outerBufferKm);
xmin = min(xb(v)) - pad;
xmax = max(xb(v)) + pad;
ymin = min(yb(v)) - pad;
ymax = max(yb(v)) + pad;
if ~(isfinite(xmin) && isfinite(xmax) && isfinite(ymin) && isfinite(ymax) && xmax > xmin && ymax > ymin)
    return;
end

dx = P.topology_grid_km;
estCells = ((xmax - xmin) / dx + 1) * ((ymax - ymin) / dx + 1);
if estCells > P.topology_max_cells
    dx = sqrt((xmax - xmin) * (ymax - ymin) / P.topology_max_cells);
    dx = max(dx, P.topology_grid_km);
end
xv = xmin:dx:xmax;
yv = ymin:dx:ymax;
if numel(xv) < 5 || numel(yv) < 5
    return;
end

[X, Y] = meshgrid(xv, yv);
domainMask = polyshape_grid_mask_local(searchPoly, X, Y);
if ~any(domainMask(:))
    return;
end

[lonQ, latQ] = local_km_to_lonlat(X(:), Y(:), lon0, lat0);
z = sample_gebco_nearest(P.gebco_path, lonQ, latQ, dx);
if isempty(z)
    return;
end
z = reshape(z, size(X));
oceanMask = domainMask & isfinite(z) & z <= P.topology_ocean_zmax_m;
landMask = domainMask & isfinite(z) & z > P.topology_ocean_zmax_m;
if nnz(oceanMask) < 20
    return;
end
majorLandMask = select_major_land_mask(landMask, dx, P);
straitBarrierMask = build_strait_barrier_mask(majorLandMask, oceanMask, dx, P);
oceanForTarget = oceanMask & ~straitBarrierMask;
if nnz(oceanForTarget) < 20
    oceanForTarget = oceanMask;
    straitBarrierMask = false(size(oceanMask));
end

seedPoly = safe_polybuffer(studyPoly, P.topology_seed_buffer_km);
seedMask = oceanForTarget & polyshape_grid_mask_local(seedPoly, X, Y);
if ~any(seedMask(:))
    seedPoly = safe_polybuffer(innerPoly, P.topology_seed_buffer_km + P.inner_refinement_buffer_km);
    seedMask = oceanForTarget & polyshape_grid_mask_local(seedPoly, X, Y);
end
if ~any(seedMask(:))
    seedPoly = safe_polybuffer(studyPoly, max(P.topology_seed_buffer_km, 0.45 * outerBufferKm));
    seedMask = oceanForTarget & polyshape_grid_mask_local(seedPoly, X, Y);
end
if ~any(seedMask(:)) && any(straitBarrierMask(:))
    seedPoly = safe_polybuffer(studyPoly, max(P.topology_seed_buffer_km, 0.45 * outerBufferKm));
    seedMask = oceanMask & polyshape_grid_mask_local(seedPoly, X, Y);
    if any(seedMask(:))
        oceanForTarget = oceanMask;
        straitBarrierMask = false(size(oceanMask));
    end
end
if ~any(seedMask(:))
    return;
end

targetOceanMask = select_target_ocean_component_mask(oceanForTarget, seedMask, X, Y, ...
    studyPoly, lon0, lat0, basinId);
if ~any(targetOceanMask(:))
    targetOceanMask = flood_fill_from_seed_mask(oceanForTarget, seedMask);
end
if ~any(targetOceanMask(:)) && ~isequal(oceanForTarget, oceanMask)
    targetOceanMask = select_target_ocean_component_mask(oceanMask, seedMask, X, Y, ...
        studyPoly, lon0, lat0, basinId);
    if ~any(targetOceanMask(:))
        targetOceanMask = flood_fill_from_seed_mask(oceanMask, seedMask);
    end
    straitBarrierMask = false(size(oceanMask));
end
if nnz(targetOceanMask) < 20
    return;
end

ctx.ok = true;
ctx.xv = xv;
ctx.yv = yv;
ctx.dx = dx;
ctx.lon0 = lon0;
ctx.lat0 = lat0;
ctx.domainMask = domainMask;
ctx.landMask = landMask;
ctx.majorLandMask = majorLandMask;
ctx.oceanMask = oceanMask;
ctx.straitBarrierMask = straitBarrierMask;
ctx.targetOceanMask = targetOceanMask;
ctx.innerMask = polyshape_grid_mask_local(innerPoly, X, Y);
ctx.z = z;
end

%% ============================================================
function majorLandMask = select_major_land_mask(landMask, dxKm, P)

majorLandMask = false(size(landMask));
components = connected_components_pixelidx(landMask);
if isempty(components)
    return;
end

cellArea = dxKm * dxKm;
minArea = 25000;
if isfield(P, 'topology_major_land_min_area_km2')
    minArea = P.topology_major_land_min_area_km2;
end

areas = zeros(numel(components), 1);
for i = 1:numel(components)
    areas(i) = numel(components{i}) * cellArea;
end

keep = areas >= minArea;
if ~any(keep)
    [~, order] = sort(areas, 'descend');
    keep(order(1:min(3, numel(order)))) = true;
end

for i = find(keep(:).')
    majorLandMask(components{i}) = true;
end
end

%% ============================================================
function barrierMask = build_strait_barrier_mask(majorLandMask, oceanMask, dxKm, P)

barrierMask = false(size(oceanMask));
if ~isfield(P, 'topology_strait_barrier_km') || P.topology_strait_barrier_km <= 0 || ...
        ~any(majorLandMask(:)) || ~any(oceanMask(:))
    return;
end

components = connected_components_pixelidx(majorLandMask);
if isempty(components)
    return;
end

cellArea = max(dxKm * dxKm, eps);
areas = zeros(numel(components), 1);
for i = 1:numel(components)
    areas(i) = numel(components{i}) * cellArea;
end
[~, order] = sort(areas, 'descend');

maxComp = 25;
if isfield(P, 'topology_strait_barrier_max_land_components')
    maxComp = max(1, round(P.topology_strait_barrier_max_land_components));
end
order = order(1:min(maxComp, numel(order)));

radiusCells = max(1, round(P.topology_strait_barrier_km / dxKm));
dilatedSum = zeros(size(oceanMask), 'uint16');
for ii = order(:).'
    tmp = false(size(oceanMask));
    tmp(components{ii}) = true;
    dilatedSum = dilatedSum + uint16(dilate_mask_local(tmp, radiusCells));
end

barrierMask = dilatedSum >= 2 & oceanMask;
smoothCells = 0;
if isfield(P, 'topology_strait_barrier_smooth_km')
    smoothCells = max(0, round(P.topology_strait_barrier_smooth_km / dxKm));
end
if smoothCells > 0 && any(barrierMask(:))
    barrierMask = close_mask_local(barrierMask, smoothCells);
    barrierMask = open_mask_local(barrierMask, smoothCells);
end
barrierMask = barrierMask & oceanMask;
end

%% ============================================================
function frac = local_fraction_mask(mask, radiusCells)

radiusCells = max(1, round(radiusCells));
[xx, yy] = meshgrid(-radiusCells:radiusCells, -radiusCells:radiusCells);
kernel = double((xx.^2 + yy.^2) <= radiusCells^2);
den = conv2(ones(size(mask)), kernel, 'same');
num = conv2(double(mask), kernel, 'same');
frac = num ./ max(den, eps);
end

%% ============================================================
function [sectorCount, oppositePairCount] = land_enclosure_sector_metrics(landMask, radiusCells)

radiusCells = max(2, round(radiusCells));
[dc, dr] = meshgrid(-radiusCells:radiusCells, -radiusCells:radiusCells);
rr = sqrt(double(dc).^2 + double(dr).^2);
within = rr > 0 & rr <= radiusCells;
ang = atan2(-double(dr), double(dc));
sectorWidth = 2 * pi / 8;
sectorId = mod(floor((ang + pi) / sectorWidth), 8) + 1;

sectorHit = false([size(landMask), 8]);
land = double(landMask);
for is = 1:8
    kernel = double(within & sectorId == is);
    if any(kernel(:))
        sectorHit(:, :, is) = conv2(land, kernel, 'same') > 0;
    end
end

sectorCount = sum(sectorHit, 3);
oppositePairCount = double(sectorHit(:, :, 1) & sectorHit(:, :, 5)) + ...
    double(sectorHit(:, :, 2) & sectorHit(:, :, 6)) + ...
    double(sectorHit(:, :, 3) & sectorHit(:, :, 7)) + ...
    double(sectorHit(:, :, 4) & sectorHit(:, :, 8));
end

%% ============================================================
function edgeMask = mask_boundary_local(mask)

if ~any(mask(:))
    edgeMask = false(size(mask));
    return;
end
n = conv2(double(mask), ones(3), 'same');
edgeMask = mask & n < 9;
end

%% ============================================================
function targetOceanMask = select_target_ocean_component_mask(oceanMask, seedMask, X, Y, studyPoly, lon0, lat0, basinId)

targetOceanMask = false(size(oceanMask));
components = connected_components_pixelidx(oceanMask);
if isempty(components)
    return;
end

[refLon, refLat, refKnown] = basin_ocean_reference_lonlat(basinId, lon0, lat0);
if refKnown
    [xRef, yRef] = lonlat_to_local_km(refLon, refLat, lon0, lat0);
else
    xRef = NaN;
    yRef = NaN;
end

[xs, ys] = boundary(studyPoly);
validStudy = isfinite(xs) & isfinite(ys);
if any(validStudy)
    xStudy = mean(xs(validStudy), 'omitnan');
    yStudy = mean(ys(validStudy), 'omitnan');
else
    xStudy = mean(X(seedMask), 'omitnan');
    yStudy = mean(Y(seedMask), 'omitnan');
end

best = 0;
bestScore = -inf;
bestRefDist = inf;
bestSeedHits = 0;

for i = 1:numel(components)
    idx = components{i};
    seedHits = nnz(seedMask(idx));
    if seedHits == 0
        continue;
    end

    step = max(1, ceil(numel(idx) / 3000));
    sampleIdx = idx(1:step:end);
    if refKnown
        refDist = min(hypot(X(sampleIdx) - xRef, Y(sampleIdx) - yRef), [], 'omitnan');
    else
        refDist = 0;
    end

    cx = mean(X(sampleIdx), 'omitnan');
    cy = mean(Y(sampleIdx), 'omitnan');
    studyDist = hypot(cx - xStudy, cy - yStudy);
    % Keep the component that the refinement seed actually touches.  The
    % basin reference is only a tie-breaker; otherwise narrow land bridges
    % can flip the final red domain to the opposite sea.
    score = double(seedHits) - 1.0e-3 * studyDist + 1.0e-6 * numel(idx);
    if refKnown
        score = score - 1.0e-4 * refDist;
    end

    if score > bestScore || ...
            (abs(score - bestScore) < 1e-9 && refDist < bestRefDist) || ...
            (abs(score - bestScore) < 1e-9 && refDist == bestRefDist && seedHits > bestSeedHits)
        best = i;
        bestScore = score;
        bestRefDist = refDist;
        bestSeedHits = seedHits;
    end
end

if best > 0
    targetOceanMask(components{best}) = true;
end
end

%% ============================================================
function waterMask = current_connected_ocean_mask(ctx, currentMask, searchMask, touchCells)

ocean = ctx.oceanMask & searchMask;
seed = ocean & dilate_mask_local(currentMask, max(1, round(touchCells)));
if any(seed(:))
    waterMask = flood_fill_from_seed_mask(ocean, seed);
else
    waterMask = ctx.targetOceanMask & searchMask;
end
waterMask = waterMask & ocean;
end

%% ============================================================
function zq = sample_gebco_nearest(gebcoPath, lonQ, latQ, dxKm)

zq = [];
persistent lonVec latVec cachedPath cacheOk
try
    if isempty(cachedPath) || ~strcmp(cachedPath, gebcoPath)
        lonVec = double(ncread(gebcoPath, 'lon'));
        latVec = double(ncread(gebcoPath, 'lat'));
        cachedPath = gebcoPath;
        cacheOk = true;
    end
catch
    cacheOk = false;
end
if isempty(cacheOk) || ~cacheOk
    return;
end

lonQ = wrapTo180_local(double(lonQ(:)));
latQ = double(latQ(:));
valid = isfinite(lonQ) & isfinite(latQ);
if ~any(valid)
    return;
end

lonlim = [min(lonQ(valid)) - 1, max(lonQ(valid)) + 1];
latlim = [max(-90, min(latQ(valid)) - 1), min(90, max(latQ(valid)) + 1)];
if diff(lonlim) > 170
    return;
end
lonlim(1) = max(min(lonVec), lonlim(1));
lonlim(2) = min(max(lonVec), lonlim(2));

iAll = find(lonVec >= lonlim(1) & lonVec <= lonlim(2));
jAll = find(latVec >= latlim(1) & latVec <= latlim(2));
if numel(iAll) < 3 || numel(jAll) < 3
    return;
end

baseDlon = median(diff(lonVec), 'omitnan');
targetDeg = max(0.025, dxKm / 111.32 / 2.2);
stride = max(1, floor(targetDeg / max(baseDlon, eps)));
i0 = iAll(1);
j0 = jAll(1);
ni = floor((iAll(end) - i0) / stride) + 1;
nj = floor((jAll(end) - j0) / stride) + 1;
try
    z = double(ncread(gebcoPath, 'elevation', [i0, j0], [ni, nj], [stride, stride]));
catch
    return;
end
lonSub = lonVec(i0 + (0:ni-1) * stride);
latSub = latVec(j0 + (0:nj-1) * stride);
try
    zq = interp2(lonSub(:).', latSub(:), z.', lonQ, latQ, 'nearest', NaN);
catch
    zq = [];
end
end

%% ============================================================
function [addMask, info] = select_topologic_shelf_additions(ctx, clippedPoly, outerPoly, P)

info = struct();
info.accepted_component_count = 0;

[X, Y] = meshgrid(ctx.xv, ctx.yv);
currentMask = polyshape_grid_mask_local(clippedPoly, X, Y);
addMask = false(size(currentMask));

if ~isfield(P, 'topology_shelf_expand') || ~P.topology_shelf_expand
    return;
end

searchMask = polyshape_grid_mask_local(safe_polybuffer(outerPoly, P.topology_shelf_extra_km), X, Y);
searchEdgeMask = mask_boundary_local(searchMask);
shelfMask = ctx.targetOceanMask & searchMask & ~currentMask & isfinite(ctx.z) & ...
    ctx.z <= P.topology_ocean_zmax_m & ctx.z >= P.topology_shelf_isobath_m;
if ~any(shelfMask(:))
    return;
end

touchCells = max(1, round(P.topology_shelf_touch_km / ctx.dx));
seedMask = shelfMask & dilate_mask_local(currentMask, touchCells);
if ~any(seedMask(:))
    return;
end

reachCells = max(touchCells, round(P.topology_shelf_reach_km / ctx.dx));
distSeed = grid_distance_steps_from_sources(seedMask, reachCells);
growthMask = shelfMask & (distSeed <= reachCells);
grownMask = flood_fill_from_seed_mask(growthMask, seedMask);
razorCells = max(0, round(P.topology_shelf_razor_km / ctx.dx));
if razorCells > 0 && any(grownMask(:))
    openedMask = open_mask_local(grownMask, razorCells);
    fillCells = max(0, round(P.topology_shelf_fill_km / ctx.dx));
    if fillCells > 0
        openedMask = close_mask_local(openedMask, fillCells);
    end
    openedMask = flood_fill_from_seed_mask(openedMask, dilate_mask_local(currentMask, touchCells + razorCells));
    grownMask = openedMask & shelfMask;
end
components = connected_components_pixelidx(grownMask);
cellArea = ctx.dx * ctx.dx;
acceptedAreaKm2 = 0;
minCells = max(topology_min_component_cells(P), ceil((0.05 * P.topology_shelf_max_area_km2) / cellArea));

for i = 1:numel(components)
    idx = components{i};
    if numel(idx) < minCells
        continue;
    end
    areaKm2 = numel(idx) * cellArea;
    if areaKm2 > P.topology_shelf_max_area_km2
        continue;
    end
    if acceptedAreaKm2 + areaKm2 > P.topology_shelf_total_max_area_km2
        continue;
    end
    searchEdgeFrac = nnz(searchEdgeMask(idx)) / max(1, numel(idx));
    if searchEdgeFrac > P.topology_shelf_search_edge_touch_max_fraction
        continue;
    end
    [majorSpanKm, minorSpanKm] = component_principal_spans_km(idx, size(grownMask), ctx.dx);
    fillFraction = areaKm2 / max(majorSpanKm * minorSpanKm, ctx.dx * ctx.dx);
    if minorSpanKm < P.topology_shelf_min_minor_width_km || ...
            fillFraction < P.topology_shelf_min_fill_fraction
        continue;
    end
    addMask(idx) = true;
    acceptedAreaKm2 = acceptedAreaKm2 + areaKm2;
    info.accepted_component_count = info.accepted_component_count + 1;
end
end

%% ============================================================
function [addMask, info] = select_topologic_shelf_corner_gap_additions(ctx, clippedPoly, outerPoly, shelfAddMask, P)

info = struct();
info.accepted_component_count = 0;

[X, Y] = meshgrid(ctx.xv, ctx.yv);
currentMask = polyshape_grid_mask_local(clippedPoly, X, Y);
baseMask = (currentMask | shelfAddMask) & ctx.domainMask;
addMask = false(size(baseMask));

searchPoly = safe_polybuffer(clippedPoly, P.topology_shelf_gap_search_km);
searchPoly = intersect(searchPoly, safe_polybuffer(outerPoly, P.topology_bay_extra_km));
searchMask = polyshape_grid_mask_local(searchPoly, X, Y);
if ~any(searchMask(:))
    return;
end
searchEdgeMask = mask_boundary_local(searchMask);
touchCells = max(1, round(P.topology_shelf_touch_km / ctx.dx));
connectedOcean = current_connected_ocean_mask(ctx, baseMask, searchMask, touchCells);

gapCloseCells = max(1, round(P.topology_shelf_gap_fill_km / ctx.dx));
closedBase = close_mask_local(baseMask, gapCloseCells);
closedBase = fill_holes_mask_local(closedBase);
gapOcean = connectedOcean & searchMask & ~baseMask & isfinite(ctx.z) & ...
    ctx.z <= P.topology_ocean_zmax_m;
if ~isfield(P, 'topology_shelf_gap_shallow_only') || P.topology_shelf_gap_shallow_only
    gapOcean = gapOcean & ctx.z >= P.topology_shelf_isobath_m;
end
candidateMask = closedBase & gapOcean;
candidateMask = candidateMask & dilate_mask_local(baseMask, gapCloseCells + 1);
if ~any(candidateMask(:))
    return;
end

components = connected_components_pixelidx(candidateMask);
cellArea = ctx.dx * ctx.dx;
minCells = topology_min_component_cells(P);
acceptedAreaKm2 = 0;
touchBand = dilate_mask_local(baseMask, touchCells);

for i = 1:numel(components)
    idx = components{i};
    if numel(idx) < minCells
        continue;
    end
    areaKm2 = numel(idx) * cellArea;
    if areaKm2 > P.topology_shelf_gap_max_area_km2 || ...
            acceptedAreaKm2 + areaKm2 > P.topology_shelf_gap_total_max_area_km2
        continue;
    end
    if nnz(touchBand(idx)) < minCells
        continue;
    end
    edgeFrac = nnz(searchEdgeMask(idx)) / max(1, numel(idx));
    if edgeFrac > P.topology_shelf_gap_edge_touch_max_fraction
        continue;
    end
    [majorSpanKm, minorSpanKm] = component_principal_spans_km(idx, size(candidateMask), ctx.dx);
    fillFraction = areaKm2 / max(majorSpanKm * minorSpanKm, ctx.dx * ctx.dx);
    if fillFraction < 0.10 || minorSpanKm < 0.75 * P.topology_shelf_min_minor_width_km
        continue;
    end
    addMask(idx) = true;
    acceptedAreaKm2 = acceptedAreaKm2 + areaKm2;
    info.accepted_component_count = info.accepted_component_count + 1;
end
end

%% ============================================================
function [addMask, info] = select_topologic_land_closed_bay_additions(ctx, clippedPoly, outerPoly, P)

info = struct();
info.core_cell_count = 0;
info.accepted_component_count = 0;

[X, Y] = meshgrid(ctx.xv, ctx.yv);
currentMask = polyshape_grid_mask_local(clippedPoly, X, Y) & ctx.domainMask;
addMask = false(size(currentMask));

searchKm = max(P.topology_land_bay_search_km, P.topology_bay_ring_recover_km);
searchPoly = safe_polybuffer(clippedPoly, searchKm);
searchPoly = intersect(searchPoly, safe_polybuffer(outerPoly, P.topology_bay_extra_km));
searchMask = polyshape_grid_mask_local(searchPoly, X, Y);
if ~any(searchMask(:))
    return;
end

majorLand = ctx.majorLandMask & searchMask;
if ~any(majorLand(:))
    return;
end

bridgeCells = max(1, round(P.topology_land_bay_bridge_km / ctx.dx));
touchCells = max(1, round(P.topology_land_bay_touch_km / ctx.dx));

% This is the bay test requested by the user: if the rounded/current domain
% boundary and major mainland nearly close a loop, fill the enclosed water.
closureBase = (currentMask | majorLand) & searchMask;
closedLoop = close_mask_local(closureBase, bridgeCells);
filledClosedLoop = fill_holes_mask_local(closedLoop);
enclosedWater = filledClosedLoop & ~closedLoop;
bridgeWater = closedLoop & ~closureBase;

nearDomain = dilate_mask_local(currentMask, bridgeCells + touchCells);
nearLand = dilate_mask_local(majorLand, bridgeCells + touchCells);
domainDist = grid_distance_steps_from_sources(currentMask, bridgeCells);
landDist = grid_distance_steps_from_sources(majorLand, bridgeCells);
candidateMask = (enclosedWater | bridgeWater) & ctx.oceanMask & searchMask & ~currentMask & ...
    isfinite(ctx.z) & ctx.z <= P.topology_ocean_zmax_m & nearDomain & nearLand & ...
    domainDist <= bridgeCells & landDist <= bridgeCells;
enclosureCells = max(2, round(P.topology_land_bay_enclosure_km / ctx.dx));
[sectorCount, oppositePairCount] = land_enclosure_sector_metrics(majorLand, enclosureCells);
candidateMask = candidateMask & ...
    (sectorCount >= P.topology_land_bay_min_land_sectors | ...
    oppositePairCount >= P.topology_land_bay_min_opposite_pairs);
if ~any(candidateMask(:))
    return;
end

searchEdgeMask = mask_boundary_local(searchMask);
domainTouch = dilate_mask_local(currentMask, touchCells);
landTouch = dilate_mask_local(majorLand, touchCells);
components = connected_components_pixelidx(candidateMask);

cellArea = ctx.dx * ctx.dx;
minCells = max(topology_min_component_cells(P), ceil(P.topology_land_bay_min_area_km2 / cellArea));
acceptedAreaKm2 = 0;

for i = 1:numel(components)
    idx = components{i};
    if numel(idx) < minCells
        continue;
    end
    areaKm2 = numel(idx) * cellArea;
    if areaKm2 > P.topology_land_bay_max_area_km2 || ...
            acceptedAreaKm2 + areaKm2 > P.topology_land_bay_total_max_area_km2
        continue;
    end
    if nnz(domainTouch(idx)) < minCells || nnz(landTouch(idx)) < minCells
        continue;
    end
    touchFrac = min(nnz(domainTouch(idx)), nnz(landTouch(idx))) / max(1, numel(idx));
    if touchFrac < P.topology_land_bay_min_touch_fraction
        continue;
    end
    edgeFrac = nnz(searchEdgeMask(idx)) / max(1, numel(idx));
    if edgeFrac > P.topology_land_bay_edge_touch_max_fraction
        continue;
    end
    [majorSpanKm, minorSpanKm] = component_principal_spans_km(idx, size(candidateMask), ctx.dx);
    if majorSpanKm > P.topology_land_bay_max_major_span_km
        continue;
    end
    fillFraction = areaKm2 / max(majorSpanKm * minorSpanKm, ctx.dx * ctx.dx);
    if minorSpanKm < P.topology_land_bay_min_minor_width_km || ...
            fillFraction < P.topology_land_bay_min_fill_fraction
        continue;
    end
    addMask(idx) = true;
    acceptedAreaKm2 = acceptedAreaKm2 + areaKm2;
    info.core_cell_count = info.core_cell_count + numel(idx);
    info.accepted_component_count = info.accepted_component_count + 1;
end
end

%% ============================================================
function [addMask, info] = select_topologic_major_embayment_additions(ctx, clippedPoly, outerPoly, P)

info = struct();
info.core_cell_count = 0;
info.accepted_component_count = 0;

[X, Y] = meshgrid(ctx.xv, ctx.yv);
currentMask = polyshape_grid_mask_local(clippedPoly, X, Y) & ctx.domainMask;
addMask = false(size(currentMask));

if isfield(P, 'topology_major_bay_enable') && ~P.topology_major_bay_enable
    return;
end

searchKm = max(P.topology_major_bay_search_km, P.topology_land_bay_search_km);
searchPoly = safe_polybuffer(clippedPoly, searchKm);
searchPoly = intersect(searchPoly, safe_polybuffer(outerPoly, P.topology_bay_extra_km));
searchMask = polyshape_grid_mask_local(searchPoly, X, Y);
if ~any(searchMask(:))
    return;
end

majorLand = ctx.majorLandMask & searchMask;
if ~any(majorLand(:))
    return;
end

closeCells = max(2, round(P.topology_major_bay_close_km / ctx.dx));
recoverCells = max(closeCells, round(P.topology_major_bay_recover_km / ctx.dx));
touchCells = max(1, round(P.topology_major_bay_touch_km / ctx.dx));
bayOcean = current_connected_ocean_mask(ctx, currentMask, searchMask, touchCells);

ringLand = dilate_mask_local(majorLand, closeCells);
bayCoreWater = bayOcean & searchMask & ~currentMask & ~ringLand & ...
    isfinite(ctx.z) & ctx.z <= P.topology_ocean_zmax_m;
if ~any(bayCoreWater(:))
    return;
end

searchEdgeMask = mask_boundary_local(searchMask);
touchBand = dilate_mask_local(currentMask, touchCells);
enclosureCells = max(2, round(P.topology_land_bay_enclosure_km / ctx.dx));
[sectorCount, oppositePairCount] = land_enclosure_sector_metrics(majorLand, enclosureCells);

components = connected_components_pixelidx(bayCoreWater);
cellArea = ctx.dx * ctx.dx;
minCells = max(topology_min_component_cells(P), ceil(P.topology_major_bay_min_area_km2 / cellArea));
acceptedAreaKm2 = 0;

for i = 1:numel(components)
    idx = components{i};
    if numel(idx) < minCells
        continue;
    end
    coreAreaKm2 = numel(idx) * cellArea;
    if coreAreaKm2 > P.topology_major_bay_max_area_km2 || ...
            acceptedAreaKm2 + coreAreaKm2 > P.topology_major_bay_total_max_area_km2
        continue;
    end
    edgeFrac = nnz(searchEdgeMask(idx)) / max(1, numel(idx));
    if edgeFrac > P.topology_major_bay_edge_touch_max_fraction
        continue;
    end
    landSectorFrac = nnz(sectorCount(idx) >= P.topology_land_bay_min_land_sectors | ...
        oppositePairCount(idx) >= P.topology_land_bay_min_opposite_pairs) / max(1, numel(idx));
    if landSectorFrac < P.topology_major_bay_min_land_sector_fraction
        continue;
    end
    touchFrac = nnz(touchBand(idx)) / max(1, numel(idx));
    if touchFrac < P.topology_major_bay_min_touch_fraction
        continue;
    end
    [touchSectorCount, touchOppositePairs] = component_touch_sector_metrics(idx, touchBand, size(bayCoreWater));
    if touchSectorCount < P.topology_major_bay_min_touch_sectors && touchOppositePairs < 1
        continue;
    end

    coreMask = false(size(currentMask));
    coreMask(idx) = true;
    distCore = grid_distance_steps_from_sources(coreMask, recoverCells);
    bayMask = bayOcean & searchMask & ~currentMask & ...
        isfinite(ctx.z) & ctx.z <= P.topology_ocean_zmax_m & distCore <= recoverCells;
    bayMask = flood_fill_from_seed_mask(bayMask | coreMask, coreMask) & ...
        bayOcean & searchMask & ~currentMask & ...
        isfinite(ctx.z) & ctx.z <= P.topology_ocean_zmax_m & distCore <= recoverCells;
    if ~any(bayMask(:))
        continue;
    end
    bayAreaKm2 = nnz(bayMask) * cellArea;
    if bayAreaKm2 > P.topology_major_bay_max_area_km2 || ...
            acceptedAreaKm2 + bayAreaKm2 > P.topology_major_bay_total_max_area_km2
        continue;
    end
    bayEdgeFrac = nnz(searchEdgeMask & bayMask) / max(1, nnz(bayMask));
    if bayEdgeFrac > P.topology_major_bay_edge_touch_max_fraction
        continue;
    end

    addMask = addMask | bayMask;
    acceptedAreaKm2 = acceptedAreaKm2 + bayAreaKm2;
    info.core_cell_count = info.core_cell_count + numel(idx);
    info.accepted_component_count = info.accepted_component_count + 1;
end
end

%% ============================================================
function [addMask, info] = select_named_embayment_additions(ctx, clippedPoly, innerPoly, outerPoly, basinId, P)

info = struct();
info.core_cell_count = 0;
info.accepted_component_count = 0;

[X, Y] = meshgrid(ctx.xv, ctx.yv);
currentMask = polyshape_grid_mask_local(clippedPoly, X, Y) & ctx.domainMask;
innerMask = polyshape_grid_mask_local(innerPoly, X, Y) & ctx.domainMask;
addMask = false(size(currentMask));

if ~isfield(P, 'topology_named_embayment_enable') || ~P.topology_named_embayment_enable
    return;
end

searchKm = topology_named_embayment_search_km(P, basinId);
searchPoly = safe_polybuffer(outerPoly, searchKm);
searchMask = polyshape_grid_mask_local(searchPoly, X, Y);
if ~any(searchMask(:))
    return;
end

touchCells = max(1, round(P.topology_named_embayment_touch_km / ctx.dx));
touchBand = dilate_mask_local(currentMask, touchCells);
minTouchCells = max(3, ceil(0.15 * topology_min_component_cells(P)));
connectedOcean = current_connected_ocean_mask(ctx, currentMask, searchMask, touchCells);
triggerKm = 180;
if isfield(P, 'topology_named_embayment_trigger_km')
    triggerKm = max(0, P.topology_named_embayment_trigger_km);
end
triggerCells = max(1, round(triggerKm / ctx.dx));
innerTriggerBand = dilate_mask_local(innerMask, triggerCells);

Specs = named_embayment_specs(basinId, P);
if isempty(Specs)
    return;
end

for ispec = 1:numel(Specs)
    S = Specs(ispec);
    win = named_embayment_window_mask(S, X, Y, ctx, P);
    triggerMode = named_embayment_trigger_mode(S);
    triggerPoly = named_embayment_trigger_poly(S, ctx.lon0, ctx.lat0);
    triggerMask = false(size(currentMask));
    if area(triggerPoly) > 0
        triggerMask = polyshape_grid_mask_local(triggerPoly, X, Y) & ctx.domainMask;
    end
    namedWater = connectedOcean & searchMask & win & isfinite(ctx.z) & ctx.z <= P.topology_ocean_zmax_m;
    if nnz(namedWater) < topology_min_component_cells(P)
        continue;
    end
    if triggerMode == "inner_only" || triggerMode == "inner_or_current"
        gateMask = innerTriggerBand;
        if any(triggerMask(:))
            gateMask = gateMask & dilate_mask_local(triggerMask, max(1, round(0.75 * triggerCells)));
        end
        if nnz(gateMask & namedWater) < max(3, minTouchCells)
            continue;
        end
    end
    if nnz(touchBand & namedWater) < minTouchCells
        continue;
    end

    components = connected_components_pixelidx(namedWater);
    if isempty(components)
        continue;
    end

    [xCtr, yCtr] = lonlat_to_local_km(S.lon_center, S.lat_center, ctx.lon0, ctx.lat0);
    best = 0;
    bestScore = -inf;
    centerSeed = namedWater & hypot(X - xCtr, Y - yCtr) <= max(80, 3 * ctx.dx);
    centerComponent = flood_fill_from_seed_mask(namedWater, centerSeed);
    if any(centerComponent(:)) && nnz(centerComponent & touchBand) >= minTouchCells
        for i = 1:numel(components)
            if any(centerComponent(components{i}))
                best = i;
                break;
            end
        end
    end
    for i = 1:numel(components)
        if best > 0
            break;
        end
        idx = components{i};
        if nnz(touchBand(idx)) < minTouchCells
            continue;
        end
        cx = mean(X(idx), 'omitnan');
        cy = mean(Y(idx), 'omitnan');
        distCtr = hypot(cx - xCtr, cy - yCtr);
        score = 0.001 * numel(idx) + 1.5 * nnz(touchBand(idx)) - 0.25 * distCtr;
        if score > bestScore
            best = i;
            bestScore = score;
        end
    end
    if best <= 0
        continue;
    end

    oneMask = false(size(addMask));
    oneMask(components{best}) = true;
    oneMask = oneMask & ~currentMask;
    smoothCells = max(1, round(P.topology_named_embayment_smooth_km / ctx.dx));
    if smoothCells > 0 && any(oneMask(:))
        oneMask = close_mask_local(oneMask, smoothCells);
        oneMask = oneMask & namedWater & ~currentMask;
    end
    if any(oneMask(:))
        addMask = addMask | oneMask;
        info.core_cell_count = info.core_cell_count + nnz(oneMask);
        info.accepted_component_count = info.accepted_component_count + 1;
    end
end

end

%% ============================================================
function Specs = named_embayment_specs(basinId, P)

Specs = struct('name', {}, 'lon_min', {}, 'lon_max', {}, 'lat_min', {}, 'lat_max', {}, ...
    'lon_center', {}, 'lat_center', {}, ...
    'trigger_lon_min', {}, 'trigger_lon_max', {}, 'trigger_lat_min', {}, 'trigger_lat_max', {}, ...
    'trigger_mode', {});
basinId = string(basinId);

if basinId == "BASIN_NATL"
    Specs(end + 1) = struct('name', "Gulf of Mexico", ...
        'lon_min', -98.8, 'lon_max', -78.0, 'lat_min', 17.0, 'lat_max', 32.5, ...
        'lon_center', -89.2, 'lat_center', 25.0, ...
        'trigger_lon_min', -97.8, 'trigger_lon_max', -81.0, ...
        'trigger_lat_min', 21.2, 'trigger_lat_max', 31.8, ...
        'trigger_mode', "inner_only");
    if ~isfield(P, 'topology_named_caribbean_embayments_enable') || ...
            P.topology_named_caribbean_embayments_enable
        Specs(end + 1) = struct('name', "Southwestern Caribbean", ...
            'lon_min', -84.8, 'lon_max', -75.0, 'lat_min', 8.0, 'lat_max', 13.3, ...
            'lon_center', -80.6, 'lat_center', 10.4, ...
            'trigger_lon_min', -84.8, 'trigger_lon_max', -75.0, ...
            'trigger_lat_min', 8.0, 'trigger_lat_max', 13.3, ...
            'trigger_mode', "current_touch");
        Specs(end + 1) = struct('name', "Gulf of Honduras", ...
            'lon_min', -89.5, 'lon_max', -83.0, 'lat_min', 15.0, 'lat_max', 19.2, ...
            'lon_center', -86.6, 'lat_center', 16.9, ...
            'trigger_lon_min', -89.5, 'trigger_lon_max', -83.0, ...
            'trigger_lat_min', 15.0, 'trigger_lat_max', 19.2, ...
            'trigger_mode', "current_touch");
    end
end

if basinId == "BASIN_WNP" && isfield(P, 'topology_named_wpac_embayments_enable') && ...
        P.topology_named_wpac_embayments_enable
    Specs(end + 1) = struct('name', "Gulf of Thailand", ...
        'lon_min', 99.0, 'lon_max', 105.8, 'lat_min', 6.5, 'lat_max', 14.4, ...
        'lon_center', 101.7, 'lat_center', 10.7, ...
        'trigger_lon_min', 99.0, 'trigger_lon_max', 105.8, ...
        'trigger_lat_min', 6.5, 'trigger_lat_max', 14.4, ...
        'trigger_mode', "current_touch");
    Specs(end + 1) = struct('name', "Beibu Gulf", ...
        'lon_min', 106.0, 'lon_max', 111.8, 'lat_min', 16.0, 'lat_max', 22.8, ...
        'lon_center', 108.8, 'lat_center', 19.6, ...
        'trigger_lon_min', 106.0, 'trigger_lon_max', 111.8, ...
        'trigger_lat_min', 16.0, 'trigger_lat_max', 22.8, ...
        'trigger_mode', "current_touch");
    Specs(end + 1) = struct('name', "Bohai Sea", ...
        'lon_min', 116.8, 'lon_max', 123.8, 'lat_min', 36.8, 'lat_max', 41.7, ...
        'lon_center', 120.3, 'lat_center', 39.1, ...
        'trigger_lon_min', 116.8, 'trigger_lon_max', 123.8, ...
        'trigger_lat_min', 36.8, 'trigger_lat_max', 41.7, ...
        'trigger_mode', "current_touch");
end
end

%% ============================================================
function [sectorCount, oppositePairCount] = component_touch_sector_metrics(idx, touchMask, maskSize)

sectorCount = 0;
oppositePairCount = 0;
idx = idx(:);
tidx = idx(touchMask(idx));
if numel(tidx) < 3
    return;
end

[rAll, cAll] = ind2sub(maskSize, idx);
[rt, ct] = ind2sub(maskSize, tidx);
cr = mean(double(rAll), 'omitnan');
cc = mean(double(cAll), 'omitnan');
ang = atan2(-(double(rt) - cr), double(ct) - cc);
sectorWidth = 2 * pi / 8;
sid = mod(floor((ang + pi) / sectorWidth), 8) + 1;
hit = false(1, 8);
hit(unique(sid)) = true;
sectorCount = nnz(hit);
oppositePairCount = double(hit(1) & hit(5)) + double(hit(2) & hit(6)) + ...
    double(hit(3) & hit(7)) + double(hit(4) & hit(8));
end

%% ============================================================
function bayDomain = bay_mask_to_smoothed_domain(bayMask, searchPoly, ctx, P)

bayDomain = polyshape();
if ~any(bayMask(:))
    return;
end

bayPoly = grid_mask_to_polyshape_local(bayMask, ctx.xv, ctx.yv);
if area(bayPoly) <= 0
    return;
end

bayDomain = safe_polybuffer(bayPoly, max(P.inner_refinement_buffer_km, 0.55 * P.outer_backshore_limit_km));
if area(bayDomain) <= 0
    return;
end
bayDomain = intersect(bayDomain, searchPoly);
if area(bayDomain) <= 0
    return;
end

baySmooth = 0;
if isfield(P, 'topology_bay_domain_smooth_km')
    baySmooth = max(0, P.topology_bay_domain_smooth_km);
end

if baySmooth > 0
    bayDomain = smooth_bay_domain_ocean_side(bayDomain, searchPoly, ctx, P, baySmooth);
end
end

%% ============================================================
function Pout = smooth_bay_domain_ocean_side(rawBayDomain, searchPoly, ctx, P, baySmooth)

Pout = rawBayDomain;
if baySmooth <= 0 || area(rawBayDomain) <= 0
    return;
end

oceanOnly = true;
if isfield(P, 'topology_bay_domain_smooth_ocean_only')
    oceanOnly = logical(P.topology_bay_domain_smooth_ocean_only);
end
if ~oceanOnly
    Pout = smooth_bay_domain_full(rawBayDomain, searchPoly, P, baySmooth);
    return;
end

protectKm = max(2 * ctx.dx, 0.9 * baySmooth);
if isfield(P, 'topology_bay_domain_land_protect_km')
    protectKm = max(0, P.topology_bay_domain_land_protect_km);
end
insetFraction = 0.70;
if isfield(P, 'topology_bay_domain_smooth_inset_fraction')
    insetFraction = max(0.35, min(0.98, P.topology_bay_domain_smooth_inset_fraction));
end
coverMin = 0.94;
if isfield(P, 'topology_bay_domain_smooth_cover_min')
    coverMin = max(0.80, min(0.995, P.topology_bay_domain_smooth_cover_min));
end
maxIncrease = 0.35;
if isfield(P, 'topology_bay_domain_smooth_max_area_increase_fraction')
    maxIncrease = max(0, min(0.80, P.topology_bay_domain_smooth_max_area_increase_fraction));
end

if protectKm <= 0
    Pout = smooth_bay_domain_full(rawBayDomain, searchPoly, P, baySmooth);
    return;
end

try
    protectMask = ctx.landMask & ctx.domainMask;
    if ~any(protectMask(:))
        Pout = smooth_bay_domain_full(rawBayDomain, searchPoly, P, baySmooth);
        return;
    end
    protectPoly = grid_mask_to_polyshape_local(protectMask, ctx.xv, ctx.yv);
    if area(protectPoly) <= 0
        Pout = smooth_bay_domain_full(rawBayDomain, searchPoly, P, baySmooth);
        return;
    end
    protectPoly = safe_polybuffer(protectPoly, protectKm);
    protectedBay = intersect(rawBayDomain, protectPoly);
    oceanBay = subtract(rawBayDomain, protectPoly);
    if area(oceanBay) <= 0
        return;
    end

    roundedOcean = safe_polybuffer(oceanBay, baySmooth);
    roundedOcean = safe_polybuffer(roundedOcean, -insetFraction * baySmooth);
    if area(roundedOcean) <= 0
        return;
    end
    roundedOcean = subtract(roundedOcean, protectPoly);
    roundedOcean = intersect(roundedOcean, searchPoly);
    if area(roundedOcean) <= 0
        return;
    end

    candidate = rmholes(union(protectedBay, roundedOcean));
    candidate = intersect(candidate, searchPoly);
    candidate = rmholes(union(candidate));
    if area(candidate) <= 0
        return;
    end

    candidate = smooth_bay_domain_curve(candidate, protectedBay, protectPoly, searchPoly, P, baySmooth);

    rawArea = area(rawBayDomain);
    newArea = area(candidate);
    coverFrac = area(intersect(candidate, rawBayDomain)) / max(rawArea, eps);
    if coverFrac < coverMin || newArea < 0.82 * rawArea || newArea > (1 + maxIncrease) * rawArea
        return;
    end
    Pout = candidate;
catch
    Pout = rawBayDomain;
end
end

%% ============================================================
function Pout = smooth_bay_domain_curve(Pin, protectedBay, protectPoly, searchPoly, P, baySmooth)

Pout = Pin;
if area(Pin) <= 0 || baySmooth <= 0
    return;
end
if ~isfield(P, 'topology_bay_domain_curve_smooth_enable') || ...
        ~P.topology_bay_domain_curve_smooth_enable
    return;
end

curveP = P;
curveP.final_open_boundary_curve_smooth_enable = true;
curveP.final_open_boundary_curve_step_km = 30;
curveP.final_open_boundary_curve_window_km = max(2.0 * baySmooth, 480);
curveP.final_open_boundary_curve_cover_pad_km = 0;
curveP.final_open_boundary_curve_cover_min = 0.86;
curveP.final_open_boundary_curve_max_area_increase_fraction = 0.55;
if isfield(P, 'topology_bay_domain_curve_step_km')
    curveP.final_open_boundary_curve_step_km = max(12, P.topology_bay_domain_curve_step_km);
end
if isfield(P, 'topology_bay_domain_curve_window_km')
    curveP.final_open_boundary_curve_window_km = max(curveP.final_open_boundary_curve_step_km, ...
        P.topology_bay_domain_curve_window_km);
end
if isfield(P, 'topology_bay_domain_curve_cover_min')
    curveP.final_open_boundary_curve_cover_min = max(0.75, min(0.98, ...
        P.topology_bay_domain_curve_cover_min));
end
if isfield(P, 'topology_bay_domain_curve_max_area_increase_fraction')
    curveP.final_open_boundary_curve_max_area_increase_fraction = max(0, min(0.80, ...
        P.topology_bay_domain_curve_max_area_increase_fraction));
end

try
    curve = smooth_outer_boundary_curve(Pin, polyshape(), protectedBay, curveP);
    if area(curve) <= 0
        return;
    end
    curveOcean = subtract(curve, protectPoly);
    curveOcean = intersect(curveOcean, searchPoly);
    candidate = rmholes(union(protectedBay, curveOcean));
    candidate = intersect(candidate, searchPoly);
    candidate = rmholes(union(candidate));
    if area(candidate) <= 0
        return;
    end
    oldArea = area(Pin);
    newArea = area(candidate);
    coverFrac = area(intersect(candidate, Pin)) / max(oldArea, eps);
    if coverFrac < curveP.final_open_boundary_curve_cover_min || ...
            newArea < 0.78 * oldArea || ...
            newArea > (1 + curveP.final_open_boundary_curve_max_area_increase_fraction) * oldArea
        return;
    end
    Pout = candidate;
catch
    Pout = Pin;
end
end

%% ============================================================
function Pout = smooth_bay_domain_final_outline(Pin, searchPoly, ctx, P)

Pout = Pin;
if area(Pin) <= 0 || ~isfield(P, 'topology_bay_domain_final_smooth_enable') || ...
        ~P.topology_bay_domain_final_smooth_enable
    return;
end

smoothKm = 0;
if isfield(P, 'topology_bay_domain_final_smooth_km')
    smoothKm = max(0, P.topology_bay_domain_final_smooth_km);
end
if smoothKm <= 0
    return;
end

protectKm = 35;
if isfield(P, 'topology_bay_domain_final_land_protect_km')
    protectKm = max(0, P.topology_bay_domain_final_land_protect_km);
end
insetFraction = 0.88;
if isfield(P, 'topology_bay_domain_final_inset_fraction')
    insetFraction = max(0.55, min(1.05, P.topology_bay_domain_final_inset_fraction));
end
coverMin = 0.72;
if isfield(P, 'topology_bay_domain_final_cover_min')
    coverMin = max(0.45, min(0.98, P.topology_bay_domain_final_cover_min));
end
minFrac = 0.60;
if isfield(P, 'topology_bay_domain_final_min_area_fraction')
    minFrac = max(0.30, min(0.98, P.topology_bay_domain_final_min_area_fraction));
end
maxIncrease = 0.85;
if isfield(P, 'topology_bay_domain_final_max_area_increase_fraction')
    maxIncrease = max(0, min(1.50, P.topology_bay_domain_final_max_area_increase_fraction));
end

try
    protectedBay = polyshape();
    protectPoly = polyshape();
    if protectKm > 0 && isfield(ctx, 'landMask') && any(ctx.landMask(:))
        protectMask = ctx.landMask & ctx.domainMask;
        protectPoly = grid_mask_to_polyshape_local(protectMask, ctx.xv, ctx.yv);
        if area(protectPoly) > 0
            protectPoly = safe_polybuffer(protectPoly, protectKm);
            protectedBay = intersect(Pin, protectPoly);
        end
    end

    oceanBay = Pin;
    if area(protectPoly) > 0
        oceanBay = subtract(Pin, protectPoly);
    end
    if area(oceanBay) <= 0
        return;
    end

    roundedOcean = safe_polybuffer(oceanBay, smoothKm);
    roundedOcean = safe_polybuffer(roundedOcean, -insetFraction * smoothKm);
    if area(protectPoly) > 0
        roundedOcean = subtract(roundedOcean, protectPoly);
    end
    roundedOcean = intersect(roundedOcean, searchPoly);
    if area(roundedOcean) <= 0
        return;
    end

    candidate = rmholes(union(protectedBay, roundedOcean));
    candidate = intersect(candidate, searchPoly);
    candidate = rmholes(union(candidate));
    if area(candidate) <= 0
        return;
    end

    if isfield(P, 'topology_bay_domain_final_curve_enable') && ...
            P.topology_bay_domain_final_curve_enable
        curveP = P;
        curveP.final_open_boundary_curve_step_km = 35;
        curveP.final_open_boundary_curve_window_km = 900;
        curveP.final_open_boundary_curve_cover_pad_km = 0;
        curveP.final_open_boundary_curve_cover_min = 0.70;
        curveP.final_open_boundary_curve_max_area_increase_fraction = 0.90;
        if isfield(P, 'topology_bay_domain_final_curve_step_km')
            curveP.final_open_boundary_curve_step_km = max(12, P.topology_bay_domain_final_curve_step_km);
        end
        if isfield(P, 'topology_bay_domain_final_curve_window_km')
            curveP.final_open_boundary_curve_window_km = max(curveP.final_open_boundary_curve_step_km, ...
                P.topology_bay_domain_final_curve_window_km);
        end
        if isfield(P, 'topology_bay_domain_final_curve_cover_min')
            curveP.final_open_boundary_curve_cover_min = max(0.45, min(0.98, ...
                P.topology_bay_domain_final_curve_cover_min));
        end
        if isfield(P, 'topology_bay_domain_final_curve_max_area_increase_fraction')
            curveP.final_open_boundary_curve_max_area_increase_fraction = max(0, min(1.50, ...
                P.topology_bay_domain_final_curve_max_area_increase_fraction));
        end

        curved = smooth_outer_boundary_curve(candidate, polyshape(), protectedBay, curveP);
        if area(curved) > 0
            curvedOcean = curved;
            if area(protectPoly) > 0
                curvedOcean = subtract(curved, protectPoly);
            end
            curvedOcean = intersect(curvedOcean, searchPoly);
            curvedCandidate = rmholes(union(protectedBay, curvedOcean));
            curvedCandidate = intersect(curvedCandidate, searchPoly);
            curvedCandidate = rmholes(union(curvedCandidate));
            if bay_domain_final_candidate_ok(curvedCandidate, Pin, coverMin, minFrac, maxIncrease)
                candidate = curvedCandidate;
            end
        end
    end

    if bay_domain_final_candidate_ok(candidate, Pin, coverMin, minFrac, maxIncrease)
        Pout = candidate;
    end
catch
    Pout = Pin;
end
end

%% ============================================================
function ok = bay_domain_final_candidate_ok(candidate, original, coverMin, minFrac, maxIncrease)

ok = false;
try
    oldArea = area(original);
    newArea = area(candidate);
    coverFrac = area(intersect(candidate, original)) / max(oldArea, eps);
catch
    return;
end
if newArea <= 0
    return;
end
if coverFrac < coverMin || newArea < minFrac * oldArea || newArea > (1 + maxIncrease) * oldArea
    return;
end
ok = true;
end

%% ============================================================
function Pout = smooth_bay_domain_full(rawBayDomain, searchPoly, P, baySmooth)

Pout = rawBayDomain;
if baySmooth <= 0 || area(rawBayDomain) <= 0
    return;
end

insetFraction = 0.70;
if isfield(P, 'topology_bay_domain_smooth_inset_fraction')
    insetFraction = max(0.35, min(0.98, P.topology_bay_domain_smooth_inset_fraction));
end
coverMin = 0.94;
if isfield(P, 'topology_bay_domain_smooth_cover_min')
    coverMin = max(0.80, min(0.995, P.topology_bay_domain_smooth_cover_min));
end
maxIncrease = 0.35;
if isfield(P, 'topology_bay_domain_smooth_max_area_increase_fraction')
    maxIncrease = max(0, min(0.80, P.topology_bay_domain_smooth_max_area_increase_fraction));
end

try
    rounded = safe_polybuffer(rawBayDomain, baySmooth);
    rounded = safe_polybuffer(rounded, -insetFraction * baySmooth);
    rounded = intersect(rounded, searchPoly);
    if area(rounded) <= 0
        return;
    end
    rawArea = area(rawBayDomain);
    newArea = area(rounded);
    coverFrac = area(intersect(rounded, rawBayDomain)) / max(rawArea, eps);
    if coverFrac < coverMin || newArea < 0.82 * rawArea || newArea > (1 + maxIncrease) * rawArea
        return;
    end
    Pout = rmholes(rounded);
catch
    Pout = rawBayDomain;
end
end

%% ============================================================
function [addMask, info] = select_topologic_bay_additions(ctx, clippedPoly, outerPoly, P)

info = struct();
info.core_cell_count = 0;
info.accepted_component_count = 0;

[X, Y] = meshgrid(ctx.xv, ctx.yv);
currentMask = polyshape_grid_mask_local(clippedPoly, X, Y);
searchMask = polyshape_grid_mask_local(safe_polybuffer(outerPoly, P.topology_bay_extra_km), X, Y);
searchEdgeMask = mask_boundary_local(searchMask);
addMask = false(size(currentMask));

closeCells = max(2, round(P.topology_bay_ring_close_km / ctx.dx));
recoverCells = max(closeCells, round(P.topology_bay_ring_recover_km / ctx.dx));
touchCells = max(1, round(P.topology_bay_ring_touch_km / ctx.dx));
bayOcean = current_connected_ocean_mask(ctx, currentMask, searchMask, touchCells);

% Ignore small islands in the bay detector.  Only major land may close a
% water body; this prevents ordinary island chains from creating fake bays.
ringLand = dilate_mask_local(ctx.majorLandMask, closeCells);
closedWater = bayOcean & searchMask & ~ringLand;
components = connected_components_pixelidx(closedWater);
if isempty(components)
    return;
end

touchBand = dilate_mask_local(currentMask, touchCells);
cellArea = ctx.dx * ctx.dx;
minCells = max(topology_min_component_cells(P), ceil(P.topology_bay_ring_min_area_km2 / cellArea));
acceptedAreaKm2 = 0;

for i = 1:numel(components)
    idx = components{i};
    if numel(idx) < minCells
        continue;
    end
    coreAreaKm2 = numel(idx) * cellArea;
    if coreAreaKm2 > P.topology_bay_ring_max_area_km2
        continue;
    end
    edgeFracCore = nnz(searchEdgeMask(idx)) / max(1, numel(idx));
    if edgeFracCore > P.topology_bay_ring_search_edge_touch_max_fraction
        continue;
    end
    if nnz(touchBand(idx)) < P.topology_bay_ring_current_touch_min_cells
        continue;
    end

    coreMask = false(size(currentMask));
    coreMask(idx) = true;
    distCore = grid_distance_steps_from_sources(coreMask, recoverCells);
    bayMask = bayOcean & searchMask & ~currentMask & (distCore <= recoverCells);
    bayMask = flood_fill_from_seed_mask(bayMask | coreMask, coreMask) & ...
        bayOcean & searchMask & ~currentMask & (distCore <= recoverCells);
    if ~any(bayMask(:))
        continue;
    end
    bayAreaKm2 = nnz(bayMask) * cellArea;
    if bayAreaKm2 < P.topology_bay_ring_min_area_km2 || bayAreaKm2 > P.topology_bay_ring_max_area_km2
        continue;
    end
    if acceptedAreaKm2 + bayAreaKm2 > P.topology_bay_ring_total_max_area_km2
        continue;
    end
    edgeFracBay = nnz(searchEdgeMask & bayMask) / max(1, nnz(bayMask));
    if edgeFracBay > P.topology_bay_ring_search_edge_touch_max_fraction
        continue;
    end
    if nnz(touchBand & bayMask) < P.topology_bay_ring_current_touch_min_cells
        continue;
    end

    addMask = addMask | bayMask;
    acceptedAreaKm2 = acceptedAreaKm2 + bayAreaKm2;
    info.core_cell_count = info.core_cell_count + numel(idx);
    info.accepted_component_count = info.accepted_component_count + 1;
end
end

%% ============================================================
function n = topology_min_component_cells(P)

n = max(6, round((70 / max(P.topology_grid_km, 1))^2));
end

%% ============================================================
function [majorSpanKm, minorSpanKm] = component_principal_spans_km(idx, maskSize, dxKm)

[r, c] = ind2sub(maskSize, idx(:));
if numel(r) < 3
    majorSpanKm = dxKm;
    minorSpanKm = dxKm;
    return;
end

x = (double(c) - mean(double(c), 'omitnan')) * dxKm;
y = (double(r) - mean(double(r), 'omitnan')) * dxKm;
pts = [x(:), y(:)];
try
    C = cov(pts);
    [V, D] = eig(C);
    [~, order] = sort(diag(D), 'descend');
    V = V(:, order);
    proj = pts * V;
    span1 = max(proj(:, 1), [], 'omitnan') - min(proj(:, 1), [], 'omitnan') + dxKm;
    span2 = max(proj(:, 2), [], 'omitnan') - min(proj(:, 2), [], 'omitnan') + dxKm;
    majorSpanKm = max(span1, span2);
    minorSpanKm = max(dxKm, min(span1, span2));
catch
    majorSpanKm = (max(c) - min(c) + 1) * dxKm;
    minorSpanKm = (max(r) - min(r) + 1) * dxKm;
    if minorSpanKm > majorSpanKm
        tmp = minorSpanKm;
        minorSpanKm = majorSpanKm;
        majorSpanKm = tmp;
    end
end
end

%% ============================================================
function dist = grid_distance_steps_from_sources(sourceMask, maxSteps)

dist = inf(size(sourceMask));
source = find(sourceMask);
if isempty(source)
    return;
end

[nrow, ncol] = size(sourceMask);
queue = zeros(numel(sourceMask), 1);
head = 1;
tail = numel(source);
queue(1:tail) = source(:);
dist(source) = 0;

while head <= tail
    q = queue(head);
    head = head + 1;
    dq = dist(q);
    if dq >= maxSteps
        continue;
    end

    [r, c] = ind2sub([nrow, ncol], q);
    for dr = -1:1
        rr = r + dr;
        if rr < 1 || rr > nrow
            continue;
        end
        for dc = -1:1
            if dr == 0 && dc == 0
                continue;
            end
            cc = c + dc;
            if cc < 1 || cc > ncol
                continue;
            end
            ii = sub2ind([nrow, ncol], rr, cc);
            if dist(ii) > dq + 1
                dist(ii) = dq + 1;
                tail = tail + 1;
                queue(tail) = ii;
            end
        end
    end
end
end

%% ============================================================
function out = flood_fill_from_seed_mask(mask, seedMask)

out = false(size(mask));
seeds = find(mask & seedMask);
if isempty(seeds)
    return;
end

nmax = nnz(mask);
queue = zeros(nmax, 1);
head = 1;
tail = numel(seeds);
queue(1:tail) = seeds(:);
out(seeds) = true;
[nrow, ncol] = size(mask);

while head <= tail
    q = queue(head);
    head = head + 1;
    [r, c] = ind2sub([nrow, ncol], q);
    for dr = -1:1
        rr = r + dr;
        if rr < 1 || rr > nrow
            continue;
        end
        for dc = -1:1
            if dr == 0 && dc == 0
                continue;
            end
            cc = c + dc;
            if cc < 1 || cc > ncol
                continue;
            end
            ii = sub2ind([nrow, ncol], rr, cc);
            if mask(ii) && ~out(ii)
                tail = tail + 1;
                if tail > numel(queue)
                    queue(end + nmax, 1) = 0; %#ok<AGROW>
                end
                queue(tail) = ii;
                out(ii) = true;
            end
        end
    end
end
end

%% ============================================================
function components = connected_components_pixelidx(mask)

components = {};
remaining = find(mask);
seen = false(size(mask));
[nrow, ncol] = size(mask);

for is = 1:numel(remaining)
    seed = remaining(is);
    if seen(seed)
        continue;
    end
    queue = zeros(nnz(mask), 1);
    head = 1;
    tail = 1;
    queue(1) = seed;
    seen(seed) = true;
    pixels = zeros(nnz(mask), 1);
    np = 0;
    while head <= tail
        q = queue(head);
        head = head + 1;
        np = np + 1;
        pixels(np) = q;
        [r, c] = ind2sub([nrow, ncol], q);
        for dr = -1:1
            rr = r + dr;
            if rr < 1 || rr > nrow
                continue;
            end
            for dc = -1:1
                if dr == 0 && dc == 0
                    continue;
                end
                cc = c + dc;
                if cc < 1 || cc > ncol
                    continue;
                end
                ii = sub2ind([nrow, ncol], rr, cc);
                if mask(ii) && ~seen(ii)
                    tail = tail + 1;
                    queue(tail) = ii;
                    seen(ii) = true;
                end
            end
        end
    end
    components{end + 1} = pixels(1:np); %#ok<AGROW>
end
end

%% ============================================================
function out = dilate_mask_local(mask, radiusCells)

radiusCells = round(radiusCells);
if radiusCells <= 0
    out = mask;
    return;
end
[xx, yy] = meshgrid(-radiusCells:radiusCells, -radiusCells:radiusCells);
kernel = (xx.^2 + yy.^2) <= radiusCells^2;
out = conv2(double(mask), double(kernel), 'same') > 0;
end

%% ============================================================
function out = erode_mask_local(mask, radiusCells)

radiusCells = round(radiusCells);
if radiusCells <= 0
    out = mask;
    return;
end
[xx, yy] = meshgrid(-radiusCells:radiusCells, -radiusCells:radiusCells);
kernel = double((xx.^2 + yy.^2) <= radiusCells^2);
need = sum(kernel(:));
out = conv2(double(mask), kernel, 'same') >= need;
end

%% ============================================================
function out = open_mask_local(mask, radiusCells)

out = dilate_mask_local(erode_mask_local(mask, radiusCells), radiusCells);
end

%% ============================================================
function out = close_mask_local(mask, radiusCells)

out = erode_mask_local(dilate_mask_local(mask, radiusCells), radiusCells);
end

%% ============================================================
function out = fill_holes_mask_local(mask)

out = mask;
if ~any(mask(:))
    return;
end
background = ~mask;
edgeSeeds = false(size(mask));
edgeSeeds(1, :) = background(1, :);
edgeSeeds(end, :) = edgeSeeds(end, :) | background(end, :);
edgeSeeds(:, 1) = edgeSeeds(:, 1) | background(:, 1);
edgeSeeds(:, end) = edgeSeeds(:, end) | background(:, end);
outside = flood_fill_from_seed_mask(background, edgeSeeds);
holes = background & ~outside;
out = mask | holes;
end

%% ============================================================
function out = areaopen_mask_local(mask, minCells)

out = false(size(mask));
if ~any(mask(:))
    return;
end
minCells = max(1, round(minCells));
components = connected_components_pixelidx(mask);
for i = 1:numel(components)
    idx = components{i};
    if numel(idx) >= minCells
        out(idx) = true;
    end
end
end

%% ============================================================
function mask = polyshape_grid_mask_local(Ps, X, Y)

mask = false(size(X));
if area(Ps) <= 0
    return;
end
try
    mask(:) = isinterior(Ps, X(:), Y(:));
catch
    mask(:) = false;
end
end

%% ============================================================
function Pout = grid_mask_to_polyshape_local(mask, xv, yv)

Pout = polyshape();
if ~any(mask(:)) || numel(xv) < 2 || numel(yv) < 2
    return;
end

dx = median(diff(xv), 'omitnan');
dy = median(diff(yv), 'omitnan');
xedge = [xv(1) - dx / 2, (xv(1:end-1) + xv(2:end)) / 2, xv(end) + dx / 2];
yedge = [yv(1) - dy / 2, (yv(1:end-1) + yv(2:end)) / 2, yv(end) + dy / 2];
[nrow, ncol] = size(mask);

starts = zeros(nnz(mask) * 4, 1);
ends = zeros(nnz(mask) * 4, 1);
ne = 0;
[rr, cc] = find(mask);
for k = 1:numel(rr)
    r = rr(k);
    c = cc(k);
    if r == 1 || ~mask(r - 1, c)
        ne = ne + 1; starts(ne) = node_id(c - 1, r - 1, ncol); ends(ne) = node_id(c, r - 1, ncol);
    end
    if c == ncol || ~mask(r, c + 1)
        ne = ne + 1; starts(ne) = node_id(c, r - 1, ncol); ends(ne) = node_id(c, r, ncol);
    end
    if r == nrow || ~mask(r + 1, c)
        ne = ne + 1; starts(ne) = node_id(c, r, ncol); ends(ne) = node_id(c - 1, r, ncol);
    end
    if c == 1 || ~mask(r, c - 1)
        ne = ne + 1; starts(ne) = node_id(c - 1, r, ncol); ends(ne) = node_id(c - 1, r - 1, ncol);
    end
end
starts = starts(1:ne);
ends = ends(1:ne);
used = false(ne, 1);

for e0 = 1:ne
    if used(e0)
        continue;
    end
    loopNodes = starts(e0);
    used(e0) = true;
    cur = ends(e0);
    guard = 0;
    while cur ~= loopNodes(1) && guard < ne + 5
        loopNodes(end + 1, 1) = cur; %#ok<AGROW>
        nxt = find(starts == cur & ~used, 1, 'first');
        if isempty(nxt)
            break;
        end
        used(nxt) = true;
        cur = ends(nxt);
        guard = guard + 1;
    end
    if cur == loopNodes(1) && numel(loopNodes) >= 4
        [x, y] = node_ids_to_xy(loopNodes, xedge, yedge, ncol);
        try
            p = polyshape(x(:), y(:), 'Simplify', true);
            if area(p) > 0
                Pout = union(Pout, p);
            end
        catch
        end
    end
end
Pout = rmholes(union(Pout));
end

%% ============================================================
function id = node_id(cNode, rNode, ncol)

id = rNode * (ncol + 1) + cNode + 1;
end

%% ============================================================
function [x, y] = node_ids_to_xy(ids, xedge, yedge, ncol)

cNode = mod(ids - 1, ncol + 1);
rNode = floor((ids - 1) / (ncol + 1));
x = xedge(cNode + 1).';
y = yedge(rNode + 1).';
end

%% ============================================================
function corridorPoly = make_backshore_clip_corridor(outerPoly, studyPoly, innerPoly, basinId, lon0, lat0, axisV, outerBufferKm, P)

corridorPoly = outerPoly;
if area(outerPoly) <= 0
    return;
end

[xs, ys] = boundary(studyPoly);
vs = isfinite(xs) & isfinite(ys);
if nnz(vs) < 3
    return;
end

axisV = axisV(:) / max(norm(axisV), eps);
tSeed = xs(vs) .* axisV(1) + ys(vs) .* axisV(2);

[refLon, refLat, sideKnown] = basin_ocean_reference_lonlat(basinId, lon0, lat0);
if ~sideKnown || isempty(tSeed)
    return;
end

[rx, ry] = lonlat_to_local_km(refLon, refLat, lon0, lat0);
tRef = rx .* axisV(1) + ry .* axisV(2);
tMid = median(tSeed, 'omitnan');
if tRef >= tMid
    oceanVec = axisV;
else
    oceanVec = -axisV;
end

backPad = max(P.inner_refinement_buffer_km + 40, P.outer_backshore_limit_km);
corridorPoly = union(safe_polybuffer(studyPoly, backPad), innerPoly);

shiftFractions = [0.35, 0.70, 1.05];
bufferFractions = [0.82, 0.78, 0.70];
for i = 1:numel(shiftFractions)
    d = shiftFractions(i) * outerBufferKm;
    shifted = translate_polyshape_local(studyPoly, oceanVec(1) * d, oceanVec(2) * d);
    if area(shifted) <= 0
        continue;
    end
    corridorPoly = union(corridorPoly, safe_polybuffer(shifted, bufferFractions(i) * outerBufferKm));
end
corridorPoly = rmholes(union(corridorPoly));
if area(corridorPoly) <= 0 || area(intersect(corridorPoly, innerPoly)) < 0.85 * area(innerPoly)
    corridorPoly = outerPoly;
end
end

%% ============================================================
function Pout = translate_polyshape_local(Pin, dx, dy)

Pout = polyshape();
try
    rr = regions(Pin);
catch
    rr = Pin;
end
for i = 1:numel(rr)
    [x, y] = boundary(rr(i));
    v = isfinite(x) & isfinite(y);
    if nnz(v) < 3
        continue;
    end
    try
        q = polyshape(x(v) + dx, y(v) + dy, 'Simplify', true);
        if area(q) > 0
            Pout = union(Pout, q);
        end
    catch
    end
end
end

function [refLon, refLat, known] = basin_ocean_reference_lonlat(basinId, lon0, lat0)

bid = upper(string(basinId));
refLon = lon0;
refLat = lat0;
known = true;

switch bid
    case "BASIN_ENP"
        refLon = lon0 - 8;
        refLat = lat0 - 2;
    case "BASIN_WNP"
        refLon = lon0 + 8;
        refLat = lat0 - 1;
    case "BASIN_NATL"
        if lon0 < -84 && lat0 <= 22
            refLon = lon0 + 5;       % Central America Caribbean side
            refLat = lat0 + 4;
        elseif lon0 < -82 && lat0 > 22
            refLon = lon0;
            refLat = lat0 - 6;       % Gulf of Mexico / Caribbean side
        elseif lon0 < -68 && lat0 >= 20
            refLon = lon0 + 7;       % North American Atlantic side
            refLat = lat0;
        else
            refLon = lon0;
            refLat = lat0 - 5;
        end
    case "BASIN_NIO"
        if lon0 < 74
            refLon = lon0 - 6;       % Arabian Sea side
            refLat = lat0 - 2;
        else
            refLon = lon0 + 3;       % Bay of Bengal / eastern Indian side
            refLat = lat0 - 6;
        end
    case "BASIN_AUSSP"
        if lon0 < 130
            refLon = lon0 - 5;
            refLat = lat0 - 4;
        else
            refLon = lon0 + 7;
            refLat = lat0;
        end
    case "BASIN_SIO"
        if lon0 < 55
            refLon = lon0 - 6;
        else
            refLon = lon0 + 6;
        end
        refLat = lat0;
    otherwise
        known = false;
end

refLon = wrapTo180_local(refLon);
refLat = max(-80, min(80, refLat));
end

%% ============================================================
function [refLon, refLat, known] = basin_tc_approach_reference_lonlat(basinId, lon0, lat0)

bid = upper(string(basinId));
refLon = lon0;
refLat = lat0;
known = true;

switch bid
    case "BASIN_ENP"
        refLon = lon0 - 9;
        refLat = lat0 - 2;
    case "BASIN_WNP"
        refLon = lon0 + 9;
        refLat = lat0 - 3;
    case "BASIN_NATL"
        if lon0 < -84 && lat0 <= 22
            refLon = lon0 + 7;       % west Caribbean: mainly from east/southeast
            refLat = lat0 + 1;
        elseif lon0 < -82 && lat0 > 22
            refLon = lon0 + 4;       % Gulf of Mexico: mainly from south/southeast
            refLat = lat0 - 6;
        elseif lon0 < -68 && lat0 >= 20
            refLon = lon0 + 8;       % western Atlantic: mainly from east/southeast
            refLat = lat0 - 2;
        else
            refLon = lon0 + 3;
            refLat = lat0 - 5;
        end
    case "BASIN_NIO"
        if lon0 < 74
            refLon = lon0 - 3;       % Arabian Sea: mainly from south/southwest
            refLat = lat0 - 7;
        else
            refLon = lon0 + 2;       % Bay of Bengal: mainly from south/southeast
            refLat = lat0 - 7;
        end
    case "BASIN_AUSSP"
        if lon0 < 130
            refLon = lon0 - 4;
            refLat = lat0 - 6;
        else
            refLon = lon0 + 6;
            refLat = lat0 - 3;
        end
    case "BASIN_SIO"
        if lon0 < 55
            refLon = lon0 - 4;
        else
            refLon = lon0 + 4;
        end
        refLat = lat0 - 6;
    otherwise
        known = false;
end

refLon = wrapTo180_local(refLon);
refLat = max(-80, min(80, refLat));
end

%% ============================================================
function Pout = keep_largest_relevant_region(Pin, innerPoly, studyPoly, P)

Pout = Pin;
if ~isfield(P, 'keep_largest_outer_region') || ~P.keep_largest_outer_region || area(Pin) <= 0
    return;
end

try
    rr = regions(Pin);
catch
    return;
end
if numel(rr) <= 1
    return;
end

if isfield(P, 'final_keep_only_inner_touching_regions') && P.final_keep_only_inner_touching_regions
    innerTouchBuffer = 0;
    if isfield(P, 'final_inner_touch_buffer_km')
        innerTouchBuffer = max(0, P.final_inner_touch_buffer_km);
    end
    innerTouchMinArea = 1;
    if isfield(P, 'final_inner_touch_min_area_km2')
        innerTouchMinArea = max(0, P.final_inner_touch_min_area_km2);
    end
    innerTouchPoly = innerPoly;
    if innerTouchBuffer > 0
        innerTouchPoly = safe_polybuffer(innerPoly, innerTouchBuffer);
    end
    keepInner = false(numel(rr), 1);
    for i = 1:numel(rr)
        try
            directTouchArea = area(intersect(rr(i), innerPoly));
            bufferedTouchArea = area(intersect(rr(i), innerTouchPoly));
            keepInner(i) = directTouchArea >= innerTouchMinArea || bufferedTouchArea >= innerTouchMinArea;
        catch
            keepInner(i) = false;
        end
    end
    if any(keepInner)
        Pkeep = polyshape();
        for k = find(keepInner(:).')
            Pkeep = union(Pkeep, rr(k));
        end
        if area(Pkeep) > 0
            Pout = rmholes(union(Pkeep));
            return;
        end
    end
end

touchArea = zeros(numel(rr), 1);
regionArea = zeros(numel(rr), 1);
for i = 1:numel(rr)
    regionArea(i) = area(rr(i));
    try
        touchArea(i) = area(intersect(rr(i), innerPoly)) + area(intersect(rr(i), studyPoly));
    catch
        touchArea(i) = 0;
    end
end

keep = touchArea > max(1, 0.01 * max(touchArea));
if ~any(keep)
    [~, best] = max(regionArea);
    Pout = rr(best);
    return;
end

kept = find(keep);
if numel(kept) == 1
    Pout = rr(kept);
    return;
end

if isfield(P, 'force_single_outer_region') && P.force_single_outer_region
    score = touchArea(kept) + 0.05 * regionArea(kept);
    [~, ibest] = max(score);
    Pout = rr(kept(ibest));
    return;
end

% Preserve separate pieces only when they all contain target refinement
% geometry.  Otherwise keep the largest relevant continuous domain.
[~, order] = sort(regionArea(kept), 'descend');
best = kept(order(1));
innerCoverage = area(intersect(rr(best), innerPoly)) / max(area(innerPoly), eps);
if innerCoverage >= 0.80
    Pout = rr(best);
else
    Pout = polyshape();
    for k = kept(:).'
        Pout = union(Pout, rr(k));
    end
end
end

%% ============================================================
function [lon, lat] = oriented_capsule_lonlat(lon0, lat0, axisU, alongHalfKm, crossHalfKm, n)

axisU = axisU(:) / norm(axisU);
axisV = [-axisU(2); axisU(1)];

r = max(1, crossHalfKm);
straightHalf = max(0, alongHalfKm - r);
nArc = max(16, round(n / 2));
ang1 = linspace(-pi/2, pi/2, nArc).';
ang2 = linspace(pi/2, 3*pi/2, nArc).';

p1 = [straightHalf + r * cos(ang1), r * sin(ang1)];
p2 = [-straightHalf + r * cos(ang2), r * sin(ang2)];
p = [p1; p2; p1(1, :)];

xy = p(:, 1) .* axisU.' + p(:, 2) .* axisV.';
[lon, lat] = local_km_to_lonlat(xy(:, 1), xy(:, 2), lon0, lat0);
end

%% ============================================================
function Pout = domain_polygons_to_local_polyshape(domain_ids, DomainPolys, lon0, lat0)

Pout = polyshape();
ids = string(domain_ids(:));
allIds = string({DomainPolys.id});
for i = 1:numel(ids)
    idx = find(allIds == ids(i), 1, 'first');
    if isempty(idx)
        continue;
    end
    lon = DomainPolys(idx).lon(:);
    lat = DomainPolys(idx).lat(:);
    v = isfinite(lon) & isfinite(lat);
    if nnz(v) < 3
        continue;
    end
    [x, y] = lonlat_to_local_km(lon(v), lat(v), lon0, lat0);
    try
        p = polyshape(x, y, 'Simplify', true);
        if ~isempty(p.Vertices) && area(p) > 0
            Pout = union(Pout, p);
        end
    catch
    end
end
end

%% ============================================================
function Pout = cells_to_local_polyshape(C, lon0, lat0)

Pout = polyshape();
if isempty(C)
    return;
end
halfDeg = 0.125;
for i = 1:height(C)
    lon = double(C.lon(i));
    lat = double(C.lat(i));
    xx = [lon-halfDeg; lon+halfDeg; lon+halfDeg; lon-halfDeg; lon-halfDeg];
    yy = [lat-halfDeg; lat-halfDeg; lat+halfDeg; lat+halfDeg; lat-halfDeg];
    [x, y] = lonlat_to_local_km(xx, yy, lon0, lat0);
    try
        Pout = union(Pout, polyshape(x, y, 'Simplify', true));
    catch
    end
end
Pout = union(Pout);
end

%% ============================================================
function Pout = lonlat_vectors_to_local_polyshape(lon, lat, lon0, lat0)

Pout = polyshape();
segments = nan_split_segments(lon(:), lat(:));
for i = 1:numel(segments)
    seg = segments{i};
    if size(seg, 1) < 3
        continue;
    end
    [x, y] = lonlat_to_local_km(seg(:, 1), seg(:, 2), lon0, lat0);
    try
        p = polyshape(x, y, 'Simplify', true);
        if area(p) > 0
            Pout = union(Pout, p);
        end
    catch
    end
end
Pout = union(Pout);
end

%% ============================================================
function [lon, lat] = local_polyshape_to_lonlat_vectors(P, lon0, lat0)

lon = [];
lat = [];
rr = regions(P);
for i = 1:numel(rr)
    [x, y] = boundary(rr(i));
    v = isfinite(x) & isfinite(y);
    if nnz(v) < 3
        continue;
    end
    [lo, la] = local_km_to_lonlat(x(v), y(v), lon0, lat0);
    if ~isempty(lon)
        lon(end + 1, 1) = NaN; %#ok<AGROW>
        lat(end + 1, 1) = NaN; %#ok<AGROW>
    end
    lon = [lon; lo(:)]; %#ok<AGROW>
    lat = [lat; la(:)]; %#ok<AGROW>
end
end

%% ============================================================
function [alongHalf, crossHalf] = half_extents_along_axis(x, y, axisU, axisV)

v = isfinite(x) & isfinite(y);
if nnz(v) < 2
    alongHalf = 0;
    crossHalf = 0;
    return;
end
x = x(v);
y = y(v);
s = x .* axisU(1) + y .* axisU(2);
t = x .* axisV(1) + y .* axisV(2);
alongHalf = 0.5 * (max(s) - min(s));
crossHalf = 0.5 * (max(t) - min(t));
end

%% ============================================================
function P2 = safe_polybuffer(P1, d)

if area(P1) == 0
    P2 = P1;
    return;
end
try
    P2 = polybuffer(P1, d, 'JointType', 'round');
catch
    P2 = P1;
end
end

%% ============================================================
function shelfWidthKm = estimate_shelf_width_km(P, studyPoly, lon0, lat0)

shelfWidthKm = P.shelf_min_width_km;
if ~P.use_gebco_shelf_width || exist(P.gebco_path, 'file') ~= 2 || area(studyPoly) == 0
    return;
end

persistent lonVec latVec cachedPath cacheOk
try
    if isempty(cachedPath) || ~strcmp(cachedPath, P.gebco_path)
        lonVec = double(ncread(P.gebco_path, 'lon'));
        latVec = double(ncread(P.gebco_path, 'lat'));
        cachedPath = P.gebco_path;
        cacheOk = true;
    end
catch
    cacheOk = false;
end
if isempty(cacheOk) || ~cacheOk
    return;
end

[xb, yb] = boundary(studyPoly);
v = isfinite(xb) & isfinite(yb);
if nnz(v) < 3
    return;
end
[lonB, latB] = local_km_to_lonlat(xb(v), yb(v), lon0, lat0);
lonUn = unwrap_lon_around(lonB, lon0);
lonlim = [min(lonUn)-4, max(lonUn)+4];
latlim = [max(-80, min(latB)-4), min(80, max(latB)+4)];
if lonlim(1) < min(lonVec) || lonlim(2) > max(lonVec) || diff(lonlim) > 80
    return;
end

iAll = find(lonVec >= lonlim(1) & lonVec <= lonlim(2));
jAll = find(latVec >= latlim(1) & latVec <= latlim(2));
if numel(iAll) < 3 || numel(jAll) < 3
    return;
end

stride = max(P.shelf_sample_step, ceil(max(numel(iAll), numel(jAll)) / 900));
i0 = iAll(1);
j0 = jAll(1);
ni = floor((iAll(end) - i0) / stride) + 1;
nj = floor((jAll(end) - j0) / stride) + 1;
try
    z = double(ncread(P.gebco_path, 'elevation', [i0, j0], [ni, nj], [stride, stride]));
catch
    return;
end
lonSub = lonVec(i0 + (0:ni-1) * stride);
latSub = latVec(j0 + (0:nj-1) * stride);
[LON, LAT] = ndgrid(lonSub, latSub);
[xg, yg] = lonlat_to_local_km(LON(:), LAT(:), lon0, lat0);
zg = z(:);

shallowShelf = isfinite(zg) & zg <= 0 & zg >= P.shelf_isobath_m;
if nnz(shallowShelf) < 20
    return;
end
xq = xg(shallowShelf);
yq = yg(shallowShelf);
if numel(xq) > 25000
    keep = round(linspace(1, numel(xq), 25000));
    xq = xq(keep);
    yq = yq(keep);
end
d = min_distance_to_vertices_km(xq, yq, xb(v), yb(v));
d = d(isfinite(d) & d > 0 & d <= P.shelf_max_width_km + 150);
if numel(d) < 20
    return;
end
shelfWidthKm = prctile(d, 85);
shelfWidthKm = max(P.shelf_min_width_km, min(P.shelf_max_width_km, shelfWidthKm));
end

%% ============================================================
function dmin = min_distance_to_vertices_km(xq, yq, xv, yv)

xq = double(xq(:));
yq = double(yq(:));
xv = double(xv(:));
yv = double(yv(:));
v = isfinite(xv) & isfinite(yv);
xv = xv(v);
yv = yv(v);
dmin = inf(numel(xq), 1);
chunk = 1000;
for i0 = 1:chunk:numel(xq)
    i1 = min(numel(xq), i0 + chunk - 1);
    dx = xq(i0:i1) - xv.';
    dy = yq(i0:i1) - yv.';
    dmin(i0:i1) = sqrt(min(dx.^2 + dy.^2, [], 2));
end
end

%% ============================================================
function make_global_block_quicklook(T_block, BlockGeom, P)

fig = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1500 820]);
cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>

useMmap = exist('m_proj', 'file') == 2 && exist('m_coast', 'file') == 2;
basins = unique(T_block.basin_id, 'stable');
colors = lines(max(1, numel(basins)));

if useMmap
    m_proj('miller', 'lon', [-180 180], 'lat', [-60 65]);
    hold on;
    try
        m_coast('patch', [0.91 0.90 0.85], 'edgecolor', [0.30 0.30 0.30]);
    catch
        m_coast('color', [0.30 0.30 0.30]);
    end
    for i = 1:height(T_block)
        ib = find(basins == T_block.basin_id(i), 1, 'first');
        c = colors(ib, :);
        plot_lonlat_line_global_mmap(BlockGeom(i).outer_lon, BlockGeom(i).outer_lat, c, 1.0);
        plot_lonlat_line_global_mmap(BlockGeom(i).inner_lon, BlockGeom(i).inner_lat, c, 1.2);
        m_text(T_block.lon(i), T_block.lat(i), erase(T_block.block_id(i), "ADC_"), ...
            'FontSize', 6.5, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'Interpreter', 'none');
    end
    m_grid('box', 'on', 'tickdir', 'in', 'fontsize', 9, 'linestyle', ':');
else
    hold on;
    axis([-180 180 -60 65]);
    box on;
    grid on;
    for i = 1:height(T_block)
        ib = find(basins == T_block.basin_id(i), 1, 'first');
        c = colors(ib, :);
        plot_lonlat_line_global_plain(BlockGeom(i).outer_lon, BlockGeom(i).outer_lat, c, 1.0);
        plot_lonlat_line_global_plain(BlockGeom(i).inner_lon, BlockGeom(i).inner_lat, c, 1.2);
        text(T_block.lon(i), T_block.lat(i), erase(T_block.block_id(i), "ADC_"), ...
            'FontSize', 6.5, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'Interpreter', 'none');
    end
    xlabel('Longitude');
    ylabel('Latitude');
end

title({sprintf('ADCIRC model blocks from TC-exposed coastal partitions (%d blocks)', height(T_block)), ...
    'thin outline: outer ADCIRC domain; thick outline: inner coastal refinement domain'}, ...
    'FontWeight', 'bold', 'Interpreter', 'none');

legendHandles = gobjects(numel(basins), 1);
legendLabels = cell(numel(basins), 1);
for ib = 1:numel(basins)
    legendHandles(ib) = plot(nan, nan, '-', 'Color', colors(ib, :), 'LineWidth', 2);
    label = T_block.basin_label(find(T_block.basin_id == basins(ib), 1, 'first'));
    legendLabels{ib} = char(label);
end
legend(legendHandles, legendLabels, 'Location', 'southoutside', ...
    'Orientation', 'horizontal', 'NumColumns', 2, 'Interpreter', 'none');

exportgraphics(fig, fullfile(P.output_dir, 'global_tc_adcirc_blocks_quicklook.png'), ...
    'Resolution', P.figure_dpi);
end

%% ============================================================
function make_basin_block_quicklooks(T_block, BlockGeom, P)

basins = unique(T_block.basin_id, 'stable');
colors = lines(max(1, height(T_block)));
useMmap = exist('m_proj', 'file') == 2 && exist('m_coast', 'file') == 2;

Index = table(strings(0, 1), strings(0, 1), nan(0, 1), strings(0, 1), ...
    'VariableNames', {'basin_id', 'basin_label', 'block_count', 'plot_path'});

for ib = 1:numel(basins)
    I = find(T_block.basin_id == basins(ib));
    if isempty(I)
        continue;
    end

    [lonlim, latlim, lon0] = limits_for_blocks(T_block(I, :), BlockGeom(I), 3.0);
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1300 850]);
    cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>

    if useMmap
        m_proj('miller', 'lon', lonlim, 'lat', latlim);
        hold on;
        try
            m_coast('patch', [0.91 0.90 0.85], 'edgecolor', [0.30 0.30 0.30]);
        catch
            m_coast('color', [0.30 0.30 0.30]);
        end
        for k = 1:numel(I)
            c = colors(I(k), :);
            lonOuter = unwrap_lon_around(BlockGeom(I(k)).outer_lon, lon0);
            lonInner = unwrap_lon_around(BlockGeom(I(k)).inner_lon, lon0);
            draw_filled_lonlat_segments_mmap(lonOuter, BlockGeom(I(k)).outer_lat, ...
                lighten_color(c, 0.60), c, 1.1, 0.25);
            m_line(lonInner, BlockGeom(I(k)).inner_lat, 'Color', c, 'LineWidth', 1.3);
            m_text(unwrap_lon_around(T_block.lon(I(k)), lon0), T_block.lat(I(k)), ...
                T_block.block_id(I(k)), 'FontSize', 7, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'Interpreter', 'none');
        end
        m_grid('box', 'on', 'tickdir', 'in', 'fontsize', 9, 'linestyle', ':');
    else
        hold on;
        axis([lonlim latlim]);
        box on;
        grid on;
        for k = 1:numel(I)
            c = colors(I(k), :);
            lonOuter = unwrap_lon_around(BlockGeom(I(k)).outer_lon, lon0);
            lonInner = unwrap_lon_around(BlockGeom(I(k)).inner_lon, lon0);
            draw_filled_lonlat_segments_plain(lonOuter, BlockGeom(I(k)).outer_lat, ...
                lighten_color(c, 0.60), c, 1.1, 0.25);
            plot(lonInner, BlockGeom(I(k)).inner_lat, '-', 'Color', c, 'LineWidth', 1.3);
            text(unwrap_lon_around(T_block.lon(I(k)), lon0), T_block.lat(I(k)), ...
                T_block.block_id(I(k)), 'FontSize', 7, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'Interpreter', 'none');
        end
    end

    label = T_block.basin_label(I(1));
    title({sprintf('%s: %d ADCIRC blocks', label, numel(I)), ...
        'transparent fill: outer ADCIRC domain; bold line: inner coastal refinement domain'}, ...
        'FontWeight', 'bold', 'Interpreter', 'none');

    out = fullfile(P.quicklook_dir, sprintf('%s_adcirc_blocks.png', sanitize_filename(basins(ib))));
    exportgraphics(fig, out, 'Resolution', P.figure_dpi);
    Index = [Index; table(basins(ib), label, numel(I), string(out), ...
        'VariableNames', Index.Properties.VariableNames)]; %#ok<AGROW>
end

writetable(Index, fullfile(P.quicklook_dir, 'adcirc_block_quicklook_index.csv'));
end

%% ============================================================
function make_block_panel_quicklook(T_block, BlockGeom, P)

n = height(T_block);
if n == 0
    return;
end

ncol = 5;
nrow = ceil(n / ncol);
figW = 520 * ncol;
figH = 380 * nrow;
fig = figure('Color', 'w', 'Visible', 'off', 'Position', [40 40 figW figH]);
cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>

tl = tiledlayout(fig, nrow, ncol, 'Padding', 'compact', 'TileSpacing', 'compact');
useMmap = exist('m_proj', 'file') == 2 && exist('m_coast', 'file') == 2;

for i = 1:n
    ax = nexttile(tl);
    axes(ax); %#ok<LAXES>
    [lonlim, latlim, lon0] = limits_for_blocks(T_block(i, :), BlockGeom(i), 1.5);
    cOuter = [0.70 0.12 0.22];
    cInner = [0.08 0.08 0.08];
    cShelf = [0.05 0.55 0.20];
    cBay = [0.00 0.45 0.85];

    lonOuter = unwrap_lon_around(BlockGeom(i).outer_lon, lon0);
    lonInner = unwrap_lon_around(BlockGeom(i).inner_lon, lon0);
    lonShelf = unwrap_lon_around(BlockGeom(i).shelf_lon, lon0);
    lonBay = unwrap_lon_around(BlockGeom(i).bay_lon, lon0);

    if useMmap
        m_proj('miller', 'lon', lonlim, 'lat', latlim);
        hold on;
        try
            m_coast('patch', [0.91 0.90 0.85], 'edgecolor', [0.45 0.45 0.42]);
        catch
            m_coast('color', [0.45 0.45 0.42]);
        end
        if any(isfinite(lonOuter))
            draw_filled_lonlat_segments_mmap(lonOuter, BlockGeom(i).outer_lat, ...
                lighten_color(cOuter, 0.72), cOuter, 1.0, 0.25);
        end
        if any(isfinite(lonShelf))
            m_line(lonShelf, BlockGeom(i).shelf_lat, 'Color', cShelf, 'LineWidth', 1.0);
        end
        if any(isfinite(lonBay))
            m_line(lonBay, BlockGeom(i).bay_lat, 'Color', cBay, 'LineWidth', 1.2);
        end
        if any(isfinite(lonInner))
            m_line(lonInner, BlockGeom(i).inner_lat, 'Color', cInner, 'LineWidth', 1.1);
        end
        m_text(unwrap_lon_around(T_block.lon(i), lon0), T_block.lat(i), T_block.block_id(i), ...
            'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'Interpreter', 'none');
        m_grid('box', 'on', 'tickdir', 'in', 'fontsize', 5.5, 'linestyle', ':');
    else
        hold on;
        draw_filled_lonlat_segments_plain(lonOuter, BlockGeom(i).outer_lat, ...
            lighten_color(cOuter, 0.72), cOuter, 1.0, 0.25);
        if any(isfinite(lonShelf))
            plot(lonShelf, BlockGeom(i).shelf_lat, '-', 'Color', cShelf, 'LineWidth', 1.0);
        end
        if any(isfinite(lonBay))
            plot(lonBay, BlockGeom(i).bay_lat, '-', 'Color', cBay, 'LineWidth', 1.2);
        end
        plot(lonInner, BlockGeom(i).inner_lat, '-', 'Color', cInner, 'LineWidth', 1.1);
        text(unwrap_lon_around(T_block.lon(i), lon0), T_block.lat(i), T_block.block_id(i), ...
            'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'Interpreter', 'none');
        axis([lonlim latlim]);
        axis equal;
        box on;
        grid on;
        set(gca, 'FontSize', 5.5);
    end

    title(sprintf('%s | shelf %.0f | bay %.0f km^2', ...
        T_block.block_id(i), T_block.shelf_add_area_km2(i), T_block.bay_add_area_km2(i)), ...
        'FontSize', 7, 'Interpreter', 'none');
end

title(tl, {'Per-block ADCIRC geometry diagnostics', ...
    'red: outer domain; black: inner refinement; green: shelf gap-fill; blue: ring-bay-add'}, ...
    'FontWeight', 'bold', 'Interpreter', 'none');
out = fullfile(P.quicklook_dir, 'adcirc_block_panels_all.png');
exportgraphics(fig, out, 'Resolution', P.figure_dpi);
end

%% ============================================================
function write_block_range_summary(path, T_block, P)

outer_geojson = repmat(string(fullfile(P.output_dir, 'global_tc_adcirc_block_outer_domains.geojson')), height(T_block), 1);
inner_geojson = repmat(string(fullfile(P.output_dir, 'global_tc_adcirc_block_inner_domains.geojson')), height(T_block), 1);
outer_ring_csv = strings(height(T_block), 1);
inner_ring_csv = strings(height(T_block), 1);
for i = 1:height(T_block)
    outer_ring_csv(i) = string(fullfile(P.block_boundary_dir, sprintf('%s_outer_ring.csv', T_block.block_id(i))));
    inner_ring_csv(i) = string(fullfile(P.block_boundary_dir, sprintf('%s_inner_ring.csv', T_block.block_id(i))));
end

Range = table(T_block.block_id, T_block.basin_id, T_block.basin_label, ...
    T_block.lon, T_block.lat, T_block.domain_count, T_block.grid_cell_count, ...
    T_block.block_span_km, T_block.inner_area_km2, T_block.outer_area_km2, ...
    T_block.inner_along_half_km, T_block.inner_cross_half_km, ...
    T_block.outer_along_half_km, T_block.outer_cross_half_km, ...
    T_block.shelf_width_km, T_block.outer_buffer_km, ...
    T_block.shelf_add_area_km2, T_block.shelf_add_cell_count, ...
    T_block.shelf_add_component_count, ...
    T_block.bay_add_area_km2, T_block.bay_add_cell_count, ...
    T_block.bay_add_core_cell_count, T_block.bay_add_component_count, ...
    T_block.auto_alongshore_radius_deg, T_block.auto_offshore_deg, ...
    T_block.auto_inner_radius, T_block.open_ocean_area_risk, ...
    outer_geojson, inner_geojson, outer_ring_csv, inner_ring_csv, ...
    'VariableNames', {'block_id', 'basin_id', 'basin_label', ...
    'center_lon', 'center_lat', 'domain_count', 'grid_cell_count', ...
    'coastal_cell_span_km', 'inner_refinement_area_km2', 'outer_adcirc_area_km2', ...
    'inner_along_half_km', 'inner_cross_half_km', ...
    'outer_along_half_km', 'outer_cross_half_km', ...
    'shelf_width_km', 'outer_buffer_km', ...
    'shelf_add_area_km2', 'shelf_add_cell_count', ...
    'shelf_add_component_count', ...
    'bay_add_area_km2', 'bay_add_cell_count', ...
    'bay_add_core_cell_count', 'bay_add_component_count', ...
    'auto_alongshore_radius_deg', 'auto_offshore_deg', ...
    'auto_inner_radius_deg', 'open_ocean_area_risk', ...
    'outer_geojson', 'inner_geojson', 'outer_ring_csv', 'inner_ring_csv'});
writetable(Range, path);
end

%% ============================================================
function write_block_ring_csvs(T_block, BlockGeom, P)

for i = 1:height(T_block)
    outerLon = unwrap_lon_around(BlockGeom(i).outer_lon(:), T_block.lon(i));
    innerLon = unwrap_lon_around(BlockGeom(i).inner_lon(:), T_block.lon(i));
    outer = table(outerLon(:), BlockGeom(i).outer_lat(:), ...
        'VariableNames', {'lon', 'lat'});
    inner = table(innerLon(:), BlockGeom(i).inner_lat(:), ...
        'VariableNames', {'lon', 'lat'});
    writetable(outer, fullfile(P.block_boundary_dir, sprintf('%s_outer_ring.csv', T_block.block_id(i))));
    writetable(inner, fullfile(P.block_boundary_dir, sprintf('%s_inner_ring.csv', T_block.block_id(i))));
end
end

function write_block_geojson(path, T_block, BlockGeom, whichRing)

fid = fopen(path, 'w');
if fid < 0
    error('Cannot open GeoJSON for writing: %s', path);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '{\n');
fprintf(fid, '  "type": "FeatureCollection",\n');
fprintf(fid, '  "name": "%s",\n', erase(string(whichRing), '"'));
fprintf(fid, '  "crs": {"type": "name", "properties": {"name": "EPSG:4326"}},\n');
fprintf(fid, '  "features": [\n');
for i = 1:height(T_block)
    if i > 1
        fprintf(fid, ',\n');
    end
    if strcmpi(whichRing, 'outer')
        lon = BlockGeom(i).outer_lon;
        lat = BlockGeom(i).outer_lat;
    else
        lon = BlockGeom(i).inner_lon;
        lat = BlockGeom(i).inner_lat;
    end
    fprintf(fid, '    {\n');
    fprintf(fid, '      "type": "Feature",\n');
    fprintf(fid, '      "properties": {');
    fprintf(fid, '"block_id": "%s", ', escape_json_string(T_block.block_id(i)));
    fprintf(fid, '"basin_id": "%s", ', escape_json_string(T_block.basin_id(i)));
    fprintf(fid, '"basin_label": "%s", ', escape_json_string(T_block.basin_label(i)));
    fprintf(fid, '"domain_count": %d, ', T_block.domain_count(i));
    fprintf(fid, '"grid_cell_count": %d, ', T_block.grid_cell_count(i));
    fprintf(fid, '"outer_area_km2": %.6f, ', T_block.outer_area_km2(i));
    fprintf(fid, '"shelf_width_km": %.6f, ', T_block.shelf_width_km(i));
    fprintf(fid, '"outer_buffer_km": %.6f, ', T_block.outer_buffer_km(i));
    fprintf(fid, '"shelf_add_area_km2": %.6f, ', T_block.shelf_add_area_km2(i));
    fprintf(fid, '"shelf_add_cell_count": %d, ', T_block.shelf_add_cell_count(i));
    fprintf(fid, '"shelf_add_component_count": %d, ', T_block.shelf_add_component_count(i));
    fprintf(fid, '"bay_add_area_km2": %.6f, ', T_block.bay_add_area_km2(i));
    fprintf(fid, '"bay_add_cell_count": %d, ', T_block.bay_add_cell_count(i));
    fprintf(fid, '"bay_add_core_cell_count": %d, ', T_block.bay_add_core_cell_count(i));
    fprintf(fid, '"bay_add_component_count": %d, ', T_block.bay_add_component_count(i));
    fprintf(fid, '"open_ocean_area_risk": %d', T_block.open_ocean_area_risk(i));
    fprintf(fid, '},\n');
    fprintf(fid, '      "geometry": ');
    write_geojson_geometry(fid, lon(:), lat(:), T_block.lon(i));
    fprintf(fid, '\n');
    fprintf(fid, '    }');
end
fprintf(fid, '\n  ]\n');
fprintf(fid, '}\n');
end

%% ============================================================
function write_geojson_geometry(fid, lon, lat, centerLon)

if nargin < 4 || ~isfinite(centerLon)
    centerLon = weighted_circular_mean_lon(lon(isfinite(lon)), ones(nnz(isfinite(lon)), 1));
end

segments = geojson_dateline_safe_segments(lon, lat, centerLon);
if numel(segments) <= 1
    fprintf(fid, '{"type": "Polygon", "coordinates": [');
    if isempty(segments)
        write_geojson_ring(fid, [NaN NaN]);
    else
        write_geojson_ring(fid, segments{1});
    end
    fprintf(fid, ']}');
    return;
end

fprintf(fid, '{"type": "MultiPolygon", "coordinates": [');
for i = 1:numel(segments)
    if i > 1
        fprintf(fid, ', ');
    end
    fprintf(fid, '[');
    write_geojson_ring(fid, segments{i});
    fprintf(fid, ']');
end
fprintf(fid, ']}');
end

%% ============================================================
function segments = geojson_dateline_safe_segments(lon, lat, centerLon)

baseSegments = nan_split_segments(lon, lat);
segments = {};
for i = 1:numel(baseSegments)
    seg = baseSegments{i};
    if size(seg, 1) < 3
        continue;
    end

    lonUn = unwrap_lon_around(seg(:, 1), centerLon);
    hasJump = any(abs(diff(seg(:, 1))) > 180) || (max(lonUn) > 180) || (min(lonUn) < -180);
    if ~hasJump
        segments{end + 1} = seg; %#ok<AGROW>
        continue;
    end

    parts = split_geojson_segment_at_dateline(seg, centerLon);
    for ip = 1:numel(parts)
        if size(parts{ip}, 1) >= 4
            segments{end + 1} = parts{ip}; %#ok<AGROW>
        end
    end
end
end

%% ============================================================
function parts = split_geojson_segment_at_dateline(seg, centerLon)

parts = {};
lonUn = unwrap_lon_around(seg(:, 1), centerLon);
lat = double(seg(:, 2));
valid = isfinite(lonUn) & isfinite(lat);
lonUn = lonUn(valid);
lat = lat(valid);
if numel(lonUn) < 3
    return;
end

try
    p = polyshape(lonUn, lat, 'Simplify', true);
catch
    p = polyshape();
end
if area(p) <= 0
    lonWrapped = wrapTo180_local(lonUn);
    if any([lonWrapped(1), lat(1)] ~= [lonWrapped(end), lat(end)])
        lonWrapped(end + 1, 1) = lonWrapped(1); %#ok<AGROW>
        lat(end + 1, 1) = lat(1); %#ok<AGROW>
    end
    parts = {[lonWrapped(:), lat(:)]};
    return;
end

bands = [-540 -180; -180 180; 180 540];
for ib = 1:size(bands, 1)
    b = bands(ib, :);
    clipBox = polyshape([b(1); b(2); b(2); b(1); b(1)], ...
        [-90; -90; 90; 90; -90], 'Simplify', true);
    try
        pc = intersect(p, clipBox);
    catch
        pc = polyshape();
    end
    if area(pc) <= 0
        continue;
    end

    shift = 0;
    if b(1) >= 180
        shift = -360;
    elseif b(2) <= -180
        shift = 360;
    end

    rr = regions(pc);
    for ir = 1:numel(rr)
        [x, y] = boundary(rr(ir));
        sub = nan_split_segments(x + shift, y);
        for is = 1:numel(sub)
            ring = sub{is};
            ring(:, 1) = max(-180, min(180, ring(:, 1)));
            if size(ring, 1) >= 4 && ~any(abs(diff(ring(:, 1))) > 180)
                parts{end + 1} = ring; %#ok<AGROW>
            end
        end
    end
end

if isempty(parts)
    lonWrapped = wrapTo180_local(lonUn);
    parts = {[lonWrapped(:), lat(:)]};
end
end

%% ============================================================
function segments = nan_split_segments(lon, lat)

lon = double(lon(:));
lat = double(lat(:));
breaks = isnan(lon) | isnan(lat);
segments = {};
startIdx = 1;
for i = 1:numel(lon)+1
    if i > numel(lon) || breaks(i)
        if i > startIdx
            seg = [lon(startIdx:i-1), lat(startIdx:i-1)];
            seg = seg(all(isfinite(seg), 2), :);
            if size(seg, 1) >= 3
                if any(seg(1, :) ~= seg(end, :))
                    seg(end + 1, :) = seg(1, :); %#ok<AGROW>
                end
                segments{end + 1} = seg; %#ok<AGROW>
            end
        end
        startIdx = i + 1;
    end
end
end

%% ============================================================
function write_oceanmesh2d_runner(path, P)

fid = fopen(path, 'w');
if fid < 0
    error('Cannot write runner: %s', path);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '%% P2_run_global_tc_adcirc_block_meshes.m\n');
fprintf(fid, '%% Optional OceanMesh2D runner generated by P1_build_adcirc_model_blocks_from_partitions.m\n');
fprintf(fid, 'clearvars; clc; close all;\n\n');
fprintf(fid, 'SCRIPT_DIR = fileparts(mfilename(''fullpath''));\n');
fprintf(fid, 'BLOCK_CSV = fullfile(SCRIPT_DIR, ''global_tc_adcirc_blocks.csv'');\n');
fprintf(fid, 'BOUNDARY_DIR = fullfile(SCRIPT_DIR, ''block_boundaries'');\n');
fprintf(fid, 'OUT_MESH_DIR = fullfile(SCRIPT_DIR, ''adcirc_meshes'');\n');
fprintf(fid, 'OUT_FIG_DIR = fullfile(SCRIPT_DIR, ''adcirc_mesh_figures'');\n');
fprintf(fid, 'if exist(OUT_MESH_DIR, ''dir'') ~= 7, mkdir(OUT_MESH_DIR); end\n\n');
fprintf(fid, 'if exist(OUT_FIG_DIR, ''dir'') ~= 7, mkdir(OUT_FIG_DIR); end\n\n');
fprintf(fid, 'repo_dir = find_parent_dir_named(SCRIPT_DIR, ''OceanMesh2D-Projection'');\n');
fprintf(fid, 'addpath(genpath(fullfile(repo_dir, ''utilities'')));\n');
fprintf(fid, 'addpath(genpath(fullfile(repo_dir, ''datasets'')));\n');
fprintf(fid, 'addpath(genpath(fullfile(repo_dir, ''m_map'')));\n');
fprintf(fid, 'addpath(genpath(fullfile(repo_dir, ''Examples'', ''Global_autofunction'')));\n\n');
fprintf(fid, 'assert(exist(''geodata'', ''class'') == 8 || exist(''geodata'', ''file'') == 2, ''OceanMesh2D geodata is not on path.'');\n');
fprintf(fid, 'assert(exist(''meshgen'', ''class'') == 8 || exist(''meshgen'', ''file'') == 2, ''OceanMesh2D meshgen is not on path.'');\n\n');
fprintf(fid, 'T = readtable(BLOCK_CSV, ''TextType'', ''string'');\n');
fprintf(fid, 'selected_blocks = strings(0,1); %% keep empty to process all blocks, or set e.g. ["ADC_WNP_01"]\n');
fprintf(fid, 'if ~isempty(selected_blocks), T = T(ismember(T.block_id, selected_blocks), :); end\n\n');
fprintf(fid, 'coastline_coarse = ''land_polygons'';\n');
fprintf(fid, 'coastline_fine = ''land_polygons'';\n');
fprintf(fid, 'dem_coarse = ''GEBCO_2025_sub_ice.nc'';\n');
fprintf(fid, 'fine_dem_dir = fullfile(repo_dir, ''datasets'', ''fine_dem'');\n\n');
fprintf(fid, 'min_el_out = 1e3; max_el_out = 20e3; max_el_ns_out = 5e3;\n');
fprintf(fid, 'min_el_in = 200; max_el_in = 5e3; max_el_ns_in = 500;\n');
fprintf(fid, 'wl = 30; g = 0.25;\n\n');
fprintf(fid, 'for i = 1:height(T)\n');
fprintf(fid, '    block_id = T.block_id(i);\n');
fprintf(fid, '    fprintf(''\\n===== [%%d/%%d] %%s =====\\n'', i, height(T), block_id);\n');
fprintf(fid, '    outer = readtable(fullfile(BOUNDARY_DIR, block_id + "_outer_ring.csv"));\n');
fprintf(fid, '    inner = readtable(fullfile(BOUNDARY_DIR, block_id + "_inner_ring.csv"));\n');
fprintf(fid, '    bbox_outer = [outer.lon, outer.lat];\n');
fprintf(fid, '    bbox_inner = [inner.lon, inner.lat];\n');
fprintf(fid, '    [dem_fine, dem_info] = auto_select_fine_dem_by_point(T.lon(i), T.lat(i), fine_dem_dir, dem_coarse); %%#ok<NASGU>\n');
fprintf(fid, '    gdat_out = geodata(''shp'', coastline_coarse, ''dem'', dem_coarse, ''bbox'', bbox_outer, ''h0'', min_el_out, ''window'', 10);\n');
fprintf(fid, '    fh_out = edgefx(''geodata'', gdat_out, ''fs'', 3, ''max_el_ns'', max_el_ns_out, ''max_el'', max_el_out, ''dt'', 5, ''g'', g, ''wl'', wl);\n');
fprintf(fid, '    gdat_in = geodata(''shp'', coastline_fine, ''dem'', dem_fine, ''bbox'', bbox_inner, ''h0'', min_el_in, ''window'', 10);\n');
fprintf(fid, '    fh_in = edgefx(''geodata'', gdat_in, ''fs'', 6, ''max_el_ns'', max_el_ns_in, ''max_el'', max_el_in, ''dt'', 5, ''g'', g);\n');
fprintf(fid, '    rng(1.23456789);\n');
fprintf(fid, '    mshopts = meshgen(''ef'', {fh_out, fh_in}, ''bou'', {gdat_out, gdat_in}, ''plot_on'', 1, ''proj'', ''lam'');\n');
fprintf(fid, '    mshopts = mshopts.build;\n');
fprintf(fid, '    m = mshopts.grd;\n');
fprintf(fid, '    m = interp(m, {gdat_out, gdat_in}, ''mindepth'', 1);\n');
fprintf(fid, '    m = lim_bathy_slope(m, 0.1, 0);\n');
fprintf(fid, '    m = bound_courant_number(m, 1.2, 0.5, 0, 10);\n');
fprintf(fid, '    if exist(''make_bc_from_capsule_boundary'', ''file'') == 2\n');
fprintf(fid, '        m = make_bc_from_capsule_boundary(m, gdat_out, ''shore_tol'', 0.1, ''min_open_edges'', 10, ''bridge_gap_edges'', 4, ''add_inner_islands'', true, ''inner_ibtype'', 21, ''plot_check'', true);\n');
fprintf(fid, '    end\n');
fprintf(fid, '    outBase = fullfile(OUT_MESH_DIR, "ADCIRC_Block_Mesh_" + block_id);\n');
fprintf(fid, '    write(m, char(outBase), ''14'');\n');
fprintf(fid, '    save(char(outBase + ".mat"), ''m'');\n');
fprintf(fid, '    fig = figure(''Color'', ''w'', ''Visible'', ''off'');\n');
fprintf(fid, '    plot(m, ''type'', ''resologmesh'');\n');
fprintf(fid, '    title("ADCIRC block mesh " + block_id, ''Interpreter'', ''none'');\n');
fprintf(fid, '    exportgraphics(fig, char(fullfile(OUT_FIG_DIR, "ADCIRC_Block_Mesh_" + block_id + ".png")), ''Resolution'', 180);\n');
fprintf(fid, '    close(fig);\n');
fprintf(fid, 'end\n\n');
fprintf(fid, 'function root = find_parent_dir_named(startDir, targetName)\n');
fprintf(fid, 'root = ''''; d = char(startDir); targetName = char(targetName);\n');
fprintf(fid, 'while true\n');
fprintf(fid, '    [parent, name] = fileparts(d);\n');
fprintf(fid, '    if strcmpi(name, targetName), root = d; return; end\n');
fprintf(fid, '    if isempty(parent) || strcmp(parent, d), return; end\n');
fprintf(fid, '    d = parent;\n');
fprintf(fid, 'end\n');
fprintf(fid, 'end\n');
end

%% ============================================================
function [lonlim, latlim, lon0] = limits_for_blocks(Tb, G, pad)

lon0 = weighted_circular_mean_lon(Tb.lon, ones(height(Tb), 1));
lonAll = [];
latAll = [];
for i = 1:numel(G)
    lonAll = [lonAll; unwrap_lon_around(G(i).outer_lon, lon0)]; %#ok<AGROW>
    latAll = [latAll; G(i).outer_lat]; %#ok<AGROW>
end
lonlim = [min(lonAll) - pad, max(lonAll) + pad];
latlim = [max(-80, min(latAll) - pad), min(80, max(latAll) + pad)];
end

%% ============================================================
function T = ensure_domain_basins(T)

if ismember('basin_id', T.Properties.VariableNames) && ismember('basin_label', T.Properties.VariableNames)
    T.basin_id = string(T.basin_id);
    T.basin_label = string(T.basin_label);
    return;
end

error('Domain table lacks basin_id/basin_label. Run partition_global_tc_coastal_model_domains.m first.');
end

%% ============================================================
function DomainPolys = read_domain_boundary_geojson(path)

assert(exist(path, 'file') == 2, 'Missing domain boundary GeoJSON: %s', path);
G = jsondecode(fileread(path));
features = G.features;
DomainPolys = struct('id', {}, 'lon', {}, 'lat', {});
for i = 1:numel(features)
    props = features(i).properties;
    if isfield(props, 'model_domain_id')
        id = string(props.model_domain_id);
    else
        id = "GTC_UNKNOWN_" + i;
    end
    geom = features(i).geometry;
    [lon, lat] = geojson_coordinates_to_nan_vectors(geom);
    DomainPolys(end + 1).id = id; %#ok<AGROW>
    DomainPolys(end).lon = lon(:);
    DomainPolys(end).lat = lat(:);
end
fprintf('  -> loaded partition boundary polygons from GeoJSON: %d\n', numel(DomainPolys));
end

%% ============================================================
function [lon, lat] = geojson_coordinates_to_nan_vectors(geom)

lon = [];
lat = [];
typ = string(geom.type);
coords = geom.coordinates;
if typ == "Polygon"
    rings = polygon_coords_to_rings(coords);
elseif typ == "MultiPolygon"
    rings = {};
    for i = 1:numel(coords)
        rr = polygon_coords_to_rings(coords{i});
        rings = [rings, rr]; %#ok<AGROW>
    end
else
    rings = {};
end

for i = 1:numel(rings)
    ring = rings{i};
    if size(ring, 1) < 3
        continue;
    end
    if ~isempty(lon)
        lon(end + 1, 1) = NaN; %#ok<AGROW>
        lat(end + 1, 1) = NaN; %#ok<AGROW>
    end
    lon = [lon; ring(:, 1)]; %#ok<AGROW>
    lat = [lat; ring(:, 2)]; %#ok<AGROW>
end
end

%% ============================================================
function rings = polygon_coords_to_rings(coords)

rings = {};
if isnumeric(coords)
    sz = size(coords);
    if numel(sz) == 3 && sz(end) == 2
        for i = 1:sz(1)
            ring = squeeze(coords(i, :, :));
            rings{end + 1} = ring; %#ok<AGROW>
        end
    elseif size(coords, 2) == 2
        rings{1} = coords;
    end
elseif iscell(coords)
    for i = 1:numel(coords)
        if isnumeric(coords{i})
            c = coords{i};
            if ndims(c) == 3 && size(c, 3) == 2
                c = squeeze(c(1, :, :));
            end
            rings{end + 1} = c; %#ok<AGROW>
        end
    end
end
end

%% ============================================================
function T = read_table_as_strings(path)

try
    T = readtable(path, 'TextType', 'string');
catch
    T = readtable(path);
end
for i = 1:numel(T.Properties.VariableNames)
    name = T.Properties.VariableNames{i};
    if isstring(T.(name)) || iscellstr(T.(name)) || ischar(T.(name))
        T.(name) = string(T.(name));
    end
end
if ismember('model_domain_id', T.Properties.VariableNames)
    T.model_domain_id = string(T.model_domain_id);
end
end

%% ============================================================
function [x, y, lon0, lat0] = lonlat_to_local_km(lon, lat, lon0, lat0)

if nargin < 3
    lon0 = weighted_circular_mean_lon(lon, ones(numel(lon), 1));
    lat0 = mean(double(lat), 'omitnan');
end
lon = double(lon(:));
lat = double(lat(:));
x = wrapTo180_local(lon - lon0) .* 111.32 .* max(cosd(lat0), 0.15);
y = (lat - lat0) .* 111.32;
end

%% ============================================================
function [lon, lat] = local_km_to_lonlat(x, y, lon0, lat0)

lon = lon0 + double(x(:)) ./ (111.32 .* max(cosd(lat0), 0.15));
lat = lat0 + double(y(:)) ./ 111.32;
lon = wrapTo180_local(lon);
end

%% ============================================================
function axisU = principal_axis_2d(x, y, w)

x = double(x(:));
y = double(y(:));
w = double(w(:));
valid = isfinite(x) & isfinite(y) & isfinite(w) & w > 0;
if nnz(valid) < 2
    axisU = [1; 0];
    return;
end
x = x(valid);
y = y(valid);
w = w(valid);
w = w ./ sum(w);
x0 = sum(w .* x);
y0 = sum(w .* y);
X = [x - x0, y - y0];
C = (X .* w).' * X;
[V, D] = eig(C);
[~, imax] = max(diag(D));
axisU = V(:, imax);
if axisU(1) < 0
    axisU = -axisU;
end
end

%% ============================================================
function r = along_range_from_indices(s, idx)

r = max(s(idx)) - min(s(idx));
end

%% ============================================================
function s = span_from_indices(x, y, idx)

s = pairwise_span_km(x(idx), y(idx));
end

%% ============================================================
function s = cell_span_for_domain_ids(T_cell, domain_ids)

C = T_cell(ismember(T_cell.model_domain_id, string(domain_ids(:))), :);
if isempty(C)
    s = 0;
    return;
end
lon0 = weighted_circular_mean_lon(C.lon, ones(height(C), 1));
lat0 = mean(double(C.lat), 'omitnan');
[x, y] = lonlat_to_local_km(C.lon, C.lat, lon0, lat0);
s = pairwise_span_km(x, y);
end

%% ============================================================
function s = pairwise_span_km(x, y)

x = double(x(:));
y = double(y(:));
if numel(x) <= 1
    s = 0;
    return;
end
s = 0;
for i = 1:numel(x)
    d = hypot(x(i) - x, y(i) - y);
    s = max(s, max(d, [], 'omitnan'));
end
end

%% ============================================================
function w = table_weights(T)

if ismember('lowlying_area_km2', T.Properties.VariableNames)
    w = double(T.lowlying_area_km2);
elseif ismember('grid_cell_count', T.Properties.VariableNames)
    w = double(T.grid_cell_count);
else
    w = ones(height(T), 1);
end
w(~isfinite(w) | w <= 0) = 1;
end

%% ============================================================
function v = table_sum_if_present(T, name)

if ismember(name, T.Properties.VariableNames)
    v = sum(double(T.(name)), 'omitnan');
else
    v = NaN;
end
end

%% ============================================================
function v = table_max_if_present(T, name)

if ismember(name, T.Properties.VariableNames)
    v = max(double(T.(name)), [], 'omitnan');
else
    v = NaN;
end
end

%% ============================================================
function areaKm2 = polygon_area_local_km2(lon, lat, lon0, lat0)

[x, y] = lonlat_to_local_km(lon, lat, lon0, lat0);
areaKm2 = abs(polyarea(x, y));
end

%% ============================================================
function lon0 = weighted_circular_mean_lon(lon, weights)

lon = double(lon(:));
weights = double(weights(:));
valid = isfinite(lon) & isfinite(weights) & weights > 0;
if ~any(valid)
    lon0 = mean(wrapTo180_local(lon(isfinite(lon))), 'omitnan');
    if ~isfinite(lon0)
        lon0 = 0;
    end
    return;
end
ang = deg2rad(wrapTo180_local(lon(valid)));
w = weights(valid);
s = sum(w .* sin(ang), 'omitnan');
c = sum(w .* cos(ang), 'omitnan');
lon0 = wrapTo180_local(rad2deg(atan2(s, c)));
end

%% ============================================================
function m = weighted_mean(x, w)

x = double(x(:));
w = double(w(:));
valid = isfinite(x) & isfinite(w) & w > 0;
if ~any(valid)
    m = mean(x, 'omitnan');
else
    m = sum(x(valid) .* w(valid), 'omitnan') / sum(w(valid), 'omitnan');
end
end

%% ============================================================
function d = local_distance_km(lon1, lat1, lon2, lat2)

R = 6371.0;
lon1 = deg2rad(double(lon1));
lat1 = deg2rad(double(lat1));
lon2 = deg2rad(double(lon2));
lat2 = deg2rad(double(lat2));
dlon = lon2 - lon1;
dlat = lat2 - lat1;
a = sin(dlat / 2).^2 + cos(lat1) .* cos(lat2) .* sin(dlon / 2).^2;
d = 2 * R * atan2(sqrt(a), sqrt(max(0, 1 - a)));
end

%% ============================================================
function short = basin_short_name(basin_id)

s = erase(string(basin_id), "BASIN_");
short = char(s);
end

%% ============================================================
function c = lighten_color(c, amount)

c = double(c(:)).';
c = c + max(0, min(1, amount)) .* (1 - c);
c = max(0, min(1, c));
end

%% ============================================================
function plot_lonlat_line_global_mmap(lon, lat, color, lineWidth)

[lonSeg, latSeg] = split_dateline_jumps(lon, lat);
m_line(lonSeg, latSeg, 'Color', color, 'LineWidth', lineWidth);
end

%% ============================================================
function plot_lonlat_line_global_plain(lon, lat, color, lineWidth)

[lonSeg, latSeg] = split_dateline_jumps(lon, lat);
plot(lonSeg, latSeg, '-', 'Color', color, 'LineWidth', lineWidth);
end

%% ============================================================
function draw_filled_lonlat_segments_mmap(lon, lat, faceColor, edgeColor, lineWidth, faceAlpha)

segments = nan_split_segments(lon, lat);
for i = 1:numel(segments)
    seg = segments{i};
    if size(seg, 1) < 4
        continue;
    end
    m_patch(seg(:, 1), seg(:, 2), faceColor, ...
        'EdgeColor', edgeColor, 'LineWidth', lineWidth, 'FaceAlpha', faceAlpha);
end
end

%% ============================================================
function draw_filled_lonlat_segments_plain(lon, lat, faceColor, edgeColor, lineWidth, faceAlpha)

segments = nan_split_segments(lon, lat);
for i = 1:numel(segments)
    seg = segments{i};
    if size(seg, 1) < 4
        continue;
    end
    patch(seg(:, 1), seg(:, 2), faceColor, ...
        'EdgeColor', edgeColor, 'LineWidth', lineWidth, 'FaceAlpha', faceAlpha);
end
end

%% ============================================================
function [lonOut, latOut] = split_dateline_jumps(lon, lat)

lon = double(lon(:));
lat = double(lat(:));
lonOut = lon;
latOut = lat;
if numel(lon) < 2
    return;
end
jump = abs(diff(lon)) > 180;
if ~any(jump)
    return;
end

lonOut = [];
latOut = [];
for i = 1:(numel(lon) - 1)
    lonOut(end + 1, 1) = lon(i); %#ok<AGROW>
    latOut(end + 1, 1) = lat(i); %#ok<AGROW>
    if jump(i)
        lonOut(end + 1, 1) = NaN; %#ok<AGROW>
        latOut(end + 1, 1) = NaN; %#ok<AGROW>
    end
end
lonOut(end + 1, 1) = lon(end);
latOut(end + 1, 1) = lat(end);
end

%% ============================================================
function write_geojson_ring(fid, coords)

fprintf(fid, '[');
for i = 1:size(coords, 1)
    if i > 1
        fprintf(fid, ', ');
    end
    fprintf(fid, '[%.8f, %.8f]', coords(i, 1), coords(i, 2));
end
fprintf(fid, ']');
end

%% ============================================================
function s = escape_json_string(s)

s = char(string(s));
s = strrep(s, '\', '\\');
s = strrep(s, '"', '\"');
end

%% ============================================================
function s = sanitize_filename(s)

s = char(string(s));
s = regexprep(s, '[^\w\-]+', '_');
end

%% ============================================================
function lon = unwrap_lon_around(lon, center_lon)

lon = center_lon + wrapTo180_local(double(lon) - center_lon);
end

%% ============================================================
function lon = wrapTo180_local(lon)

lon = mod(double(lon) + 180, 360) - 180;
lon(lon == -180) = 180;
end

%% ============================================================
function root = find_parent_dir_named(startDir, targetName)

root = '';
if isempty(startDir)
    return;
end
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

%% ============================================================
function ensure_dir(d)

if exist(d, 'dir') ~= 7
    mkdir(d);
end
end

%% ============================================================
function cleanup_previous_block_outputs(P)

delete_matching_files(P.block_boundary_dir, {'ADC_*_outer_ring.csv', 'ADC_*_inner_ring.csv'});
delete_matching_files(P.quicklook_dir, {'BASIN_*_adcirc_blocks.png', 'adcirc_block_quicklook_index.csv', ...
    'adcirc_block_panels_all.png'});
delete_matching_files(P.output_dir, { ...
    'global_tc_adcirc_blocks.csv', ...
    'global_tc_adcirc_block_members.csv', ...
    'global_tc_adcirc_block_ranges.csv', ...
    'global_tc_adcirc_block_outer_domains.geojson', ...
    'global_tc_adcirc_block_inner_domains.geojson', ...
    'global_tc_adcirc_blocks_quicklook.png', ...
    'global_tc_adcirc_model_blocks.mat', ...
    'P2_run_global_tc_adcirc_block_meshes.m', ...
    '00_run_metadata.txt'});
end

%% ============================================================
function delete_matching_files(dirPath, patterns)

if exist(dirPath, 'dir') ~= 7
    return;
end
for ip = 1:numel(patterns)
    files = dir(fullfile(dirPath, patterns{ip}));
    for i = 1:numel(files)
        if files(i).isdir
            continue;
        end
        p = fullfile(files(i).folder, files(i).name);
        if exist(p, 'file') ~= 2
            continue;
        end
        try
            warnState = warning('off', 'all');
            delete(p);
            warning(warnState);
        catch
            warning(warnState);
            warning('Could not delete old generated file: %s', p);
        end
    end
end
end

%% ============================================================
function add_existing_genpath(path_list)

for ip = 1:numel(path_list)
    p = path_list{ip};
    if exist(p, 'dir') == 7
        addpath(genpath(p));
    end
end
end

%% ============================================================
function write_run_metadata(path, P, T_block)

fid = fopen(path, 'w');
if fid < 0
    warning('Cannot write metadata: %s', path);
    return;
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'P1_build_adcirc_model_blocks_from_partitions.m\n');
fprintf(fid, 'Created by MATLAB on %s\n\n', datestr(now));
fprintf(fid, 'Purpose: prepare local ADCIRC/OceanMesh2D blocks from global TC coastal partitions.\n\n');
fprintf(fid, 'Domain CSV: %s\n', P.domain_csv);
fprintf(fid, 'Cell CSV: %s\n', P.cell_csv);
fprintf(fid, 'Output dir: %s\n', P.output_dir);
fprintf(fid, 'Block count: %d\n', height(T_block));
fprintf(fid, 'Max block alongshore km: %.3f\n', P.max_block_alongshore_km);
fprintf(fid, 'Max block span km: %.3f\n', P.max_block_span_km);
fprintf(fid, 'Max domains per block: %d\n', P.max_domains_per_block);
fprintf(fid, 'Max inner refinement area km2: %.3f\n', P.max_inner_refinement_area_km2);
fprintf(fid, 'Second-pass merge enabled: %d\n', logical(P.merge_overlapping_blocks));
fprintf(fid, 'Merge outer overlap min fraction: %.3f\n', P.merge_outer_overlap_min_fraction);
fprintf(fid, 'Merge outer union min fraction: %.3f\n', P.merge_outer_union_min_fraction);
fprintf(fid, 'Merge outer gap km: %.3f\n', P.merge_outer_gap_km);
fprintf(fid, 'Merge candidate max outer area km2: %.3f\n', P.merge_candidate_max_outer_area_km2);
fprintf(fid, 'Merge candidate max span km: %.3f\n', P.merge_candidate_max_span_km);
fprintf(fid, 'Coastal connectivity split gap km: %.3f\n', P.coastal_connect_gap_km);
fprintf(fid, 'Coastal connectivity merge gap km: %.3f\n', P.merge_coastal_gap_km);
fprintf(fid, 'Small-block absorption enabled: %d\n', logical(P.absorb_small_blocks));
fprintf(fid, 'Small-block max cells: %d\n', P.small_block_max_cells);
fprintf(fid, 'Small-block cluster max cells: %d\n', P.small_block_cluster_max_cells);
fprintf(fid, 'Small-block merge-to-major gap km: %.3f\n', P.small_merge_to_major_gap_km);
fprintf(fid, 'Small-block merge-to-small gap km: %.3f\n', P.small_merge_small_gap_km);
fprintf(fid, 'Small-block candidate max outer area km2: %.3f\n', P.small_merge_candidate_max_outer_area_km2);
fprintf(fid, 'Capsule-style topology enabled: %d\n', logical(P.topology_capsule_style));
fprintf(fid, 'Capsule inland km: %.3f\n', P.topology_capsule_inland_km);
fprintf(fid, 'Capsule bay-close km: %.3f\n', P.topology_capsule_bay_close_km);
fprintf(fid, 'Capsule razor km: %.3f\n', P.topology_capsule_razor_km);
fprintf(fid, 'Capsule fill km: %.3f\n', P.topology_capsule_fill_km);
fprintf(fid, 'Raw shelf crawling enabled: %d\n', logical(P.topology_shelf_expand));
fprintf(fid, 'Shelf isobath m: %.3f\n', P.topology_shelf_isobath_m);
fprintf(fid, 'Shelf reach km: %.3f\n', P.topology_shelf_reach_km);
fprintf(fid, 'Shelf total max area km2: %.3f\n', P.topology_shelf_total_max_area_km2);
fprintf(fid, 'Shelf razor km: %.3f\n', P.topology_shelf_razor_km);
fprintf(fid, 'Shelf min minor width km: %.3f\n', P.topology_shelf_min_minor_width_km);
fprintf(fid, 'Shelf min fill fraction: %.3f\n', P.topology_shelf_min_fill_fraction);
fprintf(fid, 'Shelf gap-fill km: %.3f\n', P.topology_shelf_gap_fill_km);
fprintf(fid, 'Shelf gap max area km2: %.3f\n', P.topology_shelf_gap_max_area_km2);
fprintf(fid, 'Ring-bay close km: %.3f\n', P.topology_bay_ring_close_km);
fprintf(fid, 'Ring-bay recover km: %.3f\n', P.topology_bay_ring_recover_km);
fprintf(fid, 'Ring-bay total max area km2: %.3f\n', P.topology_bay_ring_total_max_area_km2);
fprintf(fid, 'Ring-bay search-edge touch max fraction: %.4f\n', P.topology_bay_ring_search_edge_touch_max_fraction);
fprintf(fid, 'Named Caribbean embayments enabled: %d\n', logical(P.topology_named_caribbean_embayments_enable));
fprintf(fid, 'Named embayment extra search km: %.3f\n', P.topology_named_embayment_extra_km);
fprintf(fid, 'Named embayment window pad km: %.3f\n', P.topology_named_embayment_window_pad_km);
fprintf(fid, 'Named embayment window superellipse enabled: %d\n', logical(P.topology_named_embayment_window_use_superellipse));
fprintf(fid, 'Named embayment window scale: %.3f\n', P.topology_named_embayment_window_scale);
fprintf(fid, 'Named embayment window exponent: %.3f\n', P.topology_named_embayment_window_exponent);
fprintf(fid, 'Named embayment trigger km: %.3f\n', P.topology_named_embayment_trigger_km);
fprintf(fid, 'Bay-domain smoothing km: %.3f\n', P.topology_bay_domain_smooth_km);
fprintf(fid, 'Bay-domain ocean-only smoothing: %d\n', logical(P.topology_bay_domain_smooth_ocean_only));
fprintf(fid, 'Bay-domain land protect km: %.3f\n', P.topology_bay_domain_land_protect_km);
fprintf(fid, 'Bay-domain smoothing inset fraction: %.3f\n', P.topology_bay_domain_smooth_inset_fraction);
fprintf(fid, 'Bay-domain curve smoothing enabled: %d\n', logical(P.topology_bay_domain_curve_smooth_enable));
fprintf(fid, 'Bay-domain curve smoothing window km: %.3f\n', P.topology_bay_domain_curve_window_km);
fprintf(fid, 'Bay-domain final outline smoothing enabled: %d\n', logical(P.topology_bay_domain_final_smooth_enable));
fprintf(fid, 'Bay-domain final outline smooth km: %.3f\n', P.topology_bay_domain_final_smooth_km);
fprintf(fid, 'Bay-domain final outline land protect km: %.3f\n', P.topology_bay_domain_final_land_protect_km);
fprintf(fid, 'Bay-domain final outline curve window km: %.3f\n', P.topology_bay_domain_final_curve_window_km);
fprintf(fid, 'Final open boundary close km: %.3f\n', P.final_open_boundary_close_km);
fprintf(fid, 'Final open boundary round km: %.3f\n', P.final_open_boundary_round_km);
fprintf(fid, 'Final open boundary edge round km: %.3f\n', P.final_open_boundary_edge_round_km);
fprintf(fid, 'Final open boundary curve window km: %.3f\n', P.final_open_boundary_curve_window_km);
fprintf(fid, 'Final open boundary envelope enabled: %d\n', logical(P.final_open_boundary_envelope_enable));
fprintf(fid, 'Final open boundary max extra km: %.3f\n', P.final_open_boundary_max_extra_km);
fprintf(fid, 'Final output outer smoothing enabled: %d\n', logical(P.final_output_outer_smooth_enable));
fprintf(fid, 'Final output outer protect km: %.3f\n', P.final_output_outer_protect_km);
fprintf(fid, 'Final output outer round km: %.3f\n', P.final_output_outer_round_km);
fprintf(fid, 'Final output outer curve window km: %.3f\n', P.final_output_outer_curve_window_km);
fprintf(fid, 'Outer open-ocean area warning threshold km2: %.3f\n\n', P.outer_open_ocean_max_area_km2);
fprintf(fid, 'Optional OceanMesh2D runner: %s\n', fullfile(P.output_dir, 'P2_run_global_tc_adcirc_block_meshes.m'));
end

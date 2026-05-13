function [dem_fine, dem_info] = auto_select_fine_dem_by_point(lon0, lat0, fine_dem_dir, dem_coarse)
%AUTO_SELECT_FINE_DEM_BY_POINT
% 按点位自动选择覆盖该点的精细 DEM；若找不到，则回退到 dem_coarse
%
% 输入：
%   lon0, lat0      - 目标点经纬度
%   fine_dem_dir    - 精细 DEM 根目录（会递归搜索 .nc）
%   dem_coarse      - 找不到精细 DEM 时的默认粗 DEM（如 GEBCO）
%
% 输出：
%   dem_fine        - 选中的 DEM 文件路径（精细 nc 或 dem_coarse）
%   dem_info        - 结构体，记录选择信息
%
% 说明：
%   1) 只按“点是否落在 nc 覆盖范围内”判断
%   2) 若多个 nc 都覆盖该点，优先选择分辨率更高（网格更密）的
%   3) 若都不覆盖，则直接返回 dem_coarse

    if nargin < 4 || isempty(dem_coarse)
        error('必须提供 dem_coarse 作为默认回退 DEM。');
    end

    dem_info = struct();
    dem_info.query_lon = lon0;
    dem_info.query_lat = lat0;
    dem_info.used_fallback = false;
    dem_info.reason = '';
    dem_info.selected_file = '';
    dem_info.candidates = [];

    % 默认先回退到 coarse
    dem_fine = dem_coarse;

    if ~isfolder(fine_dem_dir)
        dem_info.used_fallback = true;
        dem_info.reason = sprintf('fine_dem_dir 不存在，回退到 coarse DEM: %s', dem_coarse);
        dem_info.selected_file = dem_fine;
        fprintf('[DEM] %s\n', dem_info.reason);
        return;
    end

    % 递归搜索所有 nc
    nc_files = dir(fullfile(fine_dem_dir, '**', '*.nc'));
    if isempty(nc_files)
        dem_info.used_fallback = true;
        dem_info.reason = sprintf('fine_dem_dir 下未找到任何 .nc，回退到 coarse DEM: %s', dem_coarse);
        dem_info.selected_file = dem_fine;
        fprintf('[DEM] %s\n', dem_info.reason);
        return;
    end

    cand = [];
    n_ok = 0;

    for i = 1:numel(nc_files)
        f = fullfile(nc_files(i).folder, nc_files(i).name);

        try
            info = inspect_nc_extent(f);

            % 记录候选信息
            n_ok = n_ok + 1;
            cand(n_ok).file = f; %#ok<AGROW>
            cand(n_ok).lon_min = info.lon_min;
            cand(n_ok).lon_max = info.lon_max;
            cand(n_ok).lat_min = info.lat_min;
            cand(n_ok).lat_max = info.lat_max;
            cand(n_ok).dx = info.dx;
            cand(n_ok).dy = info.dy;
            cand(n_ok).nx = info.nx;
            cand(n_ok).ny = info.ny;
            cand(n_ok).covers = point_in_extent(lon0, lat0, info.lon_min, info.lon_max, info.lat_min, info.lat_max);
            cand(n_ok).score = inf;  % 后面再算
            cand(n_ok).var_lon = info.var_lon;
            cand(n_ok).var_lat = info.var_lat;
        catch ME
            fprintf('[DEM] 跳过无法读取的 nc: %s\n', f);
            fprintf('      原因: %s\n', ME.message);
        end
    end

    if isempty(cand)
        dem_info.used_fallback = true;
        dem_info.reason = sprintf('未找到可解析范围的精细 nc，回退到 coarse DEM: %s', dem_coarse);
        dem_info.selected_file = dem_fine;
        fprintf('[DEM] %s\n', dem_info.reason);
        return;
    end

    dem_info.candidates = cand;

    % 只保留覆盖该点的 nc
    I = find([cand.covers]);
    if isempty(I)
        dem_info.used_fallback = true;
        dem_info.reason = sprintf('未找到覆盖点(%.4f, %.4f)的精细 nc，回退到 coarse DEM: %s', ...
                                  lon0, lat0, dem_coarse);
        dem_info.selected_file = dem_fine;
        fprintf('[DEM] %s\n', dem_info.reason);
        return;
    end

    % 多个覆盖时，优先选分辨率更高的（dx*dy 更小）
    for k = 1:numel(I)
        ii = I(k);
        if isfinite(cand(ii).dx) && isfinite(cand(ii).dy)
            cand(ii).score = abs(cand(ii).dx) * abs(cand(ii).dy);
        else
            % 没法判断分辨率时给较差分
            cand(ii).score = inf;
        end
    end

    cover_cands = cand(I);
    [~, ibest_local] = min([cover_cands.score]);
    best = cover_cands(ibest_local);

    dem_fine = best.file;
    dem_info.used_fallback = false;
    dem_info.reason = sprintf('找到覆盖点的精细 DEM: %s', best.file);
    dem_info.selected_file = dem_fine;

    fprintf('[DEM] 使用精细 DEM: %s\n', dem_fine);
    fprintf('      范围: lon=[%.4f, %.4f], lat=[%.4f, %.4f]\n', ...
        best.lon_min, best.lon_max, best.lat_min, best.lat_max);
    fprintf('      分辨率近似: dx=%.8f, dy=%.8f\n', best.dx, best.dy);
end

%% ===== 子函数：读取 nc 覆盖范围 =====
function out = inspect_nc_extent(ncfile)
    info = ncinfo(ncfile);
    varnames = {info.Variables.Name};

    lon_candidates = {'lon','longitude','x'};
    lat_candidates = {'lat','latitude','y'};

    var_lon = pick_first_existing(varnames, lon_candidates);
    var_lat = pick_first_existing(varnames, lat_candidates);

    if isempty(var_lon) || isempty(var_lat)
        error('未找到经纬度变量（lon/lat/longitude/latitude/x/y）');
    end

    lon = double(ncread(ncfile, var_lon));
    lat = double(ncread(ncfile, var_lat));

    lon = lon(:);
    lat = lat(:);

    lon = lon(isfinite(lon));
    lat = lat(isfinite(lat));

    if isempty(lon) || isempty(lat)
        error('经纬度变量为空或全为 NaN');
    end

    % 兼容 0~360，经查询点通常是 -180~180
    lon_adj = lon;
    if max(lon_adj) > 180 && any(lon_adj > 180)
        lon_adj(lon_adj > 180) = lon_adj(lon_adj > 180) - 360;
    end

    lon_min = min(lon_adj);
    lon_max = max(lon_adj);
    lat_min = min(lat);
    lat_max = max(lat);

    dx = estimate_spacing(lon_adj);
    dy = estimate_spacing(lat);

    out = struct();
    out.file = ncfile;
    out.var_lon = var_lon;
    out.var_lat = var_lat;
    out.lon_min = lon_min;
    out.lon_max = lon_max;
    out.lat_min = lat_min;
    out.lat_max = lat_max;
    out.dx = dx;
    out.dy = dy;
    out.nx = numel(unique(lon_adj));
    out.ny = numel(unique(lat));
end

%% ===== 子函数：判断点是否在范围内 =====
function tf = point_in_extent(lon0, lat0, lon_min, lon_max, lat_min, lat_max)
    tol = 1e-8;
    tf = (lon0 >= lon_min - tol) && (lon0 <= lon_max + tol) && ...
         (lat0 >= lat_min - tol) && (lat0 <= lat_max + tol);
end

%% ===== 子函数：估计分辨率 =====
function d = estimate_spacing(v)
    v = unique(v(:));
    v = v(isfinite(v));
    if numel(v) < 2
        d = inf;
        return;
    end
    dv = diff(sort(v));
    dv = dv(abs(dv) > 0);
    if isempty(dv)
        d = inf;
    else
        d = median(abs(dv));
    end
end

%% ===== 子函数：从候选名里挑第一个存在的 =====
function name = pick_first_existing(varnames, candidates)
    name = '';
    for i = 1:numel(candidates)
        idx = find(strcmpi(varnames, candidates{i}), 1);
        if ~isempty(idx)
            name = varnames{idx};
            return;
        end
    end
end
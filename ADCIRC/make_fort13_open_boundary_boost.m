function make_fort13_open_boundary_boost(auto_Name)
% ============================================================
% 自动生成 ADCIRC fort.13
%
% 功能：
% 1) 从 fort.14 读取网格节点和开边界节点
% 2) 根据“到开边界的距离”自动增大底摩擦
% 3) 根据 CF_uprange.shp 对局部易发散区域进一步增大底摩擦
% 4) 可选：浅水区采用 Manning 形式 Cf
% 5) 先生成一个参考 fort.13，再复制到 cases_<auto_Name> 下各事件目录
%
% 说明：
% - 适用于任意局地 ADCIRC 网格
% - 不再依赖固定 river/sea 分区
% - 如果找不到 CF_uprange.shp，会自动跳过局部增强，只保留开边界增强
% ============================================================

    if nargin < 1 || isempty(auto_Name)
        auto_Name = 'BoB';
    end

    % 事件目录根目录
    case_root = ['cases_' auto_Name];

    % fort.14 搜索优先级：
    % 1) ADCIRC_Capsule_Mesh_<auto_Name>.14
    % 2) fort.14
    % 3) 当前目录唯一一个 .14 文件
    fort14_file = resolve_fort14_file(auto_Name);

    % 输出参考 fort.13
    out_fort13 = 'fort.13';

    % ---------- fort.13 默认属性 ----------
    slope_limiter_default = 0.01;   % elemental_slope_limiter
    cf_default            = 0.0015; % quadratic_friction_coefficient_at_sea_floor

    % ---------- 是否启用浅水 Manning 型 Cf ----------
    use_shallow_manning = true;
    manning_n           = 0.01;
    shallow_depth_max   = 20.0;   % m，水深小于该值时允许用 Manning 型 Cf
    depth_floor         = 5;   % m，防止过浅导致 Cf 过大
    cf_cap              = 0.020;  % 全局 Cf 上限，防止过强

    % ---------- 开边界增摩擦参数 ----------
    % 距离开边界 <= inner_km：使用最大增幅 max_boost
    % 距离在 [inner_km, outer_km]：线性衰减到 1
    % 距离 > outer_km：不增强
    inner_km  = 30;     % 最强增强带
    outer_km  = 50;    % 增强衰减终点
    max_boost = 0.015/cf_default;    % 开边界处最大放大倍数

    % ---------- CF_uprange.shp 局部增强 ----------
    use_cf_uprange = true;
    cf_uprange_shp = 'H:\Global_compoundflood\ADCIRC\OceanMesh2D-Projection\datasets\CF_uprange.shp';   % 你也可以改成完整路径
    local_boost    = 10.0;                % shp 区域内再乘这个倍数
    local_cf_cap   = 0.050;              % shp 区域允许更高的上限

    % ---------- 写法 ----------
    % true  -> fort.13 中只写“非默认节点”
    % false -> 把所有节点都写进去
    write_only_nondefault = true;

    % ---------- 复制 ----------
    copy_to_case_dirs = true;

    %% =========================
    % 1) 读取 fort.14
    %% =========================
    fprintf('读取 fort.14: %s\n', fort14_file);

    G = read_fort14_nodes_and_open_boundary(fort14_file);

    fprintf('节点数 NP        = %d\n', G.NP);
    fprintf('单元数 NE        = %d\n', G.NE);
    fprintf('开边界条数 NOPE  = %d\n', G.NOPE);
    fprintf('开边界节点数     = %d\n', numel(G.open_nodes));

    if isempty(G.open_nodes)
        error('未在 fort.14 中读到开边界节点，无法构造开边界增摩擦 fort.13。');
    end

    %% =========================
    % 2) 计算每个节点到开边界的距离（m）
    %% =========================
    fprintf('计算节点到开边界的距离...\n');

    lon0 = mean(G.lon, 'omitnan');
    lat0 = mean(G.lat, 'omitnan');

    [x_all, y_all] = ll2xy_local(G.lon, G.lat, lon0, lat0);
    [x_ob,  y_ob]  = ll2xy_local(G.lon(G.open_nodes), G.lat(G.open_nodes), lon0, lat0);

    dist_to_open_m  = min_distance_to_points_chunked(x_all, y_all, x_ob, y_ob, 20000);
    dist_to_open_km = dist_to_open_m / 1000;

    %% =========================
    % 3) 计算基础 Cf
    %% =========================
    depth = G.depth(:);          % ADCIRC fort.14 通常水深为正
    depth_use = max(depth, depth_floor);

    cf_base = cf_default * ones(G.NP, 1);

    if use_shallow_manning
        I_shallow = depth > 0 & depth <= shallow_depth_max;
        cf_manning = 9.81 * manning_n^2 ./ (depth_use .^ (1/3));
        cf_manning = min(cf_manning, cf_cap);
        cf_base(I_shallow) = max(cf_base(I_shallow), cf_manning(I_shallow));
    end

    %% =========================
    % 4) 根据开边界距离放大 Cf
    %% =========================
    boost = ones(G.NP, 1);

    I1 = dist_to_open_km <= inner_km;
    I2 = dist_to_open_km > inner_km & dist_to_open_km < outer_km;

    boost(I1) = max_boost;
    boost(I2) = 1 + (max_boost - 1) * (outer_km - dist_to_open_km(I2)) / (outer_km - inner_km);

    cf_final = cf_base .* boost;
    cf_final = max(cf_final, cf_default);
    cf_final = min(cf_final, cf_cap);

    %% =========================
    % 5) 根据 CF_uprange.shp 进一步局部增加 Cf
    %% =========================
    local_mask = false(G.NP, 1);

    if use_cf_uprange
        shp_found = false;

        if exist(cf_uprange_shp, 'file')
            shp_found = true;
            shp_file_use = cf_uprange_shp;
        else
            % 自动在当前目录及 datasets 下找一下
            cand = {
                fullfile(pwd, 'CF_uprange.shp')
                fullfile('..', 'datasets', 'CF_uprange.shp')
                fullfile('..', 'datasets', 'shp', 'CF_uprange.shp')
                fullfile('..', 'datasets', 'coastalline', 'CF_uprange.shp')
                };
            shp_file_use = '';
            for ii = 1:numel(cand)
                if exist(cand{ii}, 'file')
                    shp_found = true;
                    shp_file_use = cand{ii};
                    break;
                end
            end
        end

        if shp_found
            fprintf('读取局部增摩擦 shp: %s\n', shp_file_use);
            local_mask = points_in_shapefile(G.lon, G.lat, shp_file_use);
            fprintf('CF_uprange.shp 内节点数 = %d\n', nnz(local_mask));

            if any(local_mask)
                % 在已有 cf_final 基础上进一步局部增强
                cf_final(local_mask) = cf_final(local_mask) .* local_boost;
                cf_final(local_mask) = min(cf_final(local_mask), local_cf_cap);
            end
        else
            warning('未找到 CF_uprange.shp，跳过局部增摩擦，只保留开边界增强。');
        end
    end

    fprintf('Cf 范围: [%.6f, %.6f]\n', min(cf_final), max(cf_final));

    %% =========================
    % 6) 写参考 fort.13
    %% =========================
    fprintf('写出参考 fort.13: %s\n', out_fort13);

    if write_only_nondefault
        idx_write = find(abs(cf_final - cf_default) > 1e-12);
    else
        idx_write = (1:G.NP).';
    end

    write_fort13_file(out_fort13, G.NP, slope_limiter_default, cf_default, idx_write, cf_final(idx_write));

    fprintf('fort.13 已写出。\n');
    fprintf('写入节点数 = %d / %d\n', numel(idx_write), G.NP);

    %% =========================
    % 7) 复制到各事件目录
    %% =========================
    if copy_to_case_dirs
        if ~exist(case_root, 'dir')
            warning('未找到事件目录 %s，跳过复制。', case_root);
        else
            D = dir(case_root);
            D = D([D.isdir]);
            D = D(~ismember({D.name}, {'.','..'}));

            fprintf('复制 fort.13 到 %d 个事件目录...\n', numel(D));

            n_ok = 0;
            for i = 1:numel(D)
                case_dir = fullfile(case_root, D(i).name);
                dst = fullfile(case_dir, 'fort.13');
                copyfile(out_fort13, dst);
                n_ok = n_ok + 1;
            end

            fprintf('已复制到 %d 个事件目录。\n', n_ok);
        end
    end

    %% =========================
    % 8) 检查图
    %% =========================
    figure('Color','w');
    scatter(dist_to_open_km, cf_final, 8, depth, 'filled');
    xlabel('Distance to open boundary (km)');
    ylabel('Cf');
    title('fort.13 friction vs distance to open boundary');
    colorbar;
    grid on;

    figure('Color','w');
    scatter(G.lon, G.lat, 6, cf_final, 'filled'); hold on;
    if any(local_mask)
        plot(G.lon(local_mask), G.lat(local_mask), 'k.', 'MarkerSize', 4);
    end
    clim([0.0015 0.003])
    axis equal tight;
    xlabel('Lon');
    ylabel('Lat');
    title('fort.13 spatial Cf');
    colorbar;
    grid on;

    fprintf('\n完成。\n');
end

%% ============================================================
% 解析 fort.14：节点 + 开边界
%% ============================================================
function G = read_fort14_nodes_and_open_boundary(fort14_file)

    fid = fopen(fort14_file, 'rt');
    if fid < 0
        error('无法打开 fort.14: %s', fort14_file);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    G = struct();

    G.title = strtrim(fgetl(fid));

    line2 = strtrim(fgetl(fid));
    a = sscanf(line2, '%d %d');
    if numel(a) < 2
        error('fort.14 第2行解析失败。');
    end
    G.NE = a(1);
    G.NP = a(2);

    % ---------- 节点 ----------
    node_id = zeros(G.NP,1);
    lon     = zeros(G.NP,1);
    lat     = zeros(G.NP,1);
    depth   = zeros(G.NP,1);

    for i = 1:G.NP
        s = fgetl(fid);
        v = sscanf(s, '%d %f %f %f');
        if numel(v) < 4
            error('fort.14 节点行读取失败，第 %d 个节点。', i);
        end
        node_id(i) = v(1);
        lon(i)     = v(2);
        lat(i)     = v(3);
        depth(i)   = v(4);
    end

    G.node_id = node_id;
    G.lon     = lon;
    G.lat     = lat;
    G.depth   = depth;

    % ---------- 跳过单元 ----------
    for i = 1:G.NE
        fgetl(fid);
    end

    % ---------- 开边界 ----------
    s = strtrim(fgetl(fid));
    G.NOPE = sscanf(s, '%d', 1);

    s = strtrim(fgetl(fid));
    G.NETA = sscanf(s, '%d', 1);

    open_nodes = [];

    for ib = 1:G.NOPE
        header = strtrim(fgetl(fid));
        hv = sscanf(header, '%d');
        if isempty(hv)
            error('开边界第 %d 条 header 读取失败。', ib);
        end
        nvdll = hv(1);

        tmp = zeros(nvdll,1);
        for k = 1:nvdll
            line = strtrim(fgetl(fid));
            vv = sscanf(line, '%d');
            if isempty(vv)
                error('开边界第 %d 条第 %d 个节点读取失败。', ib, k);
            end
            tmp(k) = vv(1);
        end
        open_nodes = [open_nodes; tmp]; %#ok<AGROW>
    end

    G.open_nodes = unique(open_nodes(:));
end

%% ============================================================
% 局地投影：lon/lat -> x/y (m)
%% ============================================================
function [x, y] = ll2xy_local(lon, lat, lon0, lat0)
    x = (lon - lon0) .* cosd(lat0) * 111320;
    y = (lat - lat0) * 110540;
end

%% ============================================================
% 分块计算每个点到开边界节点的最小距离
%% ============================================================
function dmin = min_distance_to_points_chunked(x, y, xb, yb, chunk_size)

    n = numel(x);
    dmin = inf(n,1);

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

        if mod(i1, chunk_size*5) == 1 || i2 == n
            fprintf('  距离计算进度: %d / %d\n', i2, n);
        end
    end
end

%% ============================================================
% 判断点是否落在 shapefile 多边形内部
%% ============================================================
function mask = points_in_shapefile(lon, lat, shpfile)

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

%% ============================================================
% 写 fort.13
%% ============================================================
function write_fort13_file(fname, np, slope_limiter_default, cf_default, idx_write, cf_vals)

    fid = fopen(fname, 'wt');
    if fid < 0
        error('无法写入 fort.13: %s', fname);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'fort13\n');
    fprintf(fid, '%d\n', np);
    fprintf(fid, '2\n');

    % 属性1：elemental_slope_limiter
    fprintf(fid, 'elemental_slope_limiter\n');
    fprintf(fid, '1\n');
    fprintf(fid, '1\n');
    fprintf(fid, '%.8f\n', slope_limiter_default);

    % 属性2：quadratic_friction_coefficient_at_sea_floor
    fprintf(fid, 'quadratic_friction_coefficient_at_sea_floor\n');
    fprintf(fid, 'm\n');
    fprintf(fid, '1\n');
    fprintf(fid, '%.8f\n', cf_default);

    % 默认值区
    fprintf(fid, 'elemental_slope_limiter\n');
    fprintf(fid, '0\n');

    fprintf(fid, 'quadratic_friction_coefficient_at_sea_floor\n');
    fprintf(fid, '%d\n', numel(idx_write));

    for i = 1:numel(idx_write)
        fprintf(fid, '%d %.8f\n', idx_write(i), cf_vals(i));

        if mod(i, 10000) == 0 || i == numel(idx_write)
            fprintf('  fort.13 写出进度: %d / %d\n', i, numel(idx_write));
        end
    end
end

%% ============================================================
% 自动寻找 fort.14
%% ============================================================
function fort14_file = resolve_fort14_file(auto_Name)

    cand1 = ['ADCIRC_Capsule_Mesh_' auto_Name '.14'];
    cand2 = 'fort.14';

    if exist(cand1, 'file')
        fort14_file = cand1;
        return;
    end

    if exist(cand2, 'file')
        fort14_file = cand2;
        return;
    end

    D = dir('*.14');
    if numel(D) == 1
        fort14_file = D(1).name;
        return;
    elseif numel(D) > 1
        names = strjoin({D.name}, '\n');
        error('当前目录有多个 .14 文件，请手动保留一个。\n%s', names);
    else
        error('当前目录未找到 fort.14 文件。');
    end
end
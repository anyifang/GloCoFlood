function auto_adcirc_domain_morphological(lon_c, lat_c, alongshore_radius, offshore_deg, inland_deg, bay_close_deg, gshhs_shp, out_shp)
    % 生成“只包含一个完整海域”的 ADCIRC 计算域
    % 核心修复：铰链式喇叭口包络！外海侧绝对垂直，内陆侧张开捕获所有端点湾区！
    
    fprintf('\n=== 启动终极形态学平滑计算域生成引擎 ===\n');
    fprintf('中心: [%.4f, %.4f]\n', lon_c, lat_c);
    
    %----------------------------
    % 0) 参数 sanity check
    %----------------------------
    mustBeFiniteScalar(alongshore_radius);
    mustBeFiniteScalar(offshore_deg);
    mustBeFiniteScalar(inland_deg);
    mustBeFiniteScalar(bay_close_deg);
    alongshore_radius = abs(alongshore_radius);
    offshore_deg      = abs(offshore_deg);
    inland_deg        = abs(inland_deg);
    bay_close_deg     = abs(bay_close_deg);
    
    %----------------------------
    % 1) 读取陆地数据（只读 bbox 内）
    %----------------------------
    bbox_size = alongshore_radius + offshore_deg + inland_deg + bay_close_deg + 5;
    bbox = [lon_c - bbox_size, lat_c - bbox_size; lon_c + bbox_size, lat_c + bbox_size];
    fprintf('1. 读取并合并真实陆地数据...\n');
    S = shaperead(gshhs_shp, 'BoundingBox', bbox);
    P_all_land = polyshape();
    for i = 1:numel(S)
        x = S(i).X(:); y = S(i).Y(:);
        v = ~isnan(x) & ~isnan(y);
        if nnz(v) < 4, continue; end
        try
            p = polyshape(x(v), y(v), 'Simplify', true);
        catch
            continue;
        end
        if isempty(p.Vertices), continue; end
        P_all_land = union(P_all_land, p);
    end
    if area(P_all_land) == 0
        error('读取到的陆地为空：请检查 gshhs_shp 路径/BoundingBox 是否正确。');
    end
    
    % 过滤太小碎片
    rgn = regions(P_all_land);
    a = area(rgn);
    P_clean_land = polyshape();
    for i = 1:numel(rgn)
        if a(i) > 0.5, P_clean_land = union(P_clean_land, rgn(i)); end
    end
    P_clean_land = rmholes(P_clean_land);
    
    %----------------------------
    % 2) 锁定主大陆并桥接宏观湾区
    %----------------------------
    fprintf('2. 锁定主大陆并桥接宏观湾区...\n');
    th = linspace(0, 2*pi, 1200);
    P_roi = polyshape(lon_c + alongshore_radius*cos(th), lat_c + alongshore_radius*sin(th));
    P_local = intersect(P_clean_land, P_roi);
    rgns_local = regions(P_local);
    if isempty(rgns_local), rgns_local = regions(P_clean_land); end
    [~, imax] = max(area(rgns_local));
    P_home = rmholes(rgns_local(imax));
    
    P_macro = safe_polybuffer(P_home,  bay_close_deg);
    P_macro = safe_polybuffer(P_macro, -bay_close_deg);
    P_macro = rmholes(P_macro);
    P_bays = subtract(P_macro, P_home);
    
    %----------------------------
    % 3) 提取目标岸段并平滑
    %----------------------------
    fprintf('3. 提取中心海岸线段并平滑...\n');
    [xb, yb] = boundary(P_macro);
    v = ~isnan(xb) & ~isnan(yb);
    xb = xb(v); yb = yb(v);
    d = hypot(xb - lon_c, yb - lat_c);
    in_circle = d <= alongshore_radius;
    
    padded = [false; in_circle(:); false];
    starts = find(diff(padded) == 1); ends = find(diff(padded) == -1) - 1;
    if isempty(starts), idx = 1:numel(xb); else
        [~, k] = max(ends - starts + 1); idx = starts(k):ends(k);
    end
    x_line = xb(idx); y_line = yb(idx);
    x_line = x_line(:); y_line = y_line(:);
    
    win = max(21, 2*round(numel(x_line)*0.02)+1); win = min(win, 301);
    x_smooth = smoothdata(x_line, 'movmean', win);
    y_smooth = smoothdata(y_line, 'movmean', win);
    
    %----------------------------
    % 4) 法向推演与端点垂直强制锚定
    %----------------------------
    fprintf('4. 执行法向推演，并强制两侧边界绝对垂直...\n');
    dx = gradient(x_smooth); dy = gradient(y_smooth);
    L  = hypot(dx, dy); L(L==0) = 1e-6;
    nx =  dy ./ L; ny = -dx ./ L;
    
    mid = round(numel(x_smooth)/2);
    if isinterior(P_macro, x_smooth(mid) + nx(mid) * 1.0, y_smooth(mid) + ny(mid) * 1.0)
        nx = -nx; ny = -ny; 
    end
    
    % 提取端点的鲁棒法向，避免边缘抖动导致侧壁歪斜
    n_pts = max(3, min(20, round(length(nx)*0.02)));
    nx_start = mean(nx(1:n_pts)); ny_start = mean(ny(1:n_pts));
    L_s = hypot(nx_start, ny_start); nx_start = nx_start/L_s; ny_start = ny_start/L_s;
    
    nx_end = mean(nx(end-n_pts+1:end)); ny_end = mean(ny(end-n_pts+1:end));
    L_e = hypot(nx_end, ny_end); nx_end = nx_end/L_e; ny_end = ny_end/L_e;
    
    x_out_raw = x_smooth + offshore_deg .* nx;
    y_out_raw = y_smooth + offshore_deg .* ny;
    
    % 强制让外海原始端点 100% 坐落在绝对法向上
    x_out_raw(1) = x_smooth(1) + offshore_deg * nx_start;
    y_out_raw(1) = y_smooth(1) + offshore_deg * ny_start;
    x_out_raw(end) = x_smooth(end) + offshore_deg * nx_end;
    y_out_raw(end) = y_smooth(end) + offshore_deg * ny_end;
    
    % 镜像翻转消除凹陷 (保留绝杀逻辑)
    P1 = [x_out_raw(1), y_out_raw(1)]; P2 = [x_out_raw(end), y_out_raw(end)];
    V = P2 - P1; L_V = norm(V);
    if L_V > 1e-3
        N_chord = [V(2), -V(1)] / L_V;
        if dot(N_chord, [mean(nx), mean(ny)]) < 0, N_chord = -N_chord; end
        x_out_flip = x_out_raw; y_out_flip = y_out_raw;
        for k = 1:length(x_out_raw)
            Pk = [x_out_raw(k), y_out_raw(k)];
            dist = dot(Pk - P1, N_chord);
            if dist < 0 
                P_flip = Pk - 2 * dist * N_chord; 
                x_out_flip(k) = P_flip(1); y_out_flip(k) = P_flip(2);
            end
        end
        sm_win_out = max(50, round(length(x_out_flip) * 0.2));
        x_out = smoothdata(x_out_flip, 'gaussian', sm_win_out);
        y_out = smoothdata(y_out_flip, 'gaussian', sm_win_out);
    else
        x_out = x_out_raw; y_out = y_out_raw;
    end
    
    % 将平滑后的外海边界强制融合锚定回绝对法向端点
    blend_len = max(10, floor(length(x_out) * 0.05));
    w = linspace(1, 0, blend_len)';
    x_out(1:blend_len) = w .* x_out_raw(1) + (1-w) .* x_out(1:blend_len);
    y_out(1:blend_len) = w .* y_out_raw(1) + (1-w) .* y_out(1:blend_len);
    
    x_out(end-blend_len+1:end) = flipud(w) .* x_out_raw(end) + (1-flipud(w)) .* x_out(end-blend_len+1:end);
    y_out(end-blend_len+1:end) = flipud(w) .* y_out_raw(end) + (1-flipud(w)) .* y_out(end-blend_len+1:end);
    
    %----------------------------
    % 5) 组装铰链式喇叭口包络框 (Hinged Flared Envelope)
    %----------------------------
    fprintf('5. 生成喇叭口包络框 (海洋段严格垂直，内陆段张开捕获湾区)...\n');
    L_deep = inland_deg + 1.5; 
    x_in  = x_smooth - L_deep .* nx;
    y_in  = y_smooth - L_deep .* ny;
    
    % 强制让内陆端点也初始坐落在绝对法向上
    x_in(1) = x_smooth(1) - L_deep * nx_start;
    y_in(1) = y_smooth(1) - L_deep * ny_start;
    x_in(end) = x_smooth(end) - L_deep * nx_end;
    y_in(end) = y_smooth(end) - L_deep * ny_end;
    
    % ★ 核心修复 3：内陆喇叭口张开！(Flaring)
    % 使得包络面在陆地内部呈巨大扇形张开，从而完美兜住所有角落的湾区！
    flare_dist = max(2.0, inland_deg * 1.5);
    
    % 起点：向“后”张开
    x_in(1) = x_in(1) + flare_dist * ny_start;
    y_in(1) = y_in(1) - flare_dist * nx_start;
    
    % 终点：向“前”张开
    x_in(end) = x_in(end) - flare_dist * ny_end;
    y_in(end) = y_in(end) + flare_dist * nx_end;
    
    % 组装多边形：必须插入 x_smooth(1) 和 x_smooth(end) 作为铰链锚点！
    % 这确保了 [x_out(1) 到 x_smooth(1)] 这段横跨海面的线依然是 100% 垂直的直尺！
    poly_x = [flipud(x_out); x_smooth(1); x_in; x_smooth(end)];
    poly_y = [flipud(y_out); y_smooth(1); y_in; y_smooth(end)];
    
    Envelope = polyshape(poly_x, poly_y, 'Simplify', true);
    Envelope = rmholes(union(Envelope)); 
    
    if area(P_bays) > 0
        bays = regions(P_bays);
        for i = 1:numel(bays)
            bay = bays(i);
            if area(bay) == 0, continue; end
            if area(intersect(Envelope, bay)) > 0.01 * area(bay)
                spill = safe_polybuffer(bay, max(inland_deg, offshore_deg)+0.2);
                Envelope = union(Envelope, spill);
            end
        end
    end
    Envelope = rmholes(union(Envelope));
    
    seed_x = x_smooth(mid) + nx(mid) * max(0.8, 0.6*offshore_deg);
    seed_y = y_smooth(mid) + ny(mid) * max(0.8, 0.6*offshore_deg);
    
    %----------------------------
    % 6) 形态学平滑陆地 + 双海域自适应向内缩
    %----------------------------
    fprintf('6. 执行自适应深层内推...\n');
    smoothing_radius = min(0.5, max(0.15, 0.4*bay_close_deg));
    P_smooth_land = safe_polybuffer(P_clean_land,  smoothing_radius);
    P_smooth_land = safe_polybuffer(P_smooth_land, -smoothing_radius);
    P_smooth_land = rmholes(P_smooth_land);
    
    [xe, ye] = boundary(Envelope);
    ve = ~isnan(xe) & ~isnan(ye);
    pad = inland_deg + max(0.2, 0.5*bay_close_deg) + 0.5;
    bbox4 = [min(xe(ve))-pad, max(xe(ve))+pad, min(ye(ve))-pad, max(ye(ve))+pad]; 
    
    dx0 = max(min(inland_deg/2, 0.02), 0.005);
    keep_margin = max(2*dx0, 0.05); 

    P_deep_land = adaptive_inland_shrink_mask(P_smooth_land, Envelope, bbox4, inland_deg, dx0, keep_margin, seed_x, seed_y);
    P_domain = subtract(Envelope, P_deep_land);
    
    %----------------------------
    % 7) 拓扑清理
    %----------------------------
    fprintf('7. 拓扑清理...\n');
    rdom = regions(P_domain);
    if isempty(rdom)
        error('P_domain 为空：请检查参数 offshore_deg / inland_deg 是否合理。');
    end
    hit = false(numel(rdom),1);
    for i = 1:numel(rdom)
        hit(i) = isinterior(rdom(i), seed_x, seed_y);
    end
    if any(hit)
        P_domain = rdom(find(hit,1,'first'));
    else
        [~, im] = max(area(rdom)); P_domain = rdom(im);
    end
    P_domain = rmholes(P_domain);
    
    % 仅保留极微小的抗栅格锯齿打磨 (0.15度)，100% 保护侧壁像刀切一样笔直！
    fillet_r = 0.15; 
    P_domain = safe_polybuffer(P_domain,  fillet_r);
    P_domain = safe_polybuffer(P_domain, -fillet_r);
    P_domain = rmholes(P_domain);
    
    %----------------------------
    % 8) 输出 shapefile
    %----------------------------
    fprintf('8. 输出 shapefile: %s\n', out_shp);
    [x_dom, y_dom] = boundary(P_domain);
    vv = ~isnan(x_dom) & ~isnan(y_dom);
    x_dom = x_dom(vv); y_dom = y_dom(vv);
    Sout = struct('Geometry', 'Polygon', ...
                  'BoundingBox', [min(x_dom), min(y_dom); max(x_dom), max(y_dom)], ...
                  'X', [x_dom', NaN], ...
                  'Y', [y_dom', NaN], ...
                  'Name', 'ADCIRC_Morph_Domain');
    out_dir = fileparts(out_shp);
    if ~isempty(out_dir) && ~exist(out_dir, 'dir'), mkdir(out_dir); end
    shapewrite(Sout, out_shp);
    fprintf('(+) 计算域已生成: %s\n', out_shp);
    fprintf('===============================================\n');
    
    % --- 可视化预览 ---
    if exist('m_proj', 'file') ~= 0
        fprintf('>>> 正在调用 m_map 绘制高质量预览图...\n');
        plot_lon_lim = [lon_c - alongshore_radius - 2, lon_c + offshore_deg + 2];
        plot_lat_lim = [lat_c - alongshore_radius - 2, lat_c + alongshore_radius + 2];
        figure('Name', 'ADCIRC Preview', 'Color', 'w', 'Position', [150, 100, 900, 750]);
        m_proj('mercator', 'lon', plot_lon_lim, 'lat', plot_lat_lim); hold on;
        for i = 1:numel(S)
            x_cst = S(i).X; y_cst = S(i).Y; v = ~isnan(x_cst) & ~isnan(y_cst);
            if nnz(v) > 3, m_patch(x_cst(v), y_cst(v), [0.9 0.85 0.75], 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.5); end
        end
        m_patch(x_dom, y_dom, 'r', 'FaceAlpha', 0.3, 'EdgeColor', [0 0.5 1], 'LineWidth', 2.5);
        m_grid('box', 'fancy', 'tickdir', 'in', 'fontsize', 10, 'linewidth', 1.5, 'linestyle', ':', 'color', 'k');
        m_line(lon_c, lat_c, 'marker', 'p', 'color', 'b', 'linewi', 2, 'markersize', 12, 'markerfacecolor', 'y');
        drawnow;
    end
end

% =====================================================================
% helper: 安全 polybuffer
% =====================================================================
function P2 = safe_polybuffer(P1, d)
    if area(P1) == 0, P2 = P1; return; end
    try P2 = polybuffer(P1, d, 'JointType', 'round'); catch, P2 = P1; end
end

% =====================================================================
% ★ 终极防线 helper: 摧毁湖泊排斥力场，真实海洋边界判定
% =====================================================================
function P_land_in = adaptive_inland_shrink_mask(P_land, Envelope, bbox4, inland_deg, dx0, keep_margin_deg, seed_x, seed_y)
    lon_min = bbox4(1); lon_max = bbox4(2); lat_min = bbox4(3); lat_max = bbox4(4);
    lon_rng = lon_max - lon_min; lat_rng = lat_max - lat_min;
    
    MAX_PIXELS = 4e6; 
    dx_need = sqrt((lon_rng * lat_rng) / MAX_PIXELS);
    dx = max(dx0, dx_need);
    nlon = max(20, ceil(lon_rng/dx) + 1); nlat = max(20, ceil(lat_rng/dx) + 1);
    
    P_rect = polyshape([lon_min lon_max lon_max lon_min], [lat_min lat_min lat_max lat_max]);
    P_local = intersect(P_land, P_rect);
    if area(P_local) == 0, P_land_in = polyshape(); return; end
    
    pix = lon_rng/(nlon-1); 
    landMask = polyshape_to_mask(P_local, lon_min, lat_min, pix, nlat, nlon);
    
    EnvMask = polyshape_to_mask(Envelope, lon_min, lat_min, pix, nlat, nlon);
    oceanMask = ~landMask & EnvMask;
    
    col_seed = round((seed_x - lon_min) / pix) + 1; row_seed = round((seed_y - lat_min) / pix) + 1;
    col_seed = max(1, min(nlon, col_seed)); row_seed = max(1, min(nlat, row_seed));
    
    TargetOceanMask = bwselect(oceanMask, col_seed, row_seed, 8);
    if ~any(TargetOceanMask(:))
        [yo, xo] = find(oceanMask);
        if ~isempty(xo)
            [~, min_idx] = min((xo - col_seed).^2 + (yo - row_seed).^2);
            TargetOceanMask = bwselect(oceanMask, xo(min_idx), yo(min_idx), 8);
        else
            TargetOceanMask = oceanMask; 
        end
    end
    
    RemainingWater = oceanMask & ~TargetOceanMask;
    
    EnvEdge = bwperim(EnvMask);
    [edge_y, edge_x] = find(EnvEdge & RemainingWater);
    if ~isempty(edge_x)
        OtherOceanMask = bwselect(RemainingWater, edge_x, edge_y, 8);
    else
        OtherOceanMask = false(size(RemainingWater));
    end
    
    D_target = bwdist(TargetOceanMask) * pix;
    if any(OtherOceanMask(:))
        D_other = bwdist(OtherOceanMask) * pix;
    else
        D_other = inf(size(D_target)); 
    end
    
    landMask_survived = landMask & ((D_target > inland_deg) | (D_target >= D_other - keep_margin_deg));
    
    P_land_in = mask_to_polyshape(landMask_survived, lon_min, lat_min, pix);
    P_land_in = rmholes(P_land_in);
end

% =====================================================================
% helper: polyshape -> mask
% =====================================================================
function mask = polyshape_to_mask(P, lon_min, lat_min, pix, nlat, nlon)
    mask = false(nlat, nlon);
    P_outer = rmholes(P); rr = regions(P_outer);
    for i = 1:numel(rr)
        [x, y] = boundary(rr(i)); v = ~isnan(x) & ~isnan(y); x = x(v); y = y(v);
        if numel(x) < 3, continue; end
        col = (x - lon_min) / pix + 1; row = (y - lat_min) / pix + 1;
        mask = mask | poly2mask(col, row, nlat, nlon);
    end
    P_holes = subtract(P_outer, P);
    if area(P_holes) > 0
        hh = regions(P_holes);
        for i = 1:numel(hh)
            [x, y] = boundary(hh(i)); v = ~isnan(x) & ~isnan(y); x = x(v); y = y(v);
            if numel(x) < 3, continue; end
            col = (x - lon_min) / pix + 1; row = (y - lat_min) / pix + 1;
            mask(poly2mask(col, row, nlat, nlon)) = false;
        end
    end
end

% =====================================================================
% helper: mask -> polyshape
% =====================================================================
function P = mask_to_polyshape(mask, lon_min, lat_min, pix)
    B = bwboundaries(mask, 8, 'noholes');
    P = polyshape();
    for k = 1:numel(B)
        rc = B{k}; r = rc(:,1); c = rc(:,2);
        x = lon_min + (c-1) * pix; y = lat_min + (r-1) * pix;
        try
            pk = polyshape(x, y, 'Simplify', true);
            if ~isempty(pk.Vertices), P = union(P, pk); end
        catch, continue; end
    end
end

function mustBeFiniteScalar(x)
    if ~isscalar(x) || ~isnumeric(x) || ~isfinite(x), error('参数必须是有限的数值标量。'); end
end
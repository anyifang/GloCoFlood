function out_shp = build_filtered_coastline_shp(src_shp, domain_shp, out_shp, ...
    min_area_km2, min_diam_m, bbox_buffer_deg, force_rebuild)
%BUILD_FILTERED_COASTLINE_SHP_V2
% 进一步过滤当前研究区附近无法被网格稳定解析的小岛，生成新的 shoreline shp。
%
% 输入：
%   src_shp         原始海岸线 shp
%   domain_shp      外圈计算域 shp（用于限定处理范围）
%   out_shp         输出过滤后的 shp
%   min_area_km2    最小保留岛面积（km^2）
%   min_diam_m      最小保留岛等效直径（m）
%   bbox_buffer_deg 读取海岸线时对 domain bbox 的外扩缓冲（deg）
%   force_rebuild   true=强制重建；false=若已存在则直接复用
%
% 规则：
%   若一个独立岛屿同时满足：
%       面积 < min_area_km2
%       等效直径 < min_diam_m
%   则删除
%
% 加速策略：
%   1) 不用 polyshape
%   2) 大多边形先用 bbox 快速放行
%   3) 只读当前 domain bbox 附近 shp
%   4) 输出若已存在则直接复用

    if nargin < 4 || isempty(min_area_km2)
        min_area_km2 = 1.0;
    end
    if nargin < 5 || isempty(min_diam_m)
        min_diam_m = 1500;
    end
    if nargin < 6 || isempty(bbox_buffer_deg)
        bbox_buffer_deg = 0.15;
    end
    if nargin < 7 || isempty(force_rebuild)
        force_rebuild = false;
    end

    % -------------------------------------------------------------
    % 若已有结果且不强制重建，直接复用
    % -------------------------------------------------------------
    if ~force_rebuild && exist(out_shp, 'file')
        fprintf('   -> 已存在过滤后海岸线，直接复用: %s\n', out_shp);
        return;
    end

    % -------------------------------------------------------------
    % 读取 domain bbox
    % -------------------------------------------------------------
    D = shaperead(domain_shp);
    if isempty(D)
        error('无法读取 domain_shp: %s', domain_shp);
    end

    xall = [];
    yall = [];
    for i = 1:numel(D)
        x = D(i).X(:);
        y = D(i).Y(:);
        v = ~isnan(x) & ~isnan(y);
        xall = [xall; x(v)]; %#ok<AGROW>
        yall = [yall; y(v)]; %#ok<AGROW>
    end

    if isempty(xall)
        error('domain_shp 中没有有效坐标: %s', domain_shp);
    end

    bbox = [ ...
        min(xall) - bbox_buffer_deg, min(yall) - bbox_buffer_deg; ...
        max(xall) + bbox_buffer_deg, max(yall) + bbox_buffer_deg];

    fprintf('   -> 快速过滤小岛...\n');
    fprintf('   -> 原始海岸线: %s\n', src_shp);
    fprintf('   -> 处理范围 bbox = [%.4f %.4f; %.4f %.4f]\n', ...
        bbox(1,1), bbox(1,2), bbox(2,1), bbox(2,2));
    fprintf('   -> 删除阈值: area < %.3f km^2 AND d_eq < %.1f m\n', ...
        min_area_km2, min_diam_m);

    % -------------------------------------------------------------
    % 只读取当前研究区附近的海岸线
    % -------------------------------------------------------------
    S = shaperead(src_shp, 'BoundingBox', bbox);

    if isempty(S)
        warning('在当前研究区附近没有读到海岸线要素：%s', src_shp);
        return;
    end

    % -------------------------------------------------------------
    % 第一遍：筛选保留的 polygon part，先放到 cell 中
    % -------------------------------------------------------------
    keepX = {};
    keepY = {};

    n_drop = 0;
    n_keep = 0;

    for i = 1:numel(S)
        x = S(i).X(:);
        y = S(i).Y(:);

        [istart, iend] = split_parts_index(x, y);

        for j = 1:numel(istart)
            xp = x(istart(j):iend(j));
            yp = y(istart(j):iend(j));

            v = ~isnan(xp) & ~isnan(yp);
            xp = xp(v);
            yp = yp(v);

            if numel(xp) < 3
                continue;
            end

            % 若首尾重复，去掉最后一个点，避免多余计算
            if xp(1) == xp(end) && yp(1) == yp(end) && numel(xp) > 3
                xp = xp(1:end-1);
                yp = yp(1:end-1);
            end

            if numel(xp) < 3
                continue;
            end

            % -----------------------------------------
            % 先用 bbox 尺寸快速判断
            % -----------------------------------------
            latm = mean(yp, 'omitnan');
            cosfac = max(cosd(latm), 0.2);

            width_m  = (max(xp) - min(xp)) * 111320 * cosfac;
            height_m = (max(yp) - min(yp)) * 111320;

            % 明显比阈值大很多的，直接保留，不算面积
            if max(width_m, height_m) >= 1.25 * min_diam_m
                keepX{end+1} = [xp(:).', xp(1), NaN]; %#ok<AGROW>
                keepY{end+1} = [yp(:).', yp(1), NaN]; %#ok<AGROW>
                n_keep = n_keep + 1;
                continue;
            end

            % -----------------------------------------
            % 对可能很小的岛，再用 polyarea 精细判定
            % -----------------------------------------
            area_deg2 = abs(polyarea(xp, yp));
            if ~isfinite(area_deg2) || area_deg2 <= 0
                continue;
            end

            area_km2 = area_deg2 * (111.32^2) * cosfac;
            d_eq_m = 2 * sqrt(area_km2 / pi) * 1000;

            is_tiny_island = (area_km2 < min_area_km2) && (d_eq_m < min_diam_m);

            if is_tiny_island
                n_drop = n_drop + 1;
                continue;
            end

            keepX{end+1} = [xp(:).', xp(1), NaN]; %#ok<AGROW>
            keepY{end+1} = [yp(:).', yp(1), NaN]; %#ok<AGROW>
            n_keep = n_keep + 1;
        end
    end

    if isempty(keepX)
        warning('过滤后没有保留任何海岸线多边形，请检查阈值设置。');
        return;
    end

    % -------------------------------------------------------------
    % 第二遍：一次性组装 struct，避免动态扩容
    % -------------------------------------------------------------
    nk = numel(keepX);
    Sout = repmat(struct( ...
        'Geometry', 'Polygon', ...
        'BoundingBox', zeros(2,2), ...
        'X', [], ...
        'Y', [], ...
        'Id', 0), nk, 1);

    for k = 1:nk
        xb = keepX{k};
        yb = keepY{k};

        v = ~isnan(xb) & ~isnan(yb);
        xv = xb(v);
        yv = yb(v);

        Sout(k).Geometry = 'Polygon';
        Sout(k).BoundingBox = [min(xv), min(yv); max(xv), max(yv)];
        Sout(k).X = xb;
        Sout(k).Y = yb;
        Sout(k).Id = k;
    end

    % -------------------------------------------------------------
    % 输出
    % -------------------------------------------------------------
    out_dir = fileparts(out_shp);
    if ~isempty(out_dir) && ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    shapewrite(Sout, out_shp);

    fprintf('   -> 保留多边形数: %d\n', n_keep);
    fprintf('   -> 删除微小岛屿数: %d\n', n_drop);
    fprintf('   -> 输出过滤后海岸线: %s\n', out_shp);
end
function obj = make_bc_from_capsule_boundary(obj, gdat, varargin)
%MAKE_BC_FROM_CAPSULE_BOUNDARY
% 完全按外部 box/capsule 几何边界来指定 open boundary。
%
% 核心思想：
%   1) 提取 mesh 最外侧 polygon
%   2) 用真实岸线（gdat.mainland + gdat.inner）判断：
%      - 离岸线较远的外边界段 => open boundary
%      - 靠近岸线的外边界段   => outer no-flux boundary
%   3) 可选：再自动补内部 islands 为 inner no-flux boundary
%
% 输入：
%   obj  - msh 对象
%   gdat - geodata 对象（至少要有 mainland / inner 中的一项）
%
% 可选参数：
%   'shore_tol'        : 判定“靠近真实岸线”的阈值（单位：度），默认 0.03
%   'min_open_edges'   : 小于该边数的 open 段会被删除，默认 10
%   'bridge_gap_edges' : open 段之间若只隔很短的 land gap，则自动桥接，默认 3
%   'add_inner_islands': 是否调用 make_bc('inner') 补内部岛屿，默认 true
%   'inner_ibtype'     : inner boundary ibtype，默认 21
%   'plot_check'       : 是否绘图检查，默认 false
%   'verbose'          : 是否打印信息，默认 true
%
% 输出：
%   obj - 已写入 op / bd 的 msh 对象
%
% 说明：
%   1) 本函数不使用 depth 判据
%   2) 只基于几何位置来指定最外侧开边界
%   3) 适合你这种“开边界由 capsule/box 外边界决定”的自动建模流程

    % -------------------------------------------------------------
    % 参数
    % -------------------------------------------------------------
    shore_tol         = 0.03;
    min_open_edges    = 10;
    bridge_gap_edges  = 3;
    add_inner_islands = true;
    inner_ibtype      = 21;
    plot_check        = false;
    verbose           = true;

    for k = 1:2:length(varargin)
        switch lower(varargin{k})
            case 'shore_tol'
                shore_tol = varargin{k+1};
            case 'min_open_edges'
                min_open_edges = varargin{k+1};
            case 'bridge_gap_edges'
                bridge_gap_edges = varargin{k+1};
            case 'add_inner_islands'
                add_inner_islands = varargin{k+1};
            case 'inner_ibtype'
                inner_ibtype = varargin{k+1};
            case 'plot_check'
                plot_check = varargin{k+1};
            case 'verbose'
                verbose = varargin{k+1};
            otherwise
                error('未知参数: %s', varargin{k});
        end
    end

    % -------------------------------------------------------------
    % 基本检查
    % -------------------------------------------------------------
    if nargin < 2 || isempty(gdat)
        error('必须提供 gdat（geodata 对象）');
    end
    if ~isa(gdat, 'geodata')
        error('第二个输入必须是 geodata 对象');
    end

    % -------------------------------------------------------------
    % 清空已有边界
    % -------------------------------------------------------------
    obj.op = [];
    obj.bd = [];

    % -------------------------------------------------------------
    % 收集真实岸线点（只用于判断“靠岸”还是“人工外边界”）
    % -------------------------------------------------------------
    shore_xy = [];

    if isprop(gdat, 'mainland') && ~isempty(gdat.mainland)
        tmp = gdat.mainland;
        tmp = tmp(~any(isnan(tmp),2), :);
        shore_xy = [shore_xy; tmp]; %#ok<AGROW>
    end

    if isprop(gdat, 'inner') && ~isempty(gdat.inner)
        tmp = gdat.inner;
        tmp = tmp(~any(isnan(tmp),2), :);
        shore_xy = [shore_xy; tmp]; %#ok<AGROW>
    end

    if isempty(shore_xy)
        error('gdat.mainland 和 gdat.inner 都为空，无法判断哪些外边界段靠近真实岸线');
    end

    % -------------------------------------------------------------
    % 提取 mesh 最外侧 polygon
    % -------------------------------------------------------------
    [etbv, ~] = extdom_edges2(obj.t, obj.p);
    [poly, poly_idx, max_ind] = extdom_polygon(etbv, obj.p, 1);

    if isempty(poly) || isempty(poly_idx)
        error('无法从 mesh 中提取外边界 polygon');
    end

    idv = poly_idx{max_ind};    % 最大 polygon 对应的节点索引（通常闭合）
    if numel(idv) < 4
        error('外边界 polygon 节点数过少');
    end

    % 确保闭合
    if idv(1) ~= idv(end)
        idv = [idv(:); idv(1)];
    else
        idv = idv(:);
    end

    nEdge = length(idv) - 1;

    % -------------------------------------------------------------
    % 对外边界每条 edge 做分类
    % edge midpoint 离真实岸线较远 => open candidate
    % -------------------------------------------------------------
    p1 = obj.p(idv(1:end-1), :);
    p2 = obj.p(idv(2:end), :);
    pmid = 0.5 * (p1 + p2);

    % 与包内部 make_bc 的距离逻辑保持一致：直接在 lon-lat 上做 KNN
    [~, dshore] = ourKNNsearch(shore_xy', pmid', 1);
    dshore = dshore(:);

    is_open_edge = dshore > shore_tol;

    % -------------------------------------------------------------
    % 1D 圆环逻辑平滑
    % 1) 填补很短的 gap
    % 2) 删除很短的 open 段
    % -------------------------------------------------------------
    if bridge_gap_edges > 0
        is_open_edge = bridge_short_false_gaps_circular(is_open_edge, bridge_gap_edges);
    end
    if min_open_edges > 1
        is_open_edge = remove_short_true_runs_circular(is_open_edge, min_open_edges);
    end

    if verbose
        fprintf('\n[make_bc_from_capsule_boundary]\n');
        fprintf('   外边界 edge 总数      : %d\n', nEdge);
        fprintf('   shore_tol            : %.4f deg\n', shore_tol);
        fprintf('   min_open_edges       : %d\n', min_open_edges);
        fprintf('   bridge_gap_edges     : %d\n', bridge_gap_edges);
        fprintf('   判定为 open 的 edge 数: %d\n', nnz(is_open_edge));
        fprintf('   判定为 land 的 edge 数: %d\n', nnz(~is_open_edge));
    end

    if ~any(is_open_edge)
        warning('没有任何外边界 edge 被判为 open boundary，请增大 shore_tol 或检查 gdat.mainland/inner');
    end

    % -------------------------------------------------------------
    % 将 circular edge 序列拆成 open / land 的连续段
    % -------------------------------------------------------------
    runs = circular_runs_logical(is_open_edge);

    open_segments = {};
    land_segments = {};

    for i = 1:numel(runs)
        idx_edge = circular_index_range(runs(i).start_idx, runs(i).end_idx, nEdge);
        node_seq = edge_run_to_node_sequence(idv, idx_edge);

        if numel(node_seq) < 2
            continue;
        end

        if runs(i).value
            open_segments{end+1} = node_seq(:); %#ok<AGROW>
        else
            land_segments{end+1} = node_seq(:); %#ok<AGROW>
        end
    end

    % -------------------------------------------------------------
    % 写 open boundary
    % ADCIRC / OceanMesh2D 里 open boundary 节点放在 obj.op
    % -------------------------------------------------------------
    if ~isempty(open_segments)
        nope = numel(open_segments);
        nvdll = zeros(1, nope);
        neta = 0;

        for i = 1:nope
            nvdll(i) = numel(open_segments{i});
            neta = neta + nvdll(i);
        end

        maxn = max(nvdll);
        nbdv = zeros(maxn, nope);

        for i = 1:nope
            nbdv(1:nvdll(i), i) = open_segments{i};
        end

        obj.op.nope   = nope;
        obj.op.neta   = neta;
        obj.op.nvdll  = nvdll;
        obj.op.ibtype = zeros(1, nope);   % open boundary
        obj.op.nbdv   = nbdv;
    else
        obj.op = [];
    end

    % -------------------------------------------------------------
    % 写 outer no-flux boundary
    % 外边界中靠岸的那部分，统一设为 ibtype = 20
    % -------------------------------------------------------------
    if ~isempty(land_segments)
        nbou = numel(land_segments);
        nvell = zeros(1, nbou);
        nvel = 0;

        for i = 1:nbou
            nvell(i) = numel(land_segments{i});
            nvel = nvel + nvell(i);
        end

        maxn = max(nvell);
        nbvv = zeros(maxn, nbou);

        for i = 1:nbou
            nbvv(1:nvell(i), i) = land_segments{i};
        end

        obj.bd.nbou   = nbou;
        obj.bd.nvel   = nvel;
        obj.bd.nvell  = nvell;
        obj.bd.ibtype = 20 * ones(1, nbou);   % outer no-flux
        obj.bd.nbvv   = nbvv;
    else
        obj.bd = [];
    end

    % -------------------------------------------------------------
    % 可选：补内部岛屿 / 内陆孔洞
    % -------------------------------------------------------------
    if add_inner_islands
        try
            obj = make_bc(obj, 'inner', inner_ibtype);
            if verbose
                fprintf('   已追加内部 inner no-flux boundaries (ibtype=%d)\n', inner_ibtype);
            end
        catch ME
            warning('调用 make_bc(obj,''inner'',...) 失败：%s', ME.message);
        end
    end

    % -------------------------------------------------------------
    % 检查绘图
    % -------------------------------------------------------------
    if plot_check
        figure('Color','w');
        plot(obj, 'type', 'bd', 'proj', 'none');
        title('Boundary-condition check');
        drawnow;
    end
end


% =========================================================================
% 将 circular logical 序列拆成 runs
% =========================================================================
function runs = circular_runs_logical(x)

    x = logical(x(:)');
    n = numel(x);

    runs = struct('value', {}, 'start_idx', {}, 'end_idx', {}, 'len', {});

    if n == 0
        return;
    end

    if all(x == x(1))
        runs(1).value = x(1);
        runs(1).start_idx = 1;
        runs(1).end_idx = n;
        runs(1).len = n;
        return;
    end

    prev = x([end, 1:end-1]);
    starts = find(x ~= prev);

    for i = 1:numel(starts)
        s = starts(i);
        if i < numel(starts)
            e = starts(i+1) - 1;
        else
            e = starts(1) - 1;
            if e <= 0
                e = e + n;
            end
        end

        idx = circular_index_range(s, e, n);

        runs(i).value = x(s);
        runs(i).start_idx = s;
        runs(i).end_idx = e;
        runs(i).len = numel(idx);
    end
end


% =========================================================================
% 圆环索引范围
% 例如 n=10, s=8, e=3 => [8 9 10 1 2 3]
% =========================================================================
function idx = circular_index_range(s, e, n)

    if s <= e
        idx = s:e;
    else
        idx = [s:n, 1:e];
    end
end


% =========================================================================
% 填补很短的 false gap（圆环）
% =========================================================================
function x = bridge_short_false_gaps_circular(x, max_gap)

    x = logical(x(:)');
    if isempty(x) || all(x == 0) || all(x == 1)
        return;
    end

    runs = circular_runs_logical(x);

    if numel(runs) < 3
        return;
    end

    for i = 1:numel(runs)
        if runs(i).value == 0 && runs(i).len <= max_gap
            il = i - 1; if il < 1, il = numel(runs); end
            ir = i + 1; if ir > numel(runs), ir = 1; end

            if runs(il).value == 1 && runs(ir).value == 1
                idx = circular_index_range(runs(i).start_idx, runs(i).end_idx, numel(x));
                x(idx) = true;
            end
        end
    end
end


% =========================================================================
% 删除很短的 true run（圆环）
% =========================================================================
function x = remove_short_true_runs_circular(x, min_len)

    x = logical(x(:)');
    if isempty(x) || all(x == 0) || all(x == 1)
        if all(x == 1) && numel(x) < min_len
            x(:) = false;
        end
        return;
    end

    runs = circular_runs_logical(x);

    for i = 1:numel(runs)
        if runs(i).value == 1 && runs(i).len < min_len
            idx = circular_index_range(runs(i).start_idx, runs(i).end_idx, numel(x));
            x(idx) = false;
        end
    end
end


% =========================================================================
% 将一段连续 edge 序列转成节点序列
% idv 为闭合 boundary node 索引，如 [1 2 3 4 1]
% edge_idx 是 edge 编号，edge i 对应 idv(i) -> idv(i+1)
% =========================================================================
function node_seq = edge_run_to_node_sequence(idv, edge_idx)

    nEdge = length(idv) - 1;
    if isempty(edge_idx)
        node_seq = [];
        return;
    end

    s = edge_idx(1);
    e = edge_idx(end);

    if s <= e && numel(edge_idx) == (e - s + 1)
        pos = s:(e+1);
    else
        % wrap around
        pos = [s:(nEdge+1), 2:(e+1)];
    end

    node_seq = idv(pos);
end
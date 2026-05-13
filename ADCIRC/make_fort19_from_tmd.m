function make_fort19_from_tmd(auto_Name)
% ============================================================
% 自动生成 ADCIRC fort.19 潮汐边界驱动文件
%
% 功能：
% 1) 从 fort.14 读取开边界节点及顺序
% 2) 对 cases_<auto_Name> 下每个事件目录：
%       - 读取 fort22_meta.txt
%       - 生成 TMD lat_lon 输入
%       - 调用 TMD_ato
%       - 读取 TMD 输出 data_*.mat
%       - 写出 fort.19
%
% 说明：
% - 适用于任意 ADCIRC 网格
% - fort.19 的时间步长直接采用 fort22_meta.txt 中的 WTIMINC
% - fort.19 第一行写入 ETIMINC（秒）
% - 后续按 ADCIRC 要求写：每个时刻依次写所有开边界节点潮位
%
% 你需要确认：
% - TMD_ato 可以在 MATLAB 中直接调用
% - TMD 目录结构与下面参数一致
% - TMD 输出为 data_1.mat, data_2.mat, ...，每个文件内含 TimeSeries
% ============================================================



    % 事件目录根目录
    case_root = ['cases_' auto_Name];

    % fort.14 自动搜索：
    % 1) ADCIRC_Capsule_Mesh_<auto_Name>.14
    % 2) fort.14
    % 3) 当前目录唯一 .14
    fort14_file = resolve_fort14_file(auto_Name);

    % TMD 根目录
    tmd_root = 'D:\MATLAB2021a\bin\m\ADCIRC\TMD2.5\TMD';

    % TMD lat_lon 输入文件
    latlon_file = fullfile(tmd_root, 'LAT_LON', 'lat_lon');

    % TMD 输出目录（按你旧脚本，主要 data_*.mat 在 tmd_root）
    tmd_out_dir = fullfile(tmd_root, 'OUT');

    % 是否每个事件运行前清理 TMD 临时输出
    clean_tmd_before_each_case = true;

    % 若 TMD 输出长度和目标步数不一致：
    % true  -> 自动裁剪/补齐
    % false -> 直接报错
    allow_trim_or_pad = true;

    %% =========================
    % 2) 基本检查
    %% =========================
    if ~exist(case_root, 'dir')
        error('未找到事件目录: %s', case_root);
    end

    if ~exist(tmd_root, 'dir')
        error('未找到 TMD 根目录: %s', tmd_root);
    end

    if ~exist(fullfile(tmd_root, 'LAT_LON'), 'dir')
        error('未找到 TMD\\LAT_LON 目录: %s', fullfile(tmd_root, 'LAT_LON'));
    end

    %% =========================
    % 3) 读取 fort.14 开边界
    %% =========================
    fprintf('读取 fort.14: %s\n', fort14_file);
    G = read_fort14_open_boundary(fort14_file);

    fprintf('总节点数 NP      = %d\n', G.NP);
    fprintf('总单元数 NE      = %d\n', G.NE);
    fprintf('开边界条数 NOPE  = %d\n', G.NOPE);
    fprintf('开边界节点总数   = %d\n', numel(G.open_nodes));

    if isempty(G.open_nodes)
        error('fort.14 中未读取到开边界节点。');
    end

    %% =========================
    % 4) 搜索事件目录
    %% =========================
    D = dir(case_root);
    D = D([D.isdir]);
    D = D(~ismember({D.name}, {'.','..'}));

    if isempty(D)
        error('在 %s 下未找到事件目录。', case_root);
    end

    fprintf('找到 %d 个事件目录。\n', numel(D));

    %% =========================
    % 5) 逐个事件生成 fort.19
    %% =========================
    n_ok = 0;
    n_skip = 0;

    for icase = 1:numel(D)
        case_dir = fullfile(case_root, D(icase).name);
        meta_file = fullfile(case_dir, 'fort22_meta.txt');
        fort19_file = fullfile(case_dir, 'fort.19');

        fprintf('\n====================================================\n');
        fprintf('处理事件目录: %s\n', case_dir);

        if ~exist(meta_file, 'file')
            fprintf('  -> 跳过：未找到 fort22_meta.txt\n');
            n_skip = n_skip + 1;
            continue;
        end

        try
            meta = read_fort22_meta(meta_file);

            WindowStart = get_required_datetime(meta, 'WindowStart');
            WindowEnd   = get_required_datetime(meta, 'WindowEnd');
            WTIMINC     = get_required_numeric(meta, 'WTIMINC');   % 秒

            dt_min = WTIMINC / 60;
            if abs(dt_min - round(dt_min)) > 1e-8
                warning('WTIMINC 不是整分钟，已四舍五入到最近整数分钟。');
            end
            dt_min = round(dt_min);

            % 与 fort22 保持一致：起止时间都包含
            nSteps = round(seconds(WindowEnd - WindowStart) / WTIMINC) + 1;
            if nSteps <= 0
                error('无效的时间窗口：WindowEnd <= WindowStart');
            end

            fprintf('  WindowStart = %s\n', datestr(WindowStart, 31));
            fprintf('  WindowEnd   = %s\n', datestr(WindowEnd, 31));
            fprintf('  WTIMINC     = %d s\n', WTIMINC);
            fprintf('  dt_min      = %d min\n', dt_min);
            fprintf('  nSteps      = %d\n', nSteps);

            % ---------- 写 lat_lon ----------
            write_tmd_latlon_file(latlon_file, G.open_lat, G.open_lon, WindowStart, dt_min, nSteps);
            fprintf('  -> 已写 TMD 输入: %s\n', latlon_file);

            % ---------- 清理旧输出 ----------
            if clean_tmd_before_each_case
                cleanup_tmd_outputs(tmd_root, tmd_out_dir);
            end

            % 某些 TMD 脚本依赖 data.out 预存在
            ensure_blank_data_out(fullfile(tmd_root, 'data.out'));

            % ---------- 调用 TMD ----------
            run_tmd_ato(tmd_root);

            % ---------- 读取 TMD 输出 ----------
            tide = read_tmd_series_matrix(tmd_root, numel(G.open_nodes), nSteps, allow_trim_or_pad);

            % ---------- 写 fort.19 ----------
            write_fort19_file(fort19_file, tide, WTIMINC);

            fprintf('  -> 已生成 fort.19: %s\n', fort19_file);

            n_ok = n_ok + 1;

        catch ME
            fprintf('  -> 失败: %s\n', ME.message);
            n_skip = n_skip + 1;
        end
    end

    fprintf('\n====================================================\n');
    fprintf('fort.19 全部处理完成\n');
    fprintf('成功数量: %d\n', n_ok);
    fprintf('失败/跳过: %d\n', n_skip);
    fprintf('====================================================\n');
end

%% ============================================================
% 读取 fort.14：节点 + 开边界顺序
%% ============================================================
function G = read_fort14_open_boundary(fort14_file)

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
        error('fort.14 第 2 行解析失败。');
    end
    G.NE = a(1);
    G.NP = a(2);

    % ---------- 节点 ----------
    lon   = zeros(G.NP,1);
    lat   = zeros(G.NP,1);
    depth = zeros(G.NP,1);

    for i = 1:G.NP
        s = fgetl(fid);
        v = sscanf(s, '%f');
        if numel(v) < 4
            error('第 %d 个节点读取失败。', i);
        end
        lon(i)   = v(2);
        lat(i)   = v(3);
        depth(i) = v(4);
    end

    G.lon   = lon;
    G.lat   = lat;
    G.depth = depth;

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

    G.open_nodes = open_nodes(:);
    G.open_lon   = lon(G.open_nodes);
    G.open_lat   = lat(G.open_nodes);
end

%% ============================================================
% 写 TMD 的 lat_lon 输入文件
% 格式沿用你旧脚本：
% lat lon yy mm dd hh mi sec dt(min) TSLength
%% ============================================================
function write_tmd_latlon_file(latlon_file, open_lat, open_lon, start_time, dt_min, nSteps)

    fid = fopen(latlon_file, 'wt');
    if fid < 0
        error('无法写入 TMD 输入文件: %s', latlon_file);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, '   lat   lon            yy   mm  dd hh mi sec dt(min) TSLength\n');

    [yy, mm, dd] = ymd(start_time);
    hh = hour(start_time);
    mi = minute(start_time);
    ss = second(start_time);

    for i = 1:numel(open_lat)
        fprintf(fid, '  %.10f  %.10f', open_lat(i), open_lon(i));
        fprintf(fid, '    %04d %02d %02d %02d %02d %02d %d %d\n', ...
            yy, mm, dd, hh, mi, ss, dt_min, nSteps);
    end
end

%% ============================================================
% 清理 TMD 旧输出
%% ============================================================
function cleanup_tmd_outputs(tmd_root, tmd_out_dir)

    delete_if_exists(fullfile(tmd_root, 'data_*.mat'));
    delete_if_exists(fullfile(tmd_root, 'data.out'));
    delete_if_exists(fullfile(tmd_root, 'data*'));

    if exist(tmd_out_dir, 'dir')
        delete_if_exists(fullfile(tmd_out_dir, 'data.out'));
        delete_if_exists(fullfile(tmd_out_dir, 'data*'));
    end
end

function delete_if_exists(pattern)
    D = dir(pattern);
    for i = 1:numel(D)
        try
            delete(fullfile(D(i).folder, D(i).name));
        catch
        end
    end
end

%% ============================================================
% 确保 data.out 是一个空白文件
%% ============================================================
function ensure_blank_data_out(data_out_file)
    fid = fopen(data_out_file, 'wt');
    if fid < 0
        error('无法写入 %s', data_out_file);
    end
    fprintf(fid, '%s\n', ' ');
    fclose(fid);
end

%% ============================================================
% 调用 TMD_ato
%% ============================================================
function run_tmd_ato(tmd_root)

    old_dir = pwd;
    cleaner = onCleanup(@() cd(old_dir)); %#ok<NASGU>
    cd(tmd_root);

    if exist('TMD_ato', 'file') ~= 2 && exist('TMD_ato', 'builtin') ~= 5
        error('未找到 TMD_ato，请确认它在 MATLAB 路径上。');
    end

    fprintf('  -> 运行 TMD_ato ...\n');
    TMD_ato;
end

%% ============================================================
% 读取 TMD 输出 data_1.mat, data_2.mat, ...
% 每个文件表示一个开边界节点的时间序列
% 输出 tide: [nOpen, nSteps]
%% ============================================================
function tide = read_tmd_series_matrix(tmd_root, nOpen, nSteps, allow_trim_or_pad)

    tide = nan(nOpen, nSteps);

    for i = 1:nOpen
        f = fullfile(tmd_root, sprintf('data_%d.mat', i));
        if ~exist(f, 'file')
            error('缺少 TMD 输出文件: %s', f);
        end

        S = load(f);
        if ~isfield(S, 'TimeSeries')
            error('%s 中未找到变量 TimeSeries', f);
        end

        ts = S.TimeSeries(:);
        ts = real(ts);

        if numel(ts) == nSteps
            tide(i,:) = ts.';
        elseif allow_trim_or_pad
            if numel(ts) > nSteps
                tide(i,:) = ts(1:nSteps).';
            else
                tmp = nan(nSteps,1);
                tmp(1:numel(ts)) = ts;
                if ~isempty(ts)
                    tmp(numel(ts)+1:end) = ts(end);
                else
                    tmp(:) = 0;
                end
                tide(i,:) = tmp.';
            end
        else
            error('TMD 输出长度与目标步数不一致：data_%d.mat, got=%d, need=%d', i, numel(ts), nSteps);
        end
    end
end

%% ============================================================
% 写 fort.19
% 第一行：ETIMINC（秒）
% 后续：每个时刻依次写所有开边界节点潮位
%% ============================================================
function write_fort19_file(fort19_file, tide, ETIMINC)

    [nOpen, nSteps] = size(tide);

    fid = fopen(fort19_file, 'wt');
    if fid < 0
        error('无法写入 fort.19: %s', fort19_file);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, '%d\n', ETIMINC);

    for k = 1:nSteps
        for j = 1:nOpen
            fprintf(fid, '%.10f\n', tide(j,k));
        end
    end
end

%% ============================================================
% 读取 fort22_meta.txt
%% ============================================================
function meta = read_fort22_meta(meta_file)
    lines = read_text_lines(meta_file);
    meta = struct();

    for i = 1:numel(lines)
        s = strtrim(lines{i});
        if isempty(s)
            continue;
        end
        if startsWith(s, 'ADCIRC', 'IgnoreCase', true) || startsWith(s, '---')
            continue;
        end

        eq_pos = strfind(s, '=');
        if isempty(eq_pos)
            continue;
        end

        k = strtrim(s(1:eq_pos(1)-1));
        v = strtrim(s(eq_pos(1)+1:end));
        meta.(k) = v;
    end
end

%% ============================================================
% 读取文本为 cell lines
%% ============================================================
function lines = read_text_lines(fname)
    fid = fopen(fname, 'rt');
    if fid < 0
        error('无法打开文件: %s', fname);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    lines = {};
    while true
        tline = fgetl(fid);
        if ~ischar(tline)
            break;
        end
        lines{end+1,1} = tline; %#ok<AGROW>
    end
end

%% ============================================================
% 必需数值字段
%% ============================================================
function val = get_required_numeric(meta, key)
    if ~isfield(meta, key)
        error('meta 中缺少字段: %s', key);
    end
    val = str2double(strtrim(meta.(key)));
    if ~isfinite(val)
        error('字段 %s 不是有效数值: %s', key, meta.(key));
    end
end

%% ============================================================
% 必需时间字段
%% ============================================================
function dt = get_required_datetime(meta, key)
    if ~isfield(meta, key)
        error('meta 中缺少时间字段: %s', key);
    end

    s = strtrim(meta.(key));
    s = strrep(s, 'T', ' ');

    try
        dt = datetime(s, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    catch
        try
            dt = datetime(s);
        catch
            error('字段 %s 不能解析为 datetime: %s', key, s);
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
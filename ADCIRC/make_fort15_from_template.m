function make_fort15_from_template(auto_Name)
% ===============================================================
% 根据模板 fort.15 + 每个事件目录下的 fort22_meta.txt
% 自动为每个事件生成 fort.15
%
% 用法：
% 1) 修改下面这几个路径/参数
% 2) 运行本脚本
%
% 说明：
% - 会在每个 case 目录下生成一个 fort.15
% - 模板中的物理参数、站点输出等保持不变
% - 仅自动改写与事件和 fort.22 强迫相关的部分
% ===============================================================

    %% =========================
    % 1) 用户参数
    %% =========================

    % 你的模板 fort.15
    template_fort15 = 'fort.15';

    % 前一步生成的事件目录根目录
    output_root = ['cases_' auto_Name];

    % 是否强制把 NWS 改成 6
    force_NWS6 = true;

    % 是否把所有输出行的结束时间 TOUTF* 改成 RNDAY
    sync_output_end_to_rnday = 0;

    % DRAMP 是否自动更新为 min(1.0, RNDAY)
    auto_update_dramp = true;

    %% =========================
    % 2) 基本检查
    %% =========================
    if ~exist(template_fort15, 'file')
        error('模板文件不存在: %s', template_fort15);
    end

    if ~exist(output_root, 'dir')
        error('事件目录不存在: %s', output_root);
    end

    case_dirs = dir(output_root);
    case_dirs = case_dirs([case_dirs.isdir]);
    case_dirs = case_dirs(~ismember({case_dirs.name}, {'.','..'}));

    if isempty(case_dirs)
        error('在 %s 下没有找到任何事件子目录。', output_root);
    end

    fprintf('找到 %d 个事件目录。\n', numel(case_dirs));

    %% =========================
    % 3) 读取模板
    %% =========================
    template_lines = read_text_lines(template_fort15);

    %% =========================
    % 4) 逐个事件生成 fort.15
    %% =========================
    n_ok = 0;
    n_skip = 0;

    for i = 1:numel(case_dirs)
        case_name = case_dirs(i).name;
        case_dir  = fullfile(output_root, case_name);

        meta_file = fullfile(case_dir, 'fort22_meta.txt');
        out_file  = fullfile(case_dir, 'fort.15');

        fprintf('\n--------------------------------------------------\n');
        fprintf('处理事件目录: %s\n', case_dir);

        if ~exist(meta_file, 'file')
            fprintf('  -> 跳过：未找到 fort22_meta.txt\n');
            n_skip = n_skip + 1;
            continue;
        end

        try
            meta = read_fort22_meta(meta_file);

            lines = template_lines;

            % ---------- 提取关键参数 ----------
            NWLAT   = get_required_numeric(meta, 'NWLAT');
            NWLON   = get_required_numeric(meta, 'NWLON');
            WLATMAX = get_required_numeric(meta, 'WLATMAX');
            WLONMIN = get_required_numeric(meta, 'WLONMIN');
            WLATINC = get_required_numeric(meta, 'WLATINC');
            WLONINC = get_required_numeric(meta, 'WLONINC');
            WTIMINC = get_required_numeric(meta, 'WTIMINC');

            WindowStart = get_required_datetime(meta, 'WindowStart');
            WindowEnd   = get_required_datetime(meta, 'WindowEnd');

            % RNDAY = 结束时间 - 开始时间（单位：天）
            RNDAY = days(WindowEnd - WindowStart);
            if RNDAY <= 0
                error('RNDAY <= 0，请检查 %s', meta_file);
            end

            if isfield(meta, 'tc_id')
                tc_id = strtrim(meta.tc_id);
            else
                tc_id = case_name;
            end

            if isfield(meta, 'Cv_used')
                cv_used = strtrim(meta.Cv_used);
            else
                cv_used = '';
            end

            % ---------- 1) 改写第1行：运行描述 ----------
            run_desc = sprintf('%s %s', auto_Name, tc_id);
            run_desc = crop_or_pad(run_desc, 32);
            lines{1} = sprintf(' %-32s ! 32 CHARACTER ALPHANUMERIC RUN DESCRIPTION', run_desc);

            % ---------- 2) 可选：改写第2行运行标识 ----------
            % 这里保留原模板第二行，不动
            % lines{2} = ...

            % ---------- 3) 强制 NWS = 6 ----------
            idx_nws = find_line_contains(lines, '! NWS - WIND STRESS AND BAROMETRIC PRESSURE OPTION PARAMETER');
            if ~isempty(idx_nws) && force_NWS6
                lines{idx_nws} = ' 6                                   ! NWS - WIND STRESS AND BAROMETRIC PRESSURE OPTION PARAMETER';
            end

            % ---------- 4) 改写气象网格参数行 ----------
            % 模板这行虽然注释写得不太标准，但数值顺序应为：
            % NWLAT NWLON WLATMAX WLONMIN WLATINC WLONINC WTIMINC
            idx_met = find_line_any(lines, {'WTIMINC', 'WITMINC', 'STIMINC'});
            if isempty(idx_met)
                error('模板中未找到气象网格参数行（含 WTIMINC/WITMINC/STIMINC 的行）');
            end

            lines{idx_met} = sprintf( ...
                ' %d %d %.2f %.2f %.2f %.2f %d    ! NWLAT NWLON WLATMAX WLONMIN WLATINC WLONINC WTIMINC', ...
                NWLAT, NWLON, WLATMAX, WLONMIN, WLATINC, WLONINC, WTIMINC);

            % ---------- 5) 改写 RNDAY ----------
            idx_rnday = find_line_contains(lines, '! RNDAY - TOTAL LENGTH OF SIMULATION (IN DAYS)');
            if isempty(idx_rnday)
                error('模板中未找到 RNDAY 行');
            end
            lines{idx_rnday} = sprintf(' %.8f                                 ! RNDAY - TOTAL LENGTH OF SIMULATION (IN DAYS)', RNDAY);

            % ---------- 6) 改写 DRAMP ----------
            if auto_update_dramp
                idx_dramp = find_line_contains(lines, '! DRAMP - DURATION OF RAMP FUNCTION (IN DAYS)');
                if ~isempty(idx_dramp)
                    DRAMP = min(1.0, RNDAY);
                    lines{idx_dramp} = sprintf(' %.8f                                 ! DRAMP - DURATION OF RAMP FUNCTION (IN DAYS)', DRAMP);
                end
            end

            % ---------- 7) 同步各输出行的 TOUTF 到 RNDAY ----------
            if sync_output_end_to_rnday
                output_tokens = { ...
                    'NOUTE,TOUTSE,TOUTFE', ...
                    'NOUTV,TOUTSV,TOUTFV', ...
                    'NOUTM,TOUTSM,TOUTFM', ...
                    'NOUTGE,TOUTSGE,TOUTFGE', ...
                    'NOUTGV,TOUTSGV,TOUTFGV', ...
                    'NOUTGW,TOUTSGW,TOUTFGW'};

                for k = 1:numel(output_tokens)
                    idx_out = find_line_contains(lines, output_tokens{k});
                    if ~isempty(idx_out)
                        lines{idx_out} = update_output_schedule_line(lines{idx_out}, RNDAY);
                    end
                end
            end

            % ---------- 8) 可选：在第2行后插入注释信息 ----------
            % 不插入，避免破坏 fort.15 结构。只在文件末尾追加一行注释。
            if ~isempty(cv_used)
                lines{end+1} = sprintf('! auto-generated from fort22_meta.txt ; tc_id=%s ; Cv_used=%s', tc_id, cv_used);
            else
                lines{end+1} = sprintf('! auto-generated from fort22_meta.txt ; tc_id=%s', tc_id);
            end

            % ---------- 9) 写出 fort.15 ----------
            write_text_lines(out_file, lines);

            fprintf('  -> 已生成: %s\n', out_file);
            fprintf('     RNDAY   = %.2f\n', RNDAY);
            fprintf('     NWLAT   = %d\n', NWLAT);
            fprintf('     NWLON   = %d\n', NWLON);
            fprintf('     WLATMAX = %.2f\n', WLATMAX);
            fprintf('     WLONMIN = %.2f\n', WLONMIN);
            fprintf('     WLATINC = %.2f\n', WLATINC);
            fprintf('     WLONINC = %.2f\n', WLONINC);
            fprintf('     WTIMINC = %d\n', WTIMINC);

            n_ok = n_ok + 1;

        catch ME
            fprintf('  -> 失败: %s\n', ME.message);
            n_skip = n_skip + 1;
        end
    end

    fprintf('\n==================================================\n');
    fprintf('完成。\n');
    fprintf('成功生成 fort.15 数量: %d\n', n_ok);
    fprintf('跳过/失败数量       : %d\n', n_skip);
    fprintf('==================================================\n');
end

%% ===============================================================
% 工具函数：读取纯文本为 cell lines
%% ===============================================================
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

%% ===============================================================
% 工具函数：写回文本
%% ===============================================================
function write_text_lines(fname, lines)
    fid = fopen(fname, 'wt');
    if fid < 0
        error('无法写入文件: %s', fname);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    for i = 1:numel(lines)
        fprintf(fid, '%s\n', lines{i});
    end
end

%% ===============================================================
% 读取 fort22_meta.txt
% 解析形如:
%   key = value
%% ===============================================================
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

        % 保存原始字符串
        meta.(k) = v;
    end
end

%% ===============================================================
% 获取必须存在的数值字段
%% ===============================================================
function val = get_required_numeric(meta, key)
    if ~isfield(meta, key)
        error('meta 中缺少字段: %s', key);
    end
    val = str2double(strtrim(meta.(key)));
    if ~isfinite(val)
        error('字段 %s 不是有效数值: %s', key, meta.(key));
    end
end

%% ===============================================================
% 获取必须存在的时间字段
% 支持格式：
%   2026-04-02 12:00:00
%   2026-04-02T12:00:00
%% ===============================================================
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

%% ===============================================================
% 查找包含某个 token 的行
%% ===============================================================
function idx = find_line_contains(lines, token)
    idx = [];
    for i = 1:numel(lines)
        if contains(lines{i}, token)
            idx = i;
            return;
        end
    end
end

%% ===============================================================
% 查找包含多个候选 token 之一的行
%% ===============================================================
function idx = find_line_any(lines, tokens)
    idx = [];
    for i = 1:numel(lines)
        for j = 1:numel(tokens)
            if contains(lines{i}, tokens{j})
                idx = i;
                return;
            end
        end
    end
end

%% ===============================================================
% 截断/补空格到固定长度
%% ===============================================================
function s = crop_or_pad(s, n)
    s = char(s);
    if numel(s) > n
        s = s(1:n);
    elseif numel(s) < n
        s = [s repmat(' ', 1, n-numel(s))];
    end
end

%% ===============================================================
% 更新输出调度行：
% 输入形如：
%   1 0.0 4.0 180  ! NOUTE,TOUTSE,TOUTFE,NSPOOLE:...
%
% 规则：
% - 保留 NOUT 和 NSPOOL
% - TOUTF 改为 RNDAY
% - 如果 TOUTS > RNDAY，则令 TOUTS = max(0, RNDAY - 1/24)
%% ===============================================================
function new_line = update_output_schedule_line(line, rnday)
    % 拆分注释
    excl = strfind(line, '!');
    if isempty(excl)
        data_part = strtrim(line);
        comment_part = '';
    else
        data_part = strtrim(line(1:excl(1)-1));
        comment_part = line(excl(1):end);
    end

    nums = sscanf(data_part, '%f');

    if numel(nums) < 4
        % 不是标准输出调度行，原样返回
        new_line = line;
        return;
    end

    NOUT   = nums(1);
    TOUTS  = nums(2);
    TOUTF  = nums(3);
    NSPOOL = nums(4);

    %#ok<NASGU>
    TOUTF = rnday;
    if TOUTS > rnday
        TOUTS = max(0, rnday - 1/24);
    end

    new_line = sprintf(' %g %.2f %.2f %g %s', NOUT, TOUTS, TOUTF, NSPOOL, comment_part);
end
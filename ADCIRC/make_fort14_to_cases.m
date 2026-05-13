function make_fort14_to_cases(auto_Name)
% ============================================================
% 把主文件夹下的 fort14 网格文件 和 sub_intel.sh
% 复制到每一个子文件夹下：
%   - fort14 统一命名为 fort.14
%   - sub_intel.sh 保持原名
%
% 用法：
%   make_fort14_to_cases('PRD')
%
% 说明：
% - 主文件夹默认是当前工作目录
% - 子文件夹默认是 cases_<auto_Name>
% - 源 fort14 文件优先级：
%   1) ADCIRC_Capsule_Mesh_<auto_Name>.14
%   2) fort.14
%   3) 当前目录下唯一一个 .14 文件
% - sub_intel.sh 默认要求位于主文件夹下
% ============================================================

    if nargin < 1 || isempty(auto_Name)
        auto_Name = 'PRD';
    end

    root_dir  = pwd;
    case_root = fullfile(root_dir, ['cases_' auto_Name]);

    fprintf('=============================================\n');
    fprintf('主目录         : %s\n', root_dir);
    fprintf('事件目录根     : %s\n', case_root);

    if ~exist(case_root, 'dir')
        error('未找到事件目录根文件夹: %s', case_root);
    end

    % ---------- 自动寻找源 fort14 ----------
    src_fort14 = resolve_source_fort14(root_dir, auto_Name);
    fprintf('源 fort14      : %s\n', src_fort14);

    % ---------- 自动寻找 sub_intel.sh ----------
    src_subsh = resolve_source_subsh(root_dir);
    fprintf('源 sub_intel.sh: %s\n', src_subsh);

    % ---------- 找到所有子文件夹 ----------
    D = dir(case_root);
    D = D([D.isdir]);
    D = D(~ismember({D.name}, {'.', '..'}));

    if isempty(D)
        warning('在 %s 下没有找到任何子文件夹。', case_root);
        return;
    end

    fprintf('待复制子文件夹数量: %d\n', numel(D));
    fprintf('=============================================\n');

    n_ok = 0;
    n_fail = 0;

    for i = 1:numel(D)
        sub_dir     = fullfile(case_root, D(i).name);
        dst_fort14  = fullfile(sub_dir, 'fort.14');
        dst_subsh   = fullfile(sub_dir, 'sub_intel.sh');

        fprintf('[%d/%d] 处理: %s\n', i, numel(D), sub_dir);

        try
            % 复制 fort14
            copyfile(src_fort14, dst_fort14, 'f');

            % 复制 sub_intel.sh
            copyfile(src_subsh, dst_subsh, 'f');

            fprintf('         已复制 fort.14      -> %s\n', dst_fort14);
            fprintf('         已复制 sub_intel.sh -> %s\n', dst_subsh);

            n_ok = n_ok + 1;

        catch ME
            fprintf('         失败 -> %s\n', sub_dir);
            fprintf('         原因: %s\n', ME.message);
            n_fail = n_fail + 1;
        end
    end

    fprintf('=============================================\n');
    fprintf('完成。\n');
    fprintf('成功: %d\n', n_ok);
    fprintf('失败: %d\n', n_fail);
    fprintf('=============================================\n');
end

%% ============================================================
% 自动寻找源 fort14
%% ============================================================
function src_fort14 = resolve_source_fort14(root_dir, auto_Name)

    cand1 = fullfile(root_dir, ['ADCIRC_Capsule_Mesh_' auto_Name '.14']);
    cand2 = fullfile(root_dir, 'fort.14');

    if exist(cand1, 'file')
        src_fort14 = cand1;
        return;
    end

    if exist(cand2, 'file')
        src_fort14 = cand2;
        return;
    end

    D = dir(fullfile(root_dir, '*.14'));

    if isempty(D)
        error('主目录下未找到任何 .14 文件。');
    elseif numel(D) == 1
        src_fort14 = fullfile(root_dir, D(1).name);
        return;
    else
        names = strjoin({D.name}, '\n');
        error(['主目录下存在多个 .14 文件，无法自动判断，请保留一个或使用标准命名。' newline '%s'], names);
    end
end

%% ============================================================
% 自动寻找源 sub_intel.sh
%% ============================================================
function src_subsh = resolve_source_subsh(root_dir)

    cand = fullfile(root_dir, 'sub_intel.sh');

    if exist(cand, 'file')
        src_subsh = cand;
    else
        error('主目录下未找到 sub_intel.sh 文件: %s', cand);
    end
end
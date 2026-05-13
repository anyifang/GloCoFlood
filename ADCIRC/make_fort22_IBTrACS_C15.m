%% ============================================================
%  auto_fort22_TCdataCv_C15_stepwrite.m
%
%  修正版：
%  1) 使用 TC_data_Cv.mat 中每个台风自己的 best_Cv
%  2) 预筛选采用：台风中心进入中心点 2° 圆后的自身最大 wind
%  3) fort.22 按 ADCIRC NWS=6 正确顺序写出：
%     k = 1..NWLAT (北->南), j = 1..NWLON (西->东), 每行 U V P
%  4) 风场计算在二维网格上完成，不再用一维向量直接写
%  5) 恢复远场平滑回背景场
%% ============================================================

clearvars -except auto_Name lon_center lat_center inner_radius
close all;
clc;

%% =========================
% 0) 只保留这三个输入
%% =========================
% auto_Name  = 'PRD';
% lon_center = 113.0;
% lat_center = 23.0;
% inner_radius = 2;
%% =========================
% 1) 参数区
%% =========================
% Winni 1407 Cv 2.32
% KATRINA 2003   Cv
%  

%auto_Name = 'YRD';lon_center = 121.5; lat_center = 31.2; year_tc = 1997;
 %auto_Name = 'PRD';lon_center = 113.4; lat_center = 22; year_tc = 2008;
%auto_Name = 'Wenzhou';lon_center = 120.9; lat_center = 27.9; alongshore_radius = 4.0;  offshore_deg = 4.0;   ellipse_ratio     =1 ;
%auto_Name = 'BoB';lon_center = 89.8; lat_center = 22.5; year_tc = 2007;
auto_Name = 'Misp';lon_center = -90; lat_center = 29.5; year_tc = 2005;
inner_radius   = 2; % 


tc_cv_mat_file = 'H:\Global_compoundflood\TC_data\wind_stationdata\global_TC_validation_C15\TC_data_Cv.mat';
load(tc_cv_mat_file)
tc_var_name    = 'TC_data';

TC_data_Cv = TC_data;
output_root  = ['cases_' auto_Name];

IMPACT_RADIUS_DEG     = inner_radius + 2;    % 候选/预筛选区域：中心点2°圆
DOMAIN_WIND_THRESHOLD = 25.0;   % 预筛选阈值：台风自身强度 wind (m/s)

PRE_DAYS  = 4;
POST_DAYS = 2;
DT_HOUR   = 1;                  % fort.22 时间步（小时）

MET_GRID_MARGIN_DEG = 0.20;
DLON = 0.10;
DLAT = 0.10;

Pc_default_hPa  = 950;
Pn_hPa          = 1013;
Re              = 6371000;

WRITE_DIAG_MAT = true;
PLOT_CHECK     = false;

if ~exist(output_root, 'dir')
    mkdir(output_root);
end

%% =========================
% 2) 读取 TC_data_Cv
%% =========================
S_tc = load(tc_cv_mat_file);
if ~isfield(S_tc, tc_var_name)
    error('在 %s 中没有找到变量 %s', tc_cv_mat_file, tc_var_name);
end
TC_data = S_tc.(tc_var_name);

fprintf('\n读取 TC_data_Cv 成功，总台风数：%d\n', numel(TC_data));

%% =========================
% 3) 预加载 WindC15Lib
%% =========================
Vmax_list = 15:1:120;
Rmax_list = 20:1:200; % km
WindC15Lib = containers.Map;

% 
% for i = 1:length(TC_data)
%     if isempty(TC_data(i).best_Cv)
%         TC_data(i).best_Cv = 1.3;
%     end
% end

for Vmax = Vmax_list
    for Rmaxkm = Rmax_list
        fname = fullfile( ...
            'H:\wind_turbine\wind_model\CLE15_windprofile_PUBLIC_2020-06-23\C15_predata\', ...
            ['Wind_C15_data_Vmax', num2str(Vmax), '_Rmax', num2str(Rmaxkm), '.mat']);
        if exist(fname, 'file')
            S = load(fname);
            key = sprintf('%d_%d', Vmax, Rmaxkm);
            WindC15Lib(key) = S.Wind_C15_data;
        end
    end
end

%% =========================
% 4) 读取 fort.14，确定 fort.22 覆盖范围
%% =========================
fort14_file = resolve_fort14_file(auto_Name);
fprintf('\n读取 fort.14: %s\n', fort14_file);

Mesh = read_fort14_nodes_only(fort14_file);

lon_min = min(Mesh.lon);
lon_max = max(Mesh.lon);
lat_min = min(Mesh.lat);
lat_max = max(Mesh.lat);

fprintf('ADCIRC 网格范围（来自 fort.14）:\n');
fprintf('  Lon: [%.4f, %.4f]\n', lon_min, lon_max);
fprintf('  Lat: [%.4f, %.4f]\n', lat_min, lat_max);

%% =========================
% 5) 构造 fort.22 网格
%% =========================
WLONMIN = lon_min - MET_GRID_MARGIN_DEG;
WLONMAX = lon_max + MET_GRID_MARGIN_DEG;
WLATMIN = lat_min - MET_GRID_MARGIN_DEG;
WLATMAX = lat_max + MET_GRID_MARGIN_DEG;

met_lon = WLONMIN : DLON : WLONMAX;   % 西 -> 东
met_lat = WLATMAX : -DLAT : WLATMIN;  % 北 -> 南

NWLON = numel(met_lon);
NWLAT = numel(met_lat);

[MET_LON, MET_LAT] = meshgrid(met_lon, met_lat);   % [NWLAT, NWLON]

fprintf('\n研究中心: %s (%.3f, %.3f)\n', auto_Name, lon_center, lat_center);
fprintf('候选/预筛选区域半径: %.2f deg\n', IMPACT_RADIUS_DEG);
fprintf('fort.22 网格覆盖 fort.14 整个区域\n');
fprintf('NWLON = %d, NWLAT = %d, WTIMINC = %d s\n', NWLON, NWLAT, DT_HOUR*3600);

%% =========================
% 6) 用中心点2°圆筛选候选台风，并要求有 best_Cv
%% =========================
ImpactList = struct([]);
nImpact = 0;

for i = 1:numel(TC_data)
    tc = TC_data(i);

    if ~isfield(tc, 'time') || isempty(tc.time) || ...
       ~isfield(tc, 'lat')  || isempty(tc.lat)  || ...
       ~isfield(tc, 'lon')  || isempty(tc.lon)  || ...
       ~isfield(tc, 'wind') || isempty(tc.wind)
        continue;
    end

    if ~isfield(tc, 'best_Cv') || isempty(tc.best_Cv) || ~isfinite(tc.best_Cv) || tc.year ~= year_tc
        continue;
    end

    dist_deg = angular_distance_deg(tc.lat, tc.lon, lat_center, lon_center);
    inbuf = dist_deg <= IMPACT_RADIUS_DEG;

    if ~any(inbuf)
        continue;
    end

    nImpact = nImpact + 1;
    ImpactList(nImpact).iTC        = i;
    ImpactList(nImpact).idx_first  = find(inbuf, 1, 'first');
    ImpactList(nImpact).idx_last   = find(inbuf, 1, 'last');
    ImpactList(nImpact).time_first = tc.time(ImpactList(nImpact).idx_first);
    ImpactList(nImpact).time_last  = tc.time(ImpactList(nImpact).idx_last);
    ImpactList(nImpact).best_Cv    = tc.best_Cv;
    ImpactList(nImpact).tc_id      = make_tc_id(tc, i);

    wind_inbuf = tc.wind(inbuf);
    wind_inbuf = wind_inbuf(isfinite(wind_inbuf));
    if isempty(wind_inbuf)
        ImpactList(nImpact).max_tc_intensity_in_region = NaN;
    else
        ImpactList(nImpact).max_tc_intensity_in_region = max(wind_inbuf);
    end

    if isfield(tc, 'rmse') && ~isempty(tc.rmse) && isfinite(tc.rmse)
        ImpactList(nImpact).rmse = tc.rmse;
    else
        ImpactList(nImpact).rmse = NaN;
    end
end

fprintf('进入中心点 %.2f°圆范围且具有 best_Cv 的候选台风数：%d\n', IMPACT_RADIUS_DEG, nImpact);

if nImpact == 0
    warning('没有找到影响该中心点且带有 best_Cv 的台风。');
    return;
end




%% =========================
% 7) 逐个台风处理
%% =========================
SummaryCell = {};

for ii = 1:nImpact
    tc = TC_data(ImpactList(ii).iTC);
    Cv_tc = ImpactList(ii).best_Cv;
    tc_id = ImpactList(ii).tc_id;

    fprintf('\n====================================================\n');
    fprintf('处理 TC %d / %d: %s (%s)\n', ii, nImpact, char(tc.name), tc_id);
    fprintf('使用 Cv = %.3f\n', Cv_tc);
    fprintf('区域内台风最大自身强度 = %.2f m/s\n', ImpactList(ii).max_tc_intensity_in_region);
    fprintf('====================================================\n');

    maxWindPre = ImpactList(ii).max_tc_intensity_in_region;
    maxWindPreTime = ImpactList(ii).time_first;

    fprintf('  预筛选强度: %.2f m/s\n', maxWindPre);

    if ~isfinite(maxWindPre) || maxWindPre < DOMAIN_WIND_THRESHOLD
        fprintf('  -> 未达到 %.2f m/s，直接跳过，不生成 fort.22\n', DOMAIN_WIND_THRESHOLD);

        SummaryCell(end+1,:) = {tc_id, char(tc.name), tc.year, Cv_tc, ...
            ImpactList(ii).time_first, ImpactList(ii).time_last, ...
            maxWindPre, maxWindPreTime, false, '', 'precheck_fail'}; %#ok<AGROW>
        continue;
    end

    impact_day = dateshift(ImpactList(ii).time_first, 'start', 'day');
    start_day  = impact_day - caldays(PRE_DAYS);
    end_day    = impact_day + caldays(POST_DAYS);
    time_data  = (start_day : hours(DT_HOUR) : end_day)';
    run_hour   = numel(time_data);

    fprintf('研究时段: %s --> %s\n', datestr(start_day), datestr(end_day));
    fprintf('时次数量: %d\n', run_hour);

    Iwin = tc.time >= start_day & tc.time <= end_day;
    tc_win = clip_tc_by_index(tc, Iwin);

    if isempty(tc_win.time)
        fprintf('  -> 该时段内没有台风记录，直接跳过\n');
        SummaryCell(end+1,:) = {tc_id, char(tc.name), tc.year, Cv_tc, start_day, end_day, ...
            maxWindPre, maxWindPreTime, false, '', 'no_track_in_window'}; %#ok<AGROW>
        continue;
    end

    TC_start_time = tc_win.time(1);
    TC_end_time   = tc_win.time(end);

    storm_tag = sprintf('%04d_%s_%s', tc.year, sanitize_name(char(tc.name)), sanitize_name(tc_id));
    case_dir = fullfile(output_root, storm_tag);
    if ~exist(case_dir, 'dir')
        mkdir(case_dir);
    end

    fort22_file = fullfile(case_dir, 'fort.22');
    meta_txt    = fullfile(case_dir, 'fort22_meta.txt');
    diag_mat    = fullfile(case_dir, 'storm_diagnostic.mat');

    fid = fopen(fort22_file, 'wt');
    if fid < 0
        error('无法写入 %s', fort22_file);
    end

    maxWindDomain = 0;
    maxWindTime   = NaT;
    maxWindSeries = nan(run_hour,1);
    RmaxSeries_km = nan(run_hour,1);

    for n = 1:run_hour
        t = time_data(n);

        [has_tc, state_now, state_next] = get_tc_state_for_step(tc_win, t, hours(DT_HOUR), Pc_default_hPa);

        if ~has_tc
            write_background_step_nws6(fid, NWLAT, NWLON);
            maxWindSeries(n) = 0;
            RmaxSeries_km(n) = NaN;
        else
            [U, V, P, Rmax_m] = calc_c15_uvp_field_grid( ...
                state_now.lat, state_now.lon, ...
                state_next.lat, state_next.lon, DT_HOUR, ...
                state_now.wind, state_now.p, ...
                MET_LAT, MET_LON, ...
                Cv_tc, Pn_hPa, Re, WindC15Lib);

            thisMax = max(hypot(U(:), V(:)));
            maxWindSeries(n) = thisMax;
            RmaxSeries_km(n) = Rmax_m / 1000;

            if thisMax > maxWindDomain
                maxWindDomain = thisMax;
                maxWindTime = t;
            end

            write_fort22_snapshot_nws6(fid, U, V, P);
        end
    end

    fclose(fid);

    write_fort22_meta(meta_txt, fort14_file, WLONMIN, WLATMAX, DLON, DLAT, ...
        NWLON, NWLAT, DT_HOUR*3600, start_day, end_day, TC_start_time, TC_end_time, Cv_tc, tc_id);

    if WRITE_DIAG_MAT
        Storm = struct();
        Storm.auto_Name      = auto_Name;
        Storm.lon_center     = lon_center;
        Storm.lat_center     = lat_center;
        Storm.impact_radius  = IMPACT_RADIUS_DEG;
        Storm.fort14_file    = fort14_file;
        Storm.fort14_bbox    = [lon_min lon_max lat_min lat_max];
        Storm.name           = char(tc.name);
        Storm.tc_id          = tc_id;
        Storm.year           = tc.year;
        Storm.Cv_used        = Cv_tc;
        Storm.max_tc_intensity_in_region = maxWindPre;
        if isfield(tc,'rmse'); Storm.rmse = tc.rmse; end
        Storm.start_day      = start_day;
        Storm.end_day        = end_day;
        Storm.TC_start_time  = TC_start_time;
        Storm.TC_end_time    = TC_end_time;
        Storm.time_data      = time_data;
        Storm.maxWindSeries  = maxWindSeries;
        Storm.maxWindDomain  = maxWindDomain;
        Storm.maxWindTime    = maxWindTime;
        Storm.RmaxSeries_km  = RmaxSeries_km;
        save(diag_mat, 'Storm', '-v7.3');
    end

    if PLOT_CHECK
        figure('Color','w');
        plot(time_data, maxWindSeries, 'k-', 'LineWidth', 1.5); hold on;
        grid on;
        title(sprintf('%s (%s), Cv=%.3f', char(tc.name), tc_id, Cv_tc), 'Interpreter', 'none');
        ylabel('Wind speed (m/s)');
        drawnow;
    end

    fprintf('  -> 已生成 fort.22: %s\n', fort22_file);
    fprintf('  -> 正式区域最大风速: %.2f m/s @ %s\n', maxWindDomain, fmt_time_safe(maxWindTime));

    SummaryCell(end+1,:) = {tc_id, char(tc.name), tc.year, Cv_tc, start_day, end_day, ...
        maxWindDomain, maxWindTime, true, fort22_file, 'kept'}; %#ok<AGROW>
end

%% =========================
% 8) 汇总输出
%% =========================
Summary = cell2table(SummaryCell, 'VariableNames', { ...
    'TCID','Name','Year','CvUsed','StartDay','EndDay', ...
    'DomainMaxWind','DomainMaxWindTime','Kept','Fort22Path','Status'});

summary_csv = fullfile(output_root, 'fort22_case_summary_withCv.csv');
writetable(Summary, summary_csv);

fprintf('\n====================================================\n');
fprintf('全部完成。\n');
fprintf('汇总表：%s\n', summary_csv);
fprintf('====================================================\n');



%% ============================================================
% 局部函数
%% ============================================================

function fort14_file = resolve_fort14_file(auto_Name)
cand1 = ['ADCIRC_Capsule_Mesh_' auto_Name '.14'];
cand2 = 'fort.14';

if exist(cand1, 'file')
    fort14_file = cand1; return;
end
if exist(cand2, 'file')
    fort14_file = cand2; return;
end

D = dir('*.14');
if numel(D) == 1
    fort14_file = D(1).name;
elseif numel(D) > 1
    names = {D.name};
    error('当前目录存在多个 .14 文件，请保留 fort.14 或 ADCIRC_Capsule_Mesh_%s.14。\n检测到：\n%s', ...
        auto_Name, strjoin(names, '\n'));
else
    error('当前目录没有找到 .14 文件。');
end
end

function Mesh = read_fort14_nodes_only(fort14_file)
fid = fopen(fort14_file, 'rt');
if fid < 0
    error('无法打开 fort.14: %s', fort14_file);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fgetl(fid);
tmp = sscanf(strtrim(fgetl(fid)), '%f');
NE = tmp(1);
NP = tmp(2);

node_id = zeros(NP,1);
lon = zeros(NP,1);
lat = zeros(NP,1);
dep = zeros(NP,1);

for i = 1:NP
    a = sscanf(fgetl(fid), '%f');
    node_id(i) = a(1);
    lon(i)     = a(2);
    lat(i)     = a(3);
    dep(i)     = a(4);
end

Mesh = struct('NE',NE,'NP',NP,'node_id',node_id,'lon',lon,'lat',lat,'dep',dep);
end

function tc_out = clip_tc_by_index(tc, I)
tc_out = tc;
tc_out.time = tc.time(I);
tc_out.lat  = tc.lat(I);
tc_out.lon  = tc.lon(I);
tc_out.wind = tc.wind(I);
tc_out.p    = tc.p(I);

fields_copy = {'name','year','best_Cv','rmse','basin','sid'};
for k = 1:numel(fields_copy)
    f = fields_copy{k};
    if isfield(tc, f)
        tc_out.(f) = tc.(f);
    end
end
end

function [has_tc, state_now, state_next] = get_tc_state_for_step(tc_win, t, dt_step, Pc_default_hPa)

if isempty(tc_win.time) || t < tc_win.time(1) || t > tc_win.time(end)
    has_tc = false;
    state_now  = empty_tc_state();
    state_next = empty_tc_state();
    return;
end

state_now = interp_tc_state(tc_win, t, Pc_default_hPa);

t2 = t + dt_step;
if t2 > tc_win.time(end)
    state_next = state_now;
else
    state_next = interp_tc_state(tc_win, t2, Pc_default_hPa);
    if ~isfinite(state_next.lat) || ~isfinite(state_next.lon) || ~isfinite(state_next.wind)
        state_next = state_now;
    end
end

if ~isfinite(state_now.lat) || ~isfinite(state_now.lon) || ~isfinite(state_now.wind)
    has_tc = false;
else
    has_tc = true;
end
end

function s = interp_tc_state(tc, t_query, Pc_default_hPa)
t_src = datenum(tc.time(:));
tq    = datenum(t_query);

[t_src_u, IA] = unique(t_src, 'stable');

lat_u  = tc.lat(IA);
lon_u  = tc.lon(IA);
wind_u = tc.wind(IA);
p_u    = tc.p(IA);

lat  = interp1(t_src_u, lat_u,  tq, 'linear', NaN);
lon  = interp1(t_src_u, lon_u,  tq, 'linear', NaN);
wind = interp1(t_src_u, wind_u, tq, 'linear', NaN);
p    = interp1(t_src_u, p_u,    tq, 'linear', NaN);

if ~isfinite(p) || p <= 0
    p = Pc_default_hPa;
end

s = struct('lat',lat,'lon',lon,'wind',wind,'p',p);
end

function s = empty_tc_state()
s = struct('lat',NaN,'lon',NaN,'wind',NaN,'p',NaN);
end

function [U, V, P, Rmax_m] = calc_c15_uvp_field_grid( ...
    lat_c, lon_c, lat_n, lon_n, dt_hr, v_c, p_c, ...
    LAT, LON, Cv, Pn_hPa, Re, WindC15Lib)

U = zeros(size(LAT));
V = zeros(size(LAT));
P = ones(size(LAT)) * 101300;
Rmax_m = NaN;

if ~isfinite(lat_c) || ~isfinite(lon_c) || ~isfinite(v_c)
    return;
end
if ~isfinite(p_c) || p_c <= 0
    p_c = 950;
end

if dt_hr <= 0 || ~isfinite(dt_hr)
    vmc = 0;
    fai = 0;
else
    [dis, alpha] = fast_track_step(lat_c, lon_c, lat_n, lon_n);
    vmc = dis / (dt_hr * 3600);
    fai = alpha;
end

vm = v_c - vmc;
if vm < 15
    return;
end
vm = max(vm, 15);
vm = min(vm, 120);

Rmax_m = Cv * 51.6 * exp(-0.0223 * vm + 0.0281 * lat_c) * 1000;
Rmax_m = max(Rmax_m, 30e3);

B = (vm^2) * 1.15 * exp(1) / max(Pn_hPa - p_c, 1) / 100;
B = max(1.0, min(B, 2.5));

dx = (LON - lon_c) .* cosd(0.5 * (LAT + lat_c)) * 111320;
dy = (LAT - lat_c) * 110540;
r  = hypot(dx, dy);
r  = max(r, 1.0);

cta = atan2d(dy, dx);
cta(cta < 0) = cta(cta < 0) + 360;

Pg_TC = (p_c + (Pn_hPa - p_c) .* exp(-(Rmax_m ./ r).^B)) * 100;

Vmax   = round(vm);
Rmaxkm = round(Rmax_m / 1000);
key    = sprintf('%d_%d', Vmax, Rmaxkm);

if isKey(WindC15Lib, key)
    C15 = WindC15Lib(key);
    rr = double(C15.rr(:));
    vg_tab = double(C15.vg(:));
    vg = interp1(rr, vg_tab, r, 'linear', 0);
else
    vg = zeros(size(r));
end

vg(~isfinite(vg) | vg < 0) = 0;

beta = zeros(size(r));
I1 = r < Rmax_m;
I2 = r >= Rmax_m & r < 1.2 * Rmax_m;
I3 = r >= 1.2 * Rmax_m;

beta(I1) = 10 .* (1 + r(I1) ./ Rmax_m);
beta(I2) = 20 + 25 .* (r(I2) ./ Rmax_m - 1);
beta(I3) = 25;

vmoc = vmc .* r .* Rmax_m ./ (r.^2 + Rmax_m^2);

Vx_TC = 0.85 .* vg .* cosd(cta + 90 + beta) + vmoc .* cosd(fai);
Vy_TC = 0.85 .* vg .* sind(cta + 90 + beta) + vmoc .* sind(fai);

Vx_TC = 0.893 .* Vx_TC;
Vy_TC = 0.893 .* Vy_TC;

% 远场平滑回背景场
R1 = 300e3;
R2 = 400e3;
lamda = zeros(size(r));
I_mid = r >= R1 & r <= R2;
I_far = r > R2;
lamda(I_mid) = (r(I_mid) - R1) ./ (R2 - R1);
lamda(I_far) = 1.0;

U = (1 - lamda) .* Vx_TC;
V = (1 - lamda) .* Vy_TC;
P = (1 - lamda) .* Pg_TC + lamda .* 101300;
end

function [dis_m, alpha_deg] = fast_track_step(lat1, lon1, lat2, lon2)
dx = (lon2 - lon1) * cosd(0.5*(lat1 + lat2)) * 111320;
dy = (lat2 - lat1) * 110540;
dis_m = hypot(dx, dy);
alpha_deg = atan2d(lat2 - lat1, lon2 - lon1);
if alpha_deg < 0
    alpha_deg = alpha_deg + 360;
end
end

function write_background_step_nws6(fid, NWLAT, NWLON)
for k = 1:NWLAT
    for j = 1:NWLON
        fprintf(fid, '%.0f %.0f %.0f\n', 0, 0, 101300);
    end
end
end

function write_fort22_snapshot_nws6(fid, U, V, P)
[nlat, nlon] = size(U);

tol_uv = 1e-8;
tol_p  = 1e-4;

for k = 1:nlat       % 北 -> 南
    for j = 1:nlon   % 西 -> 东
        u = U(k,j);
        v = V(k,j);
        p = P(k,j);

        % 背景常态值：不用小数，减小文件体积
        if abs(u) < tol_uv && abs(v) < tol_uv && abs(p - 101300) < tol_p
            fprintf(fid, '0 0 101300\n');
        else
            fprintf(fid, '%.2f %.2f %.2f\n', u, v, p);
        end
    end
end
end

function write_fort22_meta(meta_txt, fort14_file, WLONMIN, WLATMAX, DLON, DLAT, ...
    NWLON, NWLAT, WTIMINC, start_day, end_day, TC_start_time, TC_end_time, Cv_tc, tc_id)

fid = fopen(meta_txt, 'wt');
if fid < 0
    error('无法写入 %s', meta_txt);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'ADCIRC fort.22 (NWS=6) metadata\n');
fprintf(fid, '---------------------------------\n');
fprintf(fid, 'fort14       = %s\n', fort14_file);
fprintf(fid, 'tc_id        = %s\n', tc_id);
fprintf(fid, 'Cv_used      = %.2f\n', Cv_tc);
fprintf(fid, 'WindowStart  = %s\n', datestr(start_day, 31));
fprintf(fid, 'WindowEnd    = %s\n', datestr(end_day, 31));
fprintf(fid, 'TCStart      = %s\n', datestr(TC_start_time, 31));
fprintf(fid, 'TCEnd        = %s\n', datestr(TC_end_time, 31));
fprintf(fid, 'WTIMINC      = %d\n', WTIMINC);
fprintf(fid, 'NWLON        = %d\n', NWLON);
fprintf(fid, 'NWLAT        = %d\n', NWLAT);
fprintf(fid, 'WLONMIN      = %.2f\n', WLONMIN);
fprintf(fid, 'WLATMAX      = %.2f\n', WLATMAX);
fprintf(fid, 'WLONINC      = %.2f\n', DLON);
fprintf(fid, 'WLATINC      = %.2f\n', DLAT);
end

function dist_deg = angular_distance_deg(lat1, lon1, lat2, lon2)
lat1r = deg2rad(lat1);
lon1r = deg2rad(lon1);
lat2r = deg2rad(lat2);
lon2r = deg2rad(lon2);

dlat = lat1r - lat2r;
dlon = lon1r - lon2r;

a = sin(dlat/2).^2 + cos(lat1r).*cos(lat2r).*sin(dlon/2).^2;
c = 2 .* atan2(sqrt(a), sqrt(max(1-a,0)));
dist_deg = rad2deg(c);
end

function s = sanitize_name(s)
s = char(s);
s = strtrim(s);
if isempty(s), s = 'NONAME'; end
s = regexprep(s, '[^\w\-]+', '_');
end

function s = fmt_time_safe(t)
if isnat(t)
    s = 'NaT';
else
    s = datestr(t);
end
end

function tc_id = make_tc_id(tc, idx)
if isfield(tc, 'sid') && ~isempty(tc.sid)
    sid_str = char(tc.sid);
    sid_str = strtrim(sid_str);
    if ~isempty(sid_str)
        tc_id = sid_str;
        return;
    end
end

if isfield(tc, 'basin') && ~isempty(tc.basin)
    basin_str = char(tc.basin);
else
    basin_str = 'NA';
end

if isfield(tc, 'year') && ~isempty(tc.year) && isfinite(tc.year)
    year_num = tc.year;
else
    year_num = 0;
end

if isfield(tc, 'name') && ~isempty(tc.name)
    name_str = char(tc.name);
else
    name_str = 'NONAME';
end

tc_id = sprintf('%s_%04d_%s_%05d', ...
    sanitize_name(basin_str), year_num, sanitize_name(name_str), idx);
end



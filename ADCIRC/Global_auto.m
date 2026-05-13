% 1. 设置中心点 (上海为例，该算法同样完美适用于珠江口)
clear  
close all
auto_Name = 'YRD';lon_center = 121.5; lat_center = 31.2; alongshore_radius = 5.0; offshore_deg = 4.3;   ellipse_ratio     =1.2 ;
% auto_Name = 'PRD';lon_center = 113.4; lat_center = 22; alongshore_radius = 5.0;  offshore_deg = 5.0;   ellipse_ratio     =1 ;
%auto_Name = 'Wenzhou';lon_center = 120.9; lat_center = 27.9; alongshore_radius = 4.0;  offshore_deg = 4.0;   ellipse_ratio     =1 ;
%auto_Name = 'BoB';lon_center = 89.8; lat_center = 22.5; alongshore_radius = 5; offshore_deg = 5; ellipse_ratio=1 ;
%auto_Name = 'Houston';lon_center = -94.8; lat_center = 29.6; alongshore_radius = 6.0; offshore_deg = 7.0;     ellipse_ratio     =0.8 ;
%auto_Name = 'Washington';lon_center = -75.5; lat_center = 37.8; alongshore_radius = 5.0; offshore_deg = 6.0;   ellipse_ratio     =1 ;
%auto_Name = 'Misp';lon_center = -90; lat_center = 29.5; alongshore_radius = 5.0; offshore_deg = 5.2;   ellipse_ratio     =1.2 ;
%auto_Name = 'Ho Chi Minh';lon_center = 106.4; lat_center = 10.6; alongshore_radius = 5.0; offshore_deg = 6.0;   ellipse_ratio     =1 ;

% 2. 完美几何控制参数

inland_deg = 1.5;        % [吃进陆地] 
bay_close_deg = 3.0;     % [湾区识别]
 % ★ 新增：椭圆曲率(扁率)。1为正圆，0.5为压扁一半的椭圆
inner_radius   = 2; % 
% 3. 数据路径 
gshhs_shp = '..\datasets\GSHHS_shp\c\GSHHS_c_L1.shp'; 
out_dir = '..\datasets\Global_autodata';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
outer_shp = fullfile(out_dir, sprintf('%s_domain.shp', auto_Name));

% 4. 运行形态学平滑版本  auto_adcirc_domain_capsule  auto_adcirc_domain_morphological
auto_adcirc_domain_capsule(lon_center, lat_center, alongshore_radius, offshore_deg, inland_deg, ...
    bay_close_deg, ellipse_ratio ,inner_radius , gshhs_shp, outer_shp);


% =========================================================================
% 基于 auto_adcirc_domain_capsule 生成的边界，使用 OceanMesh2D 划分网格
% =========================================================================
% 添加你的 OceanMesh2D 路径
addpath(genpath('../utilities/'))
addpath(genpath('../datasets/'))
addpath(genpath('../m_map/'))

%% STEP 1: 读取我们刚刚生成的内外圈计算域边界
fprintf('1. 读取 Capsule 流体计算域边界...\n');
outer_shp_file = outer_shp;
[fpath, fname, fext] = fileparts(outer_shp);
out_inner_shp = fullfile(fpath, [fname, '_inner', fext]);
inner_shp_file = out_inner_shp;

% 提取外圈多边形坐标 (去除 NaN)
outer_shp = shaperead(outer_shp_file);
bbox_outer = [outer_shp.X', outer_shp.Y'];
bbox_outer(any(isnan(bbox_outer), 2), :) = []; 

% 提取内圈多边形坐标 (如果有局部加密)
has_inner = exist(inner_shp_file, 'file');
if has_inner
    inner_shp = shaperead(inner_shp_file);
    bbox_inner = [inner_shp.X', inner_shp.Y'];
    bbox_inner(any(isnan(bbox_inner), 2), :) = [];
end

%% STEP 2: 设置地形、海岸线与网格分辨率参数
fprintf('1.5 过滤海岸线中的微小岛屿...\n');


coastline_coarse = 'land_polygons';
coastline_fine   = 'land_polygons';


dem_coarse       = 'GEBCO_2025_sub_ice.nc'; % 大洋地形
fine_dem_dir = fullfile('..','datasets','fine_dem');

[dem_fine, dem_info] = auto_select_fine_dem_by_point( ...
    lon_center, lat_center, fine_dem_dir, dem_coarse);

fprintf('dem_coarse = %s\n', dem_coarse);
fprintf('dem_fine   = %s\n', dem_fine);
disp(dem_info);


% --- 外圈 (全局大洋区) 参数 ---
min_el_out    = 1e3;    % 外海最小分辨率 (m)
max_el_out    = 20e3;   % 外海最大分辨率 (m)
max_el_ns_out = 5e3;    % 外海近岸最大分辨率 (m)

% --- 内圈 (局部加密区) 参数 ---
min_el_in     = 200;    % 加密区最小分辨率 (m)
max_el_in     = 5e3;    % 加密区最大分辨率 (m)
max_el_ns_in  = 500;    % 加密区近岸最大分辨率 (m)
wl = 30;
g = 0.25;
%% STEP 3: 构建外海 (Outer Domain) 的 Geodata 与 Edgefx
fprintf('2. 构建外海大洋网格控制函数...\n');
gdat_out = geodata('shp', coastline_coarse, 'dem', dem_coarse, ...
                   'bbox', bbox_outer, 'h0', min_el_out, 'window', 10);

fh_out = edgefx('geodata', gdat_out, 'fs', 3, ...
                'max_el_ns', max_el_ns_out, 'max_el', max_el_out, ...
                'dt', 5, 'g', g,'wl',wl);


fprintf('3. 构建内部局部加密网格控制函数...\n');
gdat_in = geodata('shp', coastline_fine, 'dem', dem_fine, ...
    'bbox', bbox_inner, 'h0', min_el_in, 'window', 10);

fh_in = edgefx('geodata', gdat_in, 'fs', 6, ...
    'max_el_ns', max_el_ns_in, 'max_el', max_el_in, ...
    'dt',5,'g', g);

% 组合边界和边长控制函数
ef_list = {fh_out, fh_in};
bou_list = {gdat_out, gdat_in};


%% STEP 5: 融合并生成非结构网格 (Mesh Generation)
fprintf('4. 开始生成非结构三角网格 (MeshGen)...\n');
rng(1.23456789); % 固定随机种子以保证网格生成的可重复性
mshopts = meshgen('ef', ef_list, 'bou', {gdat_out gdat_in}, 'plot_on', 1, 'proj', 'lam');
mshopts = mshopts.build; 

%% STEP 6: 地形插值与水动力学限制打磨 (Post-Processing)
fprintf('5. 执行地形插值与动力学限制平滑...\n');
m = mshopts.grd;

% 1. 组合地形插值 (最小水深限制为 1m，防止干滩导致模型崩溃)
m = interp(m, bou_list, 'mindepth', 1); 

% 2. 限制深海地形坡度，防止内部斜压梯度误差 (HPG)
m = lim_bathy_slope(m, 0.1, 0);

% 3. CFL 数自动优化：通过局部调整水深/网格确保 5s 步长能稳定运行
m = bound_courant_number(m, 1.2, 0.5, 0, 10); 

% 评估建议步长
CFL = CalcCFL(m, 1);
fprintf('   -> Max CFL (1s): %.4f\n', max(CFL));
fprintf('   -> Min CFL (1s): %.4f\n', min(CFL));
fprintf('   -> 建议安全积分步长: %.2f s\n', min(CalcCFL(m)));

%% STEP 7: 自动生成开边界条件 (Open Boundaries)
fprintf('6. 搜寻并自动建立 ADCIRC 开边界条件...\n');
% 利用外部大洋边界数据 (gdat_out) 自动建立水位潮汐强迫边界
m = make_bc_from_capsule_boundary(m, gdat_out, ...
    'shore_tol', 0.1, ...
    'min_open_edges', 10, ...
    'bridge_gap_edges', 4, ...
    'add_inner_islands', true, ...
    'inner_ibtype', 21, ...
    'plot_check', true);
% m = make_bc(m,'auto',gdat_out ,'both',0.2,50);
plot(m, 'type', 'bd'); hold on; % 绘制带有水深(log scale)的网格  help msh/plot
%% STEP 8: 导出 fort.14 文件
fprintf('7. 导出 fort.14 格式网格文件...\n');
out_mesh_name = ['ADCIRC_Capsule_Mesh_',auto_Name];
write(m, out_mesh_name, '14');

% --- 绘制最终网格验证图 ---
figure('Name', 'Final ADCIRC Mesh', 'Color', 'w');
plot(m, 'type', 'resologmesh'); hold on; % 绘制带有水深(log scale)的网格  help msh/plot
save([auto_Name,'_ADCIRC_mesh'], 'm')
title('Final ADCIRC Mesh with Bathymetry');
drawnow;

fprintf('\n(+) 殿堂级计算网格生成完毕！请在当前目录下检查 %s.14 文件。\n', out_mesh_name);

make_fort22_IBTrACS_C15
make_fort15_from_template(auto_Name)
make_fort13_open_boundary_boost(auto_Name)
make_fort19_from_tmd(auto_Name)
make_fort14_to_cases(auto_Name)
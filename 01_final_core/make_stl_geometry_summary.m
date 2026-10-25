clear; clc; close all;

core_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(core_dir);

csv_dir = fullfile(project_root,'05_outputs','csv');
fig_dir = fullfile(project_root,'05_outputs','figures');

if ~exist(csv_dir,'dir'), mkdir(csv_dir); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

%% Fixed production STL path
stl_file = fullfile(project_root,'04_data','geometry','Meteor.stl');

if ~isfile(stl_file)
    error('Production STL not found: %s', stl_file);
end

fprintf('Using production STL file:\n%s\n', stl_file);

%% Read STL
[F_raw,V_raw] = read_stl_local(stl_file);

fprintf('Raw STL: vertices = %d, faces = %d\n', size(V_raw,1), size(F_raw,1));

%% Clean duplicate vertices
[F_raw,V_raw] = clean_mesh_local(F_raw,V_raw);

n_vertices = size(V_raw,1);
n_faces = size(F_raw,1);

fprintf('Clean STL: vertices = %d, faces = %d\n', n_vertices, n_faces);

%% Raw geometry
[A_raw,N_raw,C_raw] = face_geometry_local(V_raw,F_raw);
raw_surface_area = sum(A_raw);
raw_volume = abs(mesh_volume_local(V_raw,F_raw));
raw_bbox = max(V_raw,[],1) - min(V_raw,[],1);

if raw_volume > 0
    raw_equiv_diameter = (6*raw_volume/pi)^(1/3);
else
    raw_equiv_diameter = NaN;
end

%% Scale to nominal Chelyabinsk equivalent diameter
target_equiv_diameter_m = 19.8;

if isfinite(raw_equiv_diameter) && raw_equiv_diameter > 0
    scale_factor_to_m = target_equiv_diameter_m / raw_equiv_diameter;
    scale_mode = "volume_equivalent_diameter";
else
    scale_factor_to_m = target_equiv_diameter_m / max(raw_bbox);
    scale_mode = "maximum_extent_fallback";
end

V0_center = mean(V_raw,1);
V = (V_raw - V0_center) * scale_factor_to_m;
F = F_raw;

[A_face,N_face,C_face] = face_geometry_local(V,F);

surface_area_m2 = sum(A_face);
volume_m3 = abs(mesh_volume_local(V,F));
equiv_diameter_m = (6*volume_m3/pi)^(1/3);
bbox_m = max(V,[],1) - min(V,[],1);

sphere_ref_area_m2 = pi*(equiv_diameter_m^2)/4;
sphere_surface_area_m2 = pi*equiv_diameter_m^2;
sphere_volume_m3 = pi*equiv_diameter_m^3/6;

surface_area_ratio_to_sphere = surface_area_m2 / sphere_surface_area_m2;
volume_ratio_to_sphere = volume_m3 / sphere_volume_m3;

%% Projected area proxy
dirs_cardinal = [
     1  0  0
    -1  0  0
     0  1  0
     0 -1  0
     0  0  1
     0  0 -1];

Aproj_cardinal = zeros(6,1);

for i = 1:6
    Aproj_cardinal(i) = projected_area_local(A_face,N_face,dirs_cardinal(i,:));
end

Aproj_x_mean = mean(Aproj_cardinal(1:2));
Aproj_y_mean = mean(Aproj_cardinal(3:4));
Aproj_z_mean = mean(Aproj_cardinal(5:6));

rng(42);
n_orientation_samples = 500;
dirs = randn(n_orientation_samples,3);
dirs = dirs ./ sqrt(sum(dirs.^2,2));

Aproj_random = zeros(n_orientation_samples,1);

for i = 1:n_orientation_samples
    Aproj_random(i) = projected_area_local(A_face,N_face,dirs(i,:));
end

Aproj_p10 = pct_local(Aproj_random,10);
Aproj_p50 = pct_local(Aproj_random,50);
Aproj_p90 = pct_local(Aproj_random,90);
Aproj_mean = mean(Aproj_random);
Aproj_std = std(Aproj_random);

Aproj_ratio_p10 = Aproj_p10 / sphere_ref_area_m2;
Aproj_ratio_p50 = Aproj_p50 / sphere_ref_area_m2;
Aproj_ratio_p90 = Aproj_p90 / sphere_ref_area_m2;
Aproj_ratio_mean = Aproj_mean / sphere_ref_area_m2;

%% Save geometry summary
GeometrySummary = table( ...
    string(stl_file), ...
    string(scale_mode), ...
    n_vertices, ...
    n_faces, ...
    scale_factor_to_m, ...
    target_equiv_diameter_m, ...
    raw_surface_area, ...
    raw_volume, ...
    raw_equiv_diameter, ...
    bbox_m(1), bbox_m(2), bbox_m(3), ...
    surface_area_m2, ...
    volume_m3, ...
    equiv_diameter_m, ...
    sphere_ref_area_m2, ...
    surface_area_ratio_to_sphere, ...
    volume_ratio_to_sphere, ...
    Aproj_x_mean, ...
    Aproj_y_mean, ...
    Aproj_z_mean, ...
    Aproj_p10, ...
    Aproj_p50, ...
    Aproj_p90, ...
    Aproj_mean, ...
    Aproj_std, ...
    Aproj_ratio_p10, ...
    Aproj_ratio_p50, ...
    Aproj_ratio_p90, ...
    Aproj_ratio_mean, ...
    'VariableNames', { ...
    'stl_file', ...
    'scale_mode', ...
    'n_vertices', ...
    'n_faces', ...
    'scale_factor_to_m', ...
    'target_equiv_diameter_m', ...
    'raw_surface_area_units2', ...
    'raw_volume_units3', ...
    'raw_equiv_diameter_units', ...
    'bbox_x_m', ...
    'bbox_y_m', ...
    'bbox_z_m', ...
    'surface_area_m2', ...
    'volume_m3', ...
    'equiv_diameter_m', ...
    'sphere_ref_area_m2', ...
    'surface_area_ratio_to_sphere', ...
    'volume_ratio_to_sphere', ...
    'projected_area_x_mean_m2', ...
    'projected_area_y_mean_m2', ...
    'projected_area_z_mean_m2', ...
    'projected_area_random_p10_m2', ...
    'projected_area_random_p50_m2', ...
    'projected_area_random_p90_m2', ...
    'projected_area_random_mean_m2', ...
    'projected_area_random_std_m2', ...
    'projected_area_ratio_p10_to_sphere', ...
    'projected_area_ratio_p50_to_sphere', ...
    'projected_area_ratio_p90_to_sphere', ...
    'projected_area_ratio_mean_to_sphere'});

geometry_file = fullfile(csv_dir,'geometry_summary.csv');
writetable(GeometrySummary, geometry_file);

OrientationSamples = table( ...
    (1:n_orientation_samples)', ...
    dirs(:,1), dirs(:,2), dirs(:,3), ...
    Aproj_random, ...
    Aproj_random/sphere_ref_area_m2, ...
    'VariableNames', { ...
    'sample_id', ...
    'dir_x', ...
    'dir_y', ...
    'dir_z', ...
    'projected_area_m2', ...
    'projected_area_ratio_to_sphere'});

orientation_file = fullfile(csv_dir,'projected_area_orientation_samples.csv');
writetable(OrientationSamples, orientation_file);

%% Figure 11: STL geometry views
figure('Color','w','Position',[100 100 1200 900]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
trisurf(F,V(:,1),V(:,2),V(:,3),'EdgeColor','none');
axis equal; grid on;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('3D STL mesh');
view(3); camlight; lighting gouraud;

nexttile;
trisurf(F,V(:,1),V(:,2),V(:,3),'EdgeColor','none');
axis equal; grid on;
xlabel('x [m]'); ylabel('y [m]');
title('Top view');
view(2);

nexttile;
trisurf(F,V(:,1),V(:,2),V(:,3),'EdgeColor','none');
axis equal; grid on;
xlabel('x [m]'); zlabel('z [m]');
title('Side view');
view([0 0]);

nexttile;
trisurf(F,V(:,1),V(:,2),V(:,3),'EdgeColor','none');
axis equal; grid on;
ylabel('y [m]'); zlabel('z [m]');
title('Front view');
view([90 0]);

sgtitle('STL-resolved meteor geometry scaled to nominal equivalent diameter');

prepare_export_local(gcf);
exportgraphics(gcf, fullfile(fig_dir,'fig11_meteor_stl_geometry_views.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(fig_dir,'fig11_meteor_stl_geometry_views.pdf'), 'ContentType', 'vector');

%% Figure 12: projected area ratio distribution
figure('Color','w','Position',[100 100 900 650]);

histogram(Aproj_random/sphere_ref_area_m2,30);
hold on;

xline(Aproj_ratio_p10,'--','P10','LineWidth',1.2);
xline(Aproj_ratio_p50,'-','P50','LineWidth',1.4);
xline(Aproj_ratio_p90,'--','P90','LineWidth',1.2);
xline(1.0,':','Equivalent sphere','LineWidth',1.4);

xlabel('Projected area ratio to equivalent sphere [-]');
ylabel('Orientation count [-]');
title('Orientation-dependent projected-area proxy from STL mesh');
grid on;

prepare_export_local(gcf);
exportgraphics(gcf, fullfile(fig_dir,'fig12_projected_area_orientation_distribution.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(fig_dir,'fig12_projected_area_orientation_distribution.pdf'), 'ContentType', 'vector');

fprintf('\nGEOMETRY SUMMARY WRITTEN:\n');
fprintf('%s\n', geometry_file);
fprintf('%s\n', orientation_file);
fprintf('%s\n', fig_dir);

disp(GeometrySummary);

%% Local functions
function [F,V] = read_stl_local(filename)
    info = dir(filename);
    file_size = info.bytes;

    fid = fopen(filename,'rb');
    if fid < 0
        error('Could not open STL file: %s', filename);
    end

    fread(fid,80,'uint8');
    ntri = fread(fid,1,'uint32');
    fclose(fid);

    expected_binary_size = 84 + 50*double(ntri);

    if double(ntri) > 0 && expected_binary_size == file_size
        [F,V] = read_binary_stl_local(filename,double(ntri));
    else
        [F,V] = read_ascii_stl_local(filename);
    end
end

function [F,V] = read_binary_stl_local(filename,ntri)
    fid = fopen(filename,'rb');

    if fid < 0
        error('Could not open binary STL file: %s', filename);
    end

    fread(fid,80,'uint8');
    fread(fid,1,'uint32');

    V = zeros(3*ntri,3);
    F = zeros(ntri,3);

    for i = 1:ntri
        fread(fid,3,'float32');
        v1 = fread(fid,3,'float32')';
        v2 = fread(fid,3,'float32')';
        v3 = fread(fid,3,'float32')';
        fread(fid,1,'uint16');

        ii = 3*(i-1) + 1;
        V(ii,:) = v1;
        V(ii+1,:) = v2;
        V(ii+2,:) = v3;
        F(i,:) = [ii ii+1 ii+2];
    end

    fclose(fid);
end

function [F,V] = read_ascii_stl_local(filename)
    txt = fileread(filename);
    token_expr = 'vertex\s+([-+0-9EeDd\.]+)\s+([-+0-9EeDd\.]+)\s+([-+0-9EeDd\.]+)';
    tokens = regexp(txt,token_expr,'tokens');

    if isempty(tokens)
        error('ASCII STL parse failed: no vertices found.');
    end

    nv = numel(tokens);
    ntri = floor(nv/3);

    V = zeros(3*ntri,3);

    for i = 1:3*ntri
        V(i,1) = str2double(strrep(tokens{i}{1},'D','E'));
        V(i,2) = str2double(strrep(tokens{i}{2},'D','E'));
        V(i,3) = str2double(strrep(tokens{i}{3},'D','E'));
    end

    F = reshape(1:3*ntri,3,[])';
end

function [F2,V2] = clean_mesh_local(F,V)
    [V2,~,ic] = unique(V,'rows','stable');
    F2 = reshape(ic(F(:)),size(F));

    bad = F2(:,1)==F2(:,2) | F2(:,1)==F2(:,3) | F2(:,2)==F2(:,3);
    F2 = F2(~bad,:);
end

function [A,N,C] = face_geometry_local(V,F)
    v1 = V(F(:,1),:);
    v2 = V(F(:,2),:);
    v3 = V(F(:,3),:);

    e1 = v2 - v1;
    e2 = v3 - v1;

    cr = cross(e1,e2,2);
    cr_norm = sqrt(sum(cr.^2,2));

    A = 0.5*cr_norm;
    N = cr ./ max(cr_norm,eps);
    C = (v1 + v2 + v3)/3;

    mesh_center = mean(V,1);
    orient = sum(N.*(C - mesh_center),2);
    flip_idx = orient < 0;
    N(flip_idx,:) = -N(flip_idx,:);
end

function vol = mesh_volume_local(V,F)
    v1 = V(F(:,1),:);
    v2 = V(F(:,2),:);
    v3 = V(F(:,3),:);

    vol = sum(dot(v1,cross(v2,v3,2),2))/6;
end

function Aproj = projected_area_local(A_face,N_face,dir_vec)
    d = dir_vec(:)';
    d = d / norm(d);

    cosang = N_face*d';
    Aproj = sum(A_face .* max(0,cosang));
end

function y = pct_local(x,p)
    x = x(isfinite(x));
    x = sort(x(:));

    if isempty(x)
        y = NaN;
        return;
    end

    n = numel(x);

    if n == 1
        y = x(1);
        return;
    end

    pos = 1 + (n - 1)*p/100;
    lo = floor(pos);
    hi = ceil(pos);

    if lo == hi
        y = x(lo);
    else
        y = x(lo) + (pos - lo)*(x(hi) - x(lo));
    end
end

function prepare_export_local(fig)
    axs = findall(fig,'Type','axes');

    for ia = 1:numel(axs)
        try
            axtoolbar(axs(ia),'none');
        catch
        end

        try
            disableDefaultInteractivity(axs(ia));
        catch
        end

        try
            axs(ia).Toolbar.Visible = 'off';
        catch
        end
    end
end

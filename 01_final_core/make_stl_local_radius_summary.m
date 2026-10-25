clear; clc; close all;

core_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(core_dir);

csv_dir = fullfile(project_root,'05_outputs','csv');
summary_dir = fullfile(project_root,'05_outputs','summary');
fig_dir = fullfile(project_root,'05_outputs','figures');

if ~exist(csv_dir,'dir'), mkdir(csv_dir); end
if ~exist(summary_dir,'dir'), mkdir(summary_dir); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

%% Fixed production STL path
stl_file = fullfile(project_root,'04_data','geometry','Meteor.stl');

if ~isfile(stl_file)
    error('Production STL not found: %s', stl_file);
end

fprintf('Using production STL file:\n%s\n', stl_file);

%% Read and scale STL
[F_raw,V_raw] = read_stl_local(stl_file);
[F_raw,V_raw] = clean_mesh_local(F_raw,V_raw);

[A_raw,N_raw,C_raw] = face_geometry_local(V_raw,F_raw);
raw_volume = abs(mesh_volume_local(V_raw,F_raw));
raw_bbox = max(V_raw,[],1) - min(V_raw,[],1);

target_equiv_diameter_m = 19.8;

if raw_volume > 0
    raw_equiv_diameter = (6*raw_volume/pi)^(1/3);
    scale_factor_to_m = target_equiv_diameter_m / raw_equiv_diameter;
    scale_mode = "volume_equivalent_diameter";
else
    raw_equiv_diameter = NaN;
    scale_factor_to_m = target_equiv_diameter_m / max(raw_bbox);
    scale_mode = "maximum_extent_fallback";
end

V = (V_raw - mean(V_raw,1)) * scale_factor_to_m;
F = F_raw;

n_vertices = size(V,1);
n_faces = size(F,1);

[A_face,N_face,C_face] = face_geometry_local(V,F);

surface_area_m2 = sum(A_face);
volume_m3 = abs(mesh_volume_local(V,F));
equiv_diameter_m = (6*volume_m3/pi)^(1/3);

fprintf('Scaled mesh: vertices=%d, faces=%d, D_eq=%.4f m\n', n_vertices, n_faces, equiv_diameter_m);

%% Vertex dual area and area-weighted vertex normal
vertex_area = zeros(n_vertices,1);
vertex_normal_sum = zeros(n_vertices,3);

for i = 1:n_faces
    fi = F(i,:);
    Ai = A_face(i);
    ni = N_face(i,:);

    for j = 1:3
        vj = fi(j);
        vertex_area(vj) = vertex_area(vj) + Ai/3;
        vertex_normal_sum(vj,:) = vertex_normal_sum(vj,:) + Ai*ni/3;
    end
end

vertex_normal_norm = sqrt(sum(vertex_normal_sum.^2,2));
vertex_normal = vertex_normal_sum ./ max(vertex_normal_norm,eps);

%% Local edge-length scale
E = [
    F(:,[1 2])
    F(:,[2 3])
    F(:,[3 1])];

E = sort(E,2);
E = unique(E,'rows');

edge_vec = V(E(:,1),:) - V(E(:,2),:);
edge_len = sqrt(sum(edge_vec.^2,2));

edge_len_sum = zeros(n_vertices,1);
edge_len_count = zeros(n_vertices,1);

for i = 1:size(E,1)
    v1 = E(i,1);
    v2 = E(i,2);
    le = edge_len(i);

    edge_len_sum(v1) = edge_len_sum(v1) + le;
    edge_len_sum(v2) = edge_len_sum(v2) + le;

    edge_len_count(v1) = edge_len_count(v1) + 1;
    edge_len_count(v2) = edge_len_count(v2) + 1;
end

local_edge_length_m = edge_len_sum ./ max(edge_len_count,1);

%% Normal-dispersion curvature proxy
angle2_weighted_sum = zeros(n_vertices,1);
angle_weight_sum = zeros(n_vertices,1);

for i = 1:n_faces
    fi = F(i,:);
    ni = N_face(i,:);
    Ai = A_face(i);

    for j = 1:3
        vj = fi(j);
        nv = vertex_normal(vj,:);

        c = dot(ni,nv);
        c = max(-1,min(1,c));
        ang = acos(c);

        w = Ai/3;
        angle2_weighted_sum(vj) = angle2_weighted_sum(vj) + w*ang^2;
        angle_weight_sum(vj) = angle_weight_sum(vj) + w;
    end
end

normal_angle_rms_rad = sqrt(angle2_weighted_sum ./ max(angle_weight_sum,eps));

curvature_proxy_1_per_m = normal_angle_rms_rad ./ max(local_edge_length_m,eps);
local_radius_proxy_m = 1 ./ max(curvature_proxy_1_per_m,1e-9);

valid = isfinite(local_radius_proxy_m) & isfinite(curvature_proxy_1_per_m) & ...
        vertex_area > 0 & local_edge_length_m > 0;

R = local_radius_proxy_m(valid);
K = curvature_proxy_1_per_m(valid);
ang = normal_angle_rms_rad(valid);
lel = local_edge_length_m(valid);

%% Robust clipping for visualization only
R_p01 = pct_local(R,1);
R_p05 = pct_local(R,5);
R_p10 = pct_local(R,10);
R_p50 = pct_local(R,50);
R_p90 = pct_local(R,90);
R_p95 = pct_local(R,95);
R_p99 = pct_local(R,99);

K_p01 = pct_local(K,1);
K_p05 = pct_local(K,5);
K_p10 = pct_local(K,10);
K_p50 = pct_local(K,50);
K_p90 = pct_local(K,90);
K_p95 = pct_local(K,95);
K_p99 = pct_local(K,99);

angle_p50 = pct_local(ang,50);
edge_p50 = pct_local(lel,50);

R_vis = local_radius_proxy_m;
R_vis(~isfinite(R_vis)) = R_p99;
R_vis = min(max(R_vis,R_p01),R_p99);

logR_vis = log10(max(R_vis,1e-6));

%% Save vertex-level output
VertexRadius = table( ...
    (1:n_vertices)', ...
    V(:,1), V(:,2), V(:,3), ...
    vertex_area, ...
    local_edge_length_m, ...
    normal_angle_rms_rad, ...
    curvature_proxy_1_per_m, ...
    local_radius_proxy_m, ...
    log10(max(local_radius_proxy_m,1e-6)), ...
    valid, ...
    'VariableNames', { ...
    'vertex_id', ...
    'x_m', 'y_m', 'z_m', ...
    'vertex_dual_area_m2', ...
    'local_edge_length_m', ...
    'normal_angle_rms_rad', ...
    'curvature_proxy_1_per_m', ...
    'local_radius_proxy_m', ...
    'log10_local_radius_proxy_m', ...
    'valid_flag'});

vertex_file = fullfile(csv_dir,'local_radius_vertices.csv');
writetable(VertexRadius,vertex_file);

%% Save summary output
LocalRadiusSummary = table( ...
    string(stl_file), ...
    string(scale_mode), ...
    n_vertices, ...
    n_faces, ...
    scale_factor_to_m, ...
    target_equiv_diameter_m, ...
    surface_area_m2, ...
    volume_m3, ...
    equiv_diameter_m, ...
    sum(valid), ...
    mean(R,'omitnan'), ...
    std(R,'omitnan'), ...
    min(R), ...
    R_p01, R_p05, R_p10, R_p50, R_p90, R_p95, R_p99, ...
    max(R), ...
    mean(K,'omitnan'), ...
    std(K,'omitnan'), ...
    min(K), ...
    K_p01, K_p05, K_p10, K_p50, K_p90, K_p95, K_p99, ...
    max(K), ...
    angle_p50, ...
    edge_p50, ...
    'VariableNames', { ...
    'stl_file', ...
    'scale_mode', ...
    'n_vertices', ...
    'n_faces', ...
    'scale_factor_to_m', ...
    'target_equiv_diameter_m', ...
    'surface_area_m2', ...
    'volume_m3', ...
    'equiv_diameter_m', ...
    'n_valid_vertices', ...
    'local_radius_mean_m', ...
    'local_radius_std_m', ...
    'local_radius_min_m', ...
    'local_radius_p01_m', ...
    'local_radius_p05_m', ...
    'local_radius_p10_m', ...
    'local_radius_p50_m', ...
    'local_radius_p90_m', ...
    'local_radius_p95_m', ...
    'local_radius_p99_m', ...
    'local_radius_max_m', ...
    'curvature_proxy_mean_1_per_m', ...
    'curvature_proxy_std_1_per_m', ...
    'curvature_proxy_min_1_per_m', ...
    'curvature_proxy_p01_1_per_m', ...
    'curvature_proxy_p05_1_per_m', ...
    'curvature_proxy_p10_1_per_m', ...
    'curvature_proxy_p50_1_per_m', ...
    'curvature_proxy_p90_1_per_m', ...
    'curvature_proxy_p95_1_per_m', ...
    'curvature_proxy_p99_1_per_m', ...
    'curvature_proxy_max_1_per_m', ...
    'normal_angle_rms_p50_rad', ...
    'local_edge_length_p50_m'});

summary_file = fullfile(summary_dir,'local_radius_summary.csv');
writetable(LocalRadiusSummary,summary_file);

%% Figure 14: local radius map and distribution
figure('Color','w','Position',[100 100 1250 900]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
trisurf(F,V(:,1),V(:,2),V(:,3),logR_vis,'EdgeColor','none');
axis equal; grid on;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Local radius proxy map, 3D view');
view(3); camlight; lighting gouraud;
cb = colorbar;
ylabel(cb,'log_{10}(R_{local} [m])');

nexttile;
trisurf(F,V(:,1),V(:,2),V(:,3),logR_vis,'EdgeColor','none');
axis equal; grid on;
xlabel('x [m]'); ylabel('y [m]');
title('Top view');
view(2);
cb = colorbar;
ylabel(cb,'log_{10}(R_{local} [m])');

nexttile;
trisurf(F,V(:,1),V(:,2),V(:,3),logR_vis,'EdgeColor','none');
axis equal; grid on;
xlabel('x [m]'); zlabel('z [m]');
title('Side view');
view([0 0]);
cb = colorbar;
ylabel(cb,'log_{10}(R_{local} [m])');

nexttile;
histogram(log10(max(R,1e-6)),40);
grid on;
xlabel('log_{10}(local radius proxy [m])');
ylabel('Vertex count [-]');
title(sprintf('Radius distribution: P10/P50/P90 = %.2f / %.2f / %.2f m', R_p10, R_p50, R_p90));

sgtitle('STL local-radius proxy from mesh-normal variation');

prepare_export_local(gcf);
exportgraphics(gcf, fullfile(fig_dir,'fig14_local_radius_curvature_map.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(fig_dir,'fig14_local_radius_curvature_map.pdf'), 'ContentType', 'vector');

fprintf('\nLOCAL RADIUS SUMMARY WRITTEN:\n');
fprintf('%s\n', summary_file);
fprintf('%s\n', vertex_file);
fprintf('%s\n', fullfile(fig_dir,'fig14_local_radius_curvature_map.pdf'));

disp(LocalRadiusSummary);

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

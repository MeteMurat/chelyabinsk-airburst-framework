function run_r3c3b_corrected_full1000_overpressure( ...
    artifact_root,source_file,nominal_audit_file,contract_file, ...
    gate_file,status_file)
% CHELYABINSK R3C3B CORRECTED FULL-1000 OVERPRESSURE PROPAGATION
%
% Propagates the frozen corrected 1000-member airburst source matrix through
% the fixed-coefficient inverse-cube screening surrogate
%
%   DeltaP_kPa = C_overpressure * E_c_kton / s_km^3
%   s_km       = sqrt(r_ground_km^2 + h_burst_km^2)
%
% The calibration coefficient is frozen globally from the corrected nominal
% G4 anchor. It is NOT recalibrated independently for each realization.
%
% No entry repropagation, population exposure, manuscript modification, or
% canonical modification is performed.

files = string({source_file,nominal_audit_file,contract_file,gate_file,status_file});
if nargin < 6 || strlength(string(artifact_root))==0 || any(~isfile(files))
    error('R3C3B requires artifact_root and five valid frozen inputs.');
end

artifact_root = char(string(artifact_root));
output_dir = fullfile(artifact_root,'02_overpressure_outputs');
figure_dir = fullfile(output_dir,'diagnostic_figures');
report_dir = fullfile(artifact_root,'03_reports');
for p = {output_dir,figure_dir,report_dir}
    if ~exist(p{1},'dir'), mkdir(p{1}); end
end

fprintf('=== CHELYABINSK R3C3B CORRECTED FULL-1000 OVERPRESSURE PROPAGATION ===\n');
fprintf('scientific_entry_propagation_requested=False\n');
fprintf('overpressure_ensemble_propagation_requested=True\n');
fprintf('population_exposure_requested=False\n');
fprintf('candidate_manuscript_modification_requested=False\n');
fprintf('canonical_modification_requested=False\n');
fprintf('calibration_mode=GLOBAL_FIXED_NOMINAL_ANCHOR\n');
fprintf('per_realization_recalibration=False\n');
fprintf('terminal_mass_metric_publishable=False\n');
fprintf('artifact_root=%s\n',artifact_root);

SRC = read_csv_explicit(source_file);
NOM = read_csv_explicit(nominal_audit_file);
CONTRACT = read_csv_explicit(contract_file);
GATE = read_csv_explicit(gate_file);
status_text = fileread(status_file);

assert_columns(SRC,["run_id","fragmentation_seed","source_horizon_s", ...
    "outcome","h_burst_km","E_c_kton"],"full1000 pressure-source matrix");
assert_columns(NOM,["threshold_kPa","corrected_nominal_h_burst_km", ...
    "corrected_nominal_E_c_kton","C_overpressure_kPa_km3_per_kton", ...
    "R_ground_km"],"nominal threshold-radius audit");
assert_columns(CONTRACT,["contract_item","contract_symbol", ...
    "contract_value","contract_role"],"surrogate contract");
assert_columns(GATE,["gate_name","gate_pass","gate_detail"], ...
    "R3C3A gate adjudication");

upstream_status_pass = contains(status_text, ...
    'status=PASS_R3C3A_R1_THRESHOLD_RADIUS_INVERSE_ORDER_GATE_REPAIR_AND_CONTRACT_AUDIT');
upstream_all_gates_pass = all(truth_vector(GATE.gate_pass));

source_row_count_pass = height(SRC)==1000;
run_id_contract_pass = isequal(SRC.run_id,(1:1000)') && ...
    numel(unique(SRC.run_id))==1000;
seed_contract_pass = numel(unique(SRC.fragmentation_seed))==1000;
all_airburst = all(string(SRC.outcome)=="Airburst");
all_source_horizon_H360 = all(SRC.source_horizon_s==360);
source_terms_finite_positive = ...
    all(isfinite(SRC.h_burst_km)) && all(SRC.h_burst_km>0) && ...
    all(isfinite(SRC.E_c_kton)) && all(SRC.E_c_kton>0);

% Freeze nominal calibration coefficient and source.
C_values = unique(NOM.C_overpressure_kPa_km3_per_kton);
if numel(C_values)~=1
    error('R3C3B:Calibration','Nominal audit does not contain one fixed coefficient.');
end
C_overpressure = C_values(1);
h_nom_km = unique(NOM.corrected_nominal_h_burst_km);
E_nom_kton = unique(NOM.corrected_nominal_E_c_kton);
if numel(h_nom_km)~=1 || numel(E_nom_kton)~=1
    error('R3C3B:Calibration','Nominal source is not unique.');
end

threshold_kPa = [3.45;1.00;0.50;0.20];
threshold_label = ["3.45 kPa";"1.00 kPa";"0.50 kPa";"0.20 kPa"];
threshold_role = [
    "principal window-damage screening contour";
    "low-overpressure regional screening contour";
    "calibration-anchor / weak-overpressure contour";
    "far-field screening envelope"
];

% Ensure threshold contract agrees with the upstream nominal audit.
nom_thresholds = NOM.threshold_kPa;
threshold_contract_pass = isequal(nom_thresholds(:),threshold_kPa);
nominal_radius_reference = NOM.R_ground_km(:);

N = height(SRC);
nT = numel(threshold_kPa);

% Long-form 4000-row threshold-radius ensemble.
run_id = repelem(SRC.run_id,nT);
fragmentation_seed = repelem(SRC.fragmentation_seed,nT);
source_horizon_s = repelem(SRC.source_horizon_s,nT);
h_burst_km = repelem(SRC.h_burst_km,nT);
E_c_kton = repelem(SRC.E_c_kton,nT);
threshold_kPa_long = repmat(threshold_kPa,N,1);
threshold_label_long = repmat(threshold_label,N,1);
threshold_role_long = repmat(threshold_role,N,1);
C_long = repmat(C_overpressure,N*nT,1);

R_3D_km = (C_long.*E_c_kton./threshold_kPa_long).^(1/3);
R_ground_km = sqrt(max(R_3D_km.^2-h_burst_km.^2,0));
footprint_area_km2 = pi*R_ground_km.^2;
reconstructed_pressure_kPa = C_long.*E_c_kton ./ ...
    max((R_ground_km.^2+h_burst_km.^2).^(3/2),eps);
pressure_reconstruction_abs_error_kPa = ...
    abs(reconstructed_pressure_kPa-threshold_kPa_long);

RadiusLong = table(run_id,fragmentation_seed,source_horizon_s, ...
    h_burst_km,E_c_kton,threshold_label_long,threshold_kPa_long, ...
    threshold_role_long,C_long,R_3D_km,R_ground_km,footprint_area_km2, ...
    reconstructed_pressure_kPa,pressure_reconstruction_abs_error_kPa, ...
    'VariableNames',{'run_id','fragmentation_seed','source_horizon_s', ...
    'h_burst_km','E_c_kton','threshold_label','threshold_kPa', ...
    'threshold_role','C_overpressure_kPa_km3_per_kton','R_3D_km', ...
    'R_ground_km','footprint_area_km2', ...
    'reconstructed_pressure_kPa','pressure_reconstruction_abs_error_kPa'});
writetable(RadiusLong, ...
    fullfile(output_dir,'R3C3B_full1000_threshold_radius_long.csv'));

% Wide matrix for geospatial handoff.
Rmat = reshape(R_ground_km,nT,N)';
Amat = reshape(footprint_area_km2,nT,N)';
PressureWide = table(SRC.run_id,SRC.fragmentation_seed,SRC.h_burst_km, ...
    SRC.E_c_kton, ...
    Rmat(:,1),Rmat(:,2),Rmat(:,3),Rmat(:,4), ...
    Amat(:,1),Amat(:,2),Amat(:,3),Amat(:,4), ...
    'VariableNames',{'run_id','fragmentation_seed','h_burst_km','E_c_kton', ...
    'R_3_45_kPa_km','R_1_00_kPa_km','R_0_50_kPa_km','R_0_20_kPa_km', ...
    'A_3_45_kPa_km2','A_1_00_kPa_km2','A_0_50_kPa_km2','A_0_20_kPa_km2'});
writetable(PressureWide, ...
    fullfile(output_dir,'R3C3B_full1000_pressure_radius_metric_matrix.csv'));

% Per-threshold descriptive statistics.
summary_cells = cell(nT,1);
for j = 1:nT
    x = Rmat(:,j);
    a = Amat(:,j);
    summary_cells{j} = table(threshold_label(j),threshold_kPa(j), ...
        min(x),pct(x,10),median(x),mean(x),pct(x,90),max(x),std(x), ...
        min(a),pct(a,10),median(a),mean(a),pct(a,90),max(a),std(a), ...
        'VariableNames',{'threshold_label','threshold_kPa', ...
        'radius_min_km','radius_p10_km','radius_median_km', ...
        'radius_mean_km','radius_p90_km','radius_max_km','radius_std_km', ...
        'area_min_km2','area_p10_km2','area_median_km2', ...
        'area_mean_km2','area_p90_km2','area_max_km2','area_std_km2'});
end
RadiusSummary = vertcat(summary_cells{:});
writetable(RadiusSummary, ...
    fullfile(output_dir,'R3C3B_threshold_radius_summary.csv'));

% Pressure-profile percentiles from 0 to 200 km ground range.
range_km = (0:1:200)';
nR = numel(range_km);
Pmat = nan(N,nR);
for i = 1:N
    slant = sqrt(range_km.^2+SRC.h_burst_km(i)^2);
    Pmat(i,:) = (C_overpressure*SRC.E_c_kton(i)./slant.^3)';
end

Pmin = min(Pmat,[],1)';
P10 = column_percentile(Pmat,10);
P50 = column_percentile(Pmat,50);
P90 = column_percentile(Pmat,90);
Pmax = max(Pmat,[],1)';
PressureProfile = table(range_km,Pmin,P10,P50,P90,Pmax, ...
    'VariableNames',{'ground_range_km','pressure_min_kPa', ...
    'pressure_p10_kPa','pressure_median_kPa','pressure_p90_kPa', ...
    'pressure_max_kPa'});
writetable(PressureProfile, ...
    fullfile(output_dir,'R3C3B_pressure_profile_percentiles.csv'));

reference_range_km = [0;25;50;75;90;100;125;150;200];
reference_cells = cell(numel(reference_range_km),1);
for j = 1:numel(reference_range_km)
    r = reference_range_km(j);
    p = C_overpressure*SRC.E_c_kton ./ ...
        (r.^2+SRC.h_burst_km.^2).^(3/2);
    reference_cells{j} = table(r,min(p),pct(p,10),median(p), ...
        mean(p),pct(p,90),max(p),std(p), ...
        'VariableNames',{'ground_range_km','pressure_min_kPa', ...
        'pressure_p10_kPa','pressure_median_kPa','pressure_mean_kPa', ...
        'pressure_p90_kPa','pressure_max_kPa','pressure_std_kPa'});
end
ReferencePressure = vertcat(reference_cells{:});
writetable(ReferencePressure, ...
    fullfile(output_dir,'R3C3B_pressure_at_reference_ranges.csv'));

% Rank-correlation diagnostics. These are descriptive, not causal.
input_name = ["h_burst_km";"E_c_kton"];
radius_name = ["R_3_45_kPa_km";"R_1_00_kPa_km"; ...
    "R_0_50_kPa_km";"R_0_20_kPa_km"];
corr_cells = cell(8,1);
k = 0;
for i = 1:2
    if i==1, x=SRC.h_burst_km; else, x=SRC.E_c_kton; end
    for j = 1:4
        k = k+1;
        corr_cells{k}=table(input_name(i),radius_name(j), ...
            spearman_no_toolbox(x,Rmat(:,j)), ...
            'VariableNames',{'source_metric','radius_metric','spearman_rho'});
    end
end
Correlation = vertcat(corr_cells{:});
writetable(Correlation, ...
    fullfile(output_dir,'R3C3B_source_to_radius_spearman.csv'));

% Nominal reproduction against upstream audit.
R3D_nom = (C_overpressure*E_nom_kton./threshold_kPa).^(1/3);
Rground_nom = sqrt(max(R3D_nom.^2-h_nom_km.^2,0));
nominal_abs_difference_km = abs(Rground_nom-nominal_radius_reference);
NominalReproduction = table(threshold_kPa,nominal_radius_reference, ...
    Rground_nom,nominal_abs_difference_km, ...
    nominal_abs_difference_km<=1e-9, ...
    'VariableNames',{'threshold_kPa','upstream_nominal_R_ground_km', ...
    'R3C3B_nominal_R_ground_km','absolute_difference_km', ...
    'reproduction_pass'});
writetable(NominalReproduction, ...
    fullfile(output_dir,'R3C3B_nominal_radius_reproduction_audit.csv'));

% Per-run threshold ordering and reconstruction audit.
ordering_pass = false(N,1);
nonnegative_pass = false(N,1);
finite_pass = false(N,1);
pairwise_inverse_pass = false(N,1);
for i = 1:N
    rr = Rmat(i,:)';
    ordering_pass(i) = all(diff(rr)>0);
    nonnegative_pass(i) = all(rr>=0);
    finite_pass(i) = all(isfinite(rr));
    pairwise_inverse_pass(i) = pairwise_inverse_order(threshold_kPa,rr);
end
RunAudit = table(SRC.run_id,ordering_pass,nonnegative_pass, ...
    finite_pass,pairwise_inverse_pass, ...
    'VariableNames',{'run_id','radius_sequence_ascending_pass', ...
    'nonnegative_radius_pass','finite_radius_pass', ...
    'pairwise_inverse_pressure_radius_order_pass'});
writetable(RunAudit, ...
    fullfile(output_dir,'R3C3B_per_run_radius_contract_audit.csv'));

% Gates.
all_radius_rows_present = height(RadiusLong)==4000;
all_radii_finite = all(isfinite(R_ground_km));
all_radii_nonnegative = all(R_ground_km>=0);
all_radii_positive = all(R_ground_km>0);
all_run_ordering_pass = all(ordering_pass);
all_pairwise_inverse_pass = all(pairwise_inverse_pass);
pressure_reconstruction_pass = ...
    max(pressure_reconstruction_abs_error_kPa)<=1e-10;
nominal_reproduction_pass = all(NominalReproduction.reproduction_pass);
coefficient_global_fixed_pass = all(C_long==C_overpressure);
profile_percentiles_monotonic_pass = ...
    all(diff(P10)<=1e-12) && all(diff(P50)<=1e-12) && ...
    all(diff(P90)<=1e-12);
population_exposure_executed = false;
terminal_mass_metric_publishable = false;

gate_name = [
    "upstream_status_pass";
    "upstream_all_gates_pass";
    "source_row_count_pass";
    "run_id_contract_pass";
    "fragmentation_seed_contract_pass";
    "all_airburst";
    "all_source_horizon_H360";
    "source_terms_finite_positive";
    "threshold_contract_pass";
    "all_4000_threshold_radius_rows_present";
    "all_radii_finite";
    "all_radii_nonnegative";
    "all_radii_positive";
    "all_run_radius_sequences_ascending";
    "all_run_pairwise_inverse_pressure_radius_order";
    "pressure_reconstruction_pass";
    "nominal_radius_reproduction_pass";
    "global_fixed_coefficient_pass";
    "pressure_profile_percentiles_monotonic";
    "terminal_mass_metric_publishable_false";
    "population_exposure_executed_false"
];
gate_pass = [
    upstream_status_pass;
    upstream_all_gates_pass;
    source_row_count_pass;
    run_id_contract_pass;
    seed_contract_pass;
    all_airburst;
    all_source_horizon_H360;
    source_terms_finite_positive;
    threshold_contract_pass;
    all_radius_rows_present;
    all_radii_finite;
    all_radii_nonnegative;
    all_radii_positive;
    all_run_ordering_pass;
    all_pairwise_inverse_pass;
    pressure_reconstruction_pass;
    nominal_reproduction_pass;
    coefficient_global_fixed_pass;
    profile_percentiles_monotonic_pass;
    ~terminal_mass_metric_publishable;
    ~population_exposure_executed
];
gate_detail = [
    "R3C3A-R1 status is PASS";
    sprintf('upstream_gate_pass_count=%d/%d',sum(truth_vector(GATE.gate_pass)),height(GATE));
    sprintf('source_rows=%d',height(SRC));
    sprintf('unique_run_ids=%d',numel(unique(SRC.run_id)));
    sprintf('unique_fragmentation_seeds=%d',numel(unique(SRC.fragmentation_seed)));
    sprintf('airburst_count=%d',sum(string(SRC.outcome)=="Airburst"));
    sprintf('H360_count=%d',sum(SRC.source_horizon_s==360));
    sprintf('finite_positive_source_count=%d',sum(isfinite(SRC.h_burst_km)& ...
        SRC.h_burst_km>0&isfinite(SRC.E_c_kton)&SRC.E_c_kton>0));
    sprintf('thresholds=%s',strjoin(compose('%.9g',threshold_kPa),';'));
    sprintf('radius_rows=%d',height(RadiusLong));
    sprintf('finite_radius_count=%d',sum(isfinite(R_ground_km)));
    sprintf('minimum_radius_km=%.15g',min(R_ground_km));
    sprintf('strictly_positive_radius_count=%d',sum(R_ground_km>0));
    sprintf('ascending_pass_count=%d/%d',sum(ordering_pass),N);
    sprintf('pairwise_pass_count=%d/%d',sum(pairwise_inverse_pass),N);
    sprintf('maximum_pressure_reconstruction_error_kPa=%.15g', ...
        max(pressure_reconstruction_abs_error_kPa));
    sprintf('nominal_reproduction_pass_count=%d/%d', ...
        sum(NominalReproduction.reproduction_pass),height(NominalReproduction));
    sprintf('C_overpressure=%.15g; unique_C_count=1',C_overpressure);
    "P10, median, and P90 pressure curves are nonincreasing with range";
    "terminal active mass is excluded from the pressure model";
    "population exposure remains blocked until R3C4"
];
GateAudit = table(gate_name,gate_pass,gate_detail);
writetable(GateAudit, ...
    fullfile(output_dir,'R3C3B_gate_adjudication.csv'));

% Diagnostic figures.
make_radius_distribution_figure(Rmat,threshold_label,figure_dir);
make_pressure_profile_figure(PressureProfile,threshold_kPa,figure_dir);
make_source_radius_figure(SRC,Rmat,figure_dir);

all_pass = all(gate_pass);
if all_pass
    status = "PASS_R3C3B_CORRECTED_FULL1000_OVERPRESSURE_PROPAGATION";
else
    status = "FAIL_R3C3B_CORRECTED_FULL1000_OVERPRESSURE_PROPAGATION";
end

fid = fopen(fullfile(report_dir,'R3C3B_status.txt'),'w');
fprintf(fid,'=== CHELYABINSK R3C3B RESULT ===\n');
fprintf(fid,'status=%s\n',status);
fprintf(fid,'scientific_entry_propagation_executed=False\n');
fprintf(fid,'overpressure_ensemble_propagation_executed=True\n');
fprintf(fid,'population_exposure_executed=False\n');
fprintf(fid,'full1000_source_row_count=%d\n',height(SRC));
fprintf(fid,'threshold_radius_row_count=%d\n',height(RadiusLong));
fprintf(fid,'calibration_mode=GLOBAL_FIXED_NOMINAL_ANCHOR\n');
fprintf(fid,'per_realization_recalibration=False\n');
fprintf(fid,'C_overpressure_kPa_km3_per_kton=%.15g\n',C_overpressure);
fprintf(fid,'all_run_radius_sequences_ascending=%s\n',string(all_run_ordering_pass));
fprintf(fid,'all_run_pairwise_inverse_order_pass=%s\n',string(all_pairwise_inverse_pass));
fprintf(fid,'pressure_reconstruction_pass=%s\n',string(pressure_reconstruction_pass));
fprintf(fid,'nominal_radius_reproduction_pass=%s\n',string(nominal_reproduction_pass));
fprintf(fid,'radius_3_45_p10_km=%.15g\n',RadiusSummary.radius_p10_km(1));
fprintf(fid,'radius_3_45_median_km=%.15g\n',RadiusSummary.radius_median_km(1));
fprintf(fid,'radius_3_45_p90_km=%.15g\n',RadiusSummary.radius_p90_km(1));
fprintf(fid,'radius_1_00_p10_km=%.15g\n',RadiusSummary.radius_p10_km(2));
fprintf(fid,'radius_1_00_median_km=%.15g\n',RadiusSummary.radius_median_km(2));
fprintf(fid,'radius_1_00_p90_km=%.15g\n',RadiusSummary.radius_p90_km(2));
fprintf(fid,'radius_0_50_p10_km=%.15g\n',RadiusSummary.radius_p10_km(3));
fprintf(fid,'radius_0_50_median_km=%.15g\n',RadiusSummary.radius_median_km(3));
fprintf(fid,'radius_0_50_p90_km=%.15g\n',RadiusSummary.radius_p90_km(3));
fprintf(fid,'radius_0_20_p10_km=%.15g\n',RadiusSummary.radius_p10_km(4));
fprintf(fid,'radius_0_20_median_km=%.15g\n',RadiusSummary.radius_median_km(4));
fprintf(fid,'radius_0_20_p90_km=%.15g\n',RadiusSummary.radius_p90_km(4));
fprintf(fid,'terminal_mass_metric_publishable=False\n');
fprintf(fid,'pressure_statistics_publishable=%s\n',string(all_pass));
fprintf(fid,'candidate_manuscript_modified=False\n');
fprintf(fid,'canonical_modification_executed=False\n');
fprintf(fid,'next_stage=Upload R3C3B review bundle for scientific and visual adjudication; construct R3C4 only after PASS.\n');
fclose(fid);

fprintf('\n=== R3C3B RESULT ===\n');
fprintf('status=%s\n',status);
fprintf('full1000_source_row_count=%d\n',height(SRC));
fprintf('threshold_radius_row_count=%d\n',height(RadiusLong));
fprintf('all_run_radius_sequences_ascending=%s\n',string(all_run_ordering_pass));
fprintf('all_run_pairwise_inverse_order_pass=%s\n',string(all_pairwise_inverse_pass));
fprintf('pressure_reconstruction_pass=%s\n',string(pressure_reconstruction_pass));
fprintf('nominal_radius_reproduction_pass=%s\n',string(nominal_reproduction_pass));
fprintf('pressure_statistics_publishable=%s\n',string(all_pass));
fprintf('population_exposure_executed=False\n');

if ~all_pass
    error('R3C3B:GateFailure', ...
        'Corrected full1000 overpressure propagation gate failed.');
end
end

function T=read_csv_explicit(path)
T=readtable(path,'FileType','text','Delimiter',',', ...
    'ReadVariableNames',true,'TextType','string', ...
    'VariableNamingRule','preserve','Encoding','UTF-8');
names=string(T.Properties.VariableNames);
names=replace(names,char(65279),"");
names=strip(names);
names=matlab.lang.makeValidName(names,'ReplacementStyle','underscore');
names=matlab.lang.makeUniqueStrings(names);
T.Properties.VariableNames=cellstr(names);
end

function assert_columns(T,required,label)
observed=string(T.Properties.VariableNames);
missing=required(~ismember(required,observed));
if ~isempty(missing)
    error('R3C3B:Schema','Missing columns in %s: %s', ...
        char(string(label)),char(strjoin(missing,',')));
end
end

function v=truth_vector(x)
if islogical(x), v=x;
elseif isnumeric(x), v=x~=0;
else
    s=lower(strip(string(x)));
    v=s=="1"|s=="true"|s=="yes"|s=="pass";
end
v=v(:);
end

function pass=pairwise_inverse_order(P,R)
pass=true;
for i=1:numel(P)
    for j=i+1:numel(P)
        if P(i)>P(j), pass=pass&&R(i)<R(j);
        elseif P(i)<P(j), pass=pass&&R(i)>R(j);
        else, pass=pass&&abs(R(i)-R(j))<=1e-12;
        end
    end
end
end

function q=pct(x,p)
x=sort(x(:)); n=numel(x);
if n==0, q=NaN; return; end
if n==1, q=x; return; end
pos=1+(n-1)*p/100; lo=floor(pos); hi=ceil(pos);
if lo==hi, q=x(lo); else, q=x(lo)+(pos-lo)*(x(hi)-x(lo)); end
end

function q=column_percentile(X,p)
q=nan(size(X,2),1);
for j=1:size(X,2), q(j)=pct(X(:,j),p); end
end

function rho=spearman_no_toolbox(x,y)
rx=average_rank(x); ry=average_rank(y);
rx=rx-mean(rx); ry=ry-mean(ry);
den=sqrt(sum(rx.^2)*sum(ry.^2));
if den==0, rho=NaN; else, rho=sum(rx.*ry)/den; end
end

function r=average_rank(x)
[xs,ord]=sort(x(:)); n=numel(x); rs=zeros(n,1); i=1;
while i<=n
    j=i;
    while j<n && xs(j+1)==xs(i), j=j+1; end
    rs(i:j)=0.5*(i+j); i=j+1;
end
r=zeros(n,1); r(ord)=rs;
end

function make_radius_distribution_figure(R,labels,folder)
fig=figure('Visible','off','Color','w','Position',[80 80 1300 850]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for j=1:4
    nexttile(tl);
    histogram(R(:,j),25);
    grid on;
    xlabel('Ground-radius surrogate [km]');
    ylabel('Ensemble count');
    title(labels(j));
end
title(tl,'Corrected full-1000 overpressure-threshold radii');
exportgraphics(fig,fullfile(folder, ...
    'R3C3B_full1000_threshold_radius_distributions.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(folder, ...
    'R3C3B_full1000_threshold_radius_distributions.png'),'Resolution',300);
close(fig);
end

function make_pressure_profile_figure(T,thresholds,folder)
fig=figure('Visible','off','Color','w','Position',[80 80 1050 700]);
semilogy(T.ground_range_km,T.pressure_median_kPa,'LineWidth',1.8);
hold on;
semilogy(T.ground_range_km,T.pressure_p10_kPa,'--','LineWidth',1.2);
semilogy(T.ground_range_km,T.pressure_p90_kPa,'--','LineWidth',1.2);
for j=1:numel(thresholds)
    yline(thresholds(j),':',sprintf('%.2f kPa',thresholds(j)), ...
        'LineWidth',0.9);
end
grid on;
xlabel('Ground range from burst projection [km]');
ylabel('Peak-overpressure surrogate [kPa]');
title('Corrected full-1000 pressure-profile percentiles');
legend({'Median','P10','P90'},'Location','southwest');
xlim([0 200]); ylim([0.02 150]);
exportgraphics(fig,fullfile(folder, ...
    'R3C3B_pressure_profile_percentile_envelope.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(folder, ...
    'R3C3B_pressure_profile_percentile_envelope.png'),'Resolution',300);
close(fig);
end

function make_source_radius_figure(S,R,folder)
fig=figure('Visible','off','Color','w','Position',[80 80 1050 720]);
scatter(S.E_c_kton,R(:,1),28,S.h_burst_km,'filled');
grid on;
xlabel('Coupled burst energy [kt TNT]');
ylabel('3.45 kPa ground-radius surrogate [km]');
title('Corrected source-to-radius mapping');
cb=colorbar;
cb.Label.String='Main airburst altitude [km]';
exportgraphics(fig,fullfile(folder, ...
    'R3C3B_source_to_R345_mapping.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(folder, ...
    'R3C3B_source_to_R345_mapping.png'),'Resolution',300);
close(fig);
end

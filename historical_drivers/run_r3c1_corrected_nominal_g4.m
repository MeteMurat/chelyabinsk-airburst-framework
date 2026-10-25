function run_r3c1_corrected_nominal_g4(artifact_root,expected_summary_file,convergence_file)
% CHELYABINSK R3C1 CORRECTED NOMINAL DETERMINISTIC G4 REGENERATION
% Candidate-only propagation and artifact generation.
% No Monte Carlo, overpressure, exposure, manuscript, or canonical changes.

if nargin < 1 || strlength(string(artifact_root)) == 0
    error('artifact_root is required.');
end
if nargin < 2 || ~isfile(expected_summary_file)
    error('A valid frozen R3B3 factorial summary is required.');
end
if nargin < 3 || ~isfile(convergence_file)
    error('A valid frozen R3B3 convergence audit is required.');
end

artifact_root = char(string(artifact_root));
expected_summary_file = char(string(expected_summary_file));
convergence_file = char(string(convergence_file));

core_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(artifact_root,'02_nominal_outputs');
figure_dir = fullfile(output_dir,'figures');
table_dir = fullfile(output_dir,'tables');
ledger_dir = fullfile(output_dir,'ledgers');
report_dir = fullfile(artifact_root,'03_reports');

for p = {output_dir,figure_dir,table_dir,ledger_dir,report_dir}
    if ~exist(p{1},'dir'), mkdir(p{1}); end
end

addpath(core_dir);
format long g;

fprintf('=== CHELYABINSK R3C1 CORRECTED NOMINAL G4 REGENERATION ===\n');
fprintf('scientific_propagation_requested=True\n');
fprintf('scientific_model_source_repair_requested=False\n');
fprintf('fragment_capacity_contract_repair_applied=True\n');
fprintf('authoritative_source_modification_requested=False\n');
fprintf('candidate_manuscript_modification_requested=False\n');
fprintf('canonical_modification_requested=False\n');
fprintf('controlled_nominal_case=G4\n');
fprintf('fragment_capacity=128\n');
fprintf('terminal_complete_horizon_s=360\n');
fprintf('seed=42\n');
fprintf('monte_carlo_requested=False\n');
fprintf('overpressure_requested=False\n');
fprintf('population_exposure_requested=False\n');
fprintf('artifact_root=%s\n',artifact_root);

%% Frozen corrected nominal contract.
params.rho_meteor = 3300;
params.V0 = 19160;
params.entry_angle_deg = 18.3;
params.D0 = 19.8;
params.area_multiplier = 1.17161624727046;
params.h0 = 97100;
params.Yield_Stress_0 = 7.8e6;
params.min_strength_fraction = 0.45;
params.weakening_mass_coeff = 0.18;
params.Cd = 1.0;
params.Ch = 0.015;
params.Qabl = 8e6;
params.max_generations = 3;
params.max_total_fragments = 128;
params.base_fragment_count = 5;
params.dispersion_speed = 250;
params.cloud_area_gain = 2.5;
params.dt = 0.02;
params.tmax = 360;
params.energy_efficiency = 0.80;
params.energy_bin_width_m = 500;
params.seed = 42;

%% Authoritative and diagnostic parity propagation.
out_ref = meteor_single_run(params);
rng_after_ref = rng;
out_diag = meteor_single_run_r3b3_diagnostic(params);
rng_after_diag = rng;

scalar_fields = [ ...
    "first_break_h_km","first_break_time","main_airburst_h_km", ...
    "main_airburst_t","E_release_kton","E_total_lost_kton", ...
    "dE_dh_peak_kton_per_km","E_burst_raw_kton", ...
    "E_burst_coupled_kton","airburst_window_h_min_km", ...
    "airburst_window_h_max_km","final_mass","area_multiplier"];

parity_rows = cell(numel(scalar_fields)+3,1);
all_scalar_pass = true;
max_scalar_abs_diff = 0;
for i = 1:numel(scalar_fields)
    f = scalar_fields(i);
    rv = double(out_ref.(f));
    dv = double(out_diag.(f));
    abs_diff = abs(rv-dv);
    tol = max(1e-10,1e-11*max(1,abs(rv)));
    pass = (isnan(rv) && isnan(dv)) || abs_diff <= tol;
    all_scalar_pass = all_scalar_pass && pass;
    max_scalar_abs_diff = max(max_scalar_abs_diff,abs_diff);
    parity_rows{i} = table(f,rv,dv,abs_diff,tol,pass, ...
        'VariableNames',{'field','authoritative_value','diagnostic_value', ...
        'absolute_difference','tolerance','parity_pass'});
end

[ts_pass,ts_diff] = table_numeric_parity(out_ref.timeseries,out_diag.timeseries,1e-8);
[dep_pass,dep_diff] = table_numeric_parity(out_ref.deposition_profile,out_diag.deposition_profile,1e-8);
rng_pass = strcmp(rng_after_ref.Type,rng_after_diag.Type) && ...
    isequal(rng_after_ref.Seed,rng_after_diag.Seed) && ...
    isequal(rng_after_ref.State,rng_after_diag.State);

parity_rows{numel(scalar_fields)+1} = table("timeseries",NaN,NaN,ts_diff,1e-8,ts_pass, ...
    'VariableNames',{'field','authoritative_value','diagnostic_value', ...
    'absolute_difference','tolerance','parity_pass'});
parity_rows{numel(scalar_fields)+2} = table("deposition_profile",NaN,NaN,dep_diff,1e-8,dep_pass, ...
    'VariableNames',{'field','authoritative_value','diagnostic_value', ...
    'absolute_difference','tolerance','parity_pass'});
parity_rows{numel(scalar_fields)+3} = table("rng_state",NaN,NaN,0,0,rng_pass, ...
    'VariableNames',{'field','authoritative_value','diagnostic_value', ...
    'absolute_difference','tolerance','parity_pass'});
parity_audit = vertcat(parity_rows{:});
writetable(parity_audit,fullfile(report_dir,'R3C1_authoritative_diagnostic_parity.csv'));

authoritative_diagnostic_numeric_parity = all_scalar_pass && ts_pass && dep_pass;
rng_sequence_preserved = rng_pass;

%% Frozen R3B3 metric reproduction.
expected = readtable(expected_summary_file,'TextType','string');
expected = expected(expected.case_id=="G4" & expected.regime_id=="R128_T360",:);
if height(expected) ~= 1
    error('A unique frozen G4/R128_T360 row was not found.');
end

metric_name = [ ...
    "area_multiplier";"first_break_h_km";"main_airburst_h_km"; ...
    "E_total_lost_kton";"E_release_kton";"E_burst_coupled_kton"; ...
    "dE_dh_peak_kton_per_km";"q_dyn_max_MPa"];
observed_value = [ ...
    params.area_multiplier;out_ref.first_break_h_km;out_ref.main_airburst_h_km; ...
    out_ref.E_total_lost_kton;out_ref.E_release_kton;out_ref.E_burst_coupled_kton; ...
    out_ref.dE_dh_peak_kton_per_km;max(out_ref.timeseries.q_dyn_max_Pa)/1e6];
expected_value = [ ...
    expected.area_multiplier;expected.first_break_h_km;expected.main_airburst_h_km; ...
    expected.E_total_lost_kton;expected.E_release_kton;expected.E_burst_coupled_kton; ...
    expected.dE_dh_peak_kton_per_km;expected.q_dyn_max_MPa];

absolute_difference = abs(observed_value-expected_value);
tolerance = max(1e-9*ones(size(expected_value)),1e-10*max(1,abs(expected_value)));
metric_pass = absolute_difference <= tolerance;
frozen_metric_audit = table(metric_name,expected_value,observed_value, ...
    absolute_difference,tolerance,metric_pass);
writetable(frozen_metric_audit,fullfile(report_dir,'R3C1_frozen_G4_metric_reproduction.csv'));
frozen_r3b3_g4_metric_match = all(metric_pass);

conv = readtable(convergence_file,'TextType','string');
conv = conv(conv.case_id=="G4",:);
if height(conv) ~= 1
    error('A unique G4 convergence row was not found.');
end
local_airburst_horizon_convergence_reused = ...
    logical(conv.local_airburst_horizon_convergence_pass(1));

%% Scientific integrity gates.
B = out_diag.r3b.breakup_event_ledger;
C = out_diag.r3b.child_creation_ledger;
D = out_diag.r3b.deactivation_ledger;
TL = out_diag.r3b.mass_energy_time_ledger;
S = out_diag.r3b.closure_summary;
F = out_diag.r3b.terminal_fragment_ledger;

created_child_count = height(C);
if isempty(C)
    max_child_id = 1;
    all_children_inserted = false;
else
    max_child_id = max(C.child_id);
    all_children_inserted = all(C.inserted_into_active_pool);
end
breakup_event_count = height(B);
cap_hit_event_count = sum(B.cap_hit);
cap_dropped_mass_kg = sum(B.cap_dropped_mass_kg);
insertion_dropped_mass_kg = S.cumulative_insertion_dropped_mass_kg;
full_tree_realized = created_child_count==85 && max_child_id==86 && ...
    all_children_inserted && breakup_event_count==26;
terminal_complete = isempty(F) && S.final_active_mass_kg <= 1e-9;
ground_crossing_count = sum(D.reason=="post_step_ground_crossing");
mass_relative_residual = abs(S.final_mass_closure_residual_kg) / max(S.initial_mass_kg,eps);
energy_relative_residual = abs(S.final_energy_closure_residual_J) / max(S.initial_ke_J,eps);
mass_closure_pass = mass_relative_residual <= 1e-10;
energy_closure_pass = energy_relative_residual <= 1e-10;
cap_contract_pass = cap_hit_event_count==0 && abs(cap_dropped_mass_kg)<=1e-9 && ...
    abs(insertion_dropped_mass_kg)<=1e-9;
terminal_mass_metric_requested = false;

%% Derive corrected nominal diagnostics.
TS = out_ref.timeseries;
DEP = out_ref.deposition_profile;

[~,idx_break] = min(abs(TS.t_s-out_ref.first_break_time));
[fragment_count_max,idx_peak_fragment] = max(TS.fragment_count);
[q_dyn_max_Pa,idx_qmax] = max(TS.q_dyn_max_Pa);
[peak_dEdh,idx_peak_dep] = max(DEP.dE_dh_kton_per_km);

first_break_downrange_km = TS.x_km(idx_break);
first_break_speed_km_s = TS.v_m_s(idx_break)/1000;
first_break_q_MPa = B.q_Pa(1)/1e6;
peak_fragment_time_s = TS.t_s(idx_peak_fragment);
peak_fragment_altitude_km = TS.h_km(idx_peak_fragment);
qmax_time_s = TS.t_s(idx_qmax);
qmax_altitude_km = TS.h_km(idx_qmax);
peak_bin_energy_kton = DEP.dE_bin_kton(idx_peak_dep);
binned_energy_integral_kton = sum(DEP.dE_bin_kton);
peak_bin_share_pct = 100*peak_bin_energy_kton/max(binned_energy_integral_kton,eps);
burst_window_raw_kton = out_ref.E_burst_raw_kton;
terminal_completion_time_s = TL.time_s(end);
initial_ke_kton = S.initial_ke_J/4.184e12;

entry_name = [ ...
    "Initial altitude";"Initial speed";"Initial mass";"Initial kinetic energy"; ...
    "First-breakup time";"First-breakup altitude";"First-breakup downrange"; ...
    "Speed at first breakup";"Dynamic pressure at first breakup"; ...
    "Peak active fragment count";"Time of peak active fragment count"; ...
    "Altitude at peak active fragment count";"Maximum dynamic pressure"; ...
    "Time of maximum dynamic pressure";"Altitude at maximum dynamic pressure"; ...
    "Breakup-event count";"Created-child count";"Maximum fragment ID"; ...
    "Ground-crossing count";"Terminal-completion time"];
entry_symbol = [ ...
    "h_0";"V_0";"m_0";"E_{k,0}";"t_{break}";"h_{break}";"x_{break}"; ...
    "v_{break}";"q_{dyn,break}";"N_{f,max}";"t(N_{f,max})"; ...
    "h(N_{f,max})";"q_{dyn,max}";"t(q_{dyn,max})";"h(q_{dyn,max})"; ...
    "N_{break}";"N_{child}";"ID_{max}";"N_{ground}";"t_{terminal}"];
entry_value = [ ...
    params.h0/1000;params.V0/1000;S.initial_mass_kg;initial_ke_kton; ...
    out_ref.first_break_time;out_ref.first_break_h_km;first_break_downrange_km; ...
    first_break_speed_km_s;first_break_q_MPa;fragment_count_max; ...
    peak_fragment_time_s;peak_fragment_altitude_km;q_dyn_max_Pa/1e6; ...
    qmax_time_s;qmax_altitude_km;breakup_event_count;created_child_count; ...
    max_child_id;ground_crossing_count;terminal_completion_time_s];
entry_unit = [ ...
    "km";"km s^-1";"kg";"kt TNT";"s";"km";"km";"km s^-1";"MPa"; ...
    "--";"s";"km";"MPa";"s";"km";"--";"--";"--";"--";"s"];
entry_table = table(entry_name,entry_symbol,entry_value,entry_unit, ...
    'VariableNames',{'diagnostic','symbol','value','unit'});
writetable(entry_table,fullfile(table_dir,'R3C1_nominal_entry_fragmentation_diagnostics.csv'));

energy_name = [ ...
    "Main airburst altitude";"Airburst-window lower altitude"; ...
    "Airburst-window upper altitude";"Main airburst midpoint time"; ...
    "Peak-bin deposited energy";"Peak binned deposition density"; ...
    "Burst-window raw energy";"Direct kinetic-energy loss"; ...
    "Binned deposited-energy integral";"Peak-bin share of binned energy"; ...
    "Coupling factor";"Nominal coupled burst energy";"Released energy"];
energy_symbol = [ ...
    "h_{burst}";"h_{window,min}";"h_{window,max}";"t_{burst}"; ...
    "Delta E_{peak-bin}";"(dE/dh)_{max}";"E_{burst,raw}";"E_{lost}"; ...
    "sum Delta E_b";"f_{peak-bin}";"eta_c";"E_{c,nom}";"E_{release}"];
energy_value = [ ...
    out_ref.main_airburst_h_km;out_ref.airburst_window_h_min_km; ...
    out_ref.airburst_window_h_max_km;out_ref.main_airburst_t; ...
    peak_bin_energy_kton;peak_dEdh;burst_window_raw_kton; ...
    out_ref.E_total_lost_kton;binned_energy_integral_kton;peak_bin_share_pct; ...
    params.energy_efficiency;out_ref.E_burst_coupled_kton;out_ref.E_release_kton];
energy_unit = [ ...
    "km";"km";"km";"s";"kt TNT";"kt km^-1";"kt TNT";"kt TNT"; ...
    "kt TNT";"%";"--";"kt TNT";"kt TNT"];
energy_table = table(energy_name,energy_symbol,energy_value,energy_unit, ...
    'VariableNames',{'diagnostic','symbol','value','unit'});
writetable(energy_table,fullfile(table_dir,'R3C1_nominal_energy_deposition_diagnostics.csv'));

summary = table( ...
    "G4","existing_nominal_mean_area",params.area_multiplier,params.seed, ...
    params.max_total_fragments,params.tmax,out_ref.first_break_h_km, ...
    out_ref.main_airburst_h_km,out_ref.E_total_lost_kton,out_ref.E_release_kton, ...
    out_ref.E_burst_raw_kton,out_ref.E_burst_coupled_kton, ...
    out_ref.dE_dh_peak_kton_per_km,q_dyn_max_Pa/1e6,fragment_count_max, ...
    breakup_event_count,created_child_count,max_child_id,cap_hit_event_count, ...
    cap_dropped_mass_kg,insertion_dropped_mass_kg,ground_crossing_count, ...
    terminal_completion_time_s,full_tree_realized,terminal_complete, ...
    mass_relative_residual,energy_relative_residual,string(out_ref.outcome), ...
    'VariableNames',{'case_id','case_role','area_multiplier','seed', ...
    'max_total_fragments','tmax_s','first_break_h_km','main_airburst_h_km', ...
    'E_total_lost_kton','E_release_kton','E_burst_raw_kton', ...
    'E_burst_coupled_kton','dE_dh_peak_kton_per_km','q_dyn_max_MPa', ...
    'fragment_count_max','breakup_event_count','created_child_count', ...
    'max_child_id','cap_hit_event_count','cap_dropped_mass_kg', ...
    'insertion_dropped_mass_kg','ground_crossing_count', ...
    'terminal_completion_time_s','full_tree_realized','terminal_complete', ...
    'mass_relative_residual','energy_relative_residual','outcome'});
writetable(summary,fullfile(output_dir,'R3C1_corrected_nominal_G4_summary.csv'));

%% Freeze raw outputs and ledgers.
writetable(TS,fullfile(ledger_dir,'R3C1_nominal_timeseries.csv'));
writetable(DEP,fullfile(ledger_dir,'R3C1_nominal_deposition_profile.csv'));
writetable(B,fullfile(ledger_dir,'R3C1_breakup_event_ledger.csv'));
writetable(C,fullfile(ledger_dir,'R3C1_child_creation_ledger.csv'));
writetable(D,fullfile(ledger_dir,'R3C1_deactivation_ledger.csv'));
writetable(TL,fullfile(ledger_dir,'R3C1_mass_energy_time_ledger.csv'));
writetable(S,fullfile(ledger_dir,'R3C1_closure_summary.csv'));

%% Publication-ready nominal trajectory figure.
fig1 = figure('Visible','off','Color','w','Position',[100 100 1400 900]);
tl1 = tiledlayout(fig1,2,3,'TileSpacing','compact','Padding','compact');

nexttile(tl1);
plot(TS.t_s,TS.h_km,'LineWidth',1.4); grid on;
xlabel('Time [s]'); ylabel('Mean altitude [km]');
xline(out_ref.first_break_time,'--','First breakup');

nexttile(tl1);
plot(TS.x_km,TS.h_km,'LineWidth',1.4); grid on;
xlabel('Downrange [km]'); ylabel('Mean altitude [km]');

nexttile(tl1);
plot(TS.t_s,TS.v_m_s/1000,'LineWidth',1.4); grid on;
xlabel('Time [s]'); ylabel('Mean speed [km s^{-1}]');

nexttile(tl1);
plot(TS.t_s,TS.mass_kg/1e6,'LineWidth',1.4); grid on;
xlabel('Time [s]'); ylabel('Active atmospheric mass [10^6 kg]');

nexttile(tl1);
plot(TS.t_s,TS.q_dyn_max_Pa/1e6,'LineWidth',1.4); grid on;
xlabel('Time [s]'); ylabel('Maximum dynamic pressure [MPa]');
xline(TS.t_s(idx_qmax),'--','q_{max}');

nexttile(tl1);
stairs(TS.t_s,TS.fragment_count,'LineWidth',1.4); grid on;
xlabel('Time [s]'); ylabel('Active fragment count');
xline(out_ref.first_break_time,'--','First breakup');

title(tl1,'Corrected nominal G4 entry trajectory (full 86-ID fragmentation tree)');
exportgraphics(fig1,fullfile(figure_dir,'fig08_nominal_trajectory.pdf'),'ContentType','vector');
exportgraphics(fig1,fullfile(figure_dir,'fig08_nominal_trajectory.png'),'Resolution',300);
close(fig1);

%% Publication-ready energy-deposition figure.
fig2 = figure('Visible','off','Color','w','Position',[100 100 1200 720]);
tl2 = tiledlayout(fig2,1,2,'TileSpacing','compact','Padding','compact');

nexttile(tl2);
plot(DEP.dE_dh_kton_per_km,DEP.h_center_km,'LineWidth',1.5); grid on;
xlabel('Energy-deposition density [kt km^{-1}]');
ylabel('Altitude [km]');
yline(out_ref.main_airburst_h_km,'--','Main airburst');
yline(out_ref.airburst_window_h_min_km,':','Window lower bound');
yline(out_ref.airburst_window_h_max_km,':','Window upper bound');

nexttile(tl2);
plot(DEP.dE_bin_kton,DEP.h_center_km,'LineWidth',1.5); grid on;
xlabel('Deposited energy per 0.5-km bin [kt]');
ylabel('Altitude [km]');
yline(out_ref.main_airburst_h_km,'--','Main airburst');

title(tl2,'Corrected nominal G4 altitude-resolved energy deposition');
exportgraphics(fig2,fullfile(figure_dir,'fig09_nominal_energy_deposition.pdf'),'ContentType','vector');
exportgraphics(fig2,fullfile(figure_dir,'fig09_nominal_energy_deposition.png'),'Resolution',300);
close(fig2);

%% Publication-ready fragmentation-sequence figure.
fig3 = figure('Visible','off','Color','w','Position',[100 100 1200 800]);
tl3 = tiledlayout(fig3,2,1,'TileSpacing','compact','Padding','compact');

nexttile(tl3);
stairs(TS.t_s,TS.fragment_count,'LineWidth',1.5); grid on;
xlabel('Time [s]'); ylabel('Active fragment count');
xline(out_ref.first_break_time,'--','First breakup');
yline(fragment_count_max,':',sprintf('Peak = %d',fragment_count_max));

nexttile(tl3);
scatter(B.time_s,B.altitude_km,36,B.generation+1,'filled'); grid on;
xlabel('Time [s]'); ylabel('Breakup-event altitude [km]');
cb = colorbar; cb.Label.String = 'Parent generation + 1';

title(tl3,'Corrected nominal G4 fragmentation sequence');
exportgraphics(fig3,fullfile(figure_dir,'fig10_nominal_fragmentation_sequence.pdf'),'ContentType','vector');
exportgraphics(fig3,fullfile(figure_dir,'fig10_nominal_fragmentation_sequence.png'),'Resolution',300);
close(fig3);

%% Candidate-ready LaTeX snippets (not inserted into the manuscript).
write_entry_latex(fullfile(table_dir,'R3C1_nominal_entry_fragmentation_table.tex'), ...
    params,S,out_ref,first_break_downrange_km,first_break_speed_km_s, ...
    first_break_q_MPa,fragment_count_max,peak_fragment_time_s, ...
    peak_fragment_altitude_km,q_dyn_max_Pa/1e6,qmax_time_s, ...
    qmax_altitude_km,breakup_event_count,created_child_count,max_child_id, ...
    terminal_completion_time_s);

write_energy_latex(fullfile(table_dir,'R3C1_nominal_energy_deposition_table.tex'), ...
    out_ref,peak_bin_energy_kton,peak_dEdh,burst_window_raw_kton, ...
    binned_energy_integral_kton,peak_bin_share_pct,params.energy_efficiency);

captions = [
"fig08_nominal_trajectory: Corrected nominal Chelyabinsk-class G4 trajectory under the cap-free 128-slot fragmentation contract. Active atmospheric mass is shown as a time history and is not interpreted as a terminal-mass response.";
"fig09_nominal_energy_deposition: Corrected altitude-resolved energy-deposition profile for the nominal G4 case. The main airburst and half-peak burst window are obtained from the 0.5-km binned profile.";
"fig10_nominal_fragmentation_sequence: Corrected nominal fragmentation history with the complete 86-ID tree. The lower panel records individual breakup events by parent generation."
];
writelines(captions,fullfile(report_dir,'R3C1_proposed_figure_captions.txt'));

%% Final gate adjudication.
figure_files = [ ...
    "fig08_nominal_trajectory.pdf";"fig08_nominal_trajectory.png"; ...
    "fig09_nominal_energy_deposition.pdf";"fig09_nominal_energy_deposition.png"; ...
    "fig10_nominal_fragmentation_sequence.pdf";"fig10_nominal_fragmentation_sequence.png"];
figure_exists = false(numel(figure_files),1);
for i = 1:numel(figure_files)
    figure_exists(i) = isfile(fullfile(figure_dir,figure_files(i)));
end
figure_audit = table(figure_files,figure_exists, ...
    'VariableNames',{'file','exists'});
writetable(figure_audit,fullfile(report_dir,'R3C1_figure_export_audit.csv'));
figure_export_pass = all(figure_exists);

all_gates_pass = authoritative_diagnostic_numeric_parity && ...
    rng_sequence_preserved && frozen_r3b3_g4_metric_match && ...
    local_airburst_horizon_convergence_reused && full_tree_realized && ...
    cap_contract_pass && terminal_complete && mass_closure_pass && ...
    energy_closure_pass && figure_export_pass && ~terminal_mass_metric_requested;

if all_gates_pass
    status = "PASS_R3C1_CORRECTED_NOMINAL_G4_REGENERATION";
else
    status = "FAIL_R3C1_CORRECTED_NOMINAL_G4_REGENERATION";
end

status_file = fullfile(report_dir,'R3C1_status.txt');
fid = fopen(status_file,'w');
fprintf(fid,'=== CHELYABINSK R3C1 RESULT ===\n');
fprintf(fid,'status=%s\n',status);
fprintf(fid,'scientific_propagation_executed=True\n');
fprintf(fid,'scientific_model_source_repair_executed=False\n');
fprintf(fid,'fragment_capacity_contract_repair_applied=True\n');
fprintf(fid,'authoritative_diagnostic_numeric_parity=%s\n',string(authoritative_diagnostic_numeric_parity));
fprintf(fid,'rng_sequence_preserved=%s\n',string(rng_sequence_preserved));
fprintf(fid,'frozen_r3b3_g4_metric_match=%s\n',string(frozen_r3b3_g4_metric_match));
fprintf(fid,'local_airburst_horizon_convergence_reused=%s\n',string(local_airburst_horizon_convergence_reused));
fprintf(fid,'created_child_count=%d\n',created_child_count);
fprintf(fid,'max_child_id=%d\n',max_child_id);
fprintf(fid,'breakup_event_count=%d\n',breakup_event_count);
fprintf(fid,'fragment_count_max=%d\n',fragment_count_max);
fprintf(fid,'cap_hit_event_count=%d\n',cap_hit_event_count);
fprintf(fid,'cap_dropped_mass_kg=%.15g\n',cap_dropped_mass_kg);
fprintf(fid,'insertion_dropped_mass_kg=%.15g\n',insertion_dropped_mass_kg);
fprintf(fid,'full_tree_realized=%s\n',string(full_tree_realized));
fprintf(fid,'terminal_complete=%s\n',string(terminal_complete));
fprintf(fid,'terminal_completion_time_s=%.15g\n',terminal_completion_time_s);
fprintf(fid,'terminal_mass_metric_requested=False\n');
fprintf(fid,'mass_closure_pass=%s\n',string(mass_closure_pass));
fprintf(fid,'energy_closure_pass=%s\n',string(energy_closure_pass));
fprintf(fid,'mass_relative_residual=%.15g\n',mass_relative_residual);
fprintf(fid,'energy_relative_residual=%.15g\n',energy_relative_residual);
fprintf(fid,'figure_export_pass=%s\n',string(figure_export_pass));
fprintf(fid,'candidate_manuscript_modified=False\n');
fprintf(fid,'canonical_modification_executed=False\n');
fprintf(fid,'monte_carlo_executed=False\n');
fprintf(fid,'overpressure_executed=False\n');
fprintf(fid,'population_exposure_executed=False\n');
fprintf(fid,'next_stage=Visual and numerical review of corrected nominal figures/tables; then R3C2 corrected Monte Carlo source-contract campaign\n');
fclose(fid);

fprintf('\n=== R3C1 RESULT ===\n');
fprintf('status=%s\n',status);
fprintf('authoritative_diagnostic_numeric_parity=%s\n',string(authoritative_diagnostic_numeric_parity));
fprintf('rng_sequence_preserved=%s\n',string(rng_sequence_preserved));
fprintf('frozen_r3b3_g4_metric_match=%s\n',string(frozen_r3b3_g4_metric_match));
fprintf('local_airburst_horizon_convergence_reused=%s\n',string(local_airburst_horizon_convergence_reused));
fprintf('created_child_count=%d\n',created_child_count);
fprintf('max_child_id=%d\n',max_child_id);
fprintf('fragment_count_max=%d\n',fragment_count_max);
fprintf('cap_hit_event_count=%d\n',cap_hit_event_count);
fprintf('cap_dropped_mass_kg=%.15g\n',cap_dropped_mass_kg);
fprintf('full_tree_realized=%s\n',string(full_tree_realized));
fprintf('terminal_complete=%s\n',string(terminal_complete));
fprintf('mass_closure_pass=%s\n',string(mass_closure_pass));
fprintf('energy_closure_pass=%s\n',string(energy_closure_pass));
fprintf('figure_export_pass=%s\n',string(figure_export_pass));

if ~all_gates_pass
    error('R3C1:GateFailure','One or more R3C1 scientific or artifact gates failed.');
end
end

function [pass,max_diff] = table_numeric_parity(A,B,tol)
if ~isequal(A.Properties.VariableNames,B.Properties.VariableNames) || ...
        height(A)~=height(B) || width(A)~=width(B)
    pass = false;
    max_diff = Inf;
    return;
end
AA = double(A{:,:});
BB = double(B{:,:});
nan_match = isnan(AA) & isnan(BB);
finite_match = isfinite(AA) & isfinite(BB);
invalid_mismatch = xor(isnan(AA),isnan(BB)) | xor(isinf(AA),isinf(BB));
D = abs(AA-BB);
D(nan_match) = 0;
D(~finite_match & ~nan_match) = 0;
if any(invalid_mismatch,'all')
    pass = false;
    max_diff = Inf;
    return;
end
if isempty(D)
    max_diff = 0;
else
    max_diff = max(D,[],'all');
end
pass = max_diff <= tol;
end

function write_entry_latex(path,params,S,out,break_x,break_v,break_q, ...
    nfmax,t_nfmax,h_nfmax,qmax,t_qmax,h_qmax,n_break,n_child,max_id,t_terminal)
fid = fopen(path,'w');
fprintf(fid,'%% Candidate-ready R3C1 table; not inserted automatically.\n');
fprintf(fid,'\\begin{table}[H]\n\\centering\n\\small\n');
fprintf(fid,'\\caption{Corrected nominal entry and fragmentation diagnostics for the G4 mean-area case under the 128-slot full-tree contract.}\n');
fprintf(fid,'\\label{tab:nominal_entry_fragmentation_results}\n');
fprintf(fid,'\\begin{tabular}{llll}\n\\hline\nDiagnostic & Symbol & Value & Unit \\\\\n\\hline\n');
fprintf(fid,'Initial altitude & $h_0$ & %.2f & km \\\\\n',params.h0/1000);
fprintf(fid,'Initial speed & $V_0$ & %.2f & km s$^{-1}$ \\\\\n',params.V0/1000);
fprintf(fid,'Initial mass & $m_0$ & %.6e & kg \\\\\n',S.initial_mass_kg);
fprintf(fid,'Initial kinetic energy & $E_{k,0}$ & %.3f & kt TNT \\\\\n',S.initial_ke_J/4.184e12);
fprintf(fid,'First-breakup time & $t_{\\mathrm{break}}$ & %.2f & s \\\\\n',out.first_break_time);
fprintf(fid,'First-breakup altitude & $h_{\\mathrm{break}}$ & %.3f & km \\\\\n',out.first_break_h_km);
fprintf(fid,'First-breakup downrange & $x_{\\mathrm{break}}$ & %.2f & km \\\\\n',break_x);
fprintf(fid,'Speed at first breakup & $v_{\\mathrm{break}}$ & %.3f & km s$^{-1}$ \\\\\n',break_v);
fprintf(fid,'Dynamic pressure at first breakup & $q_{\\mathrm{dyn,break}}$ & %.3f & MPa \\\\\n',break_q);
fprintf(fid,'Peak active fragment count & $N_{f,\\max}$ & %d & -- \\\\\n',nfmax);
fprintf(fid,'Time of peak fragment count & $t(N_{f,\\max})$ & %.2f & s \\\\\n',t_nfmax);
fprintf(fid,'Altitude at peak fragment count & $h(N_{f,\\max})$ & %.3f & km \\\\\n',h_nfmax);
fprintf(fid,'Maximum dynamic pressure & $q_{\\mathrm{dyn,max}}$ & %.3f & MPa \\\\\n',qmax);
fprintf(fid,'Time of maximum dynamic pressure & $t(q_{\\mathrm{dyn,max}})$ & %.2f & s \\\\\n',t_qmax);
fprintf(fid,'Altitude at maximum dynamic pressure & $h(q_{\\mathrm{dyn,max}})$ & %.3f & km \\\\\n',h_qmax);
fprintf(fid,'Breakup-event count & $N_{\\mathrm{break}}$ & %d & -- \\\\\n',n_break);
fprintf(fid,'Created-child count & $N_{\\mathrm{child}}$ & %d & -- \\\\\n',n_child);
fprintf(fid,'Maximum fragment ID & $ID_{\\max}$ & %d & -- \\\\\n',max_id);
fprintf(fid,'Terminal-completion time & $t_{\\mathrm{terminal}}$ & %.2f & s \\\\\n',t_terminal);
fprintf(fid,'\\hline\n\\end{tabular}\n\\end{table}\n');
fclose(fid);
end

function write_energy_latex(path,out,peak_bin,peak_dedh,burst_raw,bin_integral,peak_share,eta)
fid = fopen(path,'w');
fprintf(fid,'%% Candidate-ready R3C1 table; not inserted automatically.\n');
fprintf(fid,'\\begin{table}[H]\n\\centering\n\\small\n');
fprintf(fid,'\\caption{Corrected nominal energy-deposition and airburst diagnostics for the G4 mean-area case.}\n');
fprintf(fid,'\\label{tab:nominal_energy_deposition_results}\n');
fprintf(fid,'\\begin{tabular}{llll}\n\\hline\nDiagnostic & Symbol & Value & Unit \\\\\n\\hline\n');
fprintf(fid,'Main airburst altitude & $h_{\\mathrm{burst}}$ & %.2f & km \\\\\n',out.main_airburst_h_km);
fprintf(fid,'Airburst window & $\\mathcal{W}_{\\mathrm{burst}}$ & %.2f--%.2f & km \\\\\n',out.airburst_window_h_min_km,out.airburst_window_h_max_km);
fprintf(fid,'Main airburst midpoint time & $t_{\\mathrm{burst}}$ & %.2f & s \\\\\n',out.main_airburst_t);
fprintf(fid,'Peak-bin deposited energy & $\\Delta E_{\\mathrm{peak-bin}}$ & %.3f & kt TNT \\\\\n',peak_bin);
fprintf(fid,'Peak binned deposition density & $(dE_{\\mathrm{dep}}/dh)_{\\max}$ & %.3f & kt km$^{-1}$ \\\\\n',peak_dedh);
fprintf(fid,'Burst-window raw energy & $E_{\\mathrm{burst,raw}}$ & %.3f & kt TNT \\\\\n',burst_raw);
fprintf(fid,'Direct kinetic-energy loss & $E_{\\mathrm{lost}}$ & %.3f & kt TNT \\\\\n',out.E_total_lost_kton);
fprintf(fid,'Binned deposited-energy integral & $\\sum_b \\Delta E_b$ & %.3f & kt TNT \\\\\n',bin_integral);
fprintf(fid,'Peak-bin share & $f_{\\mathrm{peak-bin}}$ & %.2f & \\%% \\\\\n',peak_share);
fprintf(fid,'Coupling factor & $\\eta_c$ & %.2f & -- \\\\\n',eta);
fprintf(fid,'Nominal coupled burst energy & $E_{c,\\mathrm{nom}}$ & %.3f & kt TNT \\\\\n',out.E_burst_coupled_kton);
fprintf(fid,'Released energy & $E_{\\mathrm{release}}$ & %.3f & kt TNT \\\\\n',out.E_release_kton);
fprintf(fid,'\\hline\n\\end{tabular}\n\\end{table}\n');
fclose(fid);
end

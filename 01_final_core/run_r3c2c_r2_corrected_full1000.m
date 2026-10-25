function run_r3c2c_r2_corrected_full1000( ...
    artifact_root,design_file,pilot_results_file,horizon_gate_file)
% CHELYABINSK R3C2C-R2 EXPLICIT CSV-DELIMITER REPAIR + FULL-1000 MONTE CARLO
%
% Frozen 1000-member input design. Base horizon H360. A case is extended
% to H720 and H1440 only when it remains terminal-incomplete and exceeds
% the frozen residual thresholds, or when its outcome semantics require
% adjudication. Publishable metrics use H360 for non-extended cases and
% H1440 for scientifically accepted extended cases.
%
% No overpressure, population exposure, manuscript modification, or
% canonical modification is performed.

files = string({design_file,pilot_results_file,horizon_gate_file});
if nargin < 4 || strlength(string(artifact_root))==0 || any(~isfile(files))
    error('R3C2C requires artifact_root and three valid frozen inputs.');
end

artifact_root = char(string(artifact_root));
design_file = char(string(design_file));
pilot_results_file = char(string(pilot_results_file));
horizon_gate_file = char(string(horizon_gate_file));

core_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(artifact_root,'02_full1000_outputs');
figure_dir = fullfile(output_dir,'diagnostic_figures');
flagged_dir = fullfile(output_dir,'flagged_run_ledgers');
metric_audit_dir = fullfile(output_dir,'adaptive_metric_audits');
report_dir = fullfile(artifact_root,'03_reports');
for p = {output_dir,figure_dir,flagged_dir,metric_audit_dir,report_dir}
    if ~exist(p{1},'dir'), mkdir(p{1}); end
end
addpath(core_dir);
format long g;

fprintf('=== CHELYABINSK R3C2C-R2 EXPLICIT CSV-DELIMITER REPAIR + FULL-1000 MONTE CARLO ===\n');
fprintf('prior_R3C2C_scientific_propagation_executed=False\n');
fprintf('prior_R3C2C_R1_scientific_propagation_executed=False\n');
fprintf('explicit_csv_delimiter_repair_applied=True\n');
fprintf('prior_failed_scientific_propagation_executed=False\n');
fprintf('campaign_driver_import_schema_repair_applied=True\n');
fprintf('normalized_csv_header_binding_requested=True\n');
fprintf('scientific_propagation_requested=True\n');
fprintf('full_1000_member_campaign_requested=True\n');
fprintf('base_horizon_s=360\n');
fprintf('adaptive_horizon_720_s=True\n');
fprintf('adaptive_horizon_1440_s=True\n');
fprintf('max_total_fragments=128\n');
fprintf('dt_s=0.02\n');
fprintf('area_multiplier=1.17161624727046\n');
fprintf('residual_mass_fraction_threshold=1e-3\n');
fprintf('residual_ke_fraction_threshold=1e-8\n');
fprintf('terminal_completion_hard_gate=False\n');
fprintf('terminal_mass_metric_publishable=False\n');
fprintf('scientific_model_source_modification_requested=False\n');
fprintf('candidate_manuscript_modification_requested=False\n');
fprintf('canonical_modification_requested=False\n');
fprintf('overpressure_requested=False\n');
fprintf('population_exposure_requested=False\n');
fprintf('resume_supported=True\n');
fprintf('artifact_root=%s\n',artifact_root);

DESIGN = sortrows(read_csv_normalized(design_file),'run_id');
PILOT = read_csv_normalized(pilot_results_file);
HGATE = read_csv_normalized(horizon_gate_file);

assert_columns(DESIGN,["run_id","fragmentation_seed","V0_m_s", ...
    "entry_angle_deg","D0_m","rho_meteor_kg_m3", ...
    "Yield_Stress_Pa","Ch","energy_efficiency"],"full design");
assert_columns(PILOT,["run_id","fragmentation_seed", ...
    "first_break_h_km","main_airburst_h_km", ...
    "E_burst_coupled_kton","q_dyn_max_MPa"],"pilot results");
assert_columns(HGATE,["gate_name","gate_pass","gate_detail"], ...
    "horizon gate adjudication");

schema_source = ["full_design";"pilot_results";"horizon_gate_adjudication"];
observed_variable_count = [width(DESIGN);width(PILOT);width(HGATE)];
observed_variable_names = [strjoin(string(DESIGN.Properties.VariableNames),"|"); ...
    strjoin(string(PILOT.Properties.VariableNames),"|"); ...
    strjoin(string(HGATE.Properties.VariableNames),"|")];
required_columns_present = [true;true;true];
writetable(table(schema_source,observed_variable_count, ...
    observed_variable_names,required_columns_present), ...
    fullfile(output_dir,'R3C2C_R1_import_schema_audit.csv'));

N = height(DESIGN);
design_contract_pass = N==1000 && ...
    isequal(DESIGN.run_id,(1:1000)') && ...
    numel(unique(DESIGN.fragmentation_seed))==1000;
gate_names = string(HGATE.gate_name);
gate_truth = truth_vector(HGATE.gate_pass);
horizon_contract_pass = any(gate_names== ...
    "H360_adequate_for_publishable_metrics" & gate_truth) && ...
    any(gate_names== ...
    "terminal_completion_required_for_local_airburst_metrics_false" & gate_truth);
if ~design_contract_pass || ~horizon_contract_pass
    error('R3C2C:InputContract','Frozen design or horizon contract failed.');
end

% Frozen convergence tolerances inherited from R3C2B-R1.
tol = table( ...
    ["first_break_h_km";"first_break_time_s";"main_airburst_h_km"; ...
     "main_airburst_time_s";"E_total_lost_kton";"E_release_kton"; ...
     "E_burst_raw_kton";"E_burst_coupled_kton"; ...
     "dE_dh_peak_kton_per_km";"airburst_window_h_min_km"; ...
     "airburst_window_h_max_km";"q_dyn_max_MPa"], ...
    [1e-9;1e-9;1e-9;5e-2;1e-5;1e-5;1e-5;1e-5;1e-5;1e-9;1e-9;1e-8], ...
    [0;0;0;0;1e-8;1e-8;1e-8;1e-8;1e-8;0;0;1e-8], ...
    'VariableNames',{'metric','absolute_tolerance','relative_tolerance'});
writetable(tol,fullfile(output_dir,'R3C2C_metric_convergence_tolerances.csv'));

% Initialize or resume the 1000-row checkpoint.
checkpoint_file = fullfile(output_dir,'R3C2C_full1000_results_checkpoint.csv');
if isfile(checkpoint_file)
    RESULTS = read_csv_normalized(checkpoint_file);
    if height(RESULTS)~=N || ~isequal(RESULTS.run_id,DESIGN.run_id)
        error('R3C2C:ResumeContract','Checkpoint does not match the frozen design.');
    end
    fprintf('resume_checkpoint_found=True\n');
    fprintf('resume_completed_count=%d\n',sum(logical(RESULTS.run_completed)));
else
    RESULTS = struct2table(repmat(empty_result(),N,1));
    RESULTS.run_id = DESIGN.run_id;
    RESULTS.input_design_seed = DESIGN.input_design_seed;
    RESULTS.fragmentation_seed = DESIGN.fragmentation_seed;
    RESULTS.V0_m_s = DESIGN.V0_m_s;
    RESULTS.entry_angle_deg = DESIGN.entry_angle_deg;
    RESULTS.D0_m = DESIGN.D0_m;
    RESULTS.rho_meteor_kg_m3 = DESIGN.rho_meteor_kg_m3;
    RESULTS.Yield_Stress_Pa = DESIGN.Yield_Stress_Pa;
    RESULTS.Ch = DESIGN.Ch;
    RESULTS.energy_efficiency = DESIGN.energy_efficiency;
    RESULTS.area_multiplier(:) = 1.17161624727046;
    RESULTS.max_total_fragments(:) = 128;
    RESULTS.dt_s(:) = 0.02;
    RESULTS.base_horizon_s(:) = 360;
    writetable(RESULTS,checkpoint_file);
    fprintf('resume_checkpoint_found=False\n');
end

metric_compare_cells = cell(N,2);
campaign_tic = tic;

for i = 1:N
    if logical(RESULTS.run_completed(i))
        continue;
    end

    row = DESIGN(i,:);
    fprintf('\n--- FULL1000 %d/%d | run_id=%d | seed=%d ---\n', ...
        i,N,row.run_id,row.fragmentation_seed);

    r = table2struct(RESULTS(i,:));
    run_tic = tic;

    try
        p360 = build_params(row,360);
        out360 = meteor_single_run_r3b3_diagnostic(p360);
        r.h360_runtime_s = toc(run_tic);
        r.h360_run_completed = true;

        s360 = summarize_output(out360);
        r = assign_h360(r,out360,s360);
        r.h360_semantic_output_pass = semantic_valid(out360);
        r.h360_residual_threshold_pass = ...
            s360.residual_mass_fraction<=1e-3 && ...
            s360.residual_ke_fraction<=1e-8;
        r.h360_hard_contract_pass = s360.cap_contract_pass && ...
            s360.mass_closure_pass && s360.energy_closure_pass && ...
            s360.invalid_deactivation_zero && s360.structural_bound_pass;

        r.extension_required = r.h360_hard_contract_pass && ...
            ((~s360.terminal_complete && ~r.h360_residual_threshold_pass) || ...
             ~r.h360_semantic_output_pass);

        selected_out = out360;
        selected_summary = s360;
        r.publish_metric_source_horizon_s = 360;
        r.extension_convergence_pass = true;

        if r.extension_required
            p720 = build_params(row,720);
            t720 = tic;
            out720 = meteor_single_run_r3b3_diagnostic(p720);
            r.h720_runtime_s = toc(t720);
            r.h720_run_completed = true;
            s720 = summarize_output(out720);
            r = assign_h720(r,s720);

            p1440 = build_params(row,1440);
            t1440 = tic;
            out1440 = meteor_single_run_r3b3_diagnostic(p1440);
            r.h1440_runtime_s = toc(t1440);
            r.h1440_run_completed = true;
            s1440 = summarize_output(out1440);
            r = assign_h1440(r,s1440);

            [cmp720_1440,local_pass,global_pass] = ...
                compare_outputs(out720,out1440,"H720_to_H1440",tol,row);
            metric_compare_cells{i,1} = cmp720_1440;
            r.h720_to_h1440_local_convergence_pass = local_pass;
            r.h720_to_h1440_global_convergence_pass = global_pass;
            r.extension_convergence_pass = local_pass && global_pass;

            [cmp360_1440,~,~] = ...
                compare_outputs(out360,out1440,"H360_to_H1440",tol,row);
            metric_compare_cells{i,2} = cmp360_1440;
            run_metric_audit = [cmp720_1440;cmp360_1440];
            writetable(run_metric_audit,fullfile(metric_audit_dir, ...
                sprintf('run_%04d_adaptive_metric_audit.csv',row.run_id)));

            selected_out = out1440;
            selected_summary = s1440;
            r.publish_metric_source_horizon_s = 1440;
        end

        r = assign_publish_metrics(r,selected_out,selected_summary);
        r.publish_semantic_output_pass = semantic_valid(selected_out);
        r.publish_residual_threshold_pass = ...
            selected_summary.residual_mass_fraction<=1e-3 && ...
            selected_summary.residual_ke_fraction<=1e-8;
        r.publish_hard_contract_pass = ...
            selected_summary.cap_contract_pass && ...
            selected_summary.mass_closure_pass && ...
            selected_summary.energy_closure_pass && ...
            selected_summary.invalid_deactivation_zero && ...
            selected_summary.structural_bound_pass;

        r.run_gate_pass = r.h360_run_completed && ...
            r.h360_hard_contract_pass && ...
            (~r.extension_required || ...
                (r.h720_run_completed && r.h1440_run_completed && ...
                 r.extension_convergence_pass)) && ...
            r.publish_semantic_output_pass && ...
            r.publish_residual_threshold_pass && ...
            r.publish_hard_contract_pass;

        r.run_completed = true;
        r.runtime_s = toc(run_tic);
        r.error_identifier = "";
        r.error_message = "";

        if ~r.run_gate_pass || r.extension_required
            write_flagged_ledgers(flagged_dir,row,out360,selected_out,r);
        end

        fprintf('outcome=%s\n',r.publish_outcome);
        fprintf('h360_terminal_complete=%s\n',string(r.h360_terminal_complete));
        fprintf('h360_residual_mass_fraction=%.15g\n',r.h360_residual_mass_fraction);
        fprintf('h360_residual_ke_fraction=%.15g\n',r.h360_residual_ke_fraction);
        fprintf('extension_required=%s\n',string(r.extension_required));
        fprintf('publish_metric_source_horizon_s=%d\n',r.publish_metric_source_horizon_s);
        fprintf('cap_hit_event_count=%d\n',r.publish_cap_hit_event_count);
        fprintf('mass_closure_pass=%s\n',string(r.publish_mass_closure_pass));
        fprintf('energy_closure_pass=%s\n',string(r.publish_energy_closure_pass));
        fprintf('run_gate_pass=%s\n',string(r.run_gate_pass));

        clear out360 out720 out1440 selected_out
    catch ME
        r.runtime_s = toc(run_tic);
        r.run_completed = true;
        r.run_gate_pass = false;
        r.error_identifier = string(ME.identifier);
        r.error_message = string(ME.message);
        fprintf('run_exception=True\n');
        fprintf('error_identifier=%s\n',r.error_identifier);
        fprintf('error_message=%s\n',r.error_message);
    end

    RESULTS(i,:) = struct2table(r);
    if mod(i,5)==0 || i==N || ~r.run_gate_pass || r.extension_required
        writetable(RESULTS,checkpoint_file);
    end
end

campaign_runtime_s = toc(campaign_tic);
writetable(RESULTS,checkpoint_file);
writetable(RESULTS,fullfile(output_dir,'R3C2C_full1000_results.csv'));

metric_audit = collect_metric_audits(metric_audit_dir);
writetable(metric_audit,fullfile(output_dir, ...
    'R3C2C_adaptive_horizon_metric_audit.csv'));

% Pilot reproduction at the H360 contract.
[pilot_audit,pilot_reproduction_pass] = ...
    pilot_reproduction(RESULTS,PILOT,tol);
writetable(pilot_audit,fullfile(output_dir,'R3C2C_pilot64_reproduction_audit.csv'));

% Full-design authoritative/diagnostic parity sentinels.
[sentinel_table,parity_audit] = ...
    full_campaign_parity_sentinels(DESIGN,RESULTS);
writetable(sentinel_table,fullfile(output_dir,'R3C2C_authoritative_parity_sentinels.csv'));
writetable(parity_audit,fullfile(output_dir,'R3C2C_authoritative_parity_audit.csv'));
sentinel_parity_pass = height(parity_audit)==height(sentinel_table) && ...
    all(parity_audit.parity_pass);

% Campaign gates.
all_1000_runs_completed = height(RESULTS)==1000 && all(logical(RESULTS.run_completed));
all_run_ids_exactly_once = isequal(RESULTS.run_id,(1:1000)');
all_run_gates_pass = all(logical(RESULTS.run_gate_pass));
zero_cap_hits_all_runs = all(RESULTS.publish_cap_hit_event_count==0);
zero_cap_dropped_mass_all_runs = all(abs(RESULTS.publish_cap_dropped_mass_kg)<=1e-9);
zero_insertion_dropped_mass_all_runs = ...
    all(abs(RESULTS.publish_insertion_dropped_mass_kg)<=1e-9);
mass_closure_pass_all_runs = all(logical(RESULTS.publish_mass_closure_pass));
energy_closure_pass_all_runs = all(logical(RESULTS.publish_energy_closure_pass));
invalid_deactivation_zero_all_runs = ...
    all(RESULTS.publish_invalid_deactivation_count==0);
structural_bounds_pass_all_runs = all(logical(RESULTS.publish_structural_bound_pass));
semantic_output_pass_all_runs = all(logical(RESULTS.publish_semantic_output_pass));
residual_threshold_pass_all_runs = all(logical(RESULTS.publish_residual_threshold_pass));
all_extensions_converged = all(~logical(RESULTS.extension_required) | ...
    logical(RESULTS.extension_convergence_pass));
terminal_mass_metric_publishable = false;
terminal_completion_hard_gate = false;

gate_name = [
    "design_contract_pass";
    "horizon_contract_pass";
    "all_1000_runs_completed";
    "all_run_ids_exactly_once";
    "all_run_gates_pass";
    "zero_cap_hits_all_runs";
    "zero_cap_dropped_mass_all_runs";
    "zero_insertion_dropped_mass_all_runs";
    "mass_closure_pass_all_runs";
    "energy_closure_pass_all_runs";
    "invalid_deactivation_zero_all_runs";
    "structural_bounds_pass_all_runs";
    "semantic_output_pass_all_runs";
    "residual_threshold_pass_all_runs";
    "all_adaptive_extensions_converged";
    "pilot64_reproduction_pass";
    "authoritative_diagnostic_sentinel_parity_pass";
    "terminal_completion_hard_gate_false";
    "terminal_mass_metric_publishable_false"
];
gate_pass = [
    design_contract_pass;
    horizon_contract_pass;
    all_1000_runs_completed;
    all_run_ids_exactly_once;
    all_run_gates_pass;
    zero_cap_hits_all_runs;
    zero_cap_dropped_mass_all_runs;
    zero_insertion_dropped_mass_all_runs;
    mass_closure_pass_all_runs;
    energy_closure_pass_all_runs;
    invalid_deactivation_zero_all_runs;
    structural_bounds_pass_all_runs;
    semantic_output_pass_all_runs;
    residual_threshold_pass_all_runs;
    all_extensions_converged;
    pilot_reproduction_pass;
    sentinel_parity_pass;
    ~terminal_completion_hard_gate;
    ~terminal_mass_metric_publishable
];
gate_detail = [
    sprintf('design_rows=%d',height(DESIGN));
    "R3C2B-R1 H360 adequacy gate is frozen PASS";
    sprintf('completed=%d/1000',sum(logical(RESULTS.run_completed)));
    sprintf('unique_run_ids=%d',numel(unique(RESULTS.run_id)));
    sprintf('run_gate_pass_count=%d',sum(logical(RESULTS.run_gate_pass)));
    sprintf('maximum_cap_hits=%g',max(RESULTS.publish_cap_hit_event_count));
    sprintf('maximum_abs_cap_drop_kg=%.15g',max(abs(RESULTS.publish_cap_dropped_mass_kg)));
    sprintf('maximum_abs_insertion_drop_kg=%.15g',max(abs(RESULTS.publish_insertion_dropped_mass_kg)));
    sprintf('maximum_mass_relative_residual=%.15g',max(RESULTS.publish_mass_relative_residual));
    sprintf('maximum_energy_relative_residual=%.15g',max(RESULTS.publish_energy_relative_residual));
    sprintf('maximum_invalid_deactivation_count=%g',max(RESULTS.publish_invalid_deactivation_count));
    sprintf('structural_pass_count=%d',sum(logical(RESULTS.publish_structural_bound_pass)));
    sprintf('semantic_pass_count=%d',sum(logical(RESULTS.publish_semantic_output_pass)));
    sprintf('residual_pass_count=%d',sum(logical(RESULTS.publish_residual_threshold_pass)));
    sprintf('extension_count=%d; convergence_pass_count=%d', ...
        sum(logical(RESULTS.extension_required)), ...
        sum(~logical(RESULTS.extension_required) | logical(RESULTS.extension_convergence_pass)));
    sprintf('pilot_metric_pass_rows=%d/%d',sum(pilot_audit.metric_pass),height(pilot_audit));
    sprintf('sentinel_pass_count=%d/%d',sum(parity_audit.parity_pass),height(parity_audit));
    "terminal completion is retained as a diagnostic";
    "fixed-horizon terminal active mass remains retired"
];
writetable(table(gate_name,gate_pass,gate_detail), ...
    fullfile(output_dir,'R3C2C_full1000_gate_adjudication.csv'));

% Publishable result subset and descriptive statistics.
publish_fields = [
    "run_id","fragmentation_seed","publish_metric_source_horizon_s", ...
    "publish_outcome","publish_first_break_h_km","publish_first_break_time_s", ...
    "publish_main_airburst_h_km","publish_main_airburst_time_s", ...
    "publish_E_total_lost_kton","publish_E_release_kton", ...
    "publish_E_burst_raw_kton","publish_E_burst_coupled_kton", ...
    "publish_dE_dh_peak_kton_per_km","publish_airburst_window_h_min_km", ...
    "publish_airburst_window_h_max_km","publish_q_dyn_max_MPa", ...
    "publish_fragment_count_max","publish_breakup_event_count", ...
    "publish_created_child_count","publish_terminal_complete"
];
PUB = RESULTS(:,cellstr(publish_fields));
writetable(PUB,fullfile(output_dir,'R3C2C_full1000_publishable_metric_matrix.csv'));

stats_fields = [
    "publish_first_break_h_km","publish_main_airburst_h_km", ...
    "publish_main_airburst_time_s","publish_E_total_lost_kton", ...
    "publish_E_release_kton","publish_E_burst_raw_kton", ...
    "publish_E_burst_coupled_kton","publish_dE_dh_peak_kton_per_km", ...
    "publish_q_dyn_max_MPa","publish_fragment_count_max", ...
    "publish_breakup_event_count","publish_created_child_count"
];
stats_cells = cell(numel(stats_fields),1);
for j = 1:numel(stats_fields)
    fld = stats_fields(j);
    x = RESULTS.(fld);
    x = x(isfinite(x));
    if isempty(x), v=nan(1,9);
    else
        v=[min(x),pct(x,1),pct(x,5),median(x),mean(x), ...
           pct(x,95),pct(x,99),max(x),std(x)];
    end
    stats_cells{j} = table(fld,numel(x),v(1),v(2),v(3),v(4),v(5), ...
        v(6),v(7),v(8),v(9), ...
        'VariableNames',{'metric','finite_count','minimum','p01','p05', ...
        'median','mean','p95','p99','maximum','standard_deviation'});
end
STATS = vertcat(stats_cells{:});
writetable(STATS,fullfile(output_dir,'R3C2C_full1000_summary_statistics.csv'));

[outcome_name,~,og] = unique(RESULTS.publish_outcome);
outcome_count = accumarray(og,1);
outcome_fraction = outcome_count/height(RESULTS);
writetable(table(outcome_name,outcome_count,outcome_fraction), ...
    fullfile(output_dir,'R3C2C_full1000_outcome_counts.csv'));

% Rank-correlation diagnostic between frozen inputs and principal outputs.
input_fields = ["V0_m_s","entry_angle_deg","D0_m","rho_meteor_kg_m3", ...
    "Yield_Stress_Pa","Ch","energy_efficiency"];
output_fields = [
    "publish_first_break_h_km","publish_main_airburst_h_km", ...
    "publish_E_burst_coupled_kton","publish_dE_dh_peak_kton_per_km", ...
    "publish_q_dyn_max_MPa"
];
corr_rows = cell(numel(input_fields)*numel(output_fields),1);
k=0;
for i = 1:numel(input_fields)
    for j = 1:numel(output_fields)
        k=k+1;
        x=RESULTS.(input_fields(i));
        y=RESULTS.(output_fields(j));
        mask=isfinite(x)&isfinite(y);
        rho=spearman_no_toolbox(x(mask),y(mask));
        corr_rows{k}=table(input_fields(i),output_fields(j),sum(mask),rho, ...
            'VariableNames',{'input_parameter','output_metric','pair_count','spearman_rho'});
    end
end
CORR=vertcat(corr_rows{:});
writetable(CORR,fullfile(output_dir,'R3C2C_full1000_spearman_correlations.csv'));

% Prefix-size stability is diagnostic only.
prefix_sizes=[100;250;500;750;1000];
stability_cells=cell(numel(prefix_sizes)*4,1);
metrics_stability=["publish_main_airburst_h_km","publish_E_burst_coupled_kton", ...
    "publish_dE_dh_peak_kton_per_km","publish_q_dyn_max_MPa"];
k=0;
for a=1:numel(prefix_sizes)
    nn=prefix_sizes(a);
    for b=1:numel(metrics_stability)
        k=k+1;
        x=RESULTS.(metrics_stability(b));
        x=x(1:nn); x=x(isfinite(x));
        stability_cells{k}=table(nn,metrics_stability(b),numel(x), ...
            pct(x,5),median(x),pct(x,95), ...
            'VariableNames',{'prefix_size','metric','finite_count','p05','median','p95'});
    end
end
writetable(vertcat(stability_cells{:}), ...
    fullfile(output_dir,'R3C2C_ensemble_prefix_stability.csv'));

FAILED=RESULTS(~logical(RESULTS.run_gate_pass),:);
EXTENDED=RESULTS(logical(RESULTS.extension_required),:);
writetable(FAILED,fullfile(output_dir,'R3C2C_failed_or_flagged_runs.csv'));
writetable(EXTENDED,fullfile(output_dir,'R3C2C_adaptively_extended_runs.csv'));

% Diagnostic figures. Manuscript adoption is deferred.
make_distribution_figure(RESULTS,figure_dir);
make_correlation_figure(CORR,input_fields,output_fields,figure_dir);
make_stability_figure(vertcat(stability_cells{:}),figure_dir);

all_gates_pass = all(gate_pass);
if all_gates_pass
    status = "PASS_R3C2C_R2_CORRECTED_FULL1000_MONTE_CARLO";
else
    status = "FAIL_R3C2C_R2_CORRECTED_FULL1000_MONTE_CARLO";
end

fid=fopen(fullfile(report_dir,'R3C2C_status.txt'),'w');
fprintf(fid,'=== CHELYABINSK R3C2C RESULT ===\n');
fprintf(fid,'status=%s\n',status);
fprintf(fid,'scientific_propagation_executed=True\n');
fprintf(fid,'campaign_driver_import_schema_repair_applied=True\n');
fprintf(fid,'explicit_csv_delimiter_repair_applied=True\n');
fprintf(fid,'full_1000_member_campaign_executed=True\n');
fprintf(fid,'completed_run_count=%d\n',sum(logical(RESULTS.run_completed)));
fprintf(fid,'run_gate_pass_count=%d\n',sum(logical(RESULTS.run_gate_pass)));
fprintf(fid,'adaptive_extension_count=%d\n',sum(logical(RESULTS.extension_required)));
fprintf(fid,'H360_terminal_complete_count=%d\n',sum(logical(RESULTS.h360_terminal_complete)));
fprintf(fid,'publish_terminal_complete_count=%d\n',sum(logical(RESULTS.publish_terminal_complete)));
fprintf(fid,'zero_cap_hits_all_runs=%s\n',string(zero_cap_hits_all_runs));
fprintf(fid,'zero_cap_dropped_mass_all_runs=%s\n',string(zero_cap_dropped_mass_all_runs));
fprintf(fid,'zero_insertion_dropped_mass_all_runs=%s\n',string(zero_insertion_dropped_mass_all_runs));
fprintf(fid,'mass_closure_pass_all_runs=%s\n',string(mass_closure_pass_all_runs));
fprintf(fid,'energy_closure_pass_all_runs=%s\n',string(energy_closure_pass_all_runs));
fprintf(fid,'residual_threshold_pass_all_runs=%s\n',string(residual_threshold_pass_all_runs));
fprintf(fid,'all_adaptive_extensions_converged=%s\n',string(all_extensions_converged));
fprintf(fid,'pilot64_reproduction_pass=%s\n',string(pilot_reproduction_pass));
fprintf(fid,'sentinel_count=%d\n',height(parity_audit));
fprintf(fid,'sentinel_parity_pass=%s\n',string(sentinel_parity_pass));
fprintf(fid,'maximum_publish_mass_relative_residual=%.15g\n',max(RESULTS.publish_mass_relative_residual));
fprintf(fid,'maximum_publish_energy_relative_residual=%.15g\n',max(RESULTS.publish_energy_relative_residual));
fprintf(fid,'maximum_publish_residual_mass_fraction=%.15g\n',max(RESULTS.publish_residual_mass_fraction));
fprintf(fid,'maximum_publish_residual_ke_fraction=%.15g\n',max(RESULTS.publish_residual_ke_fraction));
fprintf(fid,'terminal_completion_hard_gate=False\n');
fprintf(fid,'terminal_mass_metric_publishable=False\n');
fprintf(fid,'full1000_statistics_publishable=%s\n',string(all_gates_pass));
fprintf(fid,'overpressure_executed=False\n');
fprintf(fid,'population_exposure_executed=False\n');
fprintf(fid,'campaign_runtime_s=%.15g\n',campaign_runtime_s);
fprintf(fid,'candidate_manuscript_modified=False\n');
fprintf(fid,'canonical_modification_executed=False\n');
fprintf(fid,'next_stage=Upload R3C2C review bundle. After scientific and visual PASS, construct corrected overpressure-source handoff.\n');
fclose(fid);

fprintf('\n=== R3C2C RESULT ===\n');
fprintf('status=%s\n',status);
fprintf('completed_run_count=%d\n',sum(logical(RESULTS.run_completed)));
fprintf('run_gate_pass_count=%d\n',sum(logical(RESULTS.run_gate_pass)));
fprintf('adaptive_extension_count=%d\n',sum(logical(RESULTS.extension_required)));
fprintf('H360_terminal_complete_count=%d\n',sum(logical(RESULTS.h360_terminal_complete)));
fprintf('publish_terminal_complete_count=%d\n',sum(logical(RESULTS.publish_terminal_complete)));
fprintf('zero_cap_hits_all_runs=%s\n',string(zero_cap_hits_all_runs));
fprintf('mass_closure_pass_all_runs=%s\n',string(mass_closure_pass_all_runs));
fprintf('energy_closure_pass_all_runs=%s\n',string(energy_closure_pass_all_runs));
fprintf('residual_threshold_pass_all_runs=%s\n',string(residual_threshold_pass_all_runs));
fprintf('all_adaptive_extensions_converged=%s\n',string(all_extensions_converged));
fprintf('pilot64_reproduction_pass=%s\n',string(pilot_reproduction_pass));
fprintf('sentinel_parity_pass=%s\n',string(sentinel_parity_pass));
fprintf('full1000_statistics_publishable=%s\n',string(all_gates_pass));
fprintf('terminal_completion_hard_gate=False\n');
fprintf('terminal_mass_metric_publishable=False\n');
fprintf('campaign_runtime_s=%.15g\n',campaign_runtime_s);

if ~all_gates_pass
    error('R3C2C:GateFailure', ...
        'Full1000 scientific gates failed. Freeze outputs; do not start overpressure.');
end
end


function T=read_csv_normalized(path)
% Read a frozen comma-separated file without delimiter auto-detection.
% R3C2C-R1 showed that readtable could infer underscore as a delimiter from
% underscore-rich scientific identifiers. The comma contract is therefore
% explicit and is verified against the raw header before table creation.
raw_header = read_raw_header(path);
raw_tokens = strip(split(raw_header,","));
if numel(raw_tokens)<2
    error('R3C2C:ImportSchema', ...
        'Raw CSV header does not contain a comma delimiter: %s', ...
        char(raw_header));
end
T=readtable(path,'FileType','text','Delimiter',',', ...
    'ReadVariableNames',true,'TextType','string', ...
    'VariableNamingRule','preserve','Encoding','UTF-8');
names=string(T.Properties.VariableNames);
names=replace(names,char(65279),"");
names=strip(names);
names=matlab.lang.makeValidName(names,'ReplacementStyle','underscore');
names=matlab.lang.makeUniqueStrings(names);
T.Properties.VariableNames=cellstr(names);
normalized_header = string(T.Properties.VariableNames);
raw_normalized = replace(raw_tokens,char(65279),"");
raw_normalized = matlab.lang.makeValidName(raw_normalized, ...
    'ReplacementStyle','underscore');
raw_normalized = matlab.lang.makeUniqueStrings(raw_normalized);
if width(T)~=numel(raw_normalized) || ...
        ~isequal(normalized_header(:),string(raw_normalized(:)))
    error('R3C2C:ImportSchema', ...
        'Explicit comma import disagrees with raw header. Raw: %s; Imported: %s', ...
        char(strjoin(raw_tokens,'|')),char(strjoin(normalized_header,'|')));
end
end

function header=read_raw_header(path)
fid=fopen(path,'r','n','UTF-8');
if fid<0
    error('R3C2C:ImportSchema','Unable to open CSV: %s',path);
end
cleanup=onCleanup(@() fclose(fid));
line=fgetl(fid);
if ~ischar(line)
    error('R3C2C:ImportSchema','CSV has no readable header: %s',path);
end
header=string(line);
header=replace(header,char(65279),"");
end

function assert_columns(T,required,label)
observed=string(T.Properties.VariableNames);
missing=required(~ismember(required,observed));
if ~isempty(missing)
    error('R3C2C:ImportSchema', ...
        'Missing required columns in %s: %s. Observed: %s', ...
        char(string(label)),char(strjoin(missing,',')), ...
        char(strjoin(observed,',')));
end
end

function v=truth_vector(x)
if islogical(x)
    v=x;
elseif isnumeric(x)
    v=x~=0;
elseif isstring(x)
    s=lower(strip(x));
    v=s=="1"|s=="true"|s=="yes"|s=="pass";
elseif iscell(x)
    s=lower(strip(string(x)));
    v=s=="1"|s=="true"|s=="yes"|s=="pass";
elseif iscategorical(x)
    s=lower(strip(string(x)));
    v=s=="1"|s=="true"|s=="yes"|s=="pass";
else
    error('R3C2C:ImportSchema','Unsupported gate-pass data type: %s',class(x));
end
v=v(:);
end

function T=collect_metric_audits(folder)
files=dir(fullfile(folder,'run_*_adaptive_metric_audit.csv'));
if isempty(files)
    T=table('Size',[0 9], ...
        'VariableTypes',{'double','double','string','string', ...
        'double','double','double','double','logical'}, ...
        'VariableNames',{'run_id','fragmentation_seed','comparison','metric', ...
        'earlier_value','later_value','absolute_difference', ...
        'allowed_difference','metric_pass'});
    return;
end
cells=cell(numel(files),1);
for i=1:numel(files)
    cells{i}=read_csv_normalized(fullfile(files(i).folder,files(i).name));
end
T=vertcat(cells{:});
T=sortrows(T,{'run_id','comparison','metric'});
end

function params=build_params(row,tmax)
params.rho_meteor=row.rho_meteor_kg_m3;
params.V0=row.V0_m_s;
params.entry_angle_deg=row.entry_angle_deg;
params.D0=row.D0_m;
params.h0=97100;
params.Yield_Stress_0=row.Yield_Stress_Pa;
params.min_strength_fraction=0.45;
params.weakening_mass_coeff=0.18;
params.Cd=1.0;
params.Ch=row.Ch;
params.Qabl=8e6;
params.max_generations=3;
params.max_total_fragments=128;
params.base_fragment_count=5;
params.dispersion_speed=250;
params.cloud_area_gain=2.5;
params.area_multiplier=1.17161624727046;
params.dt=0.02;
params.tmax=tmax;
params.energy_efficiency=row.energy_efficiency;
params.energy_bin_width_m=500;
params.seed=row.fragmentation_seed;
end

function s=summarize_output(out)
S=out.r3b.closure_summary;
B=out.r3b.breakup_event_ledger;
C=out.r3b.child_creation_ledger;
D=out.r3b.deactivation_ledger;
TF=out.r3b.terminal_fragment_ledger;
TL=out.r3b.mass_energy_time_ledger;
s.terminal_active_count=height(TF);
s.terminal_complete=isempty(TF)&&abs(S.final_active_mass_kg)<=1e-6;
s.final_active_mass_kg=S.final_active_mass_kg;
s.final_active_ke_J=S.final_active_ke_J;
s.residual_mass_fraction=abs(S.final_active_mass_kg)/max(S.initial_mass_kg,eps);
s.residual_ke_fraction=abs(S.final_active_ke_J)/max(S.initial_ke_J,eps);
if isempty(TL),s.terminal_time_s=NaN;else,s.terminal_time_s=TL.time_s(end);end
s.breakup_event_count=height(B);
s.created_child_count=height(C);
if isempty(C),s.max_child_id=1;s.all_children_inserted=true;
else,s.max_child_id=max(C.child_id);s.all_children_inserted=all(C.inserted_into_active_pool);end
s.cap_hit_event_count=sum(B.cap_hit);
s.cap_dropped_mass_kg=sum(B.cap_dropped_mass_kg);
s.insertion_dropped_mass_kg=S.cumulative_insertion_dropped_mass_kg;
s.cap_contract_pass=s.cap_hit_event_count==0&&abs(s.cap_dropped_mass_kg)<=1e-9&& ...
    abs(s.insertion_dropped_mass_kg)<=1e-9;
s.mass_relative_residual=abs(S.final_mass_closure_residual_kg)/max(S.initial_mass_kg,eps);
s.energy_relative_residual=abs(S.final_energy_closure_residual_J)/max(S.initial_ke_J,eps);
s.mass_closure_pass=s.mass_relative_residual<=1e-9;
s.energy_closure_pass=s.energy_relative_residual<=1e-9;
if isempty(D),s.invalid_deactivation_count=0;
else
    valid=D.reason=="post_step_ground_crossing"|D.reason=="post_step_mass_depletion";
    s.invalid_deactivation_count=sum(~valid);
end
s.invalid_deactivation_zero=s.invalid_deactivation_count==0;
s.structural_bound_pass=s.breakup_event_count<=26&&s.created_child_count<=85&& ...
    s.max_child_id<=86&&s.all_children_inserted;
end

function r=assign_h360(r,out,s)
r.h360_outcome=string(out.outcome);
r.h360_first_break_h_km=out.first_break_h_km;
r.h360_first_break_time_s=out.first_break_time;
r.h360_main_airburst_h_km=out.main_airburst_h_km;
r.h360_main_airburst_time_s=out.main_airburst_t;
r.h360_E_total_lost_kton=out.E_total_lost_kton;
r.h360_E_release_kton=out.E_release_kton;
r.h360_E_burst_raw_kton=out.E_burst_raw_kton;
r.h360_E_burst_coupled_kton=out.E_burst_coupled_kton;
r.h360_dE_dh_peak_kton_per_km=out.dE_dh_peak_kton_per_km;
r.h360_airburst_window_h_min_km=out.airburst_window_h_min_km;
r.h360_airburst_window_h_max_km=out.airburst_window_h_max_km;
r.h360_q_dyn_max_MPa=max(out.timeseries.q_dyn_max_Pa)/1e6;
r.h360_fragment_count_max=max(out.r3b.mass_energy_time_ledger.active_count);
r.h360_breakup_event_count=s.breakup_event_count;
r.h360_created_child_count=s.created_child_count;
r.h360_max_child_id=s.max_child_id;
r.h360_cap_hit_event_count=s.cap_hit_event_count;
r.h360_cap_dropped_mass_kg=s.cap_dropped_mass_kg;
r.h360_insertion_dropped_mass_kg=s.insertion_dropped_mass_kg;
r.h360_terminal_active_count=s.terminal_active_count;
r.h360_terminal_complete=s.terminal_complete;
r.h360_terminal_time_s=s.terminal_time_s;
r.h360_residual_mass_fraction=s.residual_mass_fraction;
r.h360_residual_ke_fraction=s.residual_ke_fraction;
r.h360_mass_relative_residual=s.mass_relative_residual;
r.h360_energy_relative_residual=s.energy_relative_residual;
r.h360_mass_closure_pass=s.mass_closure_pass;
r.h360_energy_closure_pass=s.energy_closure_pass;
r.h360_invalid_deactivation_count=s.invalid_deactivation_count;
r.h360_structural_bound_pass=s.structural_bound_pass;
end

function r=assign_h720(r,s)
r.h720_terminal_complete=s.terminal_complete;
r.h720_terminal_active_count=s.terminal_active_count;
r.h720_residual_mass_fraction=s.residual_mass_fraction;
r.h720_residual_ke_fraction=s.residual_ke_fraction;
r.h720_cap_contract_pass=s.cap_contract_pass;
r.h720_mass_closure_pass=s.mass_closure_pass;
r.h720_energy_closure_pass=s.energy_closure_pass;
end

function r=assign_h1440(r,s)
r.h1440_terminal_complete=s.terminal_complete;
r.h1440_terminal_active_count=s.terminal_active_count;
r.h1440_residual_mass_fraction=s.residual_mass_fraction;
r.h1440_residual_ke_fraction=s.residual_ke_fraction;
r.h1440_cap_contract_pass=s.cap_contract_pass;
r.h1440_mass_closure_pass=s.mass_closure_pass;
r.h1440_energy_closure_pass=s.energy_closure_pass;
end

function r=assign_publish_metrics(r,out,s)
r.publish_outcome=string(out.outcome);
r.publish_first_break_h_km=out.first_break_h_km;
r.publish_first_break_time_s=out.first_break_time;
r.publish_main_airburst_h_km=out.main_airburst_h_km;
r.publish_main_airburst_time_s=out.main_airburst_t;
r.publish_E_total_lost_kton=out.E_total_lost_kton;
r.publish_E_release_kton=out.E_release_kton;
r.publish_E_burst_raw_kton=out.E_burst_raw_kton;
r.publish_E_burst_coupled_kton=out.E_burst_coupled_kton;
r.publish_dE_dh_peak_kton_per_km=out.dE_dh_peak_kton_per_km;
r.publish_airburst_window_h_min_km=out.airburst_window_h_min_km;
r.publish_airburst_window_h_max_km=out.airburst_window_h_max_km;
r.publish_q_dyn_max_MPa=max(out.timeseries.q_dyn_max_Pa)/1e6;
r.publish_fragment_count_max=max(out.r3b.mass_energy_time_ledger.active_count);
r.publish_breakup_event_count=s.breakup_event_count;
r.publish_created_child_count=s.created_child_count;
r.publish_max_child_id=s.max_child_id;
r.publish_cap_hit_event_count=s.cap_hit_event_count;
r.publish_cap_dropped_mass_kg=s.cap_dropped_mass_kg;
r.publish_insertion_dropped_mass_kg=s.insertion_dropped_mass_kg;
r.publish_terminal_active_count=s.terminal_active_count;
r.publish_terminal_complete=s.terminal_complete;
r.publish_terminal_time_s=s.terminal_time_s;
r.publish_residual_mass_fraction=s.residual_mass_fraction;
r.publish_residual_ke_fraction=s.residual_ke_fraction;
r.publish_mass_relative_residual=s.mass_relative_residual;
r.publish_energy_relative_residual=s.energy_relative_residual;
r.publish_mass_closure_pass=s.mass_closure_pass;
r.publish_energy_closure_pass=s.energy_closure_pass;
r.publish_invalid_deactivation_count=s.invalid_deactivation_count;
r.publish_structural_bound_pass=s.structural_bound_pass;
end

function pass=semantic_valid(out)
outcome=string(out.outcome);
global_values=[out.E_total_lost_kton,out.E_release_kton, ...
    max(out.timeseries.q_dyn_max_Pa)/1e6];
global_ok=all(isfinite(global_values))&&all(global_values>=0);
if outcome=="Airburst"
    local=[out.first_break_h_km,out.first_break_time,out.main_airburst_h_km, ...
        out.main_airburst_t,out.E_burst_raw_kton,out.E_burst_coupled_kton, ...
        out.dE_dh_peak_kton_per_km,out.airburst_window_h_min_km, ...
        out.airburst_window_h_max_km];
    pass=global_ok&&all(isfinite(local))&& ...
        out.main_airburst_h_km>=5&&out.E_burst_raw_kton>0&& ...
        out.E_burst_coupled_kton>0&&out.dE_dh_peak_kton_per_km>0&& ...
        out.airburst_window_h_max_km>=out.airburst_window_h_min_km;
elseif outcome=="Burn-up"||outcome=="Impact"
    pass=global_ok;
else
    pass=false;
end
end

function [T,local_pass,global_pass]=compare_outputs(A,B,label,tol,row)
a=metric_vector(A);b=metric_vector(B);
abs_diff=abs(a-b);
scale=max(abs(a),abs(b));
allowed=max(tol.absolute_tolerance,tol.relative_tolerance.*scale);
metric_pass=false(height(tol),1);
for i=1:height(tol)
    if isnan(a(i))&&isnan(b(i)),metric_pass(i)=true;abs_diff(i)=0;
    elseif isfinite(a(i))&&isfinite(b(i)),metric_pass(i)=abs_diff(i)<=allowed(i);
    end
end
T=table(repmat(row.run_id,height(tol),1), ...
    repmat(row.fragmentation_seed,height(tol),1), ...
    repmat(string(label),height(tol),1),tol.metric,a,b,abs_diff,allowed,metric_pass, ...
    'VariableNames',{'run_id','fragmentation_seed','comparison','metric', ...
    'earlier_value','later_value','absolute_difference','allowed_difference','metric_pass'});
global_mask=tol.metric=="E_total_lost_kton"|tol.metric=="E_release_kton";
local_pass=all(metric_pass(~global_mask));
global_pass=all(metric_pass(global_mask));
end

function v=metric_vector(out)
v=[out.first_break_h_km;out.first_break_time;out.main_airburst_h_km; ...
out.main_airburst_t;out.E_total_lost_kton;out.E_release_kton; ...
out.E_burst_raw_kton;out.E_burst_coupled_kton; ...
out.dE_dh_peak_kton_per_km;out.airburst_window_h_min_km; ...
out.airburst_window_h_max_km;max(out.timeseries.q_dyn_max_Pa)/1e6];
end

function [AUDIT,pass]=pilot_reproduction(R,P,tol)
fields_r=["h360_first_break_h_km","h360_first_break_time_s", ...
"h360_main_airburst_h_km","h360_main_airburst_time_s","h360_E_total_lost_kton", ...
"h360_E_release_kton","h360_E_burst_raw_kton","h360_E_burst_coupled_kton", ...
"h360_dE_dh_peak_kton_per_km","h360_airburst_window_h_min_km", ...
"h360_airburst_window_h_max_km","h360_q_dyn_max_MPa"];
fields_p=["first_break_h_km","first_break_time_s","main_airburst_h_km", ...
"main_airburst_time_s","E_total_lost_kton","E_release_kton","E_burst_raw_kton", ...
"E_burst_coupled_kton","dE_dh_peak_kton_per_km","airburst_window_h_min_km", ...
"airburst_window_h_max_km","q_dyn_max_MPa"];
cells=cell(height(P)*numel(fields_r),1);k=0;
for i=1:height(P)
    rr=R(R.run_id==P.run_id(i),:);
    for j=1:numel(fields_r)
        k=k+1;
        a=P.(fields_p(j))(i);b=rr.(fields_r(j));
        d=abs(a-b);
        allowed=max(tol.absolute_tolerance(j),tol.relative_tolerance(j)*max(abs(a),abs(b)));
        ok=(isnan(a)&&isnan(b))||(isfinite(a)&&isfinite(b)&&d<=allowed);
        cells{k}=table(P.run_id(i),fields_p(j),a,b,d,allowed,ok, ...
            'VariableNames',{'run_id','metric','frozen_pilot_value', ...
            'full_campaign_H360_value','absolute_difference','allowed_difference','metric_pass'});
    end
end
AUDIT=vertcat(cells{:});
pass=all(AUDIT.metric_pass);
end

function [SENT,PARITY]=full_campaign_parity_sentinels(D,R)
means=[19160,18.3,19.8,3300,7.8e6,0.015,0.80];
stds=[600,1.5,1.5,400,3e6,0.003,0.08];
X=[D.V0_m_s,D.entry_angle_deg,D.D0_m,D.rho_meteor_kg_m3, ...
    D.Yield_Stress_Pa,D.Ch,D.energy_efficiency];
Z=(X-means)./stds;
[~,idx_nom]=min(sum(Z.^2,2));
ids=D.run_id(idx_nom);
reasons="nominal_nearest";
names=["V0","entry_angle","D0","rho","Y0","Ch","eta"];
for j=1:7
    [~,imin]=min(X(:,j));[~,imax]=max(X(:,j));
    [ids,reasons]=append_reason(ids,reasons,D.run_id(imin),names(j)+"_minimum");
    [ids,reasons]=append_reason(ids,reasons,D.run_id(imax),names(j)+"_maximum");
end
[~,risk_mass]=maxk(R.h360_residual_mass_fraction,3);
[~,risk_ke]=maxk(R.h360_residual_ke_fraction,3);
for idx=[risk_mass(:);risk_ke(:)]'
    [ids,reasons]=append_reason(ids,reasons,R.run_id(idx),"high_residual_risk");
end
ext=find(logical(R.extension_required));
for idx=ext'
    [ids,reasons]=append_reason(ids,reasons,R.run_id(idx),"adaptive_extension_case");
end
[ids,ord]=sort(ids);reasons=reasons(ord);
SENT=table(ids,reasons,'VariableNames',{'run_id','sentinel_reason'});

P=repmat(empty_parity(),height(SENT),1);
for i=1:height(SENT)
    row=D(D.run_id==SENT.run_id(i),:);
    horizon=R.publish_metric_source_horizon_s(R.run_id==SENT.run_id(i));
    params=build_params(row,horizon);
    diag=meteor_single_run_r3b3_diagnostic(params);
    rng_diag=rng;
    ref=meteor_single_run(params);
    rng_ref=rng;
    p=empty_parity();
    p.run_id=row.run_id;
    p.horizon_s=horizon;
    p.sentinel_reason=SENT.sentinel_reason(i);
    [p.scalar_parity_pass,p.scalar_max_abs_diff]=scalar_parity(ref,diag);
    [p.timeseries_parity_pass,p.timeseries_max_abs_diff]=table_parity(ref.timeseries,diag.timeseries);
    [p.deposition_parity_pass,p.deposition_max_abs_diff]=table_parity(ref.deposition_profile,diag.deposition_profile);
    p.outcome_match=string(ref.outcome)==string(diag.outcome);
    p.rng_state_match=strcmp(rng_ref.Type,rng_diag.Type)&& ...
        isequal(rng_ref.Seed,rng_diag.Seed)&&isequal(rng_ref.State,rng_diag.State);
    p.parity_pass=p.scalar_parity_pass&&p.timeseries_parity_pass&& ...
        p.deposition_parity_pass&&p.outcome_match&&p.rng_state_match;
    P(i)=p;
    clear ref diag
end
PARITY=struct2table(P);
end

function [ids,reasons]=append_reason(ids,reasons,id,reason)
loc=find(ids==id,1);
if isempty(loc)
    ids(end+1,1)=id;reasons(end+1,1)=string(reason);
else
    reasons(loc)=reasons(loc)+";"+string(reason);
end
end

function write_flagged_ledgers(root,row,base,selected,r)
prefix=sprintf('run_%04d',row.run_id);
writetable(add_context(base.r3b.closure_summary,row,360), ...
    fullfile(root,[prefix '_H360_closure.csv']));
writetable(add_context(base.r3b.terminal_fragment_ledger,row,360), ...
    fullfile(root,[prefix '_H360_terminal_fragments.csv']));
if r.publish_metric_source_horizon_s~=360
    writetable(add_context(selected.r3b.closure_summary,row,r.publish_metric_source_horizon_s), ...
        fullfile(root,[prefix '_final_closure.csv']));
    writetable(add_context(selected.r3b.terminal_fragment_ledger,row,r.publish_metric_source_horizon_s), ...
        fullfile(root,[prefix '_final_terminal_fragments.csv']));
end
end

function T=add_context(T,row,horizon)
n=height(T);
if n==0,return;end
T=addvars(T,repmat(row.run_id,n,1),repmat(row.fragmentation_seed,n,1), ...
    repmat(horizon,n,1),'Before',1, ...
    'NewVariableNames',{'run_id','fragmentation_seed','horizon_s'});
end

function T=vertcat_nonempty(C)
keep=false(numel(C),1);
for i=1:numel(C),keep(i)=istable(C{i})&&height(C{i})>0;end
if any(keep)
    S=C(keep);names=S{1}.Properties.VariableNames;
    for i=2:numel(S)
        if ~isequal(S{i}.Properties.VariableNames,names)
            error('R3C2C:Schema','Metric-audit schema mismatch.');
        end
    end
    T=vertcat(S{:});
else,T=table();
end
end

function rho=spearman_no_toolbox(x,y)
rx=average_rank(x);ry=average_rank(y);
rx=rx-mean(rx);ry=ry-mean(ry);
den=sqrt(sum(rx.^2)*sum(ry.^2));
if den==0,rho=NaN;else,rho=sum(rx.*ry)/den;end
end

function r=average_rank(x)
[xs,ord]=sort(x(:));n=numel(x);rs=zeros(n,1);i=1;
while i<=n
    j=i;
    while j<n&&xs(j+1)==xs(i),j=j+1;end
    rs(i:j)=0.5*(i+j);i=j+1;
end
r=zeros(n,1);r(ord)=rs;
end

function make_distribution_figure(R,folder)
fig=figure('Visible','off','Color','w','Position',[80 80 1300 850]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
nexttile(tl);histogram(R.publish_main_airburst_h_km,25);grid on;
xlabel('Main airburst altitude [km]');ylabel('Count');
nexttile(tl);histogram(R.publish_E_burst_coupled_kton,25);grid on;
xlabel('Coupled burst energy [kt]');ylabel('Count');
nexttile(tl);histogram(R.publish_dE_dh_peak_kton_per_km,25);grid on;
xlabel('Peak positive-decrement density [kt km^{-1}]');ylabel('Count');
nexttile(tl);histogram(R.publish_q_dyn_max_MPa,25);grid on;
xlabel('Maximum dynamic pressure [MPa]');ylabel('Count');
title(tl,'Corrected full-1000 Monte Carlo distributions');
exportgraphics(fig,fullfile(folder,'R3C2C_full1000_distributions.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(folder,'R3C2C_full1000_distributions.png'),'Resolution',300);
close(fig);
end

function make_correlation_figure(C,inputs,outputs,folder)
M=nan(numel(inputs),numel(outputs));
for i=1:numel(inputs)
    for j=1:numel(outputs)
        row=C(C.input_parameter==inputs(i)&C.output_metric==outputs(j),:);
        if height(row)==1,M(i,j)=row.spearman_rho;end
    end
end
fig=figure('Visible','off','Color','w','Position',[80 80 1250 720]);
imagesc(M,[-1 1]);colorbar;
set(gca,'XTick',1:numel(outputs),'XTickLabel',outputs, ...
    'YTick',1:numel(inputs),'YTickLabel',inputs,'TickLabelInterpreter','none');
xtickangle(30);
for i=1:size(M,1)
    for j=1:size(M,2)
        text(j,i,sprintf('%.2f',M(i,j)),'HorizontalAlignment','center');
    end
end
title('Spearman rank correlation: frozen inputs versus corrected outputs');
exportgraphics(fig,fullfile(folder,'R3C2C_spearman_correlation_heatmap.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(folder,'R3C2C_spearman_correlation_heatmap.png'),'Resolution',300);
close(fig);
end

function make_stability_figure(T,folder)
metrics=unique(T.metric,'stable');
fig=figure('Visible','off','Color','w','Position',[80 80 1300 820]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for i=1:min(4,numel(metrics))
    nexttile(tl);
    A=T(T.metric==metrics(i),:);
    plot(A.prefix_size,A.p05,'--o',A.prefix_size,A.median,'-s', ...
        A.prefix_size,A.p95,':^','LineWidth',1.2);
    grid on;xlabel('Prefix ensemble size');ylabel('Metric value');
    title(strrep(metrics(i),'_','\_'));legend({'P05','Median','P95'},'Location','best');
end
title(tl,'Full-1000 prefix-size stability diagnostic');
exportgraphics(fig,fullfile(folder,'R3C2C_ensemble_prefix_stability.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(folder,'R3C2C_ensemble_prefix_stability.png'),'Resolution',300);
close(fig);
end

function q=pct(x,p)
x=sort(x(:));n=numel(x);
if n==0,q=NaN;return;end
if n==1,q=x;return;end
pos=1+(n-1)*p/100;lo=floor(pos);hi=ceil(pos);
if lo==hi,q=x(lo);else,q=x(lo)+(pos-lo)*(x(hi)-x(lo));end
end

function [pass,max_diff]=scalar_parity(A,B)
fields=["first_break_h_km","first_break_time","main_airburst_h_km", ...
"main_airburst_t","E_release_kton","E_total_lost_kton", ...
"dE_dh_peak_kton_per_km","E_burst_raw_kton","E_burst_coupled_kton", ...
"airburst_window_h_min_km","airburst_window_h_max_km","final_mass","area_multiplier"];
pass=true;max_diff=0;
for i=1:numel(fields)
    a=double(A.(fields(i)));b=double(B.(fields(i)));
    if isnan(a)&&isnan(b),d=0;
    elseif isfinite(a)&&isfinite(b)
        d=abs(a-b);allowed=max(1e-8,1e-12*max(1,max(abs([a,b]))));
        pass=pass&&d<=allowed;
    else,d=Inf;pass=false;
    end
    max_diff=max(max_diff,d);
end
end

function [pass,max_diff]=table_parity(A,B)
if ~isequal(A.Properties.VariableNames,B.Properties.VariableNames)|| ...
        height(A)~=height(B)||width(A)~=width(B)
    pass=false;max_diff=Inf;return;
end
AA=double(A{:,:});BB=double(B{:,:});
nan_match=isnan(AA)&isnan(BB);
bad=xor(isnan(AA),isnan(BB))|xor(isinf(AA),isinf(BB));
D=abs(AA-BB);D(nan_match)=0;D(~isfinite(D))=0;
if any(bad,'all'),pass=false;max_diff=Inf;return;end
if isempty(D),max_diff=0;else,max_diff=max(D,[],'all');end
fv=[AA(isfinite(AA));BB(isfinite(BB))];
if isempty(fv),scale=1;else,scale=max(1,max(abs(fv)));end
allowed=max(1e-8,1e-12*scale);pass=max_diff<=allowed;
end

function r=empty_result()
r=struct( ...
'run_id',NaN,'input_design_seed',NaN,'fragmentation_seed',NaN, ...
'V0_m_s',NaN,'entry_angle_deg',NaN,'D0_m',NaN, ...
'rho_meteor_kg_m3',NaN,'Yield_Stress_Pa',NaN,'Ch',NaN, ...
'energy_efficiency',NaN,'area_multiplier',NaN,'max_total_fragments',NaN, ...
'dt_s',NaN,'base_horizon_s',NaN,'run_completed',false,'runtime_s',NaN, ...
'error_identifier',"",'error_message',"", ...
'h360_run_completed',false,'h360_runtime_s',NaN,'h360_outcome',"", ...
'h360_first_break_h_km',NaN,'h360_first_break_time_s',NaN, ...
'h360_main_airburst_h_km',NaN,'h360_main_airburst_time_s',NaN, ...
'h360_E_total_lost_kton',NaN,'h360_E_release_kton',NaN, ...
'h360_E_burst_raw_kton',NaN,'h360_E_burst_coupled_kton',NaN, ...
'h360_dE_dh_peak_kton_per_km',NaN, ...
'h360_airburst_window_h_min_km',NaN,'h360_airburst_window_h_max_km',NaN, ...
'h360_q_dyn_max_MPa',NaN,'h360_fragment_count_max',NaN, ...
'h360_breakup_event_count',NaN,'h360_created_child_count',NaN, ...
'h360_max_child_id',NaN,'h360_cap_hit_event_count',NaN, ...
'h360_cap_dropped_mass_kg',NaN,'h360_insertion_dropped_mass_kg',NaN, ...
'h360_terminal_active_count',NaN,'h360_terminal_complete',false, ...
'h360_terminal_time_s',NaN,'h360_residual_mass_fraction',NaN, ...
'h360_residual_ke_fraction',NaN,'h360_mass_relative_residual',NaN, ...
'h360_energy_relative_residual',NaN,'h360_mass_closure_pass',false, ...
'h360_energy_closure_pass',false,'h360_invalid_deactivation_count',NaN, ...
'h360_structural_bound_pass',false,'h360_semantic_output_pass',false, ...
'h360_residual_threshold_pass',false,'h360_hard_contract_pass',false, ...
'extension_required',false,'h720_run_completed',false,'h720_runtime_s',NaN, ...
'h720_terminal_complete',false,'h720_terminal_active_count',NaN, ...
'h720_residual_mass_fraction',NaN,'h720_residual_ke_fraction',NaN, ...
'h720_cap_contract_pass',false,'h720_mass_closure_pass',false, ...
'h720_energy_closure_pass',false,'h1440_run_completed',false, ...
'h1440_runtime_s',NaN,'h1440_terminal_complete',false, ...
'h1440_terminal_active_count',NaN,'h1440_residual_mass_fraction',NaN, ...
'h1440_residual_ke_fraction',NaN,'h1440_cap_contract_pass',false, ...
'h1440_mass_closure_pass',false,'h1440_energy_closure_pass',false, ...
'h720_to_h1440_local_convergence_pass',false, ...
'h720_to_h1440_global_convergence_pass',false, ...
'extension_convergence_pass',false,'publish_metric_source_horizon_s',NaN, ...
'publish_outcome',"",'publish_first_break_h_km',NaN, ...
'publish_first_break_time_s',NaN,'publish_main_airburst_h_km',NaN, ...
'publish_main_airburst_time_s',NaN,'publish_E_total_lost_kton',NaN, ...
'publish_E_release_kton',NaN,'publish_E_burst_raw_kton',NaN, ...
'publish_E_burst_coupled_kton',NaN, ...
'publish_dE_dh_peak_kton_per_km',NaN, ...
'publish_airburst_window_h_min_km',NaN, ...
'publish_airburst_window_h_max_km',NaN,'publish_q_dyn_max_MPa',NaN, ...
'publish_fragment_count_max',NaN,'publish_breakup_event_count',NaN, ...
'publish_created_child_count',NaN,'publish_max_child_id',NaN, ...
'publish_cap_hit_event_count',NaN,'publish_cap_dropped_mass_kg',NaN, ...
'publish_insertion_dropped_mass_kg',NaN,'publish_terminal_active_count',NaN, ...
'publish_terminal_complete',false,'publish_terminal_time_s',NaN, ...
'publish_residual_mass_fraction',NaN,'publish_residual_ke_fraction',NaN, ...
'publish_mass_relative_residual',NaN,'publish_energy_relative_residual',NaN, ...
'publish_mass_closure_pass',false,'publish_energy_closure_pass',false, ...
'publish_invalid_deactivation_count',NaN,'publish_structural_bound_pass',false, ...
'publish_semantic_output_pass',false,'publish_residual_threshold_pass',false, ...
'publish_hard_contract_pass',false,'run_gate_pass',false);
end

function p=empty_parity()
p=struct('run_id',NaN,'horizon_s',NaN,'sentinel_reason',"", ...
'scalar_parity_pass',false,'scalar_max_abs_diff',NaN, ...
'timeseries_parity_pass',false,'timeseries_max_abs_diff',NaN, ...
'deposition_parity_pass',false,'deposition_max_abs_diff',NaN, ...
'outcome_match',false,'rng_state_match',false,'parity_pass',false);
end

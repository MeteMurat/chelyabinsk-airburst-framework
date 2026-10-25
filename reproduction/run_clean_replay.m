function result = run_clean_replay(repositoryRoot, mode)
% Portable R1-C6 replay harness. Historical scientific kernels are unchanged.
% mode: 'smoke' (5 shape cases + 8 MC rows) or 'full1000' (plus exact original campaign).
% Inputs under expected/ are read ONLY for comparison, not as simulation outputs.
if nargin < 1, repositoryRoot = fileparts(fileparts(mfilename('fullpath'))); end
if nargin < 2, mode = 'smoke'; end
repositoryRoot = char(repositoryRoot); mode = lower(char(mode));
assert(any(strcmp(mode, {'smoke','full1000'})), 'Unknown replay mode.');
root = repositoryRoot;
outdir = fullfile(root,'replay_outputs',mode);
if ~isfolder(outdir), mkdir(outdir); end
validationDir = fullfile(root,'validation');
if ~isfolder(validationDir), mkdir(validationDir); end
checks = struct('name',{},'passed',{},'max_absolute_error',{});
result = struct('schema_version',1,'mode',mode,'status','RUNNING', ...
    'matlab_forward_replay_executed',false,'mc_rows_replayed',0, ...
    'population_rows_recomputed',0,'scientific_model_modified',false, ...
    'full_reproduction_verified',false);
result.started_utc = char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
try
    % Remove user/project path entries for this batch session only. Never savepath.
    oldPath = path;
    cleanupPath = onCleanup(@() path(oldPath)); %#ok<NASGU>
    entries = strsplit(path,pathsep); mr = [matlabroot filesep];
    keep = cellfun(@(p) startsWith([p filesep],mr,'IgnoreCase',ispc),entries);
    path(strjoin(entries(keep),pathsep));
    addpath(fullfile(root,'01_final_core'),fullfile(root,'reproduction'));
    cd(root);
    rng(42,'twister');
    set(groot,'defaultFigureVisible','off');
    result.matlab_release = version('-release');
    result.matlab_version = version;
    result.architecture = computer('arch');
    result.replay_execution = 'CPU sequential; no parpool or gpuArray requested by replay harness';
    result.rng_generator = 'twister';
    check_flag('mapping_readgeoraster_available', exist('readgeoraster','file')~=0);
    manifest = jsondecode(fileread(fullfile(root,'SOURCE_MANIFEST.json')));
    for k=1:numel(manifest.files)
        f=manifest.files(k);
        check_flag(['source_identity:' f.path],strcmpi(sha256file(fullfile(root,f.path)),f.sha256));
    end
    bins=jsondecode(fileread(fullfile(root,'inputs','required_binary_inputs.json')));
    for k=1:numel(bins)
        check_flag(['binary_identity:' bins(k).filename], ...
            strcmpi(sha256file(fullfile(root,bins(k).destination)),bins(k).sha256));
    end
    for name={'meteor_single_run','meteor_single_run_r3b3_diagnostic','run_r3c2c_r2_corrected_full1000'}
        actual=which(name{1}); intended=fullfile(root,'01_final_core',[name{1} '.m']);
        check_flag(['path_resolution:' name{1}],strcmpi(strrep(actual,'\','/'),strrep(intended,'\','/')));
    end
    C=jsondecode(fileread(fullfile(root,'inputs','replay_config.json')));
    validate_population_table_shapes(C);
    expected=@(name) read_csv(fullfile(root,'expected',name));
    fprintf('\n=== CLEAN MATLAB REPLAY: GEOMETRY AND LOCAL RADIUS ===\n');
    invoke_geometry_summary();
    rng(42,'twister');
    invoke_local_radius_summary();
    geo=read_csv(fullfile(root,'05_outputs','csv','geometry_summary.csv'));
    rad=read_csv(fullfile(root,'05_outputs','summary','local_radius_summary.csv'));
    check_table('geometry',geo,expected('geometry.csv'),1e-9,1e-9);
    check_table('local_radius',rad,expected('local_radius.csv'),1e-9,1e-9);
    % Derive the area inputs from newly computed geometry, then compare to the contract.
    derivedChi=[1;geo.projected_area_ratio_p10_to_sphere; ...
        geo.projected_area_ratio_p50_to_sphere;geo.projected_area_ratio_p90_to_sphere; ...
        geo.projected_area_ratio_mean_to_sphere];
    check_num('geometry_to_G0_G4_contract',derivedChi,[C.g0g4.area_multiplier]',1e-12,1e-10);

    fprintf('\n=== CLEAN MATLAB REPLAY: NOMINAL AND G0-G4 ===\n');
    gRows=cell(5,1); nominal=[]; nominalRng=[];
    for k=1:5
        p=C.nominal; p.area_multiplier=C.g0g4(k).area_multiplier;
        rng(42,'twister');
        fprintf('Replaying %s with chi_A=%.15g\n',C.g0g4(k).case_id,p.area_multiplier);
        o=meteor_single_run_r3b3_diagnostic(p);
        check_kernel_gates(o,false);
        s=summarize_case(o,p,C.g0g4(k).case_id);
        gr=struct('case_id',string(C.g0g4(k).case_id),'regime_id',"R128_T360", ...
            'chi_A',p.area_multiplier,'first_breakup_altitude_km',o.first_break_h_km, ...
            'coupled_burst_energy_kt',o.E_burst_coupled_kton, ...
            'peak_deposition_density_kt_per_km',o.dE_dh_peak_kton_per_km, ...
            'max_dynamic_pressure_MPa',s.q_dyn_max_MPa, ...
            'max_total_fragments',p.max_total_fragments,'tmax_s',p.tmax, ...
            'theoretical_full_tree_id_count',86,'cap_hit_event_count',s.cap_hit_event_count, ...
            'cap_dropped_mass_kg',s.cap_dropped_mass_kg, ...
            'insertion_dropped_mass_kg',s.insertion_dropped_mass_kg, ...
            'full_tree_realized',s.full_tree_realized,'terminal_complete',s.terminal_complete, ...
            'closure_pass',s.mass_relative_residual<=1e-9 && s.energy_relative_residual<=1e-9);
        gRows{k}=struct2table(gr);
        if strcmp(C.g0g4(k).case_id,'G4')
            nominal=o; nominalSummary=s; nominalRng=rng;
        end
    end
    G=vertcat(gRows{:}); check_table('G0-G4',G,expected('g0g4.csv'),1e-8,1e-10);
    writetable(G,fullfile(outdir,'g0g4_replayed.csv'));
    check_flag('nominal_selected',~isempty(nominal));
    check_table('nominal',struct2table(nominalSummary),expected('nominal.csv'),1e-8,1e-10);
    check_table('nominal_timeseries',nominal.timeseries,expected('nominal_timeseries.csv'),1e-8,1e-10);
    check_table('nominal_deposition',nominal.deposition_profile,expected('nominal_deposition.csv'),1e-8,1e-10);
    rng(42,'twister'); ref=meteor_single_run(C.nominal); referenceRng=rng;
    fields={'first_break_h_km','first_break_time','main_airburst_h_km','main_airburst_t', ...
        'E_release_kton','E_total_lost_kton','E_burst_raw_kton','E_burst_coupled_kton', ...
        'dE_dh_peak_kton_per_km','airburst_window_h_min_km','airburst_window_h_max_km'};
    for k=1:numel(fields), check_num(['reference_diagnostic:' fields{k}], ...
            ref.(fields{k}),nominal.(fields{k}),1e-8,1e-12); end
    check_flag('reference_diagnostic_rng_state',isequal(nominalRng,referenceRng));
    check_table('reference_diagnostic_timeseries',ref.timeseries,nominal.timeseries,1e-8,1e-12);
    check_table('reference_diagnostic_deposition',ref.deposition_profile,nominal.deposition_profile,1e-8,1e-12);
    writetable(struct2table(nominalSummary),fullfile(outdir,'nominal_replayed.csv'));
    writetable(nominal.timeseries,fullfile(outdir,'nominal_timeseries_replayed.csv'));
    writetable(nominal.deposition_profile,fullfile(outdir,'nominal_deposition_replayed.csv'));
    save(fullfile(outdir,'nominal_diagnostic.mat'),'nominal','-v7');
    result.matlab_forward_replay_executed=true;
    ref=[]; %#ok<NASGU>

    fprintf('\n=== CLEAN MATLAB REPLAY: EIGHT FROZEN MC RUNS ===\n');
    design=read_csv(fullfile(root,'inputs','mc_design.csv'));
    check_flag('design_is_1000_unique_runs',height(design)==1000 && isequal(design.run_id,(1:1000)'));
    smokeRows=cell(numel(C.smoke_run_ids),1);
    for k=1:numel(C.smoke_run_ids)
        id=C.smoke_run_ids(k); row=design(design.run_id==id,:);
        check_flag('smoke_input_identity',height(row)==1);
        p=params_from_design(row,C.nominal);
        rng(42,'twister');
        fprintf('Smoke MC run_id=%d, fragmentation_seed=%d\n',row.run_id,row.fragmentation_seed);
        o=meteor_single_run_r3b3_diagnostic(p);
        check_kernel_gates(o,true);
        smokeRows{k}=struct2table(as_publishable(o,row,p));
    end
    smoke=vertcat(smokeRows{:}); archived=expected('mc_metrics.csv');
    target=archived(ismember(archived.run_id,C.smoke_run_ids),:);
    smoke=sortrows(smoke,'run_id'); target=sortrows(target,'run_id');
    check_table('smoke_MC',smoke,target,1e-8,1e-10);
    writetable(smoke,fullfile(outdir,'smoke_MC_replayed.csv'));
    result.mc_rows_replayed=height(smoke);
    Psmoke=compute_pressure(smoke,C);
    pressureExpected=expected('pressure.csv');
    check_table('smoke_pressure',Psmoke, ...
        pressureExpected(ismember(pressureExpected.run_id,smoke.run_id),:),1e-8,1e-10);

    fprintf('\n=== CLEAN MATLAB REPLAY: NATIVE RASTER EXPOSURE ===\n');
    raster=prepare_population(fullfile(root,'04_data','population','turkiye.tif'),C);
    popExpected=expected('population.csv');
    popSmoke=compute_population(Psmoke,raster,C);
    compare_population('smoke_population',popSmoke, ...
        popExpected(ismember(popExpected.run_id,smoke.run_id),:));
    result.population_rows_recomputed=height(popSmoke);
    writetable(popSmoke,fullfile(outdir,'smoke_population_replayed.csv'));
    nrow=table(0,C.nominal.seed,nominal.main_airburst_h_km,nominal.E_burst_coupled_kton, ...
        'VariableNames',{'run_id','fragmentation_seed','publish_main_airburst_h_km','publish_E_burst_coupled_kton'});
    Pnom=compute_pressure(nrow,C); NomPressure=expected('nominal_pressure.csv');
    rr=[Pnom.R_3_45_kPa_km;Pnom.R_1_00_kPa_km;Pnom.R_0_50_kPa_km;Pnom.R_0_20_kPa_km];
    check_num('nominal_pressure_from_new_G4',rr,NomPressure.R_ground_km,1e-8,1e-10);
    NP=compute_population(Pnom,raster,C); NE=expected('nominal_population.csv');
    % The old population wrapper rounded nominal E_c to 331.6302 kt.
    % New chained radius uses unrounded rerun E_c; compare population/cell outputs,
    % and record the small historical radius difference separately (no recalibration).
    check_num('nominal_population_counts',NP.included_cell_count,NE.included_cell_count,0,0);
    check_num('nominal_population_people',NP.exposed_population,NE.exposed_population,1e-4,1e-11);
    check_num('nominal_population_support_area',NP.population_grid_support_area_km2,NE.affected_area_km2,1e-6,1e-11);
    result.nominal_historical_radius_delta_km=max(abs(NP.radius_km-NE.radius_km));
    check_flag('historical_nominal_rounding_difference_documented',result.nominal_historical_radius_delta_km<1e-5);
    writetable(NP,fullfile(outdir,'nominal_population_chained_replayed.csv'));
    result.status='PASS_SMOKE_REPLAY_FULL1000_PENDING';
    result.check_count=numel(checks); result.smoke_pass=true; result.all_checks_pass=all([checks.passed]);
    write_json(fullfile(validationDir,'smoke_result.json'),result);
    writetable(struct2table(checks),fullfile(outdir,'checks.csv'));

    if strcmp(mode,'full1000')
        fprintf('\n=== CLEAN MATLAB REPLAY: ORIGINAL FULL-1000 DRIVER ===\n');
        campaignRoot=fullfile(root,'replay_outputs','campaign_full1000');
        % Original driver supports checkpoints. The Python launcher binds workspace
        % reuse to the same source+input manifest; expected matrices are never copied here.
        rng(42,'twister');
        run_r3c2c_r2_corrected_full1000(campaignRoot, ...
            fullfile(root,'inputs','mc_design.csv'),fullfile(root,'inputs','pilot64.csv'), ...
            fullfile(root,'inputs','historical_horizon_gate.csv'));
        campaignOut=fullfile(campaignRoot,'02_full1000_outputs');
        PUB=read_csv(fullfile(campaignOut,'R3C2C_full1000_publishable_metric_matrix.csv'));
        PUB=sortrows(PUB,'run_id');
        check_flag('full1000_has_exact_run_ids',height(PUB)==1000 && isequal(PUB.run_id,(1:1000)'));
        check_table('full1000_publishable_metrics',PUB,archived,1e-8,1e-10);
        gates=read_csv(fullfile(campaignOut,'R3C2C_full1000_gate_adjudication.csv'));
        check_flag('original_full1000_all_gates_pass',all(truth(gates.gate_pass)));
        raw=read_csv(fullfile(campaignOut,'R3C2C_full1000_results.csv'));
        old=expected('mc_diagnostics.csv');
        check_flag('full1000_all_runs_completed',all(truth(raw.run_completed)) && all(truth(raw.run_gate_pass)));
        for f={'publish_cap_hit_event_count','publish_cap_dropped_mass_kg', ...
                'publish_insertion_dropped_mass_kg','publish_max_child_id', ...
                'publish_terminal_active_count','publish_terminal_complete','extension_required'}
            check_num(['full1000_diagnostic:' f{1}],raw.(f{1}),old.(f{1}),0,0);
        end
        result.terminal_incomplete_count=sum(~truth(raw.publish_terminal_complete));
        result.terminal_incomplete_is_not_automatic_failure=true;
        P=compute_pressure(PUB,C);
        check_table('full1000_pressure',P,pressureExpected,1e-8,1e-10);
        E=compute_population(P,raster,C);
        compare_population('full1000_population',E,popExpected);
        writetable(P,fullfile(outdir,'full1000_pressure_replayed.csv'));
        writetable(E,fullfile(outdir,'full1000_population_replayed.csv'));
        result.mc_rows_replayed=height(PUB);
        result.population_rows_recomputed=height(E);
        result.mc_burst_altitude_quantiles_km=linear_quantile(PUB.publish_main_airburst_h_km,[.1,.5,.9]);
        result.mc_coupled_energy_quantiles_kt=linear_quantile(PUB.publish_E_burst_coupled_kton,[.1,.5,.9]);
        check_num('full1000_altitude_quantiles',result.mc_burst_altitude_quantiles_km, ...
            linear_quantile(archived.publish_main_airburst_h_km,[.1,.5,.9]),1e-9,1e-10);
        check_num('full1000_energy_quantiles',result.mc_coupled_energy_quantiles_kt, ...
            linear_quantile(archived.publish_E_burst_coupled_kton,[.1,.5,.9]),1e-8,1e-10);
        result.full_reproduction_verified=true;
        result.status='PASS_CLEAN_MATLAB_FORWARD_REPLAY_FULL1000_AND_POPULATION';
    end
    % Recheck immutable public source/input/reference files after all computations.
    for k=1:numel(manifest.files)
        f=manifest.files(k);
        check_flag(['source_unchanged:' f.path],strcmpi(sha256file(fullfile(root,f.path)),f.sha256));
    end
    result.finished_utc=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
    result.check_count=numel(checks); result.all_checks_pass=all([checks.passed]);
    result.nominal_burst_altitude_km=nominal.main_airburst_h_km;
    result.nominal_coupled_energy_kt=nominal.E_burst_coupled_kton;
    writetable(struct2table(checks),fullfile(outdir,'checks.csv'));
    write_json(fullfile(validationDir,[mode '_result.json']),result);
    write_json(fullfile(validationDir,'replay_status.json'),result);
    fprintf('\nmatlab_forward_replay_executed=True\nmc_rows_replayed=%d\npopulation_rows_recomputed=%d\nstatus=%s\n', ...
        result.mc_rows_replayed,result.population_rows_recomputed,result.status);
catch ME
    result.status='FAIL_CLEAN_MATLAB_REPLAY'; result.full_reproduction_verified=false;
    result.error_identifier=ME.identifier;
    % Full stack/log remains local. Do not publish machine paths from exceptions.
    result.error_summary='A replay check failed. Inspect local matlab_console.log and checks.csv.';
    result.check_count=numel(checks);
    if ~isempty(checks), writetable(struct2table(checks),fullfile(outdir,'checks.csv')); end
    write_json(fullfile(validationDir,[mode '_result.json']),result);
    write_json(fullfile(validationDir,'replay_status.json'),result);
    fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off'));
    rethrow(ME);
end

    function check_flag(name,ok)
        checks(end+1)=struct('name',char(name),'passed',logical(ok),'max_absolute_error',0);
        assert(ok,'Replay:CheckFailed','Check failed: %s',name);
    end
    function check_num(name,a,b,atol,rtol)
        a=double(a); b=double(b);
        if ~isequal(size(a),size(b))
            checks(end+1)=struct('name',char(name),'passed',false,'max_absolute_error',Inf);
            error('Replay:ShapeMismatch','Shape mismatch in %s',name);
        end
        bothNan=isnan(a)&isnan(b); finite=isfinite(a)&isfinite(b);
        d=abs(a-b); d(bothNan)=0;
        % Same elementwise tolerances; split the expression to avoid delimiter errors.
        scale = max(abs(a), abs(b));
        tolerance = max(atol, rtol * scale);
        withinTolerance = finite & (d <= tolerance);
        ok = all(bothNan(:) | withinTolerance(:));
        if isempty(d),maxd=0;else,maxd=max(d(:));end
        checks(end+1)=struct('name',char(name),'passed',logical(ok),'max_absolute_error',maxd);
        assert(ok,'Replay:NumericMismatch','Numeric mismatch in %s; max_abs_error=%.16g',name,maxd);
    end
    function check_table(label,A,B,atol,rtol)
        check_flag([label ':row_count'],height(A)==height(B));
        for j=1:width(B)
            f=B.Properties.VariableNames{j};
            check_flag([label ':column:' f],any(strcmp(f,A.Properties.VariableNames)));
            if isnumeric(B.(f)) || islogical(B.(f))
                check_num([label ':' f],A.(f),B.(f),atol,rtol);
            else
                check_flag([label ':' f],isequal(string(A.(f)),string(B.(f))));
            end
        end
    end
    function compare_population(label,A,B)
        A=sortrows(A,{'run_id','threshold_kpa'},{'ascend','descend'});
        B=sortrows(B,{'run_id','threshold_kpa'},{'ascend','descend'});
        check_num([label ':run_id'],A.run_id,B.run_id,0,0);
        check_num([label ':threshold'],A.threshold_kpa,B.threshold_kpa,0,0);
        check_num([label ':radius'],A.radius_km,B.radius_km,1e-8,1e-10);
        check_num([label ':included_cells'],A.included_cell_count,B.included_cell_count,0,0);
        check_num([label ':population'],A.exposed_population,B.exposed_population,1e-4,1e-11);
        check_num([label ':support_area'],A.population_grid_support_area_km2,B.affected_area_km2,1e-6,1e-11);
        check_num([label ':fraction'],A.exposure_fraction,B.exposure_fraction,1e-12,1e-10);
    end
    function check_kernel_gates(o,residualRequired)
        S=o.r3b.closure_summary; B=o.r3b.breakup_event_ledger; Ch=o.r3b.child_creation_ledger;
        check_flag('kernel_zero_cap_hits',sum(B.cap_hit)==0);
        check_flag('kernel_no_cap_mass_loss',abs(sum(B.cap_dropped_mass_kg))<=1e-9);
        check_flag('kernel_no_insertion_mass_loss',abs(S.cumulative_insertion_dropped_mass_kg)<=1e-9);
        closureTol=1e-9; if ~residualRequired,closureTol=1e-10;end
        check_flag('kernel_mass_closure',abs(S.final_mass_closure_residual_kg)/max(S.initial_mass_kg,eps)<=closureTol);
        check_flag('kernel_energy_closure',abs(S.final_energy_closure_residual_J)/max(S.initial_ke_J,eps)<=closureTol);
        check_flag('kernel_tree_bound',height(B)<=26 && height(Ch)<=85 && ...
            (isempty(Ch)||(max(Ch.child_id)<=86 && all(Ch.inserted_into_active_pool))));
        de=o.r3b.deactivation_ledger;
        check_flag('kernel_valid_deactivation_reasons',isempty(de)||all( ...
            de.reason=="post_step_ground_crossing"|de.reason=="post_step_mass_depletion"));
        check_flag('kernel_airburst_semantics',string(o.outcome)=="Airburst" && ...
            isfinite(o.main_airburst_h_km) && o.main_airburst_h_km>=5 && ...
            isfinite(o.E_burst_coupled_kton) && o.E_burst_coupled_kton>0);
        if residualRequired
            check_flag('kernel_residual_thresholds',abs(S.final_active_mass_kg)/max(S.initial_mass_kg,eps)<=1e-3 && ...
                abs(S.final_active_ke_J)/max(S.initial_ke_J,eps)<=1e-8);
        end
    end
end

function invoke_geometry_summary()
% A separate function workspace contains the historical script's clear command.
make_stl_geometry_summary;
end
function invoke_local_radius_summary()
make_stl_local_radius_summary;
end
function T=read_csv(p)
T=readtable(p,'Delimiter',',','TextType','string','VariableNamingRule','preserve');
end
function p=params_from_design(r,p)
p.rho_meteor=r.rho_meteor_kg_m3; p.V0=r.V0_m_s;
p.entry_angle_deg=r.entry_angle_deg; p.D0=r.D0_m;
p.Yield_Stress_0=r.Yield_Stress_Pa; p.Ch=r.Ch;
p.energy_efficiency=r.energy_efficiency; p.seed=r.fragmentation_seed;
end
function s=summarize_case(o,p,id)
S=o.r3b.closure_summary; B=o.r3b.breakup_event_ledger; Ch=o.r3b.child_creation_ledger;
De=o.r3b.deactivation_ledger; TL=o.r3b.mass_energy_time_ledger;
s.case_id=string(id); s.case_role="existing_nominal_mean_area";
s.area_multiplier=p.area_multiplier; s.seed=p.seed;
s.max_total_fragments=p.max_total_fragments; s.tmax_s=p.tmax;
for f={'first_break_h_km','main_airburst_h_km','E_total_lost_kton','E_release_kton', ...
        'E_burst_raw_kton','E_burst_coupled_kton','dE_dh_peak_kton_per_km'}
    s.(f{1})=o.(f{1});
end
s.q_dyn_max_MPa=max(o.timeseries.q_dyn_max_Pa)/1e6;
s.fragment_count_max=max(o.timeseries.fragment_count);
s.breakup_event_count=height(B);s.created_child_count=height(Ch);
if isempty(Ch),s.max_child_id=1;inserted=false;else,s.max_child_id=max(Ch.child_id);inserted=all(Ch.inserted_into_active_pool);end
s.cap_hit_event_count=sum(B.cap_hit);s.cap_dropped_mass_kg=sum(B.cap_dropped_mass_kg);
s.insertion_dropped_mass_kg=S.cumulative_insertion_dropped_mass_kg;
s.ground_crossing_count=sum(De.reason=="post_step_ground_crossing");
s.terminal_completion_time_s=TL.time_s(end);
s.full_tree_realized=height(Ch)==85 && s.max_child_id==86 && inserted && height(B)==26;
s.terminal_complete=isempty(o.r3b.terminal_fragment_ledger)&&S.final_active_mass_kg<=1e-9;
s.mass_relative_residual=abs(S.final_mass_closure_residual_kg)/max(S.initial_mass_kg,eps);
s.energy_relative_residual=abs(S.final_energy_closure_residual_J)/max(S.initial_ke_J,eps);
s.outcome=string(o.outcome);
end
function s=as_publishable(o,r,p)
s.run_id=r.run_id;s.fragmentation_seed=r.fragmentation_seed;
s.publish_metric_source_horizon_s=p.tmax;s.publish_outcome=string(o.outcome);
names={'first_break_h_km','first_break_time_s','main_airburst_h_km','main_airburst_time_s', ...
    'E_total_lost_kton','E_release_kton','E_burst_raw_kton','E_burst_coupled_kton', ...
    'dE_dh_peak_kton_per_km','airburst_window_h_min_km','airburst_window_h_max_km'};
source={'first_break_h_km','first_break_time','main_airburst_h_km','main_airburst_t', ...
    'E_total_lost_kton','E_release_kton','E_burst_raw_kton','E_burst_coupled_kton', ...
    'dE_dh_peak_kton_per_km','airburst_window_h_min_km','airburst_window_h_max_km'};
for i=1:numel(names),s.(['publish_' names{i}])=o.(source{i});end
s.publish_q_dyn_max_MPa=max(o.timeseries.q_dyn_max_Pa)/1e6;
s.publish_fragment_count_max=max(o.r3b.mass_energy_time_ledger.active_count);
s.publish_breakup_event_count=height(o.r3b.breakup_event_ledger);
s.publish_created_child_count=height(o.r3b.child_creation_ledger);
s.publish_terminal_complete=isempty(o.r3b.terminal_fragment_ledger) && ...
    abs(o.r3b.closure_summary.final_active_mass_kg)<=1e-6;
end
function P=compute_pressure(PUB,C)
h=PUB.publish_main_airburst_h_km;E=PUB.publish_E_burst_coupled_kton;
% Match original two-operation cube root then square evaluation.
ss=(C.pressure_coefficient_kPa_km3_per_kt.*E./C.thresholds_kPa(:)').^(1/3);
rr=sqrt(max(ss.^2-h.^2,0));aa=pi*rr.^2;
P=table(PUB.run_id,PUB.fragmentation_seed,h,E,rr(:,1),rr(:,2),rr(:,3),rr(:,4), ...
    aa(:,1),aa(:,2),aa(:,3),aa(:,4),'VariableNames', ...
    {'run_id','fragmentation_seed','h_burst_km','E_c_kton','R_3_45_kPa_km', ...
    'R_1_00_kPa_km','R_0_50_kPa_km','R_0_20_kPa_km', ...
    'A_3_45_kPa_km2','A_1_00_kPa_km2','A_0_50_kPa_km2','A_0_20_kPa_km2'});
end
function D=prepare_population(file,C)
[A,R]=readgeoraster(file,'OutputType','double','StandardizeMissing',true);
if ndims(A)>2,A=A(:,:,1);end
assert(isequal(size(A),C.raster_size(:)'),'Replay:RasterSize','Unexpected raster size.');
valid=isfinite(A)&A>=0;
assert(nnz(valid)==C.valid_cells,'Replay:RasterMask','Unexpected valid-cell count.');
assert(nnz(valid&A>0)==C.positive_cells,'Replay:RasterMask','Unexpected positive-cell count.');
[ri,ci]=find(valid);[lat,lon]=intrinsicToGeographic(R,ci,ri);
pop=A(valid);RE=C.earth_radius_km;
dlat=deg2rad(lat)-deg2rad(C.scenario_center_lat_deg);dlon=deg2rad(lon-C.scenario_center_lon_deg);
a=sin(dlat/2).^2+cos(deg2rad(C.scenario_center_lat_deg)).*cos(deg2rad(lat)).*sin(dlon/2).^2;
a=min(1,max(0,a));distance=2*RE*asin(sqrt(a));
area=RE^2*deg2rad(R.CellExtentInLongitude).* ...
    (sin(deg2rad(min(90,lat+R.CellExtentInLatitude/2)))-sin(deg2rad(max(-90,lat-R.CellExtentInLatitude/2))));
D.total=sum(pop,'double');assert(abs(D.total-C.population_domain)<1e-3,'Replay:PopulationDomain','Population-domain mismatch.');
[D.distance,ord]=sort(distance,'ascend');
D.cumPop=cumsum(double(pop(ord)));D.cumArea=cumsum(double(area(ord)));
end
function E=compute_population(P,D,C)
rr=P{:,{'R_3_45_kPa_km','R_1_00_kPa_km','R_0_50_kPa_km','R_0_20_kPa_km'}};
radii=reshape(rr',[],1);ids=repelem(P.run_id,4);ids=ids(:);thresholds=repmat(C.thresholds_kPa(:),height(P),1);
pop=zeros(numel(radii),1);areas=pop;counts=pop;
for i=1:numel(radii)
    lo=0;hi=numel(D.distance);
    while lo<hi
        mid=floor((lo+hi+1)/2);
        if D.distance(mid)<=radii(i),lo=mid;else,hi=mid-1;end
    end
    counts(i)=lo;
    if lo>0,pop(i)=D.cumPop(lo);areas(i)=D.cumArea(lo);end
end
E=table(ids,thresholds,radii,pi*radii.^2,areas,pop,pop/D.total,counts, ...
    'VariableNames',{'run_id','threshold_kpa','radius_km','physical_footprint_area_km2', ...
    'population_grid_support_area_km2','exposed_population','exposure_fraction','included_cell_count'});
end
function validate_population_table_shapes(C)
% Synthetic 1/8/1000-row regression of the ACTUAL local compute_population.
% No scientific kernel, raster, archived result, or random stream is touched.
D.distance = [0.5;1.5;2.5;3.5;4.5;5.5];
D.cumPop = [10;30;60;100;150;210];
D.cumArea = [1;3;6;10;15;21];
D.total = 210;
for n = [1 8 1000]
    id = (1:n)';
    one = ones(n,1);
    P = table(id,2*one,3*one,4*one,5*one, ...
        'VariableNames',{'run_id','R_3_45_kPa_km','R_1_00_kPa_km', ...
        'R_0_50_kPa_km','R_0_20_kPa_km'});
    E = compute_population(P,D,C);
    idExpected = repelem(id,4);
    idExpected = idExpected(:);
    thresholdsExpected = repmat(C.thresholds_kPa(:),n,1);
    radiusExpected = repmat([2;3;4;5],n,1);
    popExpected = repmat([30;60;100;150],n,1);
    countExpected = repmat([2;3;4;5],n,1);
    assert(height(E)==4*n && width(E)==8, ...
        'Replay:PopulationTableShape','Wrong population table shape for n=%d.',n);
    assert(iscolumn(E.run_id) && isequal(E.run_id,idExpected), ...
        'Replay:PopulationRunIDShape','Wrong repeated run IDs for n=%d.',n);
    assert(isequal(E.threshold_kpa,thresholdsExpected) && ...
        isequal(E.radius_km,radiusExpected), ...
        'Replay:PopulationOrdering','Pressure threshold or radius order changed for n=%d.',n);
    assert(isequal(E.exposed_population,popExpected) && ...
        isequal(E.included_cell_count,countExpected), ...
        'Replay:PopulationCounts','Unexpected synthetic exposure results for n=%d.',n);
    assert(iscolumn(E.population_grid_support_area_km2) && ...
        iscolumn(E.physical_footprint_area_km2), ...
        'Replay:PopulationAreas','Area output orientation incorrect for n=%d.',n);
end
fprintf('native_population_shape_regression=PASS_1_8_1000\n');
end

function q=linear_quantile(x,p)
x=sort(x(:));pos=1+(numel(x)-1)*p;lo=floor(pos);hi=ceil(pos);
q=zeros(size(p));for i=1:numel(p),q(i)=x(lo(i))+(pos(i)-lo(i))*(x(hi(i))-x(lo(i)));end
end
function t=truth(x)
if isnumeric(x)||islogical(x),t=logical(x);else,t=ismember(lower(strtrim(string(x))),["1","true","yes"]);end
end
function h=sha256file(file)
fid=fopen(file,'rb');assert(fid>=0,'Replay:MissingFile','Missing file: %s',file);
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid)
    b=fread(fid,1048576,'*uint8');if ~isempty(b),md.update(b);end
end
h=upper(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function write_json(file,s)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0,'Cannot write result JSON.');
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(s,'PrettyPrint',true));
end

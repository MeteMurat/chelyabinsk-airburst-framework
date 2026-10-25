function out = meteor_single_run_r3b3_diagnostic(params)
% =========================================================
% TEK METEOR KOSUSU
% Cikti:
%   out.first_break_h_km
%   out.main_airburst_h_km
%   out.main_airburst_t
%   out.E_release_kton
%   out.final_mass
%   out.outcome
% =========================================================

if isfield(params,'seed') && ~isempty(params.seed)
    rng(params.seed);
else
    rng(42);
end

%% 1) SABITLER
rho_sea = 1.225;
H_scale = 7200;
g = 9.81;
R_earth = 6371000;

rho_meteor = params.rho_meteor;
V0 = params.V0;
entry_angle_deg = params.entry_angle_deg;
theta0 = deg2rad(entry_angle_deg);

D0 = params.D0;
r0 = D0/2;
m0 = rho_meteor*(4/3)*pi*r0^3;

h0 = params.h0;
x0 = 0;

Yield_Stress_0 = params.Yield_Stress_0;
min_strength_fraction = params.min_strength_fraction;
weakening_mass_coeff = params.weakening_mass_coeff;

Cd = params.Cd;
Ch = params.Ch;
Qabl = params.Qabl;

max_generations = params.max_generations;
max_total_fragments = params.max_total_fragments;
base_fragment_count = params.base_fragment_count;
dispersion_speed = params.dispersion_speed;
cloud_area_gain = params.cloud_area_gain;

if isfield(params,'area_multiplier') && ~isempty(params.area_multiplier)
    area_multiplier = params.area_multiplier;
else
    area_multiplier = 1.0;
end

if ~isfinite(area_multiplier) || area_multiplier <= 0
    error('params.area_multiplier must be a positive finite scalar.');
end

dt = params.dt;
tmax = params.tmax;
Nmax = ceil(tmax/dt);

%% 2) FRAGMANLAR
fragment_template = struct( ...
    'active', true, ...
    'gen', 0, ...
    'm', 0, ...
    'r', 0, ...
    'A', 0, ...
    'h', 0, ...
    'x', 0, ...
    'v', 0, ...
    'theta', 0, ...
    'Y0', 0, ...
    'Yeff', 0, ...
    'broken', false, ...
    'id', 0);

fragments = repmat(fragment_template, max_total_fragments, 1);

fragments(1).active = true;
fragments(1).gen    = 0;
fragments(1).m      = m0;
fragments(1).r      = r0;
fragments(1).A      = area_multiplier * pi*r0^2;
fragments(1).h      = h0;
fragments(1).x      = x0;
fragments(1).v      = V0;
fragments(1).theta  = theta0;
fragments(1).Y0     = Yield_Stress_0;
fragments(1).Yeff   = Yield_Stress_0;
fragments(1).broken = false;
fragments(1).id     = 1;


next_fragment_id = 2;

%% R3B INSTRUMENTATION STATE (diagnostic only; no random draws)
birth_mass_by_id = nan(max_total_fragments,1);
birth_mass_by_id(1) = m0;
child_inserted_by_id = false(max_total_fragments,1);
child_inserted_by_id(1) = true;
deactivation_recorded_by_id = false(max_total_fragments,1);

r3b_event_counter = 0;
r3b_child_counter = 0;
r3b_deactivation_counter = 0;

cumulative_ablation_mass_kg = 0;
cumulative_cap_dropped_mass_kg = 0;
cumulative_insertion_dropped_mass_kg = 0;
cumulative_deactivated_mass_kg = 0;

cumulative_propagation_ke_delta_J = 0;
cumulative_breakup_ke_delta_J = 0;
cumulative_dispersion_ke_delta_J = 0;
cumulative_cap_reference_ke_J = 0;
cumulative_insertion_dropped_ke_J = 0;
cumulative_deactivated_ke_J = 0;
cumulative_drag_work_estimate_J = 0;
cumulative_gravity_work_estimate_J = 0;

% Breakup-event ledger.
b_event_id = zeros(0,1);
b_it = zeros(0,1);
b_time_s = zeros(0,1);
b_altitude_m = zeros(0,1);
b_q_Pa = zeros(0,1);
b_Y0_Pa = zeros(0,1);
b_Yeff_Pa = zeros(0,1);
b_strength_factor = zeros(0,1);
b_strength_clamped = false(0,1);
b_generation = zeros(0,1);
b_parent_slot = zeros(0,1);
b_parent_id = zeros(0,1);
b_parent_mass_kg = zeros(0,1);
b_parent_birth_mass_kg = zeros(0,1);
b_parent_mass_over_global_m0 = zeros(0,1);
b_parent_mass_over_birth_mass = zeros(0,1);
b_intended_nfrag = zeros(0,1);
b_created_count = zeros(0,1);
b_planned_child_mass_kg = zeros(0,1);
b_created_child_mass_kg = zeros(0,1);
b_cap_dropped_mass_kg = zeros(0,1);
b_cap_hit = false(0,1);
b_parent_ke_J = zeros(0,1);
b_created_child_ke_J = zeros(0,1);
b_created_child_ke_at_parent_speed_J = zeros(0,1);
b_dispersion_ke_delta_J = zeros(0,1);
b_cap_reference_ke_J = zeros(0,1);
b_total_breakup_ke_delta_J = zeros(0,1);
b_active_count_before_break = zeros(0,1);
b_next_fragment_id_before = zeros(0,1);
b_next_fragment_id_after = zeros(0,1);

% Child-creation ledger.
c_event_id = zeros(0,1);
c_child_order = zeros(0,1);
c_parent_id = zeros(0,1);
c_child_id = zeros(0,1);
c_generation = zeros(0,1);
c_weight = zeros(0,1);
c_mass_kg = zeros(0,1);
c_birth_strength_Pa = zeros(0,1);
c_dv_m_s = zeros(0,1);
c_dtheta_rad = zeros(0,1);
c_parent_v_m_s = zeros(0,1);
c_child_v_m_s = zeros(0,1);
c_parent_ke_share_J = zeros(0,1);
c_child_ke_J = zeros(0,1);

% Deactivation ledger.
d_event_id = zeros(0,1);
d_it = zeros(0,1);
d_time_s = zeros(0,1);
d_fragment_id = zeros(0,1);
d_generation = zeros(0,1);
d_reason = strings(0,1);
d_mass_kg = zeros(0,1);
d_ke_J = zeros(0,1);
d_altitude_m = zeros(0,1);
d_velocity_m_s = zeros(0,1);

% Diagnostic time ledger.
diag_time_s = zeros(0,1);
diag_it = zeros(0,1);
diag_active_mass_kg = zeros(0,1);
diag_active_ke_J = zeros(0,1);
diag_active_count = zeros(0,1);
diag_min_active_mass_kg = zeros(0,1);
diag_next_fragment_id = zeros(0,1);
diag_breakup_event_count = zeros(0,1);
diag_cum_ablation_mass_kg = zeros(0,1);
diag_cum_cap_dropped_mass_kg = zeros(0,1);
diag_cum_insertion_dropped_mass_kg = zeros(0,1);
diag_cum_deactivated_mass_kg = zeros(0,1);
diag_mass_closure_residual_kg = zeros(0,1);
diag_cum_propagation_ke_delta_J = zeros(0,1);
diag_cum_breakup_ke_delta_J = zeros(0,1);
diag_cum_dispersion_ke_delta_J = zeros(0,1);
diag_cum_cap_reference_ke_J = zeros(0,1);
diag_cum_insertion_dropped_ke_J = zeros(0,1);
diag_energy_closure_residual_J = zeros(0,1);

%% 3) KAYIT
time_arr       = nan(Nmax,1);
total_ke_arr   = nan(Nmax,1);
mean_h_arr     = nan(Nmax,1);
dE_dh_arr      = nan(Nmax,1);
total_mass_arr = nan(Nmax,1);
mean_v_arr     = nan(Nmax,1);
mean_x_arr     = nan(Nmax,1);
q_dyn_max_arr  = nan(Nmax,1);
active_count_arr = nan(Nmax,1);

first_break_time = NaN;
first_break_h    = NaN;
first_break_done = false;

prev_total_ke = 0;
prev_mean_h = h0;

E0_system = 0.5 * m0 * V0^2;
final_ke_system = NaN;
final_mass_system = NaN;

%% 4) ANA DONGU
for it = 1:Nmax
    t = (it-1)*dt;

    total_mass = 0;
    total_ke = 0;
    total_h_mass_weighted = 0;
    total_v_mass_weighted = 0;
    total_x_mass_weighted = 0;
    q_dyn_max = 0;
    active_count = 0;

    new_fragments_to_add = [];

    for k = 1:max_total_fragments
        if ~fragments(k).active
            continue;
        end

        m = fragments(k).m;
        r = fragments(k).r;
        A = fragments(k).A;
        h = fragments(k).h;
        x = fragments(k).x;
        v = fragments(k).v;
        theta = fragments(k).theta;
        Y0 = fragments(k).Y0;
        gen = fragments(k).gen;

        if m <= 0 || r <= 0 || h <= 0 || v <= 0
            frag_id_now = fragments(k).id;
            if frag_id_now >= 1 && frag_id_now <= max_total_fragments && ...
                    ~deactivation_recorded_by_id(frag_id_now)
                deleted_mass_now = max(m,0);
                deleted_ke_now = 0.5*deleted_mass_now*max(v,0)^2;
                if deleted_mass_now > 0
                    cumulative_deactivated_mass_kg = cumulative_deactivated_mass_kg + deleted_mass_now;
                    cumulative_deactivated_ke_J = cumulative_deactivated_ke_J + deleted_ke_now;
                end
                r3b_deactivation_counter = r3b_deactivation_counter + 1;
                d_event_id(end+1,1) = r3b_deactivation_counter;
                d_it(end+1,1) = it;
                d_time_s(end+1,1) = t;
                d_fragment_id(end+1,1) = frag_id_now;
                d_generation(end+1,1) = gen;
                d_reason(end+1,1) = "pre_step_invalid_state";
                d_mass_kg(end+1,1) = deleted_mass_now;
                d_ke_J(end+1,1) = deleted_ke_now;
                d_altitude_m(end+1,1) = h;
                d_velocity_m_s(end+1,1) = v;
                deactivation_recorded_by_id(frag_id_now) = true;
            end
            fragments(k).active = false;
            continue;
        end

        rho_air = rho_sea * exp(-h / H_scale);
        q = 0.5 * rho_air * v^2;
        q_dyn_max = max(q_dyn_max, q);

        mass_loss_fraction = max(0, min(1, 1 - m/m0));
        strength_factor = 1 - weakening_mass_coeff * mass_loss_fraction;
        strength_factor = max(min_strength_fraction, strength_factor);
        Yeff = Y0 * strength_factor;
        fragments(k).Yeff = Yeff;

        do_break = false;
        if (q >= Yeff) && (~fragments(k).broken) && (gen < max_generations)
            do_break = true;
        end

        if do_break
            fragments(k).broken = true;

            if ~first_break_done
                first_break_done = true;
                first_break_time = t;
                first_break_h = h;
            end

            nfrag = max(2, round(base_fragment_count - gen));
            weights = rand(1, nfrag).^1.5;
            weights = weights / sum(weights);

            % R3B event-local diagnostics. These statements do not draw random numbers
            % and do not alter the original fragmentation state transition.
            r3b_event_counter = r3b_event_counter + 1;
            event_id_now = r3b_event_counter;
            parent_id_now = fragments(k).id;
            parent_birth_mass_now = NaN;
            if parent_id_now >= 1 && parent_id_now <= max_total_fragments
                parent_birth_mass_now = birth_mass_by_id(parent_id_now);
            end
            if ~isfinite(parent_birth_mass_now) || parent_birth_mass_now <= 0
                parent_birth_mass_now = m;
            end
            parent_ke_now = 0.5*m*v^2;
            created_count_now = 0;
            created_mass_now = 0;
            created_ke_now = 0;
            created_ke_parent_speed_now = 0;
            next_fragment_id_before_now = next_fragment_id;
            active_count_before_break_now = sum([fragments.active]);

            fragments(k).active = false;

            for j = 1:nfrag
                if next_fragment_id > max_total_fragments
                    break;
                end

                m_child = m * weights(j);
                r_child = ((3*(m_child/rho_meteor))/(4*pi))^(1/3);
                A_geom = area_multiplier * pi*r_child^2;
                A_child = A_geom * (1 + cloud_area_gain*(gen+1));
                Y_child = Y0 * (0.85 + 0.30*rand);

                dv = dispersion_speed * (2*rand - 1) * 0.25;
                dtheta = deg2rad((2*rand - 1) * 1.5);

                frag = fragment_template;
                frag.active = true;
                frag.gen = gen + 1;
                frag.m = m_child;
                frag.r = r_child;
                frag.A = A_child;
                frag.h = h;
                frag.x = x;
                frag.v = max(1, v + dv);
                frag.theta = max(deg2rad(5), min(deg2rad(85), theta + dtheta));
                frag.Y0 = Y_child;
                frag.Yeff = Y_child;
                frag.broken = false;
                frag.id = next_fragment_id;

                child_id_now = next_fragment_id;
                birth_mass_by_id(child_id_now) = m_child;

                created_count_now = created_count_now + 1;
                created_mass_now = created_mass_now + m_child;
                created_ke_now = created_ke_now + 0.5*m_child*frag.v^2;
                created_ke_parent_speed_now = created_ke_parent_speed_now + 0.5*m_child*v^2;

                r3b_child_counter = r3b_child_counter + 1;
                c_event_id(end+1,1) = event_id_now;
                c_child_order(end+1,1) = j;
                c_parent_id(end+1,1) = parent_id_now;
                c_child_id(end+1,1) = child_id_now;
                c_generation(end+1,1) = gen + 1;
                c_weight(end+1,1) = weights(j);
                c_mass_kg(end+1,1) = m_child;
                c_birth_strength_Pa(end+1,1) = Y_child;
                c_dv_m_s(end+1,1) = dv;
                c_dtheta_rad(end+1,1) = dtheta;
                c_parent_v_m_s(end+1,1) = v;
                c_child_v_m_s(end+1,1) = frag.v;
                c_parent_ke_share_J(end+1,1) = 0.5*m_child*v^2;
                c_child_ke_J(end+1,1) = 0.5*m_child*frag.v^2;

                new_fragments_to_add = [new_fragments_to_add; frag]; %#ok<AGROW>
                next_fragment_id = next_fragment_id + 1;
            end

            cap_hit_now = created_count_now < nfrag;
            if cap_hit_now
                cap_dropped_mass_now = sum(m * weights((created_count_now+1):end));
            else
                cap_dropped_mass_now = 0;
            end
            dispersion_ke_delta_now = created_ke_now - created_ke_parent_speed_now;
            cap_reference_ke_now = 0.5*cap_dropped_mass_now*v^2;
            total_breakup_ke_delta_now = created_ke_now - parent_ke_now;

            cumulative_cap_dropped_mass_kg = cumulative_cap_dropped_mass_kg + cap_dropped_mass_now;
            cumulative_breakup_ke_delta_J = cumulative_breakup_ke_delta_J + total_breakup_ke_delta_now;
            cumulative_dispersion_ke_delta_J = cumulative_dispersion_ke_delta_J + dispersion_ke_delta_now;
            cumulative_cap_reference_ke_J = cumulative_cap_reference_ke_J + cap_reference_ke_now;

            b_event_id(end+1,1) = event_id_now;
            b_it(end+1,1) = it;
            b_time_s(end+1,1) = t;
            b_altitude_m(end+1,1) = h;
            b_q_Pa(end+1,1) = q;
            b_Y0_Pa(end+1,1) = Y0;
            b_Yeff_Pa(end+1,1) = Yeff;
            b_strength_factor(end+1,1) = strength_factor;
            b_strength_clamped(end+1,1) = strength_factor <= (min_strength_fraction + 10*eps);
            b_generation(end+1,1) = gen;
            b_parent_slot(end+1,1) = k;
            b_parent_id(end+1,1) = parent_id_now;
            b_parent_mass_kg(end+1,1) = m;
            b_parent_birth_mass_kg(end+1,1) = parent_birth_mass_now;
            b_parent_mass_over_global_m0(end+1,1) = m/max(m0,eps);
            b_parent_mass_over_birth_mass(end+1,1) = m/max(parent_birth_mass_now,eps);
            b_intended_nfrag(end+1,1) = nfrag;
            b_created_count(end+1,1) = created_count_now;
            b_planned_child_mass_kg(end+1,1) = m;
            b_created_child_mass_kg(end+1,1) = created_mass_now;
            b_cap_dropped_mass_kg(end+1,1) = cap_dropped_mass_now;
            b_cap_hit(end+1,1) = cap_hit_now;
            b_parent_ke_J(end+1,1) = parent_ke_now;
            b_created_child_ke_J(end+1,1) = created_ke_now;
            b_created_child_ke_at_parent_speed_J(end+1,1) = created_ke_parent_speed_now;
            b_dispersion_ke_delta_J(end+1,1) = dispersion_ke_delta_now;
            b_cap_reference_ke_J(end+1,1) = cap_reference_ke_now;
            b_total_breakup_ke_delta_J(end+1,1) = total_breakup_ke_delta_now;
            b_active_count_before_break(end+1,1) = active_count_before_break_now;
            b_next_fragment_id_before(end+1,1) = next_fragment_id_before_now;
            b_next_fragment_id_after(end+1,1) = next_fragment_id;

            continue;
        end

        Drag = 0.5 * Cd * rho_air * A * v^2;
        a_drag = Drag / m;

        dv_dt = g*sin(theta) - a_drag;
        dm_dt = -(Ch * rho_air * A * v^3) / (2*Qabl);

        dh_dt = -v * sin(theta);
        dx_dt =  v * cos(theta);
        dtheta_dt = (g*cos(theta)/max(v,1)) - (v*cos(theta)/R_earth);

        v_new = max(1, v + dv_dt*dt);
        m_new = max(0, m + dm_dt*dt);
        h_new = h + dh_dt*dt;
        x_new = x + dx_dt*dt;
        theta_new = theta + dtheta_dt*dt;

        % R3B exact discrete propagation ledger.
        ablation_loss_now = max(0, m - m_new);
        cumulative_ablation_mass_kg = cumulative_ablation_mass_kg + ablation_loss_now;
        cumulative_drag_work_estimate_J = cumulative_drag_work_estimate_J + Drag*v*dt;
        cumulative_gravity_work_estimate_J = cumulative_gravity_work_estimate_J + m*g*sin(theta)*v*dt;

        pre_step_ke_now = 0.5*m*v^2;
        post_step_ke_candidate_now = 0.5*m_new*v_new^2;
        post_step_ke_effective_now = post_step_ke_candidate_now;

        if m_new <= 0 || h_new <= 0 || v_new <= 0
            post_step_ke_effective_now = 0;
            frag_id_now = fragments(k).id;
            if frag_id_now >= 1 && frag_id_now <= max_total_fragments && ...
                    ~deactivation_recorded_by_id(frag_id_now)
                deleted_mass_now = max(m_new,0);
                deleted_ke_now = post_step_ke_candidate_now;
                if deleted_mass_now > 0
                    cumulative_deactivated_mass_kg = cumulative_deactivated_mass_kg + deleted_mass_now;
                    cumulative_deactivated_ke_J = cumulative_deactivated_ke_J + deleted_ke_now;
                end
                r3b_deactivation_counter = r3b_deactivation_counter + 1;
                d_event_id(end+1,1) = r3b_deactivation_counter;
                d_it(end+1,1) = it;
                d_time_s(end+1,1) = t;
                d_fragment_id(end+1,1) = frag_id_now;
                d_generation(end+1,1) = gen;
                if m_new <= 0
                    d_reason(end+1,1) = "post_step_mass_depletion";
                elseif h_new <= 0
                    d_reason(end+1,1) = "post_step_ground_crossing";
                else
                    d_reason(end+1,1) = "post_step_velocity_invalid";
                end
                d_mass_kg(end+1,1) = deleted_mass_now;
                d_ke_J(end+1,1) = deleted_ke_now;
                d_altitude_m(end+1,1) = h_new;
                d_velocity_m_s(end+1,1) = v_new;
                deactivation_recorded_by_id(frag_id_now) = true;
            end
        end

        cumulative_propagation_ke_delta_J = cumulative_propagation_ke_delta_J + ...
            (post_step_ke_effective_now - pre_step_ke_now);

        if m_new <= 0 || h_new <= 0 || v_new <= 0
            fragments(k).active = false;
            continue;
        end

        V_new = m_new / rho_meteor;
        r_new = ((3*V_new)/(4*pi))^(1/3);
        A_new = area_multiplier * pi*r_new^2 * (1 + cloud_area_gain*gen);

        fragments(k).m = m_new;
        fragments(k).r = r_new;
        fragments(k).A = A_new;
        fragments(k).h = h_new;
        fragments(k).x = x_new;
        fragments(k).v = v_new;
        fragments(k).theta = max(deg2rad(3), min(deg2rad(85), theta_new));

        total_mass = total_mass + m_new;
        total_ke = total_ke + 0.5*m_new*v_new^2;
        total_h_mass_weighted = total_h_mass_weighted + m_new*h_new;
        total_v_mass_weighted = total_v_mass_weighted + m_new*v_new;
        total_x_mass_weighted = total_x_mass_weighted + m_new*x_new;
        active_count = active_count + 1;
    end

    if ~isempty(new_fragments_to_add)
        for n = 1:length(new_fragments_to_add)
            idx_empty = find(~[fragments.active], 1, 'first');
            if isempty(idx_empty)
                for nn = n:length(new_fragments_to_add)
                    dropped_frag_now = new_fragments_to_add(nn);
                    cumulative_insertion_dropped_mass_kg = cumulative_insertion_dropped_mass_kg + dropped_frag_now.m;
                    cumulative_insertion_dropped_ke_J = cumulative_insertion_dropped_ke_J + ...
                        0.5*dropped_frag_now.m*dropped_frag_now.v^2;
                end
                break;
            end
            fragments(idx_empty) = new_fragments_to_add(n);
            inserted_id_now = new_fragments_to_add(n).id;
            if inserted_id_now >= 1 && inserted_id_now <= max_total_fragments
                child_inserted_by_id(inserted_id_now) = true;
            end
        end
    end

    % FORCE_RECOUNT_ACTIVE_FRAGMENTS_BEFORE_TERMINATION
    % Recompute aggregate state from the currently active fragment list.
    % This is required when a parent fragment breaks and child fragments are
    % inserted after the propagation loop of the same time step.
    active_idx_now = find([fragments.active]);

    if ~isempty(active_idx_now)
        total_mass = 0;
        total_ke = 0;
        total_h_mass_weighted = 0;
        total_v_mass_weighted = 0;
        total_x_mass_weighted = 0;
        q_dyn_max = 0;
        active_count = 0;

        for kk = active_idx_now
            m_now = fragments(kk).m;
            h_now = fragments(kk).h;
            x_now = fragments(kk).x;
            v_now = fragments(kk).v;

            if m_now <= 0 || h_now <= 0 || v_now <= 0
                frag_id_now = fragments(kk).id;
                if frag_id_now >= 1 && frag_id_now <= max_total_fragments && ...
                        ~deactivation_recorded_by_id(frag_id_now)
                    deleted_mass_now = max(m_now,0);
                    deleted_ke_now = 0.5*deleted_mass_now*max(v_now,0)^2;
                    if deleted_mass_now > 0
                        cumulative_deactivated_mass_kg = cumulative_deactivated_mass_kg + deleted_mass_now;
                        cumulative_deactivated_ke_J = cumulative_deactivated_ke_J + deleted_ke_now;
                    end
                    r3b_deactivation_counter = r3b_deactivation_counter + 1;
                    d_event_id(end+1,1) = r3b_deactivation_counter;
                    d_it(end+1,1) = it;
                    d_time_s(end+1,1) = t;
                    d_fragment_id(end+1,1) = frag_id_now;
                    d_generation(end+1,1) = fragments(kk).gen;
                    d_reason(end+1,1) = "recount_invalid_state";
                    d_mass_kg(end+1,1) = deleted_mass_now;
                    d_ke_J(end+1,1) = deleted_ke_now;
                    d_altitude_m(end+1,1) = h_now;
                    d_velocity_m_s(end+1,1) = v_now;
                    deactivation_recorded_by_id(frag_id_now) = true;
                end
                fragments(kk).active = false;
                continue;
            end

            rho_air_now = rho_sea * exp(-h_now / H_scale);
            q_now = 0.5 * rho_air_now * v_now^2;

            total_mass = total_mass + m_now;
            total_ke = total_ke + 0.5*m_now*v_now^2;
            total_h_mass_weighted = total_h_mass_weighted + m_now*h_now;
            total_v_mass_weighted = total_v_mass_weighted + m_now*v_now;
            total_x_mass_weighted = total_x_mass_weighted + m_now*x_now;
            q_dyn_max = max(q_dyn_max, q_now);
            active_count = active_count + 1;
        end
    end
    % R3B time-ledger snapshot after insertion and authoritative recount.
    if isempty(active_idx_now)
        min_active_mass_now = NaN;
    else
        active_idx_valid_now = find([fragments.active]);
        if isempty(active_idx_valid_now)
            min_active_mass_now = NaN;
        else
            min_active_mass_now = min([fragments(active_idx_valid_now).m]);
        end
    end

    mass_closure_residual_now = m0 - (total_mass + cumulative_ablation_mass_kg + ...
        cumulative_cap_dropped_mass_kg + cumulative_insertion_dropped_mass_kg + ...
        cumulative_deactivated_mass_kg);
    energy_closure_residual_now = E0_system + cumulative_propagation_ke_delta_J + ...
        cumulative_breakup_ke_delta_J - cumulative_insertion_dropped_ke_J - total_ke;

    diag_it(end+1,1) = it;
    diag_time_s(end+1,1) = t;
    diag_active_mass_kg(end+1,1) = total_mass;
    diag_active_ke_J(end+1,1) = total_ke;
    diag_active_count(end+1,1) = active_count;
    diag_min_active_mass_kg(end+1,1) = min_active_mass_now;
    diag_next_fragment_id(end+1,1) = next_fragment_id;
    diag_breakup_event_count(end+1,1) = r3b_event_counter;
    diag_cum_ablation_mass_kg(end+1,1) = cumulative_ablation_mass_kg;
    diag_cum_cap_dropped_mass_kg(end+1,1) = cumulative_cap_dropped_mass_kg;
    diag_cum_insertion_dropped_mass_kg(end+1,1) = cumulative_insertion_dropped_mass_kg;
    diag_cum_deactivated_mass_kg(end+1,1) = cumulative_deactivated_mass_kg;
    diag_mass_closure_residual_kg(end+1,1) = mass_closure_residual_now;
    diag_cum_propagation_ke_delta_J(end+1,1) = cumulative_propagation_ke_delta_J;
    diag_cum_breakup_ke_delta_J(end+1,1) = cumulative_breakup_ke_delta_J;
    diag_cum_dispersion_ke_delta_J(end+1,1) = cumulative_dispersion_ke_delta_J;
    diag_cum_cap_reference_ke_J(end+1,1) = cumulative_cap_reference_ke_J;
    diag_cum_insertion_dropped_ke_J(end+1,1) = cumulative_insertion_dropped_ke_J;
    diag_energy_closure_residual_J(end+1,1) = energy_closure_residual_now;

    if active_count == 0
        final_ke_system = 0;
        final_mass_system = 0;

        time_arr = time_arr(1:it-1);
        total_ke_arr = total_ke_arr(1:it-1);
        mean_h_arr = mean_h_arr(1:it-1);
        dE_dh_arr = dE_dh_arr(1:it-1);
        break;
    end

    mean_h = total_h_mass_weighted / max(total_mass, eps);

    mean_v = total_v_mass_weighted / max(total_mass, eps);

    mean_x = total_x_mass_weighted / max(total_mass, eps);mean_v = total_v_mass_weighted / max(total_mass, eps);
    mean_x = total_x_mass_weighted / max(total_mass, eps);

    if it == 1
        dE_dh = 0;
    else
        dE = max(0, prev_total_ke - total_ke);
        dh = abs(prev_mean_h - mean_h);
        dE_dh = dE / max(dh, 1e-6);
    end

    time_arr(it) = t;
    total_ke_arr(it) = total_ke;
    mean_h_arr(it) = mean_h;
    dE_dh_arr(it) = dE_dh;
    total_mass_arr(it) = total_mass;
    mean_v_arr(it) = mean_v;
    mean_x_arr(it) = mean_x;
    q_dyn_max_arr(it) = q_dyn_max;
    active_count_arr(it) = active_count;total_mass_arr(it) = total_mass;
    mean_v_arr(it) = mean_v;
    mean_x_arr(it) = mean_x;
    q_dyn_max_arr(it) = q_dyn_max;
    active_count_arr(it) = active_count;

    prev_total_ke = total_ke;
    prev_mean_h = mean_h;
    final_ke_system = total_ke;
    final_mass_system = total_mass;
end

valid = ~isnan(time_arr);
time_arr = time_arr(valid);
total_ke_arr = total_ke_arr(valid);
mean_h_arr = mean_h_arr(valid);
dE_dh_arr = dE_dh_arr(valid);

%% 5) CIKTILAR

% dE_dh_arr is a raw stepwise energy-deposition proxy with units of J/m.
% For manuscript-grade airburst identification, use an altitude-binned
% deposition profile. This avoids selecting a single time-step spike as the
% whole burst window after discrete fragmentation.
if isfield(params,'energy_bin_width_m') && ~isempty(params.energy_bin_width_m)
    energy_bin_width_m = params.energy_bin_width_m;
else
    energy_bin_width_m = 500; % 0.5 km default bin width
end

h_prof = mean_h_arr(:);
ke_prof = total_ke_arr(:);
t_prof = time_arr(:);

valid_prof = isfinite(h_prof) & isfinite(ke_prof) & isfinite(t_prof);

h_prof = h_prof(valid_prof);
ke_prof = ke_prof(valid_prof);
t_prof = t_prof(valid_prof);

if numel(h_prof) < 2
    idx_airburst = [];
    dE_dh_peak_J_per_m = NaN;
    dE_dh_peak_kton_per_km = NaN;
    main_airburst_h_km = NaN;
    main_airburst_t = NaN;
    burst_idx = [];
    airburst_window_h_min_km = NaN;
    airburst_window_h_max_km = NaN;
    E_burst_raw_J = NaN;

    dep_h_center_m = [];
    dep_t_mid_s = [];
    dep_dE_bin_J = [];
    dep_dE_dh_J_per_m = [];
else
    h_mid = 0.5*(h_prof(1:end-1) + h_prof(2:end));
    t_mid = 0.5*(t_prof(1:end-1) + t_prof(2:end));
    dE_step_J = max(0, ke_prof(1:end-1) - ke_prof(2:end));

    h_min_edge = floor(min(h_mid)/energy_bin_width_m)*energy_bin_width_m;
    h_max_edge = ceil(max(h_mid)/energy_bin_width_m)*energy_bin_width_m;

    if h_max_edge <= h_min_edge
        h_max_edge = h_min_edge + energy_bin_width_m;
    end

    edges = h_min_edge:energy_bin_width_m:h_max_edge;

    if numel(edges) < 2
        edges = [h_min_edge, h_min_edge + energy_bin_width_m];
    end

    [~,~,bin_idx] = histcounts(h_mid, edges);
    nbins = numel(edges) - 1;

    valid_bin = bin_idx > 0 & bin_idx <= nbins & isfinite(dE_step_J);

    dep_dE_bin_J = accumarray(bin_idx(valid_bin), dE_step_J(valid_bin), [nbins 1], @sum, 0);
    dep_h_center_m = 0.5*(edges(1:end-1) + edges(2:end))';
    dep_dE_dh_J_per_m = dep_dE_bin_J / energy_bin_width_m;

    % Energy-weighted time for each altitude bin; fall back to NaN if no energy.
    weighted_t_num = accumarray(bin_idx(valid_bin), t_mid(valid_bin).*dE_step_J(valid_bin), [nbins 1], @sum, 0);
    dep_t_mid_s = nan(nbins,1);
    nz = dep_dE_bin_J > 0;
    dep_t_mid_s(nz) = weighted_t_num(nz) ./ dep_dE_bin_J(nz);

    if isempty(dep_dE_dh_J_per_m) || max(dep_dE_dh_J_per_m) <= 0
        idx_airburst = [];
        dE_dh_peak_J_per_m = NaN;
        dE_dh_peak_kton_per_km = NaN;
        main_airburst_h_km = NaN;
        main_airburst_t = NaN;
        burst_idx = [];
        airburst_window_h_min_km = NaN;
        airburst_window_h_max_km = NaN;
        E_burst_raw_J = NaN;
    else
        [dE_dh_peak_J_per_m, idx_airburst] = max(dep_dE_dh_J_per_m);
        dE_dh_peak_kton_per_km = dE_dh_peak_J_per_m * 1000 / 4.184e12;

        main_airburst_h_km = dep_h_center_m(idx_airburst)/1000;
        main_airburst_t = dep_t_mid_s(idx_airburst);

        half_peak = 0.5*dE_dh_peak_J_per_m;

        left_idx = idx_airburst;
        while left_idx > 1 && dep_dE_dh_J_per_m(left_idx-1) >= half_peak
            left_idx = left_idx - 1;
        end

        right_idx = idx_airburst;
        while right_idx < numel(dep_dE_dh_J_per_m) && dep_dE_dh_J_per_m(right_idx+1) >= half_peak
            right_idx = right_idx + 1;
        end

        burst_idx = left_idx:right_idx;

        airburst_window_h_min_km = (min(dep_h_center_m(burst_idx)) - 0.5*energy_bin_width_m)/1000;
        airburst_window_h_max_km = (max(dep_h_center_m(burst_idx)) + 0.5*energy_bin_width_m)/1000;

        E_burst_raw_J = sum(dep_dE_bin_J(burst_idx));
    end
end

if isnan(final_ke_system)
    final_ke_system = 0;
end
if isnan(final_mass_system)
    final_mass_system = 0;
end

energy_efficiency = params.energy_efficiency;

E_total_lost_J = max(0, E0_system - final_ke_system);
E_total_lost_kton = E_total_lost_J / 4.184e12;

% Backward-compatible total effective released energy.
E_release_J = energy_efficiency * E_total_lost_J;
E_release_kton = E_release_J / 4.184e12;

% Local airburst-window energy.
if isfinite(E_burst_raw_J)
    E_burst_raw_J = min(max(E_burst_raw_J,0), E_total_lost_J);
end

E_burst_raw_kton = E_burst_raw_J / 4.184e12;
E_burst_coupled_J = energy_efficiency * E_burst_raw_J;
E_burst_coupled_kton = E_burst_coupled_J / 4.184e12;

% Sonuc siniflandirma
% Event-based manuscript label:
% A physically resolved airburst has priority over terminal mass depletion.
valid_airburst = isfinite(main_airburst_h_km) && main_airburst_h_km >= 5 && ...
                 isfinite(E_burst_coupled_kton) && E_burst_coupled_kton > 0 && ...
                 isfinite(dE_dh_peak_kton_per_km) && dE_dh_peak_kton_per_km > 0;

if valid_airburst
    outcome = "Airburst";
elseif final_mass_system <= 0.01*m0
    outcome = "Burn-up";
else
    outcome = "Impact";
end

out.first_break_h_km = first_break_h / 1000;
out.first_break_time = first_break_time;
out.main_airburst_h_km = main_airburst_h_km;
out.main_airburst_t = main_airburst_t;
out.E_release_kton = E_release_kton;
out.E_release_J = E_release_J;

out.E_total_lost_J = E_total_lost_J;
out.E_total_lost_kton = E_total_lost_kton;

out.dE_dh_peak_J_per_m = dE_dh_peak_J_per_m;
out.dE_dh_peak_kton_per_km = dE_dh_peak_kton_per_km;

out.E_burst_raw_J = E_burst_raw_J;
out.E_burst_raw_kton = E_burst_raw_kton;
out.E_burst_coupled_J = E_burst_coupled_J;
out.E_burst_coupled_kton = E_burst_coupled_kton;

out.airburst_window_h_min_km = airburst_window_h_min_km;
out.airburst_window_h_max_km = airburst_window_h_max_km;
out.final_mass = final_mass_system;
out.area_multiplier = area_multiplier;
out.deposition_profile = table( ...
    dep_h_center_m/1000, ...
    dep_t_mid_s, ...
    dep_dE_bin_J/4.184e12, ...
    dep_dE_dh_J_per_m, ...
    dep_dE_dh_J_per_m*1000/4.184e12, ...
    'VariableNames', { ...
    'h_center_km', ...
    't_mid_s', ...
    'dE_bin_kton', ...
    'dE_dh_J_per_m', ...
    'dE_dh_kton_per_km'});

out.outcome = outcome;

% Time-series output for nominal-case manuscript figures.
% Force all stored histories to a common length before building the table.
ts_lengths = [ ...
    numel(time_arr), ...
    numel(mean_h_arr), ...
    numel(mean_x_arr), ...
    numel(mean_v_arr), ...
    numel(total_mass_arr), ...
    numel(q_dyn_max_arr), ...
    numel(active_count_arr), ...
    numel(total_ke_arr), ...
    numel(dE_dh_arr)];

n_ts = min(ts_lengths);

time_arr = time_arr(1:n_ts);
mean_h_arr = mean_h_arr(1:n_ts);
mean_x_arr = mean_x_arr(1:n_ts);
mean_v_arr = mean_v_arr(1:n_ts);
total_mass_arr = total_mass_arr(1:n_ts);
q_dyn_max_arr = q_dyn_max_arr(1:n_ts);
active_count_arr = active_count_arr(1:n_ts);
total_ke_arr = total_ke_arr(1:n_ts);
dE_dh_arr = dE_dh_arr(1:n_ts);

% Ensure column-vector orientation.
time_arr = time_arr(:);
mean_h_arr = mean_h_arr(:);
mean_x_arr = mean_x_arr(:);
mean_v_arr = mean_v_arr(:);
total_mass_arr = total_mass_arr(:);
q_dyn_max_arr = q_dyn_max_arr(:);
active_count_arr = active_count_arr(:);
total_ke_arr = total_ke_arr(:);
dE_dh_arr = dE_dh_arr(:);

dE_dh_kton_per_km_arr = dE_dh_arr * 1000 / 4.184e12;

out.timeseries = table( ...
    time_arr, ...
    mean_h_arr/1000, ...
    mean_x_arr/1000, ...
    mean_v_arr, ...
    total_mass_arr, ...
    q_dyn_max_arr, ...
    active_count_arr, ...
    total_ke_arr, ...
    dE_dh_arr, ...
    dE_dh_kton_per_km_arr, ...
    'VariableNames', { ...
    't_s', ...
    'h_km', ...
    'x_km', ...
    'v_m_s', ...
    'mass_kg', ...
    'q_dyn_max_Pa', ...
    'fragment_count', ...
    'E_kin_J', ...
    'dE_dh_J_per_m', ...
    'dE_dh_kton_per_km'});

%% R3B DIAGNOSTIC OUTPUTS
child_inserted = false(numel(c_child_id),1);
for rr = 1:numel(c_child_id)
    cid = c_child_id(rr);
    if cid >= 1 && cid <= max_total_fragments
        child_inserted(rr) = child_inserted_by_id(cid);
    end
end

out.r3b.breakup_event_ledger = table( ...
    b_event_id, b_it, b_time_s, b_altitude_m/1000, b_q_Pa, b_Y0_Pa, b_Yeff_Pa, ...
    b_strength_factor, b_strength_clamped, b_generation, b_parent_slot, b_parent_id, ...
    b_parent_mass_kg, b_parent_birth_mass_kg, b_parent_mass_over_global_m0, ...
    b_parent_mass_over_birth_mass, b_intended_nfrag, b_created_count, ...
    b_planned_child_mass_kg, b_created_child_mass_kg, b_cap_dropped_mass_kg, b_cap_hit, ...
    b_parent_ke_J, b_created_child_ke_J, b_created_child_ke_at_parent_speed_J, ...
    b_dispersion_ke_delta_J, b_cap_reference_ke_J, b_total_breakup_ke_delta_J, ...
    b_active_count_before_break, b_next_fragment_id_before, b_next_fragment_id_after, ...
    'VariableNames', { ...
    'event_id','iteration','time_s','altitude_km','q_Pa','Y0_Pa','Yeff_Pa', ...
    'strength_factor','strength_clamped','generation','parent_slot','parent_id', ...
    'parent_mass_kg','parent_birth_mass_kg','parent_mass_over_global_m0', ...
    'parent_mass_over_birth_mass','intended_nfrag','created_count', ...
    'planned_child_mass_kg','created_child_mass_kg','cap_dropped_mass_kg','cap_hit', ...
    'parent_ke_J','created_child_ke_J','created_child_ke_at_parent_speed_J', ...
    'dispersion_ke_delta_J','cap_reference_ke_J','total_breakup_ke_delta_J', ...
    'active_count_before_break','next_fragment_id_before','next_fragment_id_after'});

out.r3b.child_creation_ledger = table( ...
    c_event_id, c_child_order, c_parent_id, c_child_id, c_generation, c_weight, ...
    c_mass_kg, c_birth_strength_Pa, c_dv_m_s, c_dtheta_rad, c_parent_v_m_s, ...
    c_child_v_m_s, c_parent_ke_share_J, c_child_ke_J, child_inserted, ...
    'VariableNames', { ...
    'event_id','child_order','parent_id','child_id','generation','normalized_weight', ...
    'mass_kg','birth_strength_Pa','dv_m_s','dtheta_rad','parent_v_m_s', ...
    'child_v_m_s','parent_ke_share_J','child_ke_J','inserted_into_active_pool'});

out.r3b.deactivation_ledger = table( ...
    d_event_id, d_it, d_time_s, d_fragment_id, d_generation, d_reason, ...
    d_mass_kg, d_ke_J, d_altitude_m/1000, d_velocity_m_s, ...
    'VariableNames', { ...
    'deactivation_event_id','iteration','time_s','fragment_id','generation','reason', ...
    'remaining_mass_kg','remaining_ke_J','altitude_km','velocity_m_s'});

out.r3b.mass_energy_time_ledger = table( ...
    diag_it, diag_time_s, diag_active_mass_kg, diag_active_ke_J, diag_active_count, ...
    diag_min_active_mass_kg, diag_next_fragment_id, diag_breakup_event_count, ...
    diag_cum_ablation_mass_kg, diag_cum_cap_dropped_mass_kg, ...
    diag_cum_insertion_dropped_mass_kg, diag_cum_deactivated_mass_kg, ...
    diag_mass_closure_residual_kg, diag_cum_propagation_ke_delta_J, ...
    diag_cum_breakup_ke_delta_J, diag_cum_dispersion_ke_delta_J, ...
    diag_cum_cap_reference_ke_J, diag_cum_insertion_dropped_ke_J, ...
    diag_energy_closure_residual_J, ...
    'VariableNames', { ...
    'iteration','time_s','active_mass_kg','active_ke_J','active_count', ...
    'minimum_active_fragment_mass_kg','next_fragment_id','breakup_event_count', ...
    'cumulative_ablation_mass_kg','cumulative_cap_dropped_mass_kg', ...
    'cumulative_insertion_dropped_mass_kg','cumulative_deactivated_mass_kg', ...
    'mass_closure_residual_kg','cumulative_propagation_ke_delta_J', ...
    'cumulative_breakup_ke_delta_J','cumulative_dispersion_ke_delta_J', ...
    'cumulative_cap_reference_ke_J','cumulative_insertion_dropped_ke_J', ...
    'energy_closure_residual_J'});

final_mass_closure_residual_kg = m0 - (final_mass_system + cumulative_ablation_mass_kg + ...
    cumulative_cap_dropped_mass_kg + cumulative_insertion_dropped_mass_kg + ...
    cumulative_deactivated_mass_kg);
final_energy_closure_residual_J = E0_system + cumulative_propagation_ke_delta_J + ...
    cumulative_breakup_ke_delta_J - cumulative_insertion_dropped_ke_J - final_ke_system;

out.r3b.closure_summary = table( ...
    m0, E0_system, final_mass_system, final_ke_system, ...
    cumulative_ablation_mass_kg, cumulative_cap_dropped_mass_kg, ...
    cumulative_insertion_dropped_mass_kg, cumulative_deactivated_mass_kg, ...
    final_mass_closure_residual_kg, ...
    cumulative_propagation_ke_delta_J, cumulative_breakup_ke_delta_J, ...
    cumulative_dispersion_ke_delta_J, cumulative_cap_reference_ke_J, ...
    cumulative_insertion_dropped_ke_J, cumulative_deactivated_ke_J, ...
    cumulative_drag_work_estimate_J, cumulative_gravity_work_estimate_J, ...
    final_energy_closure_residual_J, r3b_event_counter, r3b_child_counter, ...
    'VariableNames', { ...
    'initial_mass_kg','initial_ke_J','final_active_mass_kg','final_active_ke_J', ...
    'cumulative_ablation_mass_kg','cumulative_cap_dropped_mass_kg', ...
    'cumulative_insertion_dropped_mass_kg','cumulative_deactivated_mass_kg', ...
    'final_mass_closure_residual_kg', ...
    'cumulative_propagation_ke_delta_J','cumulative_breakup_ke_delta_J', ...
    'cumulative_dispersion_ke_delta_J','cumulative_cap_reference_ke_J', ...
    'cumulative_insertion_dropped_ke_J','cumulative_deactivated_ke_J', ...
    'cumulative_drag_work_estimate_J','cumulative_gravity_work_estimate_J', ...
    'final_energy_closure_residual_J','breakup_event_count','created_child_count'});

%% R3B3 TERMINAL FRAGMENT SNAPSHOT (diagnostic only)
terminal_slots = find([fragments.active])';
if isempty(terminal_slots)
    out.r3b.terminal_fragment_ledger = table( ...
        zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
        zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
        zeros(0,1), zeros(0,1), false(0,1), ...
        'VariableNames', {'slot','fragment_id','generation','mass_kg','radius_m', ...
        'area_m2','altitude_km','downrange_km','velocity_m_s','theta_deg', ...
        'Y0_Pa','Yeff_Pa','broken'});
else
    terminal_id = zeros(numel(terminal_slots),1);
    terminal_gen = zeros(numel(terminal_slots),1);
    terminal_mass = zeros(numel(terminal_slots),1);
    terminal_radius = zeros(numel(terminal_slots),1);
    terminal_area = zeros(numel(terminal_slots),1);
    terminal_h_km = zeros(numel(terminal_slots),1);
    terminal_x_km = zeros(numel(terminal_slots),1);
    terminal_v = zeros(numel(terminal_slots),1);
    terminal_theta_deg = zeros(numel(terminal_slots),1);
    terminal_Y0 = zeros(numel(terminal_slots),1);
    terminal_Yeff = zeros(numel(terminal_slots),1);
    terminal_broken = false(numel(terminal_slots),1);
    for zz = 1:numel(terminal_slots)
        ss = terminal_slots(zz);
        terminal_id(zz) = fragments(ss).id;
        terminal_gen(zz) = fragments(ss).gen;
        terminal_mass(zz) = fragments(ss).m;
        terminal_radius(zz) = fragments(ss).r;
        terminal_area(zz) = fragments(ss).A;
        terminal_h_km(zz) = fragments(ss).h/1000;
        terminal_x_km(zz) = fragments(ss).x/1000;
        terminal_v(zz) = fragments(ss).v;
        terminal_theta_deg(zz) = rad2deg(fragments(ss).theta);
        terminal_Y0(zz) = fragments(ss).Y0;
        terminal_Yeff(zz) = fragments(ss).Yeff;
        terminal_broken(zz) = fragments(ss).broken;
    end
    out.r3b.terminal_fragment_ledger = table( ...
        terminal_slots, terminal_id, terminal_gen, terminal_mass, terminal_radius, ...
        terminal_area, terminal_h_km, terminal_x_km, terminal_v, terminal_theta_deg, ...
        terminal_Y0, terminal_Yeff, terminal_broken, ...
        'VariableNames', {'slot','fragment_id','generation','mass_kg','radius_m', ...
        'area_m2','altitude_km','downrange_km','velocity_m_s','theta_deg', ...
        'Y0_Pa','Yeff_Pa','broken'});
end

out.r3b.instrumentation_version = "R3B3-CAP-HORIZON-DIAGNOSTIC";
out.r3b.additional_random_draw_count = 0;

end






function out = meteor_single_run(params)
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

                new_fragments_to_add = [new_fragments_to_add; frag]; %#ok<AGROW>
                next_fragment_id = next_fragment_id + 1;
            end
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
            if isempty(idx_empty), break; end
            fragments(idx_empty) = new_fragments_to_add(n);
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

end






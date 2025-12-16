function [prdData, info] = predict_Danio_rerio(par, data, auxData)
info = 1;
% unpack par, data, auxData
cPar = parscomp_st(par); vars_pull(par);
vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

%% compute temperature correction factors
TC_ab = tempcorr(temp.ab, T_ref, T_A);
TC_ap = tempcorr(temp.ap, T_ref, T_A);
TC_am = tempcorr(temp.am, T_ref, T_A);
TC_BestAdat2010 = tempcorr(temp.tL_BestAdat2010, T_ref, T_A);
TC_28 = tempcorr(temp.tJO, T_ref, T_A);

%% zero-variate data
% initial
pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
E_0 = U_E0 * p_Am;          % J, energy in egg

% life cycle
[a_b, a_j, ~, L_b, L_j, ~, info] = getAgeAndLengthAtTransitions(par, f, TC_ab, E_0);
if ~info; prdData = []; return; end

% birth
Lw_b = L_b/ del_Mt;                % cm, total length at birth

% metamorphosis
s_M_f = L_j/ L_b;                   % -, acceleration factor at f

% ultimate
L_i = f * L_m * s_M_f;                  % cm, ultimate structural length at f
Lw_i = L_i/ del_Mt;                % cm, ultimate total length at f

% life span
[tau_m, ~, ~] = get_tm_mod('abj', par, f);
aT_m = tau_m/ k_M/ TC_am;               % d, mean life span at T

% puberty at f_EatoFarl1974b
U_E0_p = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
E_0_p = U_E0_p * p_Am;          % J, energy in egg
[~, ~, a_p, ~, ~, L_p, info] = getAgeAndLengthAtTransitions(par, f, TC_ap, E_0_p);
if ~info; prdData = []; return; end

Lw_p = L_p/ del_Mt;                % cm, total length at puberty at F

% pack to output
prdData.ab = a_b;
prdData.ap = a_p;
prdData.am = aT_m;
prdData.Lb = Lw_b;
prdData.Lp = Lw_p;
prdData.Li = Lw_i;

%% validation
V = (L_m * s_M_f)^3; E = E_m * V; E_H = E_Hp; E_R = 0; 
[~, ~, ~, ~, ~, p_R, ~] = compute_powers(V, E, E_H, E_R, s_M_f, TC_28, f, par);
E_R_max = p_R / (kap_R * v / L_m); 
p_C2_max = E_R_max * v / L_m * kap_R;
RT_i_28 =  kap_R * p_C2_max / E_0; % Temperature correction is done inside the compute_powers function
prdData.Ri_28 = RT_i_28;

%% uni-variate data

%% tL_BestAdat2010:
init_cond = [V_0; E_0; 0; 0; 1; 0];
TC = TC_BestAdat2010;  F = f_BestAdat2010; 
a = [0; tL_BestAdat2010(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order

[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); L = V.^(1/3);

L = L(it_sort); L = L(2:end); % reconstuct L
prdData.tL_BestAdat2010 = L / del_Mt;

%% validation
init_cond = [V_0; E_0; 0; 0; 1; 0];
F = f;
[~, VEHRsMG] = ode45(@ode_VEHRsMG, [0; data.tL_BarrFern2010(:,1)], init_cond, [], par, F, TC_28);
V = VEHRsMG(2:end, 1); E = VEHRsMG(2:end, 2); E_H = VEHRsMG(2:end, 3); E_R = VEHRsMG(2:end, 4); s_M = VEHRsMG(2:end, 5); E_egg = VEHRsMG(2:end, 6);

[p_A, ~, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R, s_M, TC_28, F, par);
p_D = compute_dissipation_power(p_S, p_J, p_R, p_C2, E_H, E_Hp, kap_R);

eta_M = -inv(n_M) * n_O * eta_O;
J_O = [p_A, p_D, p_G] * eta_M(3, :)'; % mol/d
J_O = -J_O * 32 * 1e6 / 24; % mug/hr
W = 1 * V + (E + E_R) * w_E / mu_E / d_E; % g

prdData.tL_BarrFern2010 = V.^(1/3) / del_Mt * 10; % mm
prdData.tWw_BarrFern2010 = W * 1e3; % mg
prdData.tJO = J_O./W; % mumol/g/h

end

%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p_A, p_C, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R, s_M, TC, f, p)

p_Am = TC * p.z * p.p_M / p.kap; v = TC * p.v; p_M = TC * p.p_M; k_J = TC * p.k_J;
kap = p.kap; E_G = p.E_G; E_Hb = p.E_Hb;

% Powers
p_A = f * p_Am .* s_M .* V.^(2/3) .* (E_H >= E_Hb);
p_C = E .* (E_G * v ./ V.^(1/3) .* s_M + p_M) ./ (kap * E ./ V + E_G);
p_S = p_M * V;
p_G = kap * p_C - p_S;
p_J = k_J * E_H;
p_R = (1 - kap) * p_C - p_J;
p_C2 = v .* s_M ./ V.^(1/3) .* E_R;

end

function p_D = compute_dissipation_power(p_S, p_J, p_R, p_C2, E_H, E_Hp, kap_R)
p_D = p_S + p_J + p_R .* (E_H < E_Hp) + (1 - kap_R) * p_C2 .* (E_H >= E_Hp);
end

function dVEHRsMG = ode_VEHRsMG(t, VEHRsM, p, f, TC)

V  = VEHRsM(1); % cm, volumetric structural length
E  = VEHRsM(2); % J,   energy in reserve
E_H = VEHRsM(3); % J, E_H maturity
E_R = VEHRsM(4); % J, E_R reproduction buffer
s_M = VEHRsM(5);

E_G = p.E_G; E_Hb = p.E_Hb; E_Hj = p.E_Hj; E_Hp = p.E_Hp; kap_R = p.kap_R;

[p_A, p_C, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R, s_M, TC, f, p);

% State changes
dE = p_A - p_C;
dV = p_G / E_G;
dE_H = p_R * (E_H < E_Hp);
dE_R = (p_R - p_C2) * (E_H >= E_Hp);
ds_M = s_M / 3 / V * dV * ((E_Hb < E_H) && (E_H < E_Hj));
dEgg = kap_R * p_C2 * (E_H >= E_Hp);

dVEHRsMG = [dV; dE; dE_H; dE_R; ds_M; dEgg];

end

function [value, isterminal, direction] = lifeStageTransitions(t, VEHRsMG, par)
% Event function that triggers when y reaches 0
E_H = VEHRsMG(3);
value =  [E_H - par.E_Hb; E_H - par.E_Hj; E_H - par.E_Hp]; % The condition to trigger the event
isterminal = [0; 0; 1]; % Stop the integration when puberty occurs
direction = [0; 0; 0]; % All directions
end

function [a_b, a_j, a_p, L_b, L_j, L_p, info] = getAgeAndLengthAtTransitions(par, f, TC, E_0)
% Define initial state
stateAtFertilization = [par.V_0, E_0, 0, 0, 1, 0];
% Create event to detect life stage transitions
lifeStageEvents = @(t, VEHRsMG) lifeStageTransitions(t, VEHRsMG, par);
options = odeset('Events', lifeStageEvents);
% Solve ODE
ode = @(t, VEHRsMG) ode_VEHRsMG(t, VEHRsMG, par, f, TC);
sol = ode45(ode, [0 1000], stateAtFertilization, options);
if length(sol.xe) ~= 3
    a_b = NaN; a_j = NaN; a_p = NaN; L_b = NaN; L_j = NaN; L_p = NaN;
    info = 0; return;
end
info = 1;
a_b = sol.xe(1);
a_j = sol.xe(2);
a_p = sol.xe(3);
L_b = sol.ye(1, 1)^(1/3);
L_j = sol.ye(1, 2)^(1/3);
L_p = sol.ye(1, 3)^(1/3);

end




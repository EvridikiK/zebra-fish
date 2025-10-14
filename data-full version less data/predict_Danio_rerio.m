function [prdData, info] = predict_Danio_rerio(par, data, auxData)

% unpack par, data, auxData
cPar = parscomp_st(par); vars_pull(par);
vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

% customized filter
filterChecks =   s_G < 0 || (kap_P + kap_X) > 1;
%  f_EatoFarl1974 >1 || f_ValKwa2022 >1 || f_LawrEber2002_high > 1 || ...

if filterChecks
    info = 0;
    prdData = {};
    return;
end

%% compute temperature correction factors
TC_ab = tempcorr(temp.ab, T_ref, T_A);
TC_aj = tempcorr(temp.aj, T_ref, T_A);
TC_ap = tempcorr(temp.ap, T_ref, T_A);
TC_am = tempcorr(temp.am, T_ref, T_A);
TC_Ri = tempcorr(temp.Ri, T_ref, T_A);
TC_GSI = tempcorr(temp.GSI, T_ref, T_A);
% TC_egg = tempcorr(temp.tMC, T_ref, T_A);
TC_Schi2002 = tempcorr(temp.tL_Schi2002, T_ref, T_A);
TC_EatoFarl1974 = tempcorr(temp.tL_EatoFarl1974, T_ref, T_A);
TC_BagaPels2001 = tempcorr(temp.tL_BagaPels2001, T_ref, T_A);
TC_BestAdat2010 = tempcorr(temp.tL_BestAdat2010, T_ref, T_A);
TC_LawrEber2002 = tempcorr(temp.tL_LawrEber2002_high, T_ref, T_A);
TC_28 = tempcorr(temp.tL_YangYama2019, T_ref, T_A);
TC_26_72 = tempcorr(temp.tL_ValKwa2022, T_ref, T_A);

TC_BeauGous2015 = tempcorr(temp.tLf1_BeauGous2015, T_ref, T_A); % juvenile growth trials
TC_tN = tempcorr(temp.tN, T_ref, T_A); % for the reproduction trials in BeauGous2015

% standard length is 80% of the total length:
del_Ms = del_Mt * 1.25; % -,

%% zero-variate data

% life cycle
pars_tj = [g k l_T v_Hb v_Hj v_Hp];
pars_lj = [g; k; l_T; v_Hb; v_Hj];

[t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);

% initial
pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
E_0 = U_E0 * p_Am ;          % J, energy in egg
Wd_0 = E_0 * w_E/ mu_E;      % g, egg dry weight
V0 = Wd_0/ d_E;             % cm^3, egg volume
Lw_0 = (6 * V0/ pi)^(1/3);  % cm, egg diameter

% birth
L_b = L_m * l_b;                  % cm, structural length at birth at f
Lw_b = L_b/ del_Mt;                % cm, total length at birth
aT_b = t_b/ k_M/ TC_ab;           % d, age at birth

% metamorphosis
L_j = l_j * L_m;                  % cm, structural length at metam at f
Lw_j= L_j/ del_Mt;                 % cm, total length at metam at f
s_M_f = L_j/ L_b;                   % -, acceleration factor at f
aT_j = t_j/ k_M/ TC_aj;           % d, age at metam

% puberty is here at end of 0-variate data because of deviating f
Lp_tN = l_p * L_m / del_Mt;
% [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
% L_p = L_m * l_p;                  % cm, structural length at puberty at F
% Lw_p = L_p/ del_Mt;                % cm, total length at puberty at F
% aT_p = t_p/ k_M/ TC_ap;           % d, time since birth at puberty at F and T

% ultimate
L_i = L_m * l_i;                  % cm, ultimate structural length at f
Lw_i = L_i/ del_Mt;                % cm, ultimate total length at f
Ww_i = L_i^3 * (1 + f * w);       % g, ultimate wet weight

% reproduction
% pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector
% RT_i = TC_Ri * reprod_rate_j(L_i, f, pars_R);                 % ultimate reproduction rate
% t_R  = 3; % d, period of accumulaton of reprod buffer at T
% GSIT = (t_R * TC_GSI * k_M * g/ f^3)/ (f + kap * g * y_V_E);
% GSIT = GSIT * ((1 - kap) * f^3 - k_J * U_Hp/ L_m^2/ s_M_f^3);    % -, GSI

% Maximum reproduction rate
V = (L_m * s_M_f)^3; E = E_m * V; E_H = E_Hp; E_R = 0; 
[~, ~, ~, ~, ~, p_R, ~] = compute_powers(V, E, E_H, E_R, s_M_f, TC_Ri, f, par);
E_R_max = p_R / (kap_R * v / L_m); 
p_C2_max = E_R_max * v / L_m * kap_R;
RT_i =  kap_R * p_C2_max / E_0; % Temperature correction is done inside the compute_powers function
% Gonado-somatic index
stateUltimateSizeMaxReprodBuffer = [V, E, E_H, E_R_max, s_M_f, 0];
[~, VEHRsMG] = ode45(@ode_VEHRsMG, [0; init.GSI], stateUltimateSizeMaxReprodBuffer, [], par, f, TC_GSI); % Simulate evolution for 3 days
MV = VEHRsMG(end, 1) * M_V * w_V;
ME = VEHRsMG(end, 2) / mu_E * w_E;
ME_R = VEHRsMG(end, 4) / mu_E * w_E;
ME_e = VEHRsMG(end, 6) / mu_E * w_E;
GSIT = (ME_e) / (ME + MV + ME_R + ME_e); % Temperature correction is done inside the compute_powers function

% life span
% pars_tm = [g; k; v_Hb; v_Hj; v_Hp; h_a/k_M^2; s_G];
[tau_m, ~, ~] = get_tm_mod('abj', par, f);
% if ~info; prdData=[]; return; end
aT_m = tau_m/ k_M/ TC_am;               % d, mean life span at T

% puberty at f_EatoFarl1974b
F = f_EatoFarl1974;
[t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, F);
L_p = L_m * l_p;                  % cm, structural length at puberty at F
Lw_p = L_p/ del_Mt;                % cm, total length at puberty at F
aT_p = t_p/ k_M/ TC_ap;           % d, time since birth at puberty at F and T

% pack to output
prdData.ab = aT_b;
prdData.aj = aT_j;
prdData.ap = aT_p;
prdData.am = aT_m;
prdData.L0 = Lw_0;
prdData.Lb = Lw_b;
prdData.Lj = Lw_j;
prdData.Lp = Lw_p;
prdData.Li = Lw_i;
prdData.Wd0 = Wd_0;
prdData.Wwi = Ww_i;
prdData.Ri = RT_i;
prdData.GSI = GSIT;

%% uni-variate data

%% Feeding data (Valentine and Kwasek 2022)

init_cond = [V_0; E_0; 0; 0; 1; 0];
F = f_ValKwa2022;
[~, VEHRsMG] = ode45(@ode_VEHRsMG, [0; data.tL_ValKwa2022(:,1)], init_cond, [], par, F, TC_26_72);
V = VEHRsMG(2:end, 1); E = VEHRsMG(2:end, 2); E_H = VEHRsMG(2:end, 3); E_R = VEHRsMG(2:end, 4); s_M = VEHRsMG(2:end, 5); E_egg = VEHRsMG(2:end, 6);
[p_A, ~, ~, ~, ~, ~, ~] = compute_powers(V, E, E_H, E_R, s_M, TC_26_72, F, par);

L = V.^(1/3);
W = 1 * V + (E + E_R) * w_E / mu_E / d_E;
JX = w_X / kap_X / mu_X .* p_A;

prdData.tL_ValKwa2022 = L / del_Mt * 10; % mm
prdData.tWw_ValKwa2022 = W * 1000; % mg
prdData.tJX_ValKwa2022 = JX ./ W * 100 / auxData.init.tJX_ValKwa2022; % (%bw/d)

%% Data from Yang et al. 2019
init_cond = [V_0; E_0; 0; 0; 1; 0];
F = f_YangYama2019;
[~, VEHRsMG] = ode45(@ode_VEHRsMG, [0; data.tL_YangYama2019(:,1)], init_cond, [], par, F, TC_28);
V = VEHRsMG(2:end, 1); E = VEHRsMG(2:end, 2); E_H = VEHRsMG(2:end, 3); E_R = VEHRsMG(2:end, 4); s_M = VEHRsMG(2:end, 5); E_egg = VEHRsMG(2:end, 6);

[p_A, ~, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R, s_M, TC_28, F, par);
p_D = compute_dissipation_power(p_S, p_J, p_R, p_C2, E_H, E_Hp, kap_R);

eta_M = -inv(n_M) * n_O * eta_O;
J_O = [p_A, p_D, p_G] * eta_M(3, :)'; % mol/d
J_O = -J_O * 32 * 1e6 / 24; % mug/hr
W = 1 * V + (E + E_R) * w_E / mu_E / d_E; % g

prdData.tL_YangYama2019 = V.^(1/3) / del_Mt; % cm
prdData.tWw_YangYama2019 = W; % g
prdData.tJO_YangYama2019 = J_O; % mumol/h

%% BangGron2004: tMC and tMN, T = 25 �C
% % F = f; TC = TC_egg;
% % init_cond = [V_0; E_0; 0; 0; 1; 0];
% % [~, VEHRsMG] = ode45(@ode_VEHRsMG, [0; data.tMC(:,1)], init_cond, [], par, F, TC);
% % V = VEHRsMG(2:end, 1); E = VEHRsMG(2:end, 2);
% % MV = V * M_V; % mol of V
% % ME = E / mu_E; % mol of E
% % prdData.tMC = 1e6 * (MV + ME); % mumol, carbon mass
% % prdData.tMN = 1e6 * (n_O(4,2) * MV + n_O(4,3) * ME); % mumol, nitroen mass

% compilation of literature growth curves (see previous version and AuguGagn2011 for more discussion on this):

% initial conditions at start:
% hardly any structure, initial enegy of the egg (J), no a maturiaty at start

%% Length data from BeauGous2015
init_cond = [V_0; E_0; 0; 0; 1; 0];
TC = TC_BeauGous2015;

% --------- high food level:
F = f_BeauGous2015L; 
a = [0; tLf1_BeauGous2015(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order
[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); L = V.^(1/3);
L = L(it_sort); L = L(2:end); % reconstuct L
prdData.tLf1_BeauGous2015 = L / del_Mt * 10;

% ----------- medium food level:
F = f_BeauGous2015L * 0.9; 
a = [0; tLf2_BeauGous2015(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order
[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); L = V.^(1/3);
L = L(it_sort); L = L(2:end); % reconstuct L
prdData.tLf2_BeauGous2015 = L / del_Mt * 10;

% --------------- low food level:
F = f_BeauGous2015L * 0.8; 
a = [0; tLf3_BeauGous2015(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order
[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); L = V.^(1/3);
L = L(it_sort); L = L(2:end); % reconstuct L
prdData.tLf3_BeauGous2015 = L / del_Mt * 10;

%% tL_LawrEber2002 at low and high
init_cond = [V_0; E_0; 0; 0; 1; 0];
TC = TC_LawrEber2002;

% --- Low:
F = f_LawrEber2002_low;  
a = [0; tL_LawrEber2002_low(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order

[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); L = V.^(1/3);

L = L(it_sort); L = L(2:end); % reconstuct L
prdData.tL_LawrEber2002_low = L / del_Mt;

% --- High:
F = f_LawrEber2002_high; 
a = [0; tL_LawrEber2002_high(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order

[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); L = V.^(1/3);

L = L(it_sort); L = L(2:end); % reconstuct L
prdData.tL_LawrEber2002_high = L / del_Mt;


%% tL_BestAdat2010:
init_cond = [V_0; E_0; 0; 0; 1; 0];
TC = TC_BestAdat2010;  F = f_BestAdat2010; 
a = [0; tL_BestAdat2010(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order

[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); L = V.^(1/3);

L = L(it_sort); L = L(2:end); % reconstuct L
prdData.tL_BestAdat2010 = L / del_Mt;


%% tL of Schi2002
init_cond = [V_0; E_0; 0; 0; 1; 0];
TC = TC_Schi2002; F = f_Schi2002; 
a = [0; tL_Schi2002(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order

[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); L = V.^(1/3);

L = L(it_sort); L = L(2:end); % reconstuct L
prdData.tL_Schi2002 = L / del_Mt;


%% tL of EatoFarl1974
init_cond = [V_0; E_0; 0; 0; 1; 0];
TC = TC_EatoFarl1974; F = f_EatoFarl1974; 
a = [0; tL_EatoFarl1974(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order

[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); L = V.^(1/3);

L = L(it_sort); L = L(2:end); % reconstuct L
prdData.tL_EatoFarl1974 = L / del_Mt;


%% tL, tWw and tWd of BagaPels2001
init_cond = [V_0; E_0; 0; 0; 1; 0];
F = f_BagaPels2001; TC = TC_BagaPels2001;

[~, VEHRsMG] = ode45(@ode_VEHRsMG, [0; tL_BagaPels2001(:,1)], init_cond, [], par, F, TC);
V = VEHRsMG(2:end, 1); E = VEHRsMG(2:end, 2); E_R = VEHRsMG(2:end, 4); E_egg = VEHRsMG(2:end, 6);

prdData.tWw_BagaPels2001  = 1 * V + w_E / mu_E / d_E * (E + E_R); % g, wet weight
prdData.tWd_BagaPels2001  = d_V * V + w_E / mu_E * (E + E_R); % g, dry weight
prdData.tL_BagaPels2001   = V.^(1/3) / del_Mt; % cm, total length

%% reproduction trials with individual females BeauGous2015
F = f_BeauGous2015R; TC = TC_tN;
pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
U_E0 = initial_scaled_reserve(F, pars_UE0); % d.cm^2, initial scaled reserve
E_0 = U_E0 * p_Am;  
[~, ~, ~, l_j, l_p, l_b, ~, ~, ~, info] = get_tj(pars_tj, F);
if ~info; prdData = []; return; end
s_M_F = l_j/l_b;
V_init =  (max(l_p*L_m, init.tN * del_Mt))^3; % Either initial length or puberty length, whichever is largest
E_init   = F * E_m * V_init;                  % J, inital energy in reserve
max_E_R = get_max_E_R(V_init, par, E_m, F, s_M_F);

init_cond = [V_init; E_init; E_Hp; max_E_R; s_M_F; 0];
[~, VEHRsMG] = ode45(@ode_VEHRsMG, [0; data.tN(:,1)], init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); E = VEHRsMG(:, 2); E_H = VEHRsMG(:, 3); E_R = VEHRsMG(:, 4); s_M = VEHRsMG(:, 5); E_egg = VEHRsMG(:, 6);

prdData.tN = E_egg(2:end) / E_0;
prdData.tL1 = V([1 end]).^(1/3) / del_Ms;
prdData.Wwt = 1 * V(end) + w_E / mu_E / d_E * (E(end) + E_R(end));

end

%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p_A, p_C, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R, s_M, TC, f, p)

p_Am = TC * p.z * p.p_M / p.kap; v = TC * p.v; p_M = TC * p.p_M; k_J = TC * p.k_J;
kap = p.kap; E_G = p.kap; E_Hb = p.E_Hb;

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

E_G = p.kap; E_Hb = p.E_Hb; E_Hj = p.E_Hj; E_Hp = p.E_Hp; kap_R = p.kap_R;

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

function max_E_R = get_max_E_R(V, par, E_m, f, s_M_max)
E = f * V * E_m;
E_H = par.E_Hp;
E_R = 0;
TC = 1;
[~, ~, ~, ~, ~, p_R, ~] = compute_powers(V, E, E_H, E_R, s_M_max, TC, f, par);
max_E_R = p_R ./ (par.v * s_M_max ./ V.^(1/3));
end







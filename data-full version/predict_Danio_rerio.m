function [prdData, info] = predict_Danio_rerio(par, data, auxData)

% unpack par, data, auxData
cPar = parscomp_st(par); vars_pull(par);
vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

% customized filter
filterChecks =   E_R_init_DrewRodn2008 < 0 || E_R_init_DrewRodn2008 > 2000 || ...
    ~reach_birth(g, k, v_Hb, f_DrewRodn2008) || ...
    f_DrewRodn2008 > 1 || f_EatoFarl1974 >1 || f_ValKwa2022 >1 || f_LawrEber2002_high > 1 || ...
    s_shrink < 0 || s_G < 0 || del_X < 0 || (kap_P + kap_X) > 1;

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
TC_tS = tempcorr(temp.tS, T_ref, T_A);
TC_BestAdat2010 = tempcorr(temp.tL_BestAdat2010, T_ref, T_A);
TC_LawrEber2002 = tempcorr(temp.tL_LawrEber2002_high, T_ref, T_A);
TC_28 = tempcorr(temp.tL_YangYama2019, T_ref, T_A);
TC_26_72 = tempcorr(temp.tL_ValKwa2022, T_ref, T_A);

TC_starv = tempcorr(temp.tW, T_ref, T_A);

TC_BeauGous2015 = tempcorr(temp.tLf1_BeauGous2015, T_ref, T_A); % juvenile growth trials
TC_tN = tempcorr(temp.tN, T_ref, T_A); % for the reproduction trials in BeauGous2015
TC_tSstarv = tempcorr(temp.tS_starv, T_ref, T_A);

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
s_M_f = L_j/ L_b;                   % -, acceleration factor
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
pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector
RT_i = TC_Ri * reprod_rate_j(L_i, f, pars_R);                 % ultimate reproduction rate
t_R  = 3; % d, period of accumulaton of reprod buffer at T
GSIT = (t_R * TC_GSI * k_M * g/ f^3)/ (f + kap * g * y_V_E);
GSIT = GSIT * ((1 - kap) * f^3 - k_J * U_Hp/ L_m^2/ s_M_f^3);    % -, GSI

% init_cond = [1e-10; E_0; 0; 0; 1; 0];
% opts = odeset('Events',@ultimate_size);
% [t, VEHRsMG] = ode45(@ode_VEHRsMG, [0; 10*365], init_cond, [], par, f, TC_Ri);
% E_Ri = VEHRsMG(end, 4);
% RT_i = kap_R*v*s_M_f*L_i*E_Ri;                 % ultimate reproduction rate

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

init_cond = [1e-10; E_0; 0; 0; 1; 0];
F = f_ValKwa2022;
[~, VEHRsMG] = ode45(@ode_VEHRsMG, [0; data.tL_ValKwa2022(:,1)], init_cond, [], par, F, TC_26_72);
V = VEHRsMG(2:end, 1); E = VEHRsMG(2:end, 2); E_H = VEHRsMG(2:end, 3); E_R = VEHRsMG(2:end, 4); s_M = VEHRsMG(2:end, 5); E_egg = VEHRsMG(2:end, 5);
[p_A, ~, ~, ~, ~, ~, ~] = compute_powers(V, E, E_H, E_R, s_M, TC_26_72, F, par);

L = V.^(1/3);
W = 1 * V + (E + E_R + E_egg) * w_E / mu_E / d_E;
JX = w_X / kap_X / mu_X .* p_A;

prdData.tL_ValKwa2022 = L / del_Mt * 10; % mm
prdData.tWw_ValKwa2022 = W * 1000; % mg
prdData.tJX_ValKwa2022 = JX ./ W * 100 / auxData.init.tJX_ValKwa2022; % (%bw/d)


%% Oxygen consumption
% % init_cond = [1e-10; E_0; 0; 0; 1; 0];
% % F = f_BarrFern2010;
% % [~, VEHRsMG] = ode45(@ode_VEHRsMG, [0; data.tW_BarrFern2010(:,1)], init_cond, [], par, F, TC_28);
% % V = VEHRsMG(2:end, 1); E = VEHRsMG(2:end, 2); E_H = VEHRsMG(2:end, 3); E_R = VEHRsMG(2:end, 4); s_M = VEHRsMG(2:end, 5); E_egg = VEHRsMG(2:end, 5);
% % 
% % [p_A, ~, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R, s_M, TC_28, F, par);
% % p_D = compute_dissipation_power(p_S, p_J, p_R, p_C2, E_H, E_Hp, kap_R);
% % 
% % eta_M = -inv(n_M) * n_O * eta_O;
% % J_O = [p_A, p_D, p_G] * eta_M(3, :)' .* (L_m^2 * p_Am) * TC_28; % mol/d
% % J_O = -J_O * 1e6 / 24; % mumol/g/d
% % W = 1 * V + (E + E_R + E_egg) * w_E / mu_E / d_E; % g
% % 
% % [~, idx_tL] = ismember(data.tL_BarrFern2010(:,1) , data.tW_BarrFern2010(:,1));
% % prdData.tL_BarrFern2010 = V(idx_tL).^(1/3) / del_Mt * 10; % mm
% % prdData.tW_BarrFern2010 = W * 1e3; % mg
% % [~, idx_tJO] = ismember(data.tJO_BarrFern2010(:,1) , data.tW_BarrFern2010(:,1));
% % prdData.tJO_BarrFern2010 = J_O(idx_tJO) ./ W(idx_tJO); % mumol/g/h


%% Oxygen consumption
% % TC_JO = tempcorr(C2K(temp.tTJO_BarrBurg1999), T_ref, T_A);
% % init_cond = [1e-10; E_0; 0; 0; 1; 0];
% % F = f_BarrBurg1999;
% % for i=1:length(TC_JO)
% %     [~, VEHRsMG] = ode45(@ode_VEHRsMG, [0; data.tTL_BarrBurg1999(:,1)], init_cond, [], par, F, TC_JO(i));
% %     TV(:,i) = VEHRsMG(2:end, 1); TE(:,i) = VEHRsMG(2:end, 2); TE_H(:,i) = VEHRsMG(2:end, 3); TE_R(:,i) = VEHRsMG(2:end, 4); Ts_M(:, i) = VEHRsMG(2:end, 5); TE_egg(:,i) = VEHRsMG(2:end, 5);
% % 
% %     [p_A, ~, p_S, p_G, p_J, p_R, p_C2] = compute_powers(TV(:,i), TE(:,i), TE_H(:,i), TE_R(:,i), Ts_M(:,i), TC_JO(i), F, par);
% %     p_D = compute_dissipation_power(p_S, p_J, p_R, p_C2, TE_H(:,i), E_Hp, kap_R);
% % 
% %     eta_M = -inv(n_M) * n_O * eta_O;
% %     J_O = [p_A, p_D, p_G] * eta_M(3, :)' * TC_JO(i); % mol/d
% %     TJO(:,i) = -J_O * 1e6 / 24; % mumol/g/d
% %     TW(:,i) = 1 * TV(:,i) + (TE(:,i) + TE_R(:,i) + TE_egg(:,i)) * w_E / mu_E / d_E; % g
% % end
% % 
% % prdData.tTL_BarrBurg1999 = TV.^(1/3) / del_Mt * 10; % mm
% % prdData.tTWw_BarrBurg1999 = TW * 1e3; % mg
% % prdData.tTJO_BarrBurg1999 = TJO ./ TW; % mumol/g/h

%% Data from Yang et al. 2019
init_cond = [1e-10; E_0; 0; 0; 1; 0];
F = f_YangYama2019;
[~, VEHRsMG] = ode45(@ode_VEHRsMG, [0; data.tL_YangYama2019(:,1)], init_cond, [], par, F, TC_28);
V = VEHRsMG(2:end, 1); E = VEHRsMG(2:end, 2); E_H = VEHRsMG(2:end, 3); E_R = VEHRsMG(2:end, 4); s_M = VEHRsMG(2:end, 5); E_egg = VEHRsMG(2:end, 5);

[p_A, ~, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R, s_M, TC_28, F, par);
p_D = compute_dissipation_power(p_S, p_J, p_R, p_C2, E_H, E_Hp, kap_R);

eta_M = -inv(n_M) * n_O * eta_O;
J_O = [p_A, p_D, p_G] * eta_M(3, :)'; % mol/d
J_O = -J_O * 32 * 1e6 / 24; % mug/hr
W = 1 * V + (E + E_R + E_egg) * w_E / mu_E / d_E; % g

prdData.tL_YangYama2019 = V.^(1/3) / del_Mt; % cm
prdData.tWw_YangYama2019 = W; % g
prdData.tJO_YangYama2019 = J_O; % mumol/h

%% t-S data for larvae post hatch (GeffSimo2013)
% assumption: birth is hatch; no metabolic acceleration during starvation

% temperature correct parameters:
pT_Am = p_Am *  TC_tSstarv; vT = v * TC_tSstarv; pT_M = p_M * TC_tSstarv;
kT_J  = k_J * TC_tSstarv; kT_J1 = kT_J;
% solve ODEs:
[tt, LEHS] =  ode45(@dget_LEHS, tS_starv(:,1), [L_b; E_m; E_Hb; 1; L_b; E_Hb] ,[],...
    pT_Am, vT, pT_M, E_G, kap, kap_G, kT_J, kT_J1, s_rejuv, s_shrink, del_X, 0);

prdData.tS_starv = LEHS(:,4);

%% tS of GerhKauf2002
F = f;%f_EatoFarl1974;
a     = tS(:,1); % d, time since birth
TC    = TC_tS;      % -, temp corr factor
s_M   = l_j/ l_b;             % -, acceleration factor
L_i   = l_i * L_m;            % cm, assumption that no more growth occurs
h3_W  = TC^3 * h_a * F * v * s_M/ (6 * L_i);
h_G   = TC * s_G * F^4 * v * s_M/ L_i;
prdData.tS = exp(6 * h3_W /h_G^3 * (1 - exp(h_G * a) + h_G * a + h_G^2 * a.^2 /2));

%% BangGron2004: tMC and tMN, T = 25 �C
% % F = f; TC = TC_egg;
% % init_cond = [1e-8; E_0; 0; 0; 1; 0];
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
init_cond = [1e-10; E_0; 0; 0; 1; 0];
TC = TC_BeauGous2015;

% --------- high food level:
F = f_BeauGous2015; 
a = [0; tLf1_BeauGous2015(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order
[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); L = V.^(1/3);
L = L(it_sort); L = L(2:end); % reconstuct L
prdData.tLf1_BeauGous2015 = L / del_Mt * 10;

% ----------- medium food level:
F = f_BeauGous2015 * 0.9; 
a = [0; tLf2_BeauGous2015(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order
[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); L = V.^(1/3);
L = L(it_sort); L = L(2:end); % reconstuct L
prdData.tLf2_BeauGous2015 = L / del_Mt * 10;

% --------------- low food level:
F = f_BeauGous2015 * 0.8; 
a = [0; tLf3_BeauGous2015(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order
[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); L = V.^(1/3);
L = L(it_sort); L = L(2:end); % reconstuct L
prdData.tLf3_BeauGous2015 = L / del_Mt * 10;

%% tL_LawrEber2002 at low and high
init_cond = [1e-10; E_0; 0; 0; 1; 0];
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
init_cond = [1e-10; E_0; 0; 0; 1; 0];
TC = TC_BestAdat2010;  F = f_BestAdat2010; 
a = [0; tL_BestAdat2010(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order

[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); L = V.^(1/3);

L = L(it_sort); L = L(2:end); % reconstuct L
prdData.tL_BestAdat2010 = L / del_Mt;


%% tL of Schi2002
init_cond = [1e-10; E_0; 0; 0; 1; 0];
TC = TC_Schi2002; F = f_Schi2002; 
a = [0; tL_Schi2002(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order

[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); L = V.^(1/3);

L = L(it_sort); L = L(2:end); % reconstuct L
prdData.tL_Schi2002 = L / del_Mt;


%% tL of EatoFarl1974
init_cond = [1e-10; E_0; 0; 0; 1; 0];
TC = TC_EatoFarl1974; F = f_EatoFarl1974; 
a = [0; tL_EatoFarl1974(:,1)];
[t_sort, ~, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order

[~, VEHRsMG] = ode45(@ode_VEHRsMG, t_sort, init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); L = V.^(1/3);

L = L(it_sort); L = L(2:end); % reconstuct L
prdData.tL_EatoFarl1974 = L / del_Mt;


%% tL, tWw and tWd of BagaPels2001
init_cond = [1e-10; E_0; 0; 0; 1; 0];
F = f_BagaPels2001; TC = TC_BagaPels2001;

[~, VEHRsMG] = ode45(@ode_VEHRsMG, [0; tL_BagaPels2001(:,1)], init_cond, [], par, F, TC);
V = VEHRsMG(2:end, 1); E = VEHRsMG(2:end, 2); E_R = VEHRsMG(2:end, 4); E_egg = VEHRsMG(2:end, 5);

prdData.tWw_BagaPels2001  = 1 * V + w_E / mu_E / d_E * (E + E_R + E_egg); % g, wet weight
prdData.tWd_BagaPels2001  = d_V * V + w_E / mu_E * (E + E_R + E_egg); % g, dry weight
prdData.tL_BagaPels2001   = V.^(1/3) / del_Mt; % cm, total length

%% reproduction trials with individual females BeauGous2015
F = f; TC = TC_tN;
V_init =  max(Lp_tN*del_Mt, L_i)^3;
E_init   = f * E_m * V_init;                  % J, inital energy in reserve
init_cond = [V_init; E_init; E_Hp; E_R_init_BeauGous2015; s_M_f; 0];
% init_cond = [V_init; E_init; E_Hp; 0; s_M_f; 0];

[~, VEHRsMG] = ode45(@ode_VEHRsMG, data.tN(:,1), init_cond, [], par, F, TC);
V = VEHRsMG(:, 1); E = VEHRsMG(:, 2); E_H = VEHRsMG(:, 3); E_R = VEHRsMG(:, 4); s_M = VEHRsMG(:, 5); E_egg = VEHRsMG(:, 6);

prdData.tN = E_egg / E_0;
prdData.tL1 = V([1 end]).^(1/3) / del_Ms;
prdData.Wwt = 1 * V(end) + w_E / mu_E / d_E * (E(end) + E_R(end) + E_egg(end));

%% Starvation data: adults
TC = TC_starv; F = f_DrewRodn2008;
[lj, ~, lb, info] = get_lj(pars_lj, F);
if ~info; prdData=[]; return; end

% initial conditions for the ODE simulations
L_init  = 2.49 * del_Ms; % cm, structural length at start
E_init = F * E_m * L_init^3; % J, inital energy in reserve
InitCond = [L_init; E_init; E_Hp; E_R_init_DrewRodn2008]; % concatenate initial conditions
s_M = lj / lb;

% growth during the first 11 days in the fed condition:
t = tW(:,1);
[tt, LEHR] = ode45(@ode_LEHR, [t(t<=11); 12], InitCond,[], par, cPar, F, s_M, TC);
LEHR(end,:) = [];
% unpack vars and calculate output:
L   = LEHR(:,1); % cm,structural length
E   = LEHR(:,2); % J, energy in reserve
E_R = LEHR(:,4); % J,  reproduction buffer
EWw1 = w_E/ mu_E/ d_E * (E + E_R) + L.^3; % g, total wet weight
ELw  = [L(1); L(end)]/ del_Ms * 10;             % mm, standard length
% growth after day 11 in the fed condition:
[tt, LEHR] = ode45(@ode_LEHR, t(t>11), LEHR(end,:), [], par, cPar, F*0.7, s_M, TC);
% unpack pars and calculate output:
L   = LEHR(:,1); % cm,structural length
E   = LEHR(:,2); % J, energy in reserve
E_R = LEHR(:,4); % J,  reproduction buffer
EWw2 = L.^3 + w_E/ mu_E/ d_E * (E + E_R); % g, total wet weight
EWw = [EWw1; EWw2]; % g, concatenate total wet weight over full experiment
% growth during the first 11 when individual were fed:
t = tWs(:,1); % d, time
[tt, LEHR] = ode45(@ode_LEHR, [t(t<=11);12], InitCond,[], par, cPar, F, s_M, TC);
LEHR(end,:) = []; % remove the last dummy time
% unpack vars:
L    = LEHR(:,1); % cm,structural length
E    = LEHR(:,2); % J, energy in reserve
E_R  = LEHR(:,4); % J,  unripe buffer
EWw1 = L.^3 + w_E/mu_E/ d_E * (E + E_R); % g, total wet weight
[tt, LEHR] = ode45(@ode_LEHR, t(t>11), LEHR(end,:),[], par, cPar, 0, s_M, TC);
L   = LEHR(:,1); % cm,structural length
E   = LEHR(:,2); % J, energy in reserve
E_R = LEHR(:,4); % J,  unripe buffer
EWw2 = L.^3 + w_E/ mu_E/ d_E * (E + E_R); % g, total wet weight
EWws = [EWw1; EWw2]; % g, concatenate total wet weight over full experiment

prdData.tL = ELw;
prdData.tW = EWw;
prdData.tWs = EWws;


end

%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p_A, p_C, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R, s_M, TC, f, p)

p_Am = TC * p.z * p.p_M / p.kap; v = TC * p.v; p_M = TC * p.p_M; k_J = TC * p.k_J;
kap = p.kap; E_G = p.kap; E_Hb = p.E_Hb; kap_R = p.kap_R;

% Powers
p_A = f * p_Am .* s_M .* V.^(2/3) .* (E_H >= E_Hb);
p_C = E .* (E_G * v ./ V.^(1/3) .* s_M + p_M) ./ (kap * E ./ V + E_G);
p_S = p_M * V;
p_G = kap * p_C - p_S;
p_J = k_J * E_H;
p_R = (1 - kap) * p_C - p_J;
p_C2 = kap_R * v .* s_M ./ V.^(1/3) .* E_R;

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

E_G = p.kap; E_Hb = p.E_Hb; E_Hj = p.E_Hj; E_Hp = p.E_Hp;

[p_A, p_C, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R, s_M, TC, f, p);

% State changes
dE = p_A - p_C;
dV = p_G / E_G;
dE_H = p_R * (E_H < E_Hp);
dE_R = (p_R - p_C2) * (E_H >= E_Hp);
ds_M = s_M / 3 / V * dV * ((E_Hb < E_H) && (E_H < E_Hj));
dEgg = p_C2 * (E_H >= E_Hp);

dVEHRsMG = [dV; dE; dE_H; dE_R; ds_M; dEgg];

end

% --------------------------------------------------
% ODE only with unripe buffer
% --------------------------------------------------
function dLEHR = ode_LEHR(t, LEHR, p, c, f, s_M, TC)
%
% Input:
% p: structure 'par'
% c: structure 'Cpar' obtained by cPar = parscomp_st(par)
% f: scaled, scaled functional response,
% s_M: scalar, -, acceleration factor post metamorphosis
% TC, scalar, -, temperature correction factor

% --------------- unpack LEHR ------------------------------------------
L   =  max(0,LEHR(1)); % cm, volumetric structural length
E   =  max(0,LEHR(2)); % J,   energy in reserve
EH  =  min(p.E_Hp,LEHR(3)); % J, E_H maturity
ER  =  max(0,LEHR(4)); % J, E_R reproduction buffer

% Temperature correct the relevant paramters
vT = p.v * TC * s_M; pT_M = p.p_M * TC; kT_J = p.k_J * TC; pT_Am = c.p_Am * TC * s_M;
%
pA   = f * pT_Am * L^2 * (EH >= p.E_Hb);           % J/d, assimilation
%
if EH < p.E_Hp % juveniles cannot cover somatic maintenance with the buffer
    r  = (E * vT/ L - pT_M * L^3/ p.kap)/ (E + p.E_G * L^3/ p.kap) * ...
        (E >= pT_M * L^4/ (p.kap * vT * s_M)) + ...
        (E * vT/ L - pT_M * L^3/ p.kap)/ (E + c.kap_G * p.E_G * L^3/ p.kap) ...
        * (E < pT_M * L^4/ (p.kap * vT * s_M));
    pC   = E * (vT/ L - r); % J/d, mobilisation
    dE   = pA - pC;                                          % J/d, change in energy in reserve
    dL   = r/ 3 * L;                                         % cm/d, change in structural length
    dEH  = max(0, (1 - p.kap) * pC - kT_J * EH) * (EH < p.E_Hp);    % J/d, change in cum energy invested in maturation (it is implied here that no rejuvenation occurs).
    dER  = 0;
else % EH = EHp: adults
    pC = E * (p.E_G * vT * L^2 + pT_M * L^3)/ (p.kap * E + p.E_G * L^3);
    if p.kap * pC >= pT_M * L^3   % enough energy in reserve to cover somatic maintenance and enough to make a batch
        r    = (E * vT/ L^4 - pT_M/ p.kap)/ (E/ L^3 + p.E_G/ p.kap); % d^-1, specific growth rate
        dE   = pA - pC;                                          % J/d, change in energy in reserve
        dL   = r/ 3 * L;                                         % cm/d, change in structural length
        dEH  = 0;    % J/d, change in cum energy invested in maturation (it is implied here that no rejuvenation occurs).
        dER  = ((1 - p.kap) * pC - kT_J * p.E_Hp) ;       % J, change in cumulated energy invested in the unripe reproduction buffer
    else  % not enough energy in reserve to cover somatic maintenance
        if ER > 0
            r = 0;
        else
            r    =  (E * vT/ L - pT_M * L^3/ p.kap)/ ...
                (E + c.kap_G * p.E_G  * L^3/ p.kap); % d^-1, specific growth rate
        end
        dE   = pA - pC;                                         % J/d, change in energy in reserve
        dL   = r/ 3 * L;                                        % cm/d, change in structural length
        dEH  = 0;                                               % J/d, change in cum energy invested in maturation (it is implied here that no rejuvenation occurs).
        dER  = (1 - p.kap) * pC - kT_J * p.E_Hp;
        dER  = (dER  - (pT_M * L^3 - p.kap * pC)) * (ER > 0) ;
    end
end
% pack dLEHR
dLEHR = [dL; dE; dEH; dER];
end


function dLEHS = dget_LEHS(t,LEHS,...
    p_Am, v, p_M, E_G, kap, kap_G, k_J, k_J1, s_rejuv, s_shrink, del_X, f)

% assumption: no metabolic acceleration
% the function will not handle re-growth

% unpack LEHRU
L   =  LEHS(1); % cm, volumetric structural length
E   =  LEHS(2); % J cm^{-3}, [E] reserve density
EH  =  LEHS(3); % J, E_H maturity
S   =  min(1,LEHS(4)); % survival probability
maxL = LEHS(5); % maximum length
maxEH = LEHS(6); % maximum maturity

L_m = kap * p_Am/ p_M; % cm, ultimate length
k_M = p_M/ E_G;  % 1/d, maturity maintenance rate coefficient
E_m = p_Am/ v;   % J/cm^3, max reserve density
l = L/ L_m; e = E/E_m;  % -, scaled structural length and scaled res. dens.

p_C = E * (E_G * v/ L + p_M)/ (kap * E + E_G );   % J/cm^3 (2.12, Kooy2010) specific mobilisation flux
dE =  (f * p_Am - E * v)/ L; % J day^{-1} cm^{-3} (2.10, Kooy2010)

% pp.42 comments DEB3 equ. 4.2
if e < l
    r = (E * v/ L - p_M/ kap)/ (E + E_G * kap_G/ kap);
else
    r = (E * v/ L - p_M/ kap)/ (E + E_G/ kap);
end
dL  = L * r/ 3;

dmaxL = max(0,dL); % cm, max stuctural length
dEH = (1 - kap) * p_C * L^3 - k_J * EH; % J/d

if dEH < 0 % rejuvination
    dEH = -k_J1 * (EH - (1- kap) * p_C * L^3/ k_J);   % J/d,
end

dmaxEH = max(0,dEH); % J, max maturity

h_J  =  k_M * s_rejuv * kap/ (1 - kap)/ (E_G * L_m^3) * (maxEH - EH); % 1/d, hazard from rejuvenation
h_sh = s_shrink * k_M * max(0, maxL - L)/maxL * (L < del_X * maxL);

dS = - S * (h_J + h_sh); % 1/d, survival fraction

% pack state variables
dLEHS = [dL; dE; dEH; dS; dmaxL; dmaxEH];

end

function [value,isterminal,direction] = ultimate_size(t, VEHRsM, p, f, TC)
  L_m = p.z; 
  if VEHRsM(3) >= p.E_Hp
      a = (VEHRsM(5)* L_m)^3;
  else
      a = 0;
  end
  value = VEHRsM(1) - a;
  isterminal = 1;
  direction = 0;
end







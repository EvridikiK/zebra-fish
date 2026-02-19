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

%% zero-variate data
% initial
E_0 = getE0(f, par, cPar); % J, energy in egg

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

% puberty
F = f;
E_0_p = getE0(F, par, cPar); % J, energy in egg
[~, ~, a_p, ~, ~, L_p, info] = getAgeAndLengthAtTransitions(par, F, TC_ap, E_0_p);
if ~info; prdData = []; return; end

Lw_p = L_p/ del_Mt;                % cm, total length at puberty at F

% pack to output
prdData.ab = a_b;
prdData.ap = a_p;
prdData.am = aT_m;
prdData.Lb = Lw_b;
prdData.Lp = Lw_p;
prdData.Li = Lw_i;

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
end




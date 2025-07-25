function [prdData, info] = predict_Danio_rerio(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
      
  % customized filter  
  filterChecks = 0;
  
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
  TC_BagaPels2001 = tempcorr(temp.tL_BagaPels2001, T_ref, T_A);
  TC_BestAdat2010 = tempcorr(temp.tL_BestAdat2010, T_ref, T_A);
  TC_LawrEber2002 = tempcorr(temp.tL_LawrEber2002_high, T_ref, T_A);


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

  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Lw_b = L_b/ del_Mt;                % cm, total length at birth
  aT_b = t_b/ k_M/ TC_ab;           % d, age at birth 

  % metamorphosis
  L_j = l_j * L_m;                  % cm, structural length at metam at f               
  Lw_j= L_j/ del_Mt;                 % cm, total length at metam at f
  s_M = L_j/ L_b;                   % -, acceleration factor
  aT_j = t_j/ k_M/ TC_aj;           % d, age at metam

  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at F
  Lw_p = L_p/ del_Mt;                % cm, total length at puberty at F
  aT_p = t_p/ k_M/ TC_ap;           % d, time since birth at puberty at F and T

  % ultimate
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_Mt;                % cm, ultimate total length at f
  Ww_i = L_i^3 * (1 + f * w);       % g, ultimate wet weight 

  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector
  RT_i = TC_Ri * reprod_rate_j(L_i, f, pars_R);                 % ultimate reproduction rate
  t_R  = 3; % d, period of accumulaton of reprod buffer at T
  GSIT = (t_R * TC_GSI * k_M * g/ f^3)/ (f + kap * g * y_V_E);
  GSIT = GSIT * ((1 - kap) * f^3 - k_J * U_Hp/ L_m^2/ s_M^3);    % -, GSI

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T

  % puberty at f_EatoFarl1974b
  % F = f_EatoFarl1974;
  % [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, F);
  % L_p = L_m * l_p;                  % cm, structural length at puberty at F
  % Lw_p = L_p/ del_Mt;                % cm, total length at puberty at F
  % aT_p = t_p/ k_M/ TC_ap;           % d, time since birth at puberty at F and T

  % pack to output
  prdData.ab = aT_b;
  prdData.aj = aT_j;
  prdData.ap = aT_p;
  prdData.am = aT_m;
  prdData.Lb = Lw_b;
  prdData.Lj = Lw_j;
  prdData.Lp = Lw_p;
  prdData.Li = Lw_i;
  prdData.Wwi = Ww_i;
  prdData.Ri = RT_i;
  prdData.GSI = GSIT;

  %% uni-variate data
  
    
% compilation of literature growth curves (see previous version and AuguGagn2011 for more discussion on this):
 
 % initial conditions at start:
 init_cond = [1e-10; E_0; 0]; % hardly any structure, initial enegy of the egg (J), no a maturiaty at start

% tL_LawrEber2002 at low and high
  TC = TC_LawrEber2002;
  % Low: 
  F = f_LawrEber2002_low;  a = [0; tL_LawrEber2002_low(:,1)];
  [t_sort, it, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order
  L_j = get_lj(pars_lj, F) * L_m; % cm, structural length at metamorphosis at F
  [tt, LEH] = ode45(@ode_LEH, t_sort, init_cond, [], par, F, TC, L_b, L_j);
  L = LEH(:,1); % cm, structural length
  Lw = L/ del_Mt; Lw = Lw(it_sort); % reconstruct L
  ELw_LawrEber2002_low = Lw(2:end); % cm, total length
  prdData.tL_LawrEber2002_low = ELw_LawrEber2002_low;

  % High:
  F = f_LawrEber2002_high; a = [0; tL_LawrEber2002_high(:,1)];
  [t_sort, it, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order
  L_j = get_lj(pars_lj, F) * L_m; % cm, structural length at metamorphosis at F
  [tt, LEH] = ode45(@ode_LEH, t_sort, init_cond, [], par, F, TC, L_b, L_j);
  L = LEH(:,1); % cm, structural length
  Lw = L/ del_Mt; Lw = Lw(it_sort); % reconstruct L
  ELw_LawrEber2002_high = Lw(2:end); % cm, total length
  prdData.tL_LawrEber2002_high = ELw_LawrEber2002_high;

  % tL_BestAdat2010:
  TC = TC_BestAdat2010;  F = f_BestAdat2010; a = [0; tL_BestAdat2010(:,1)];
  L_j = get_lj(pars_lj, F) * L_m; % cm, structural length at metamorphosis at F
  [t_sort, it, it_sort] = unique(a,'sorted'); % returns the unique values in t in sorted order
  [tt, LEH] = ode45(@ode_LEH, t_sort, init_cond, [], par, F, TC, L_b, L_j);
  L    = LEH(:,1); % cm, structural length
  Lw   = L/ del_Mt; Lw   = Lw(it_sort); % reconstruct L
  ELw_BestAdat2010 = Lw(2:end); % cm, total length
  prdData.tL_BestAdat2010 = ELw_BestAdat2010;

  
  % tL, tWw and tWd of BagaPels2001
  F = f_BagaPels2001; TC = TC_BagaPels2001; 
  a = [0; tL_BagaPels2001(:,1)];
  l_j = get_lj(pars_lj, F);  L_j = l_j * L_m;
  [tt, LEH] = ode45(@ode_LEH, a, init_cond, [], par, F, TC, L_b, L_j);
  LEH(1,:) = []; L    = LEH(:,1);  E    = LEH(:,2); % cm, stuct length; J, energy in resrve
  prdData.tWw_BagaPels2001  = L.^3 + w_E/ mu_E/ d_E * E; % g, wet weight
  prdData.tL_BagaPels2001   = L/ del_Mt; % cm, total length
  

  
  
  
%      
end

%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% used for growth:
function dLEH = ode_LEH(t, LEH, p, f, TC, L_b, L_j)
% Input: 
% p: structure 'par' 
% c: structure 'Cpar' obtained by cPar = parscomp_st(par)
% f: scaled, scaled functional response, 
% s_M: scalar, -, acceleration factor post metamorphosis
% TC, scalar, -, temperature correction factor
% L_b, scaler, cm, structural length at birth at f
% L_j, scaler, cm, structural length at metamorphosis at f


% --------------- unpack LEHR ------------------------------------------
L   =  max(0,LEH(1)); % cm, volumetric structural length
E   =  max(0,LEH(2)); % J,   energy in reserve 
EH  =  max(0,LEH(3)); % J, E_H maturity

% shape correction function:
if EH < p.E_Hb
    s_M = 1;
elseif EH >= p.E_Hb && EH < p.E_Hj
    s_M = L/L_b;
else
    s_M = L_j/L_b;
end
% Temperature and shape correct the relevant paramters
vT    = s_M * p.v * TC; 
pT_Am = s_M * p.z * p.p_M/ p.kap * TC;
pT_M  = p.p_M * TC; 
kT_J  = p.k_J * TC; 
%
pA  = f * pT_Am * L^2 * (EH >= p.E_Hb);           % J/d, assimilation
r   = (E * vT/ L - pT_M * L^3/ p.kap)/ (E + p.E_G * L^3/ p.kap);
pC  = E * (vT/ L - r); % J/d, mobilisation 
dE  = pA - pC;               % J/d, change in energy in reserve
dL  = r/ 3 * L;              % cm/d, change in structural length
dEH = ((1 - p.kap) * pC - kT_J * EH) * (EH < p.E_Hp);    % J/d, change in cum energy invested in maturation (it is implied here that no rejuvenation occurs).

% pack dLEHR
dLEH = [dL; dE; dEH];    
end


function [prdData, info] = predict_Danio_rerio(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
  
%% compute temperature correction factors
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC_ap = tempcorr(temp.ap, T_ref, T_A);
  TC_am = tempcorr(temp.am, T_ref, T_A);
  TC_BestAdat2010 = tempcorr(temp.tL_BestAdat2010, T_ref, T_A);  
  

 %% zero-variate data
  
  % life cycle
  pars_tj = [g k l_T v_Hb v_Hj v_Hp];
  pars_lj = [g; k; l_T; v_Hb; v_Hj];
  
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
  
  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Lw_b = L_b/ del_Mt;                % cm, total length at birth
  aT_b = t_b/ k_M/ TC_ab;           % d, age at birth 

  % puberty (AmP estimation assumes different f, trying with f=1)
  L_p = L_m * l_p;                  % cm, structural length at puberty at F
  Lw_p = L_p/ del_Mt;                % cm, total length at puberty at F
  aT_p = t_p/ k_M/ TC_ap;           % d, time since birth at puberty at F and T

  % ultimate
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_Mt;                % cm, ultimate total length at f

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T

  % pack to output
  prdData.ab = aT_b;
  prdData.ap = aT_p;
  prdData.am = aT_m;
  prdData.Lb = Lw_b;
  prdData.Lp = Lw_p;
  prdData.Li = Lw_i;

  %% uni-variate data
  
  % L-Ww

  EWw = (LWw(:,1) * del_Mt).^3 * (1 + w * f); % g, wet weight

  % tL_BestAdat2010:
  % initial reserve
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  [U_E0, ~, info] = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
  if ~info; prdData=[]; return; end
  E_0 = U_E0 * p_Am ;          % J, energy in egg
  init_cond = [1e-10; E_0; 0]; % hardly any structure, initial enegy of the egg (J), no a maturiaty at start
  % length at metamorphosis 
  L_j = L_m * l_j;

  TC = TC_BestAdat2010;   a = [0; tL_BestAdat2010(:,1)];
  [~, LEH] = ode45(@ode_LEH, a, init_cond, [], par, f, TC, L_b, L_j);
  L    = LEH(:,1); % cm, structural length
  Lw   = L/ del_Mt; 
  ELw_BestAdat2010 = Lw(2:end); % cm, total length
  
  % pack to output
  prdData.LWw = EWw;
  prdData.tL_BestAdat2010 = ELw_BestAdat2010;

  
  
 
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







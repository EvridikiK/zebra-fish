clear
close all; 
global pets 

pets = {'Danio_rerio'}; 

load(['results_' pets{1} '.mat']);
cPar = parscomp_st(par); vars_pull(par); vars_pull(cPar);


%% Functions

function [p_A, p_C, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R, s_M, TC, f, p)

p_Am = TC * p.z * p.p_M / p.kap; v = TC * p.v; p_M = TC * p.p_M; k_J = TC * p.k_J;
kap = p.kap; E_G = p.kap; E_Hb = p.E_Hb; 

% Powers
p_A = f * p_Am .* s_M .* V.^(2/3) .* (E_H >= E_Hb);
p_S = p_M * V;
p_C = E .* (E_G * v ./ V.^(1/3) .* s_M + p_M) ./ (kap * E ./ V + E_G);
% p_C = (E ./ V) .* (E_G .* s_M .* v .* (V.^(2/3)) + p_S) / (kap .* (E ./ V) + E_G);
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

%% Reproduction at ultimate size (empty reproduction buffer)
% Theoretical maximum reproduction buffer
% Get s_M
pars_tj = [g k l_T v_Hb v_Hj v_Hp];
[t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
s_M_max = l_j/l_b;
TC = 1;

V = (L_m*s_M_max)^3;
E = E_m * V;
E_H = E_Hp;
E_R = 0;
E_e = 0;
stateUltimateSizeEmptyReprodBuffer = [V, E, E_H, E_R, s_M_max, E_e];
[p_A, p_C, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R, s_M_max, TC, f, par);
E_R_max_theoretical = p_R / (v / L_m) 

% Simulation maximum reproduction buffer
[simul_t, VEHRsMG] = ode45(@ode_VEHRsMG, [0; 1000], stateUltimateSizeEmptyReprodBuffer, [], par, f, TC);

E_R_max_simul = VEHRsMG(end, 4)

% Maximum reproduction rate
pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
E_0 = U_E0 * p_Am ;          % J, energy in egg
p_Ce_max = E_R_max_theoretical * v / L_m;
R_i_theoretical =  p_Ce_max / E_0

[p_A, p_C, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R_max_theoretical, s_M_max, TC, f, par);
R_i_simul =  p_C2 / E_0

%% Gonado-somatic index (starting with empty reproduction buffer)
t_GSI = 28;
% Simulate evolution for t_GSI days
[~, VEHRsMG] = ode45(@ode_VEHRsMG, [0; t_GSI], stateUltimateSizeEmptyReprodBuffer, [], par, f, TC);
MV = VEHRsMG(end, 1) * M_V * w_V;
ME = VEHRsMG(end, 2) / mu_E * w_E;
ME_R = VEHRsMG(end, 4) / mu_E * w_E;
ME_e = VEHRsMG(end, 6) / mu_E * w_E;

GSIemptyER = (ME_R + ME_e) / (ME + MV + ME_R + ME_e)
GSIemptyERonlyEe = (ME_e) / (ME + MV + ME_R + ME_e)
% GSIemptyER = (ME_R + ME_e) / (ME + MV)
% GSIemptyERonlyEe = (ME_e) / (ME + MV)

%% Gonado-somatic index (starting with maximum reproduction buffer)
% State for GSI computation
V = (L_m*s_M_max)^3;
E = E_m * V;
E_H = E_Hp;
E_R = E_R_max_theoretical;
E_e = 0;
stateUltimateSizeMaxReprodBuffer = [V, E, E_H, E_R, s_M_max, E_e];
% Simulate evolution for t_GSI days
[~, VEHRsMG] = ode45(@ode_VEHRsMG, [0; t_GSI], stateUltimateSizeMaxReprodBuffer, [], par, f, TC);
MV = VEHRsMG(end, 1) * M_V * w_V;
ME = VEHRsMG(end, 2) / mu_E * w_E;
ME_R = VEHRsMG(end, 4) / mu_E * w_E;
ME_e = VEHRsMG(end, 6) / mu_E * w_E;

GSImaxER = (ME_R + ME_e) / (ME + MV + ME_R + ME_e)
GSImaxERonlyEe = (ME_e) / (ME + MV + ME_R + ME_e)
% GSImaxER = (ME_R + ME_e) / (ME + M_V)
% GSImaxERonlyEe = (ME_e) / (ME + M_V)

%% Reproduction dynamics
% From ultimate size, empty reproduction buffer
[tt, VEHRsMG] = ode45(@ode_VEHRsMG, [0; 400], stateUltimateSizeEmptyReprodBuffer, [], par, f, TC);
V = VEHRsMG(:, 1); E = VEHRsMG(:, 2); E_H = VEHRsMG(:, 3); E_R = VEHRsMG(:, 4); s_M = VEHRsMG(:, 5); E_egg = VEHRsMG(:, 6);
[p_A, p_C, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R, s_M, TC, f, par);
figure
plot(tt, p_R, 'b', tt, p_C2, 'r')
title('Simulation from ultimate size, empty reproduction buffer')

% From puberty
V = (l_p * L_m)^3;
E = f * E_m * V;
E_H = E_Hp;
E_R = 0;
E_e = 0;
statePuberty = [V, E, E_Hp, E_R, s_M_max, 0];
[tt, VEHRsMG] = ode45(@ode_VEHRsMG, [0; 1000], statePuberty, [], par, f, TC);
V = VEHRsMG(:, 1); E = VEHRsMG(:, 2); E_H = VEHRsMG(:, 3); E_R = VEHRsMG(:, 4); s_M = VEHRsMG(:, 5); E_egg = VEHRsMG(:, 6);
[p_A, p_C, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R, s_M, TC, f, par);
figure
plot(tt, p_R, 'b', tt, p_C2, 'r')
title('Simulation from puberty')

% Maximum size of reproduction buffer overtime
max_E_R = p_R ./ (v * s_M_max ./ V.^(1/3));
figure
plot(tt, E_R, 'b', tt, max_E_R, 'r')
title('Maximum reproduction buffer and reproduction buffer since puberty')


%% Simulate birth as a function of initial structure
function [value, isterminal, direction] = maturityTransition(t, VEHRsMG, E_H, endSimulation)
    % Event function that triggers when y reaches 0
    value = VEHRsMG(3) - E_H; % The condition to trigger the event
    isterminal = endSimulation; % Stop the integration when the event occurs
    direction = 0; % All directions
end
birthEvent = @(t, VEHRsMG) maturityTransition(t, VEHRsMG, E_Hb, 1);
options = odeset('Events', birthEvent);

% From fertilization
initialStructures = 10.^(-(1:10));
nVs = length(initialStructures);
structureAtBirth = zeros(nVs, 1);
weightAtBirth = zeros(nVs, 1);
ageAtBirth = zeros(nVs, 1);
for i=1:nVs
    V = initialStructures(i);
    stateFertilization = [V, E_0, 0, 0, 1, 0];
    ode = @(t, VEHRsMG) ode_VEHRsMG(t, VEHRsMG, par, f, TC);
    [tt, VEHRsMG, te, VEHRsMG_b, ie ] = ode45(ode, [0; 100], stateFertilization, options);
    V = VEHRsMG_b(1); E = VEHRsMG_b(2); E_H = VEHRsMG_b(3); E_R = VEHRsMG_b(4); s_M = VEHRsMG_b(5); E_egg = VEHRsMG_b(6);
    structureAtBirth(i) = V(end);
    weightAtBirth(i) = 1 * V(end) + (E(end) + E_R(end)) * w_E / mu_E / d_E;
    ageAtBirth(i) = te;
end
figure
semilogx(initialStructures, structureAtBirth, 'b', initialStructures, weightAtBirth, 'r')
title('Effect of initial structure on state at birth')


%% Compare get_tj with differential equation
function [value, isterminal, direction] = allMaturityTransitions(t, VEHRsMG, par)
    % Event function that triggers when y reaches 0
    E_H = VEHRsMG(3);
    value =  [E_H - par.E_Hb; E_H - par.E_Hj; E_H - par.E_Hp]; % The condition to trigger the event
    isterminal = [0; 0; 1]; % Stop the integration when puberty occurs
    direction = [0; 0; 0]; % All directions
end
birthEvent = @(t, VEHRsMG) allMaturityTransitions(t, VEHRsMG, par);
options = odeset('Events', birthEvent);

f = 1;
TC = tempcorr(C2K(28.5), T_ref, T_A)
% From fertilization
V = 1e-10;
stateFertilization = [V, E_0, 0, 0, 1, 0];
ode = @(t, VEHRsMG) ode_VEHRsMG(t, VEHRsMG, par, f, TC);
sol = ode45(ode, [0 1000], stateFertilization, options);

[t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
fprintf('get_tj lengths: L_b %.4e L_j %.4e L_p %.4e s_M %.4e \n', l_b*L_m, l_j*L_m, l_p*L_m, l_j/l_b)
L_b_simul = sol.ye(1, 1)^(1/3);
L_j_simul = sol.ye(1, 2)^(1/3);
L_p_simul = sol.ye(1, 3)^(1/3);
s_M_simul = L_j_simul / L_b_simul;
fprintf('simul  lengths: L_b %.4e L_j %.4e L_p %.4e s_M %.4e \n', L_b_simul, L_j_simul, L_p_simul, s_M_simul)
% fprintf('simul  lengths: L_b %.4e \n', L_b_simul)

a_b_get_tj = t_b / k_M / TC;
a_j_get_tj = t_j / k_M / TC;
a_p_get_tj = t_p / k_M / TC;

fprintf('get_tj ages: a_b %.4e a_j %.4e a_p %.4e \n', a_b_get_tj, a_j_get_tj, a_p_get_tj)
fprintf('simul  ages: a_b %.4e a_j %.4e a_p %.4e \n', sol.xe(1), sol.xe(2), sol.xe(3))

%% Full simulation with get_tj
[tau_b, l_b] = get_tb(pars_tj([1 2 4]), f);
tel_b = [tau_b, f, l_b];
[tvel, t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f, tel_b, (0:1000*k_M)');


%% Initial scaled reserve as a function of temperature
T = 50
TC = tempcorr(C2K(T), T_ref, T_A);

pars_UE0 = [V_Hb / TC; g; k_J * TC; k_M * TC; v * TC]; % compose parameter vector
U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
E_0 = U_E0 * p_Am * TC          % J, energy in egg

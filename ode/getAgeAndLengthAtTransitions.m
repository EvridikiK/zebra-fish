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

function [value, isterminal, direction] = lifeStageTransitions(t, VEHRsMG, par)
% Event function that triggers when y reaches 0
E_H = VEHRsMG(3);
value =  [E_H - par.E_Hb; E_H - par.E_Hj; E_H - par.E_Hp]; % The condition to trigger the event
isterminal = [0; 0; 1]; % Stop the integration when puberty occurs
direction = [0; 0; 0]; % All directions
end

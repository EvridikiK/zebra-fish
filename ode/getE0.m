function [E_0] = getE0(f, par, cpar)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    f (1, 1) double
    par (1, 1) struct
    cpar (1, 1) struct
end

arguments (Output)
    E_0 (1, 1) double
end

pars_UE0 = [cpar.V_Hb; cpar.g; par.k_J; cpar.k_M; par.v]; % compose parameter vector
U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
E_0 = U_E0 * cpar.p_Am;          % J, energy in egg

end
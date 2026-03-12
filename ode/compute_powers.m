function [p_A, p_C, p_S, p_G, p_J, p_R, p_C2] = compute_powers(V, E, E_H, E_R, s_M, TC, f, p)

p_Am = p.z * p.p_M / p.kap; v = p.v; p_M = p.p_M; k_J = p.k_J;
kap = p.kap; E_G = p.E_G; E_Hb = p.E_Hb;

% Powers
p_A = TC .* f * p_Am .* s_M .* V.^(2/3) .* (E_H >= E_Hb);
p_C = TC .* E .* (E_G * v ./ V.^(1/3) .* s_M + p_M) ./ (kap .* E ./ V + E_G);
p_S = TC .* p_M .* V;
p_G = kap .* p_C - p_S;
p_J = TC .* k_J .* E_H;
p_R = (1 - kap) .* p_C - p_J;
p_C2 = TC .* v .* s_M ./ V.^(1/3) .* E_R;

end
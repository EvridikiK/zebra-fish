function max_E_R = get_max_E_R(V, par, E_m, f, s_M_max)
E = f * V * E_m;
E_H = par.E_Hp;
E_R = 0;
TC = 1;
[~, ~, ~, ~, ~, p_R, ~] = compute_powers(V, E, E_H, E_R, s_M_max, TC, f, par);
max_E_R = p_R ./ (par.v * s_M_max ./ V.^(1/3));
end

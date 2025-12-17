function [par, metaPar, txtPar] = pars_init_Danio_rerio(metaData)

metaPar.model = 'abj'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 
p = 1; % estimate parameters
k = 1; % estimate efficiencies and h_a
e = 1; % Estimate food levels?

%% core primary parameters 
par.T_A = 8000;       free.T_A   = 0;   units.T_A = 'K';          label.T_A = 'Arrhenius temp'; 
par.z = 0.3264;       free.z     = p;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 6.5;        free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.3163;   free.kap_X = k;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1000;   free.kap_P = k;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.06714;      free.v     = p;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.4156;     free.kap   = p;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.p_M = 55.5502;    free.p_M   = p;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 5127.2141;  free.E_G   = k;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 2.826e-01; free.E_Hb  = p;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hj = 4.036e+00; free.E_Hj  = p;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; 
par.E_Hp = 9.444e+02; free.E_Hp  = p;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.h_a = 1.936e-08;  free.h_a   = p;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.del_Mt = 0.1857;  free.del_Mt = p;   units.del_Mt = '-';       label.del_Mt = 'shape coefficient for adult'; 
par.V_0 = 1e-10;      free.V_0 = 0;      units.V_0 = 'cm^3';       label.V_0 = 'initial structure of egg';

par.f = 1;                        free.f     = 0;   units.f = '-';            label.f = 'scaled functional response'; 
par.f_BagaPels2001 = 1.0000;      free.f_BagaPels2001 = e;   units.f_BagaPels2001 = '-';  label.f_BagaPels2001 = 'scaled functional response of BagaPels2001 data'; 
par.f_BeauGous2015L = 1.0000;     free.f_BeauGous2015L = e;   units.f_BeauGous2015L = '-';  label.f_BeauGous2015L = 'sc. func. resp., BeauGous2015, length data'; 
par.f_BeauGous2015R = 1.0000;     free.f_BeauGous2015R = e;   units.f_BeauGous2015R = '-';  label.f_BeauGous2015R = 'sc. func. resp., BeauGous2015, reproduction data'; 
par.f_BestAdat2010 = 1.0000;      free.f_BestAdat2010 = e;   units.f_BestAdat2010 = '-';  label.f_BestAdat2010 = 'scaled functional response of BestAdat2010 data'; 
par.f_EatoFarl1974 = 1.0000;      free.f_EatoFarl1974 = e;   units.f_EatoFarl1974 = '-';  label.f_EatoFarl1974 = 'scaled functional response of EatoFarl1974b data'; 
par.f_LawrEber2002_high = 1.0000; free.f_LawrEber2002_high = e;   units.f_LawrEber2002_high = '-';  label.f_LawrEber2002_high = 'scaled functional response of LawrEber2002_high data'; 
par.f_LawrEber2002_low = 1.0000;  free.f_LawrEber2002_low = e;   units.f_LawrEber2002_low = '-';  label.f_LawrEber2002_low = 'scaled functional response of LawrEber2002_low data'; 
par.f_Schi2002 = 1.0000;          free.f_Schi2002 = e;   units.f_Schi2002 = '-';   label.f_Schi2002 = 'scaled functional response of Schi2002 data'; 
par.f_ValKwa2022 = 1.0000;        free.f_ValKwa2022 = e;   units.f_ValKwa2022 = '-';  label.f_ValKwa2022 = 'scaled functional response of ValKwa2022 data'; 
par.f_YangYama2019 = 1.0000;      free.f_YangYama2019 = e;   units.f_YangYama2019 = '-';  label.f_YangYama2019 = 'scaled functional response of YangYama2019 data'; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 

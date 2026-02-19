function [par, metaPar, txtPar] = pars_init_Danio_rerio(metaData)

metaPar.model = 'abj';

%% reference parameter (not to be changed)
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature';

%% core primary parameters
par.T_A = 8000;       free.T_A   = 0;   units.T_A = 'K';          label.T_A = 'Arrhenius temp';
par.z = 0.26412;      free.z     = 1;   units.z = '-';            label.z = 'zoom factor';
par.F_m = 6.5;        free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate';
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve';
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces';
par.v = 0.024217;     free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance';
par.kap = 0.45827;    free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma';
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency';
par.p_M = 117.8926;   free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint';
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint';
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient';
par.E_G = 5204.2857;  free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure';
par.E_Hb = 2.347e-01; free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth';
par.E_Hj = 8.037e+00; free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam';
par.E_Hp = 5.852e+02; free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty';
par.h_a = 2.082e-09;  free.h_a   = 1;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration';
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient';

%% other parameters
par.V_0 = 1e-10;      free.V_0   = 0;   units.V_0 = 'cm^3';       label.V_0 = 'initial structure of egg';
par.del_Mt = 0.12077;  free.del_Mt = 1;   units.del_Mt = '-';       label.del_Mt = 'shape coefficient for adult';
par.f = 1;            free.f     = 0;   units.f = '-';            label.f = 'scaled functional response';
par.f_BagaPels2001 = 0.76793;  free.f_BagaPels2001 = 1;   units.f_BagaPels2001 = '-';  label.f_BagaPels2001 = 'scaled functional response of BagaPels2001 data';
par.f_BestAdat2010 = 1;         free.f_BestAdat2010 = 1;   units.f_BestAdat2010 = '-';  label.f_BestAdat2010 = 'scaled functional response of BestAdat2010 data';
par.f_LawrEber2002_high = 0.80977;  free.f_LawrEber2002_high = 1;   units.f_LawrEber2002_high = '-';  label.f_LawrEber2002_high = 'scaled functional response of LawrEber2002_high data';
par.f_LawrEber2002_low = 0.54302;  free.f_LawrEber2002_low = 1;   units.f_LawrEber2002_low = '-';  label.f_LawrEber2002_low = 'scaled functional response of LawrEber2002_low data';

%% set chemical parameters from Kooy2010
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class);

%% Pack output:
txtPar.units = units; txtPar.label = label; par.free = free;

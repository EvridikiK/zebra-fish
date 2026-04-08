function [data, auxData, metaData, txtData, weights] = mydata_Danio_rerio

%% set metadata
metaData.phylum     = 'Chordata';
metaData.class      = 'Actinopterygii';
metaData.order      = 'Cypriniformes';
metaData.family     = 'Danionidae';
metaData.species    = 'Danio_rerio';
metaData.species_en = 'Zebra fish';

metaData.ecoCode.climate = {'Am'};
metaData.ecoCode.ecozone = {'TPi'};
metaData.ecoCode.habitat = {'0iFp', '0iFm'};
metaData.ecoCode.embryo  = {'Fh'};
metaData.ecoCode.migrate = {};
metaData.ecoCode.food    = {'biCi'};
metaData.ecoCode.gender  = {'D'};
metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(25); % K, body temp
metaData.data_0     = {'ab'; 'aj'; 'ap'; 'am'; 'L0'; 'Lb'; 'Lj'; 'Lp'; 'Li'; 'Wd0'; 'Wwi'; 'Ri'; 'GSI'};
metaData.data_1     = {'t-Le'; 't-Wwe'; 't-Wde'; 't-MCe'; 't-MNe'; 't-L_fT'; 't-Ww_f'; 't-N'; 't-S'; 'L-Ww'};

metaData.COMPLETE = 5; % using criteria of LikaKear2011

metaData.author   = {'Evridiki Klagkou', 'Tjui Yuew Tan', 'María-José Lagunes', 'Diogo F. Oliveira'};
metaData.email    = {'evridiki13@hotmail.gr', 'tan.tjuiyeuw@wur.nl', 'mj.laguneslopez@gmail.com', 'diogo.miguel.oliveira@tecnico.ulisboa.pt'};

%% load externalized calibration data
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'data_pipeline'));

%% set data
% zero-variate data
zeroVariateSpecs = { ...
    'ab', ...   % age at birth, BestAdat2010
    'aj', ...   % age at metam, BestAdat2010
    'ap', ...   % age at puberty, EatoFarl1974a
    'am', ...   % life span, GerhKauf2002
    'L0', ...   % egg diameter, ForbPres2010
    'Lb', ...   % total length at birth, Schi2002
    'Lj', ...   % total length at metam, Schi2002
    'Lp', ...   % standard length at puberty, EatoFarl1974a
    'Li', ...   % ultimate total length, SpenGerl2008 and Schi2002
    'Wd0', ...  % egg dry weight, Augu2011
    'Wwi', ...  % ultimate wet weight, Augu2011
    'Ri', ...   % max reproduction, EatoFarl1974b
    'GSI' ...   % Gonado Somatic Index, ForbPres2010
    };

% tim-length, larval growth curve T = 28.5 + 273 K, total length
schi2002Specs = { ...
    'tL_Schi2002' ... % time-length data, Schi2002
    };

% time-length at T = 25.5 + 273 K, standard length
eatoFarl1974Specs = { ...
    'tL_EatoFarl1974' ... % time-length data, EatoFarl1974b
    };

% time-length, wet-weight, and dry-weight data at T = 25 C
bagaPels2001Specs = { ...
    'tL_BagaPels2001', ... % age vs standard length, BagaPels2001
    'tWw_BagaPels2001', ... % age vs wet weight, BagaPels2001
    'tWd_BagaPels2001' ... % age vs dry weight, BagaPels2001
    };

% BestAdat2010 T = 25 C, larval time-length data digitized from the study
bestAdat2010Specs = { ...
    'tL_BestAdat2010' ... % larval time-length data, BestAdat2010
    };

% time-length from LawrEber2008 at 28.5 C under high and low food conditions
lawrEber2008Specs = { ...
    'tL_LawrEber2008_high', ... % high-food growth curve, LawrEber2008
    'tL_LawrEber2008_low' ... % low-food growth curve, LawrEber2008
    };

% BeauGous2015 fasting treatments and reproduction summaries
beauGous2015Specs = { ...
    'tLf1_BeauGous2015', ... % no fasting length trajectory, BeauGous2015
    'tLf2_BeauGous2015', ... % fasting every 3 days length trajectory, BeauGous2015
    'tLf3_BeauGous2015', ... % fasting every other day length trajectory, BeauGous2015
    'tL1', ... % juvenile standard-length series from courtesy raw data
    'Wwt', ... % final wet weight summary from courtesy raw data
    'tN' ... % cumulative egg production from courtesy raw data
    };

% Valentine and Kwasek feeding-rate experiment
valKwa2022Specs = { ...
    'tWw_ValKwa2022', ... % age vs wet weight, ValKwa2022
    'tL_ValKwa2022', ... % age vs total length, ValKwa2022
    'tJX_ValKwa2022' ... % age vs feed intake, ValKwa2022
    };

% Yang and Yamamoto growth, body mass, and oxygen-consumption data
yangYama2019Specs = { ...
    'tL_YangYama2019', ... % age vs total length, YangYama2019
    'tWw_YangYama2019', ... % age vs wet weight, YangYama2019
    'tJO_YangYama2019' ... % age vs oxygen consumption, YangYama2019
    };

datasetSpecs = [ ...
    zeroVariateSpecs, ...
    schi2002Specs, ...
    eatoFarl1974Specs, ...
    bagaPels2001Specs, ...
    bestAdat2010Specs, ...
    lawrEber2008Specs, ...
    beauGous2015Specs, ...
    valKwa2022Specs, ...
    yangYama2019Specs];

[data, auxData, externalTxtData] = load_calibration_data(datasetSpecs);

units = externalTxtData.units;
label = externalTxtData.label;
bibkey = externalTxtData.bibkey;
comment = externalTxtData.comment;
temp = auxData.temp;
treat = auxData.treat;
init = auxData.init;

%% set weights for all real data
weights = setweights(data, []);

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% Auto add titles to univariate data
title = struct();
datanames = fieldnames(data);
for i=1:numel(datanames)
  datum = datanames{i};
  % Skip zero-variate data
  if isscalar(data.(datum))
    continue
  end
  if ~isfield(title, datum)
    dataTitle = [strjoin(label.(datum), ' vs ') ', ' strrep(datum, '_', '\_');];
    if isfield(bibkey, datum)
      dataTitle = [dataTitle ', ' bibkey.(datum){:}];
    end
    title.(datum) = dataTitle;
  end
end

%% pack auxData and txtData for output
auxData.temp = temp;
auxData.treat = treat;
auxData.init = init;

txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;
txtData.title = title;

%% Group plots
set3 = {'tL_LawrEber2008_high','tL_LawrEber2008_low'}; subtitle3 = {'LawrEber2008 data at high, low food'};
set5 = {'tLf1_BeauGous2015','tLf2_BeauGous2015','tLf3_BeauGous2015'}; subtitle5 = {'no fasting, fasting every 3 d, fasting every other day'};
metaData.grp.sets = {set3,  set5};
metaData.grp.subtitle = {subtitle3,  subtitle5};

%% Discussion points
D1 = '';
D2 = '';
metaData.discussion = struct('D1', D1, 'D2', D2);

%% Links
metaData.links.id_CoL = '3443J'; % Cat of Life
metaData.links.id_ITIS = '163668'; % ITIS
metaData.links.id_EoL = '204011'; % Ency of Life
metaData.links.id_Wiki = 'Danio_rerio'; % Wikipedia
metaData.links.id_ADW = 'Danio_rerio'; % ADW
metaData.links.id_Taxo = '172875'; % Taxonomicon
metaData.links.id_WoRMS = '1026595'; % WoRMS
metaData.links.id_fishbase = 'Danio-rerio'; % fishbase

%% References
metaData.biblist = build_biblist_from_bibkeys(txtData.bibkey);

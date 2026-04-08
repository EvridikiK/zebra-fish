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

metaData.COMPLETE = 2.5; % using criteria of LikaKear2011

metaData.author   = {'Evridiki Klagkou', 'Tjui Yuew Tan', 'María-José Lagunes', 'Diogo F. Oliveira'};
metaData.email    = {'evridiki13@hotmail.gr', 'tan.tjuiyeuw@wur.nl', 'mj.laguneslopez@gmail.com', 'diogo.miguel.oliveira@tecnico.ulisboa.pt'};


%% load externalized calibration data
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'data_pipeline'));

%% set data
% zero-variate data
zeroVariateIds = { ...
    'ab', ...   % age at birth, BestAdat2010
    'ap', ...   % age at puberty, EatoFarl1974a
    'am', ...   % life span, GerhKauf2002
    'Lb', ...   % total length at birth, Schi2002
    'Lp', ...   % standard length at puberty, EatoFarl1974a
    'Li', ...   % ultimate total length, SpenGerl2008 and Schi2002
    'Wwi', ...  % ultimate wet weight, Augu2011
    'Ri' ...    % max reproduction, EatoFarl1974b
    };

  
% uni-variate data
% BestAdat2010 T = 25 C, rotifers from day 5 till 9, then transition to
% regular feeding from 9 till 12 was originally given in personnal
% communication with C. Lawrence
% I digitalized all the points which in fact
% correspond to 2 different strains (AB and nacre)
% time-length
bestAdat2010Ids = { ...
    'tL_BestAdat2010' ... % larval time-length data, BestAdat2010
    };

datasetIds = [zeroVariateIds, bestAdat2010Ids];
[data, auxData, externalTxtData] = load_calibration_data(datasetIds);

units = externalTxtData.units;
label = externalTxtData.label;
comment = externalTxtData.comment;
bibkey = externalTxtData.bibkey;
temp = auxData.temp;

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
        dataTitle = [strjoin(label.(datum), ' vs ') ', ' datum];
        if isfield(bibkey, datum)
            dataTitle = [dataTitle ', ' bibkey.(datum){:}];
        end
        title.(datum) = dataTitle;
    end
end

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;
txtData.title = title;

%% Group plots

%% Discussion points
D1 = 'version 2018 08 08: KimBall95 and Augu2011 no longer included';
D2 = 'version 2018 08 08: studies by GomeConc2010, SchaRyan2006, BarrFern2010, BarrBurg1999, are not longer included. Previous version included them but gave them zero weight because growth deviates from what is assumed to be the normal patterns of growth at constant food (see EatoFarly1974,LawrEber2008,  Schi2002, BestAdat2010 and discussion in AuguGagn2011)';     
D3 = 'version 2018 08 08: it is no longer possible to implement different Arrhennius temperatures for embryo and adult';
D4 = 'version 2018 08 08: Egg respiration data from BangGron2004 is not longer included, the data show that the respiration stops quite early on while development and growth is very fast. We think this is an artefact. You can find back the data in the previous version.';
D5 = 'version 2018 08 08: inclusion of juvenile growth and adult reproduction data from BeauGous2015.';
D6 = 'version 2018 08 08: Buffer handling rules from AuguGagn2012 are used for modelling starvation response. We  estimate at adult females have  200 and 800 J in the reproduction buffer at the start of the trials. This is in line with finding from AuguGang2012'; 
D7 = 'version 2018 08 08: standard length is 80% of the total length:';
metaData.discussion = struct('D1',D1,'D2',D2,'D3',D3,'D4',D4,'D5',D5, 'D6', D6, 'D7', D7);
                                 
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


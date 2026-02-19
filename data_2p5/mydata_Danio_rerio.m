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

metaData.author   = {'Starrlight Augustine'};        
metaData.date_subm = [2011 05 02];                           
metaData.email    = {'sta@akvaplan.niva.no'};                 
metaData.address  = {'Akvaplan-niva AS, Fram Centre, P.O. Box 6606 Langnes, 9296 Tromso, Norway'}; 

metaData.author_mod_1   = {'Starrlight Augustine'};        
metaData.date_mod_1 = [2018 08 08];                           
metaData.email_mod_1    = {'starrlight.augustine@akvaplan.niva.no'};                 
metaData.address_mod_1  = {'Akvaplan-niva AS, Fram Centre, P.O. Box 6606 Langnes, 9296 Tromso, Norway'}; 

metaData.curator     = {'Bas Kooijman'};
metaData.email_cur   = {'bas.kooijman@vu.nl'}; 
metaData.date_acc    = [2018 08 09]; 

%% set data
% zero-variate data
data.ab = 5;      units.ab = 'd';    label.ab = 'age at birth';           bibkey.ab = 'BestAdat2010'; 
  temp.ab = C2K(25);  units.temp.ab = 'K'; label.temp.ab = 'temperature';
data.ap = 75;     units.ap = 'd';    label.ap = 'age at puberty';         bibkey.ap = 'EatoFarl1974a'; 
  temp.ap = C2K(25.5);  units.temp.ap = 'K'; label.temp.ap = 'temperature';
  comment.ap = '(74 - 75 day of development)';
data.am = 4.5 * 365; units.am = 'd'; label.am = 'life span';              bibkey.am = 'GerhKauf2002';   
  temp.am = C2K(26); units.temp.am = 'K'; label.temp.am = 'temperature';

data.Lb = 4e-1;   units.Lb = 'cm';   label.Lb = 'total length at birth';  bibkey.Lb = 'Schi2002';
  comment.Lb = 'Swim bladder inflates, active feeding, pronephric tubules';
data.Lp = 2.6*1.25; units.Lp = 'cm'; label.Lp = 'total length at puberty';bibkey.Lp = 'EatoFarl1974a';
  comment.Lp = 'Female are between 2.4 - 2.6 cm SL in the study at first egg laying';
data.Li = 5;      units.Li = 'cm';   label.Li = 'ultimate total length';  bibkey.Li = {'SpenGerl2008','Schi2002'}; 
  comment.Li = 'also Lawr pers. comm';

data.Wwi = 1;      units.Wwi = 'g';    label.Wwi = 'ultimate wet weight';    bibkey.Wwi = 'Augu2009';

data.Ri = 240;    units.Ri = '#/d';  label.Ri = 'max reproduction';       bibkey.Ri = {'EatoFarl1974b'}; 
  temp.Ri = C2K(26);  units.temp.Ri = 'K'; label.temp.Ri = 'temperature';

  
% uni-variate data

% BestAdat2010 T = 25 C, rotifers from day 5 till 9, then transition to
% regular feeding from 9 till 12 was originally given in personnal
% communication with C. Lawrence
% I digitalized all the points which in fact
% correspond to 2 different strains (AB and nacre)
% time-length
data.tL_BestAdat2010 = [... age (dpf), total length (cm)
5.00741992	0.30707071;
5.01034883	0.42828283;
5.97532950	0.36363636;
5.97698922	0.43232323;
6.99107806	0.40000000;
6.99215200	0.44444444;
7.96045210	0.51717172;
8.00731478	0.45656566;
8.97571251	0.53333333;
9.94401261	0.60606061;
10.95819908	0.57777778;
10.95946828	0.63030303;
11.97404528	0.61818182;
11.97648604	0.71919192;
16.03869925	0.83232323;
16.03996845	0.88484848;
23.10393150	1.22424242;
23.10481018	1.26060606;
30.17170214	1.72121212;
90 3.19;
90 3.41];
units.tL_BestAdat2010 = {'d', 'cm'}; label.tL_BestAdat2010 = {'time since fertilization', 'total length'}; 
temp.tL_BestAdat2010 = C2K(25); units.temp.tL_BestAdat2010 = 'K'; label.temp.tL_BestAdat2010 = 'temperature';
bibkey.tL_BestAdat2010 = {'BestAdat2010'};

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

%
bibkey = 'Augu2009'; type = 'Misc'; bib = [ ...
'year = {2009}, ' ...
'author = {S. Augustine}, ' ...
'note = {unpublished 2009, f = 1}'];  
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

%
bibkey = 'BestAdat2010'; type = 'Article'; bib = [ ...    
'author = {Best, J. and Adatto, I. and Cockington, J. and James, A. and Lawrence, C.}, ' ...
'year  = {2010}, ' ...
'title = {A novel method for rearing first feeding larval zebrafish: polyculture with {T}ype {L} saltwater rotifers (\emph{Brachionus plicatilis}).}, ' ...  
'journal = {Zebrafish}, ' ...
'volume = {352 7(3)}, ' ...
'pages = {89-295}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

%
bibkey = 'EatoFarl1974a'; type = 'Article'; bib = [ ...    
'author = {Eaton, R. C. and Farley, R. D.}, ' ...
'year  = {1974}, ' ...
'title = {Growth and the reduction of depensation of zebrafish, \emph{Brachydanio rerio}, reared in the laboratory}, ' ...  
'journal = {Copeia}, ' ...
'volume = {1}, ' ...
'pages = {204-209}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

%
bibkey = 'EatoFarl1974b'; type = 'Article'; bib = [ ...    
'author = {Eaton, R. C. and Farley, R. D.}, ' ...
'year  = {1974}, ' ...
'title = {Spawning cycle and egg production of zebrafish, \emph{Brachydanio rerio} in the laboratory}, ' ...  
'journal = {Copeia}, ' ...
'volume = {1}, ' ...
'pages = {195-204}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

%
bibkey = 'GerhKauf2002'; type = 'Article'; bib = [ ...    
'author = {Gerhard, G. S. and Kauffman, E. J. and Wang, X. and Stewart, R. and Moore, J. L. and Kasales, C. J. and Demidenko, E. and Cheng, K. C.}, ' ...
'year  = {2002}, ' ...
'title = {Life spans and senescent phenotypes in two strains of Zebrafish (\emph{Danio rerio})}, ' ...  
'journal = {Experimental Gerontology}, ' ...
'volume = {37(8-9)}, ' ...
'pages = {1055-1068}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

%
bibkey = 'Schi2002'; type = 'Incollection'; bib = [ ...    
'author = {Schilling, T. F.}, ' ...
'year  = {2002}, ' ...
'booktitle = {Zebrafish: a practical guide}, ' ...
'editor = {Nusslein-Volhard, C. and Dahm, R.}, ' ...
'title = {The morphology of larval and adult zebrafish}, ' ...  
'publisher = {Oxford University Press Inc.}, ' ...
'address = {New York}, ' ...
'pages = {59-83}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

%
bibkey = 'SpenGerl2008'; type = 'Article'; bib = [ ...  
'author = {Spence, R. and Gerlach, G. and Lawrence, C. and Smith, C.}, ' ...
'year = {2008}, ' ...
'title = {The behaviour and ecology of the zebrafish, \emph{Danio rerio}}, ' ... 
'journal = {Biol. Rev.}, ' ...
'volume = {83}, '...
'pages = {13-34}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];


function [msolosyls, mduetsyls, fsolosyls, fduetsyls, spon] = wData
% This function only returns the identities of syllables for use in the
% analysis of the wren chronic data.  Data ChronicCompleat2018d.mat


%% List of Chronic singing data with syllable indices and locations for spontaneous activity

% 1-2: m17

    %msolosyls{1} = [1 2 3 4]; % Orig
        msolosyls{1} = [1 2 3]; % Solitary
    mduetsyls{1} = [5 7 9 11 13];
    fsolosyls{1} = []; 
    fduetsyls{1} = [6 8 10 12 14];
    spon(:,1) = [-5.5, -0.5]; % This is a mess. 
    
% 3-4: j160806

    %msolosyls{2} = 2; % Orig
        msolosyls{2} = []; % Solitary
    mduetsyls{2} = [4 6 8 10 12];
    fsolosyls{2} = 1; 
    fduetsyls{2} = [3 5 7 9 11 13];    
    spon(:,2) = [-5.0, 0.0]; 
    
% 5-6: j160807
    
    msolosyls{3} = []; 
    mduetsyls{3} = [2 4 6 8 10 12 14 16 18];
    fsolosyls{3} = 1; % Solitary
    fduetsyls{3} = [3 4 7 9 11 13 15 17 19];
    spon(:,3) = [-5.0, 0.0];
    
% 7-8: j160815
    
    msolosyls{4} = []; 
    mduetsyls{4} = [3 5 7];
    %fsolosyls{4} = [1 2]; % Orig
        fsolosyls{4} = 1; % Solitary
    fduetsyls{4} = [4 6 8];
    spon(:,4) = [-2.5, -1.0];

% 9-10: j161009
    
    msolosyls{5} = []; 
    mduetsyls{5} = [3 5 7 9 11 13 15 17 19];
    %fsolosyls{5} = [1 2]; % Orig
        fsolosyls{5} = 1; % Solitary
    fduetsyls{5} = [4 6 8 10 12 14 16 18];
    spon(:,5) = [-5.0, 0.0];

% 11-12: j161022
    
    msolosyls{6} = []; 
    mduetsyls{6} = [4 6 8 10 12 14];
    %fsolosyls{6} = [1 2 3]; % Orig
        fsolosyls{6} = [1 2]; % Solitary
    fduetsyls{6} = [5 7 9 11 13];
    spon(:,6) = [-3.5, 0.0];

% 13-14: j17060848
    
    %msolosyls{7} = [1 2 3 4 5]; % Orig
        msolosyls{7} = [1 2 3 4]; % Solitary
    mduetsyls{7} = [7 9 11 13 15 17];
    fsolosyls{7} = []; 
    fduetsyls{7} = [6 8 10 12 14 16];
    spon(:,7) = [-3.0, 0.0];

% 15-16: j170081733 

% This is a special case as the birds sang at the same time for two
% syllables.  Syllables 6 and 11 are overlapping.  They are omitted
% from the list below because I don't have a separate category for 
% them. This is good for general analysis of autogenous and heterogenous
% contributions. Use a custom analysis to example what happened with 
% the overlapping syllables and the associated intervals.

    msolosyls{8} = []; 
    mduetsyls{8} = [3 8];
    fsolosyls{8} = 1; 
    fduetsyls{8} = [2 4 5 7 9 10];
    spon(:,8) = [-5, 0];
    
% 17-18: j160734 

% This a sequence with just one male syllable and varying female amplitudes.

    msolosyls{9} = []; 
    mduetsyls{9} = 2;
    fsolosyls{9} = 1; 
    fduetsyls{9} = [3 4];
%    spon(:,9) = [0, 1.5]; % Both of the spontaneous ranges are not ideal
    spon(:,9) = [0, 6.5]; % Both of the spontaneous ranges are not ideal
    
% 19-20: j160800a 

% This just one female syllable.

    msolosyls{10} = []; 
    mduetsyls{10} = [];
    fsolosyls{10} = 1; 
    fduetsyls{10} = [];
    spon(:,10) = [0, 2]; 

% 21-22: j160800b 

% Another solo female syllable.

    msolosyls{11} = []; 
    mduetsyls{11} = [];
    fsolosyls{11} = 1; 
    fduetsyls{11} = [];
    spon(:,11) = [0, 1.5]; 
    
% 23-24: j160928a 

% Female solo syllable.

    msolosyls{12} = []; 
    mduetsyls{12} = [];
    fsolosyls{12} = 1; 
    fduetsyls{12} = [];
    spon(:,12) = [0, 2]; 
    
% 25-26: j160928b 

% Female solo syllable.

    msolosyls{13} = []; 
    mduetsyls{13} = [];
    fsolosyls{13} = 1; 
    fduetsyls{13} = [];
    spon(:,13) = [0, 2]; 
    
% 27-28: j160934a 

% Female solo syllable.

    msolosyls{14} = []; 
    mduetsyls{14} = [];
    fsolosyls{14} = 1; 
    fduetsyls{14} = [];
    spon(:,14) = [0, 2]; 
    
    
    
    
    

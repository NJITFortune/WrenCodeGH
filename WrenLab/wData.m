function [msolosyls, mduetsyls, fsolosyls, fduetsyls, Cspon, Aspon] = wData
% Usage [msolosyls, mduetsyls, fsolosyls, fduetsyls, Cspon, Aspon] = wData
% This function only returns the identities of syllables for use in the
% analysis of the wren chronic data.  Data ChronicCompleat2018d.mat


%% List of Chronic singing data with syllable indices and locations for spontaneous activity

% 1-2: m17 SLOPE

    msolosyls{1} = [1 2 3 4 5]; % Low Amplitude
    mduetsyls{1} = [7 9 11 13];
    %   msolosyls{1} = [1 2 3 4]; % Solitary
    %   mduetsyls{1} = [5 7 9 11 13];
    fsolosyls{1} = []; 
    fduetsyls{1} = [6 8 10 12 14];
%    Cspon(:,1) = [-5.5, -0.5]; % This is a mess. 
    Cspon(:,1) = [-5, -1]; % This is a mess. 
    Aspon(:,1) = [-5, -1];
    
% 3-4: j160806

    msolosyls{2} = 2; % Low Amplitude
    mduetsyls{2} = [4 6 8 10 12];
    %   msolosyls{2} = []; % Solitary
    %   mduetsyls{2} = [2 4 6 8 10 12];
    mduetsyls{2} = [4 6 8 10 12];
    fsolosyls{2} = 1; 
    fduetsyls{2} = [3 5 7 9 11 13];    
    Cspon(:,2) = [-5.0, -1]; 
    Aspon(:,2) = [-5.0, -1];
    
% 5-6: j160807
    
    msolosyls{3} = 2; 
    mduetsyls{3} = [4 6 8 10 12 14 16 18];
    fsolosyls{3} = 1; % Solitary
    fduetsyls{3} = [3 5 7 9 11 13 15 17 19];
    Cspon(:,3) = [-5.0, -1];
    Aspon(:,3) = [-5.0, -1];
    
% 7-8: j160815
    
    msolosyls{4} = []; 
    mduetsyls{4} = [3 5 7];
    % fsolosyls{4} = [1 2]; % Orig
    % fduetsyls{4} = [4 6 8];
    fsolosyls{4} = 1; % Solitary
    fduetsyls{4} = [2 4 6 8];
    Cspon(:,4) = [4.5, 7.0];
    Aspon(:,4) = [-3.0, 0.0];

% 9-10: j161009
    
    msolosyls{5} = []; 
    mduetsyls{5} = [3 5 7 9 11 13 15 17 19];
    fsolosyls{5} = [1 2]; % Orig
    %    fsolosyls{5} = 1; % Solitary
    fduetsyls{5} = [4 6 8 10 12 14 16 18];
    Cspon(:,5) = [-4.0, 0.0];
    Aspon(:,5) = [-4.0, 0.0];

% 11-12: j161022
    
    msolosyls{6} = []; 
    mduetsyls{6} = [4 6 8 10 12 14];
    %   fsolosyls{6} = [1 2 3]; % Orig
    %   fduetsyls{6} = [5 7 9 11 13];
    fsolosyls{6} = [1 2]; % Solitary
    fduetsyls{6} = [3 5 7 9 11 13];
    Cspon(:,6) = [8, 12.0];
    Aspon(:,6) = [-4.0, 0.0];

% 13-14: j17060848
    
    %   msolosyls{7} = [1 2 3 4 5]; % Orig
    %   mduetsyls{7} = [7 9 11 13 15 17];
    msolosyls{7} = [1 2 3 4]; % Solitary
    mduetsyls{7} = [5 7 9 11 13 15 17];
    fsolosyls{7} = []; 
    fduetsyls{7} = [6 8 10 12 14 16];
    Cspon(:,7) = [11.0, 15.0];

% 15-16: j170081733 

% This is a special case as the birds sang at the same time for two
% syllables. Syllables 6 and 11 are overlapping.  They are omitted
% from the list below because I don't have a separate category for 
% them. This is good for general analysis of autogenous and heterogenous
% contributions. Use a custom analysis to example what happened with 
% the overlapping syllables and the associated intervals.

    msolosyls{8} = []; 
    mduetsyls{8} = [3 8];
    fsolosyls{8} = 1; 
    fduetsyls{8} = [2 4 5 7 9 10];
    Cspon(:,8) = [-5, -1];
    
% 17-18: j160734 

% This a sequence with just one male syllable and varying female amplitudes.

    msolosyls{9} = 2; % Male sang syllable 2 but I don't think that they were quite duetting yet
    mduetsyls{9} = []; % Very low amplitude
    fsolosyls{9} = 1; 
    fduetsyls{9} = [3 4]; % Female sang loudly after male, syllables 3 and 4
    Cspon(:,9) = [0, 1.5]; % Both of the spontaneous ranges are not ideal
%    spon(:,9) = [0, 6.5]; % Both of the spontaneous ranges are not ideal


% DATA FROM HERE ARE FROM ONE PAIR. WE WERE WORRIED ABOUT A STRONG PSUEDO-REPLICATION
% ISSUE, AND SO WE DID THE ANALYSES WITHOUT THESE DATA AND GET THE SAME
% MAJOR RESULTS.

% 19-20: j160800a 

% This just one female syllable.

    msolosyls{10} = []; 
    mduetsyls{10} = [];
    fsolosyls{10} = 1; 
    fduetsyls{10} = [];
    Cspon(:,10) = [0, 2]; 

% 21-22: j160800b 

% Another solo female syllable.

    msolosyls{11} = []; 
    mduetsyls{11} = [];
    fsolosyls{11} = 1; 
    fduetsyls{11} = [];
    Cspon(:,11) = [0, 1.5]; 
    
% 23-24: j160928a 

% Female solo syllable.

    msolosyls{12} = []; 
    mduetsyls{12} = [];
    fsolosyls{12} = 1; 
    fduetsyls{12} = [];
    Cspon(:,12) = [0, 2]; 
    
% 25-26: j160928b 

% Female solo syllable.

    msolosyls{13} = []; 
    mduetsyls{13} = [];
    fsolosyls{13} = 1; 
    fduetsyls{13} = [];
    Cspon(:,13) = [0, 2]; 
    
% 27-28: j160934a 

% Female solo syllable.

    msolosyls{14} = []; 
    mduetsyls{14} = [];
    fsolosyls{14} = 1; 
    fduetsyls{14} = [];
    Cspon(:,14) = [0, 2]; 
    
% 29-30: j160934b 

% Female solo syllable.

    msolosyls{15} = []; 
    mduetsyls{15} = [];
    fsolosyls{15} = 1; 
    fduetsyls{15} = [];
    Cspon(:,15) = [0, 2]; 
    
% 31-32: j160934c 

% Female solo syllable.

    msolosyls{16} = []; 
    mduetsyls{16} = [];
    fsolosyls{16} = 1; 
    fduetsyls{16} = [];
    Cspon(:,16) = [3, 5]; 
    
% 33-34: j160934d 

% Female solo syllable.

    msolosyls{17} = []; 
    mduetsyls{17} = [];
    fsolosyls{17} = 1; 
    fduetsyls{17} = [];
    Cspon(:,17) = [0, 2]; 
    
% 35-36: j160951a 

% Female solo syllable.

    msolosyls{18} = []; 
    mduetsyls{18} = [];
    fsolosyls{18} = 1; 
    fduetsyls{18} = [];
    Cspon(:,18) = [2.5, 4.0]; 

 % 37-38: j160951b 

% Female solo syllable.

    msolosyls{19} = []; 
    mduetsyls{19} = [];
    fsolosyls{19} = 1; 
    fduetsyls{19} = [];
    Cspon(:,19) = [0, 2]; 
   
% 39-40: j160951c 

% Female solo syllable.

    msolosyls{20} = []; 
    mduetsyls{20} = [];
    fsolosyls{20} = 1; 
    fduetsyls{20} = [];
    Cspon(:,20) = [0, 2]; 
    
% 41-42: j160951d 

% Female solo syllable.

    msolosyls{21} = []; 
    mduetsyls{21} = [];
    fsolosyls{21} = 1; 
    fduetsyls{21} = [];
    Cspon(:,21) = [3.5, 5.5]; 

    
    


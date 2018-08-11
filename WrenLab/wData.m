function [msolosyls, mduetsyls, fsolosyls, spon] = wData
% This function only returns the identities of syllables for use in the
% analysis of the wren chronic data.


%% List of Chronic singing data with syllable indices and locations for spontaneous activity

% 1-2: m17

    msolosyls{1} = [1 2 3 4]; 
    mduetsyls{1} = [5 7 9 11 13];
    fsolosyls{1} = []; 
    fduetsyls{1} = [6 8 10 12 14];
    spon(:,1) = [-5.5, -0.5]; % This is a mess. 
    
% 3-4: j160806

    msolosyls{2} = [2]; 
    mduetsyls{2} = [4 6 8 10 12];
    fsolosyls{2} = [1]; 
    fduetsyls{2} = [3 5 7 9 11 13];    
    spon(:,2) = [-5.0, 0.0]; 
    
% 5-6: j160807
    
    msolosyls{3} = []; 
    mduetsyls{3} = [3 6 8 11 13 16 18 21 23];
    fsolosyls{3} = [1 2]; 
    fduetsyls{3} = [4 5 7 9 10 13 14 15 17 19 20 22 24];
    spon(:,3) = [-5.0, 0.0];
    
% 7-8: j160815
    
    msolosyls{4} = []; 
    mduetsyls{4} = [3 5 7];
    fsolosyls{4} = [1 2]; 
    fduetsyls{4} = [4 6 8];
    spon(:,4) = [-2.5, -1.0];

% 9-10: j161009
    
    msolosyls{5} = []; 
    mduetsyls{5} = [3 5 7 9 11 14 16 19 21];
    fsolosyls{5} = [1 2]; 
    fduetsyls{5} = [4 6 8 10 12 13 15 17 18 20];
    spon(:,5) = [-5.0, 0.0];

% 11-12: j161022
    
    msolosyls{6} = []; 
    mduetsyls{6} = [4 6 8 10 12 14];
    fsolosyls{6} = [1 2 3]; 
    fduetsyls{6} = [5 7 9 11 13];
    spon(:,6) = [-3.5, 0.0];

% 13-14: j17060848
    
    msolosyls{7} = [1 2 3 4 5]; 
    mduetsyls{7} = [8 10 13 15 18 20];
    fsolosyls{7} = []; 
    fduetsyls{7} = [6 7 9 11 12 14 17 19];
    spon(:,7) = [-3.0, 0.0];

% 15-16: j170081733 

%     msolosyls{8} = []; 
%     mduetsyls{8} = [3 8];
%     fsolosyls{8} = 1; 
%     fduetsyls{8} = [2 4 5 7 9 10];
%     spon(:,8) = [-5, 0];

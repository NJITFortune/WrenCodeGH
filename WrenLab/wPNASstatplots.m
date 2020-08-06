function [out, sumdat, stts] = wPNASstatplots(in, whichlist, uc, padding)
% Usage: [out, sumdat, stts] = wRS_Chronic(w, padding)
% Calculates response strength to solo and duet syllables.
% Load the Chronic data structure first:
% load ChronicCompleat2019f.mat (Current as of 8-Jan-2020)
% in is the structure from that file with all of the data (w)
%
% uc is 1 for urethane and 2 for chronic
%
% 'padding' is a critical variable - it is the shift in the window around the
% clicked boundaries of the syllables for the calculation of RS. 
% The boundaries of each syllable are used for this analysis. This is
% inherently problematic as we expect pre-motor activity in awake animals 
% to occur PRIOR to the sound and auditory feedback to occur AFTER the sound.
% We are comfortable with using a value of '0' it is a rather unbiased, which
% is the default. Change the value of padding in seconds (e.g. 0.020 or 0.150) 
% to move the premotor earlier and the auditory feedback windows later. 
%
% This relies on the function wData.m. 

%% Preparations

pad = 0.050; 

% The user can specify the padding via an argin for convenience.
if nargin == 4; pad = padding; end

% List of Chronic singing data with syllable indices and locations for spontaneous activity

[msolosyls, mduetsyls, fsolosyls, fduetsyls, Cspon, ~, birdlist] = wData;

% Initializing variables

% Normalized and Raw response strengths for DUET syllables: Female and Male, Autogenous and Heterogenous 
    sumdat.fDuetAuto.rsNorm = []; sumdat.fDuetAuto.rsRaw = [];
    sumdat.mDuetAuto.rsNorm = []; sumdat.mDuetAuto.rsRaw = [];
    sumdat.fDuetHetero.rsNorm = []; sumdat.fDuetHetero.rsRaw = [];
    sumdat.mDuetHetero.rsNorm = []; sumdat.mDuetHetero.rsRaw = [];

% Normalized and Raw response strengths for SOLO syllables: Female and Male, Autogenous and Heterogenous 
    sumdat.fSoloAuto.rsNorm = []; sumdat.fSoloAuto.rsRaw = [];
    sumdat.mSoloAuto.rsNorm = []; sumdat.mSoloAuto.rsRaw = [];
    sumdat.fSoloHetero.rsNorm = []; sumdat.fSoloHetero.rsRaw = [];
    sumdat.mSoloHetero.rsNorm = []; sumdat.mSoloHetero.rsRaw = [];

% Spike rates (spikes per second) for solo and duet syllables, Female and Male, Autogenous and Heterogenous 
    sumdat.mSoloHetero.SPS = []; sumdat.fSoloHetero.SPS = [];
    sumdat.mDuetHetero.SPS = []; sumdat.fDuetHetero.SPS = [];
    sumdat.mSoloAuto.SPS = []; sumdat.fSoloAuto.SPS = [];
    sumdat.mDuetAuto.SPS = []; sumdat.fDuetAuto.SPS = [];

% LEGACY: Extract appropriate wData references from birdlist

    PairList = birdlist{whichlist}(2:2:end)/2; % User chooses which data to analyze

%% Loop to calculate SPS and RS values for each pair of wrens   

    for curpair = PairList
        
    % Solo syllables MALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(msolosyls{curpair}) % Male sang solo syllables
        
        % Calculate RS values
        out(curpair).fSoloHetero = rs(in(curpair*2), msolosyls{curpair}, Cspon(:,curpair), pad, uc);
        out(curpair).mSoloAuto = rs(in((curpair*2)-1), msolosyls{curpair}, Cspon(:,curpair), -pad, uc);
        
        for kk = 1:length(msolosyls{curpair})
            sumdat.fSoloHetero.rsNorm(end+1) = out(curpair).fSoloHetero(kk).rsNorm;
            sumdat.fSoloHetero.rsRaw(end+1) = out(curpair).fSoloHetero(kk).rsRaw;
            sumdat.fSoloHetero.SPS(end+1) = out(curpair).fSoloHetero(kk).spikerate;
            sumdat.mSoloAuto.rsNorm(end+1) = out(curpair).mSoloAuto(kk).rsNorm;
            sumdat.mSoloAuto.rsRaw(end+1) = out(curpair).mSoloAuto(kk).rsRaw;
            sumdat.mSoloAuto.SPS(end+1) = out(curpair).mSoloAuto(kk).spikerate;
        end
    end
    
    % Solo syllables FEMALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(fsolosyls{curpair}) % Female sang solo syllables
        
        out(curpair).mSoloHetero = rs(in((curpair*2)-1), fsolosyls{curpair}, Cspon(:,curpair), pad, uc);
        out(curpair).fSoloAuto = rs(in(curpair*2), fsolosyls{curpair}, Cspon(:,curpair), -pad, uc);
        
        for kk = 1:length(fsolosyls{curpair})
            sumdat.mSoloHetero.rsNorm(end+1) = out(curpair).mSoloHetero(kk).rsNorm;
            sumdat.mSoloHetero.rsRaw(end+1) = out(curpair).mSoloHetero(kk).rsRaw;
            sumdat.mSoloHetero.SPS(end+1) = out(curpair).mSoloHetero(kk).spikerate;
            sumdat.fSoloAuto.rsNorm(end+1) = out(curpair).fSoloAuto(kk).rsNorm;
            sumdat.fSoloAuto.rsRaw(end+1) = out(curpair).fSoloAuto(kk).rsRaw;
            sumdat.fSoloAuto.SPS(end+1) = out(curpair).fSoloAuto(kk).spikerate;
        end
    end
    
    %% Duet syllables MALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(mduetsyls{curpair}) % Male sang duet syllables
        
        out(curpair).mDuetAuto = rs(in((curpair*2)-1), mduetsyls{curpair}, Cspon(:,curpair), -pad, uc);
        out(curpair).fDuetHetero = rs(in(curpair*2), mduetsyls{curpair}, Cspon(:,curpair), pad, uc);

        for kk = 1:length(mduetsyls{curpair})
            sumdat.mDuetAuto.rsNorm(end+1) = out(curpair).mDuetAuto(kk).rsNorm; 
            sumdat.mDuetAuto.rsRaw(end+1) = out(curpair).mDuetAuto(kk).rsRaw;
            sumdat.mDuetAuto.SPS(end+1) = out(curpair).mDuetAuto(kk).spikerate;
            sumdat.fDuetHetero.rsNorm(end+1) = out(curpair).fDuetHetero(kk).rsNorm;
            sumdat.fDuetHetero.rsRaw(end+1) = out(curpair).fDuetHetero(kk).rsRaw;
            sumdat.fDuetHetero.SPS(end+1) = out(curpair).fDuetHetero(kk).spikerate;
        end
    end 
    
    %% Duet syllables FEMALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(fduetsyls{curpair}) % Female sang solo syllables
        
        out(curpair).mDuetHetero = rs(in((curpair*2)-1), fduetsyls{curpair}, Cspon(:,curpair), pad, uc);
        out(curpair).fDuetAuto = rs(in(curpair*2), fduetsyls{curpair}, Cspon(:,curpair), -pad, uc);

        for kk = 1:length(fduetsyls{curpair})
            sumdat.mDuetHetero.rsNorm(end+1) = out(curpair).mDuetHetero(kk).rsNorm;
            sumdat.mDuetHetero.rsRaw(end+1) = out(curpair).mDuetHetero(kk).rsRaw;
            sumdat.mDuetHetero.SPS(end+1) = out(curpair).mDuetHetero(kk).spikerate;
            sumdat.fDuetAuto.rsNorm(end+1) = out(curpair).fDuetAuto(kk).rsNorm;
            sumdat.fDuetAuto.rsRaw(end+1) = out(curpair).fDuetAuto(kk).rsRaw;
            sumdat.fDuetAuto.SPS(end+1) = out(curpair).fDuetAuto(kk).spikerate;
        end
    end

    end % End of calculations

%% Plot MOTOR

% For Normalized RS data

    meanNorm(1) = mean(sumdat.mDuetAuto.rsNorm); s(1) = std(sumdat.mDuetAuto.rsNorm);
    meanNorm(2) = mean(sumdat.mSoloAuto.rsNorm); s(2) = std(sumdat.mSoloAuto.rsNorm);
    meanNorm(3) = mean(sumdat.fDuetAuto.rsNorm); s(3) = std(sumdat.fDuetAuto.rsNorm);
    meanNorm(4) = mean(sumdat.fSoloAuto.rsNorm); s(4) = std(sumdat.fSoloAuto.rsNorm);

% For raw RS data

    meanRaw(1) = mean(sumdat.mDuetAuto.rsRaw); sraw(1) = std(sumdat.mDuetAuto.rsRaw);
    meanRaw(2) = mean(sumdat.mSoloAuto.rsRaw); sraw(2) = std(sumdat.mSoloAuto.rsRaw);
    meanRaw(3) = mean(sumdat.fDuetAuto.rsRaw); sraw(3) = std(sumdat.fDuetAuto.rsRaw);
    meanRaw(4) = mean(sumdat.fSoloAuto.rsRaw); sraw(4) = std(sumdat.fSoloAuto.rsRaw);

% For SPS data
    
    meanSPS(1) = mean(sumdat.mDuetAuto.SPS); sps(1) = std(sumdat.mDuetAuto.SPS);
    meanSPS(2) = mean(sumdat.mSoloAuto.SPS); sps(2) = std(sumdat.mSoloAuto.SPS);
    meanSPS(3) = mean(sumdat.fDuetAuto.SPS); sps(3) = std(sumdat.fDuetAuto.SPS);
    meanSPS(4) = mean(sumdat.fSoloAuto.SPS); sps(4) = std(sumdat.fSoloAuto.SPS);
    
figure(1); % Spikes Per Second PLOTS    
hold on; title('Auto Spikes/Second'); 
    if whichlist == 4 % MALE SOLO/DUET DATA
    plot([2 1], meanSPS(1:2), 'b.', 'MarkerSize', 16); 
    errorbar([2 1], meanSPS(1:2), sps(1:2), 'b');
        for p=1:length(sumdat.mDuetAuto.SPS); plot(2.1, sumdat.mDuetAuto.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.mSoloAuto.SPS); plot(1.1, sumdat.mSoloAuto.SPS(p), 'k.', 'MarkerSize', 8); end
    end
    if whichlist == 2 % FEMALE SOLO/DUET DATA
    plot([4 3], meanSPS(3:4), 'm.', 'MarkerSize', 16); 
    errorbar([4 3], meanSPS(3:4), sps(3:4), 'm' );
        for p=1:length(sumdat.fDuetAuto.SPS); plot(4.1, sumdat.fDuetAuto.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.fSoloAuto.SPS); plot(3.1, sumdat.fSoloAuto.SPS(p), 'k.', 'MarkerSize', 8); end
    %ylim([-5 40]); 
    xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');
    xticklabels({' ','S',' ','D',' ','S',' ','D',' '})
    end

% figure(2); % RAW RS AUTOGENOUS PLOTS 
% hold on; title('Auto Raw RS'); 
%    if whichlist == 4 % MALE SOLO/DUET DATA
%     plot([2 1], meanRaw(1:2), 'b.', 'MarkerSize', 16); 
%     errorbar([2 1], meanRaw(1:2), sraw(1:2), 'b');
%         for p=1:length(sumdat.mDuetAuto.rsRaw); plot(2.1, sumdat.mDuetAuto.rsRaw(p), 'k.', 'MarkerSize', 8); end
%         for p=1:length(sumdat.mSoloAuto.rsRaw); plot(1.1, sumdat.mSoloAuto.rsRaw(p), 'k.', 'MarkerSize', 8); end
%    end
%    if whichlist == 2 % FEMALE SOLO/DUET DATA
%     plot([4 3], meanRaw(3:4), 'm.', 'MarkerSize', 16); 
%     errorbar([4 3], meanRaw(3:4), sraw(3:4), 'm' );
%         for p=1:length(sumdat.fDuetAuto.rsRaw); plot(4.1, sumdat.fDuetAuto.rsRaw(p), 'k.', 'MarkerSize', 8); end
%         for p=1:length(sumdat.fSoloAuto.rsRaw); plot(3.1, sumdat.fSoloAuto.rsRaw(p), 'k.', 'MarkerSize', 8); end
%     ylim([-10 65]); xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');
%     xticklabels({' ','S',' ','D',' ','S',' ','D',' '})
%    end

% figure(3); % NORM AUTOGENOUS PLOTS    
% hold on; title('Auto Norm RS'); 
%    if whichlist == 4 % MALE SOLO/DUET DATA
%     plot([2 1], meanNorm(1:2), 'b.', 'MarkerSize', 16); 
%     errorbar([2 1], meanNorm(1:2), s(1:2), 'b');
%         for p=1:length(sumdat.mDuetAuto.rsNorm); plot(2.1, sumdat.mDuetAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
%         for p=1:length(sumdat.mSoloAuto.rsNorm); plot(1.1, sumdat.mSoloAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
%    end
%    if whichlist == 2 % FEMALE SOLO/DUET DATA
%     plot([4 3], meanNorm(3:4), 'm.', 'MarkerSize', 16); 
%     errorbar([4 3], meanNorm(3:4), s(3:4), 'm' );
%         for p=1:length(sumdat.fDuetAuto.rsNorm); plot(4.1, sumdat.fDuetAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
%         for p=1:length(sumdat.fSoloAuto.rsNorm); plot(3.1, sumdat.fSoloAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
%     ylim([-5 40]); xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');
%     xticklabels({' ','S',' ','D',' ','S',' ','D',' '})
%    end
    
%% Plot Sensory

% Recall that PairList{1} = birdlist{2}(2:2:end)/2; % Index 1 is the FEMALE SOLO/DUET data
% Recall that PairList{2} = birdlist{4}(2:2:end)/2; % Index 2 is the MALE SOLO/DUET data

% Calculate

% For Normalized RS data

    meanNorm(1) = mean(sumdat.mDuetHetero.rsNorm); s(1) = std(sumdat.mDuetHetero.rsNorm);
    meanNorm(2) = mean(sumdat.mSoloHetero.rsNorm); s(2) = std(sumdat.mSoloHetero.rsNorm);
    meanNorm(3) = mean(sumdat.fDuetHetero.rsNorm); s(3) = std(sumdat.fDuetHetero.rsNorm);
    meanNorm(4) = mean(sumdat.fSoloHetero.rsNorm); s(4) = std(sumdat.fSoloHetero.rsNorm);

% For raw RS data

    meanRaw(1) = mean(sumdat.mDuetHetero.rsRaw); sraw(1) = std(sumdat.mDuetHetero.rsRaw);
    meanRaw(2) = mean(sumdat.mSoloHetero.rsRaw); sraw(2) = std(sumdat.mSoloHetero.rsRaw);
    meanRaw(3) = mean(sumdat.fDuetHetero.rsRaw); sraw(3) = std(sumdat.fDuetHetero.rsRaw);
    meanRaw(4) = mean(sumdat.fSoloHetero.rsRaw); sraw(4) = std(sumdat.fSoloHetero.rsRaw);

% For SPS data
    
    meanSPS(1) = mean(sumdat.mDuetHetero.SPS); sps(1) = std(sumdat.mDuetHetero.SPS);
    meanSPS(2) = mean(sumdat.mSoloHetero.SPS); sps(2) = std(sumdat.mSoloHetero.SPS);
    meanSPS(3) = mean(sumdat.fDuetHetero.SPS); sps(3) = std(sumdat.fDuetHetero.SPS);
    meanSPS(4) = mean(sumdat.fSoloHetero.SPS); sps(4) = std(sumdat.fSoloHetero.SPS);
    
figure(4); % RAW HETEROGENOUS PLOTS
hold on; title('Hetero SPS');
    if whichlist == 2 % FEMALE SOLO/DUET DATA (hetero requires we plot OTHER bird)
    plot([2 1], meanSPS(1:2), 'b.', 'MarkerSize', 16); 
    errorbar([2 1], meanSPS(1:2), sps(1:2), 'b' );
        for p=1:length(sumdat.mDuetHetero.SPS); plot(2.1, sumdat.mDuetHetero.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.mSoloHetero.SPS); plot(1.1, sumdat.mSoloHetero.SPS(p), 'k.', 'MarkerSize', 8); end
    end
    if whichlist == 4 % MALE SOLO/DUET DATA (hetero requires we plot OTHER bird)
    plot([4 3], meanSPS(3:4), 'm.', 'MarkerSize', 16); 
    errorbar([4 3], meanSPS(3:4), sps(3:4), 'm' );
        for p=1:length(sumdat.fDuetHetero.SPS); plot(4.1, sumdat.fDuetHetero.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.fSoloHetero.SPS); plot(3.1, sumdat.fSoloHetero.SPS(p), 'k.', 'MarkerSize', 8); end
    % ylim([-10 65]); 
    xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');
    xticklabels({' ','S',' ','D',' ','S',' ','D',' '})
    end
    
% figure(5); clf; % RAW HETEROGENOUS PLOTS
% hold on; title('Hetero Raw RS');
%     plot([2 1], meanRaw(1:2), 'b.', 'MarkerSize', 16); 
%     errorbar([2 1], meanRaw(1:2), sraw(1:2), 'b' );
%         for p=1:length(sumdat.mDuetHetero.rsRaw); plot(2.1, sumdat.mDuetHetero.rsRaw(p), 'k.', 'MarkerSize', 8); end
%         for p=1:length(sumdat.mSoloHetero.rsRaw); plot(1.1, sumdat.mSoloHetero.rsRaw(p), 'k.', 'MarkerSize', 8); end
%     plot([4 3], meanRaw(3:4), 'm.', 'MarkerSize', 16); 
%     errorbar([4 3], meanRaw(3:4), sraw(3:4), 'm' );
%         for p=1:length(sumdat.fDuetHetero.rsRaw); plot(4.1, sumdat.fDuetHetero.rsRaw(p), 'k.', 'MarkerSize', 8); end
%         for p=1:length(sumdat.fSoloHetero.rsRaw); plot(3.1, sumdat.fSoloHetero.rsRaw(p), 'k.', 'MarkerSize', 8); end
%     ylim([-10 65]); xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');
%     xticklabels({' ','S',' ','D',' ','S',' ','D',' '})
% 
% figure(6); clf; % NORM HETEROGENOUS PLOTS
% hold on; title('Hetero Norm RS');
%     plot([2 1], meanNorm(1:2), 'b.', 'MarkerSize', 16); 
%     errorbar([2 1], meanNorm(1:2), s(1:2), 'b' );
%         for p=1:length(sumdat.mDuetHetero.rsNorm); plot(2.1, sumdat.mDuetHetero.rsNorm(p), 'k.', 'MarkerSize', 8); end
%         for p=1:length(sumdat.mSoloHetero.rsNorm); plot(1.1, sumdat.mSoloHetero.rsNorm(p), 'k.', 'MarkerSize', 8); end
%     plot([4 3], meanNorm(3:4), 'm.', 'MarkerSize', 16); 
%     errorbar([4 3], meanNorm(3:4), s(3:4), 'm' );
%         for p=1:length(sumdat.fDuetHetero.rsNorm); plot(4.1, sumdat.fDuetHetero.rsNorm(p), 'k.', 'MarkerSize', 8); end
%         for p=1:length(sumdat.fSoloHetero.rsNorm); plot(3.1, sumdat.fSoloHetero.rsNorm(p), 'k.', 'MarkerSize', 8); end
%     ylim([-5 40]); xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');
%     xticklabels({' ','S',' ','D',' ','S',' ','D',' '})
    
    
%% Compute stats     
    
% Autogenous duet RS significant from 0?
    stts.m.dNAuto.mean = mean(sumdat.mDuetAuto.rsNorm);
    [stts.m.dNAuto.H, stts.m.dNAuto.P, stts.m.dNAuto.CI, stts.m.dNAuto.stats]  = ttest(sumdat.mDuetAuto.rsNorm);
    stts.m.dRAuto.mean = mean(sumdat.mDuetAuto.rsRaw);
    [stts.m.dRAuto.H, stts.m.dRAuto.P, stts.m.dRAuto.CI, stts.m.dRAuto.stats]  = ttest(sumdat.mDuetAuto.rsRaw);
%    stts.m.dSAuto.mean = mean(sumdat.mDuetAuto.SPS);
%    [stts.m.dSAuto.H, stts.m.dSAuto.P, stts.m.dSAuto.CI, stts.m.dSAuto.stats]  = ttest(sumdat.mDuetAuto.SPS);
    
    stts.f.dNAuto.mean = mean(sumdat.fDuetAuto.rsNorm);
    [stts.f.dNAuto.H, stts.f.dNAuto.P, stts.f.dNAuto.CI, stts.f.dNAuto.stats]  = ttest(sumdat.fDuetAuto.rsNorm);
    stts.f.dRAuto.mean = mean(sumdat.fDuetAuto.rsRaw);
    [stts.f.dRAuto.H, stts.f.dRAuto.P, stts.f.dRAuto.CI, stts.f.dRAuto.stats]  = ttest(sumdat.fDuetAuto.rsRaw);
%    stts.f.dSAuto.mean = mean(sumdat.fDuetAuto.SPS);
%    [stts.f.dSAuto.H, stts.f.dSAuto.P, stts.f.dSAuto.CI, stts.f.dSAuto.stats]  = ttest(sumdat.mDuetAuto.SPS);
    
% Autogenous Solo RS significant from 0?

    stts.m.sNAuto.mean = mean(sumdat.mSoloAuto.rsNorm);
    [stts.m.sNAuto.H, stts.m.sNAuto.P, stts.m.sNAuto.CI, stts.m.sNAuto.stats]  = ttest(sumdat.mSoloAuto.rsNorm);
    stts.m.sRAuto.mean = mean(sumdat.mSoloAuto.rsRaw);
    [stts.m.sRAuto.H, stts.m.sRAuto.P, stts.m.sRAuto.CI, stts.m.sRAuto.stats]  = ttest(sumdat.mSoloAuto.rsRaw);
    
    stts.f.sNAuto.mean = mean(sumdat.fSoloAuto.rsNorm);
    [stts.f.sNAuto.H, stts.f.sNAuto.P, stts.f.sNAuto.CI, stts.f.sNAuto.stats]  = ttest(sumdat.fSoloAuto.rsNorm);
    stts.m.sRAuto.mean = mean(sumdat.fSoloAuto.rsRaw);
    [stts.f.sRAuto.H, stts.f.sRAuto.P, stts.f.sRAuto.CI, stts.f.sRAuto.stats]  = ttest(sumdat.fSoloAuto.rsRaw);

% Difference between Autogenous Duet and Solo ?

[stts.m.SvsDNAuto.P, stts.m.SvsDNAuto.H, stts.m.SvsDNAuto.stats]  = ranksum(sumdat.mSoloAuto.rsNorm, sumdat.mDuetAuto.rsNorm);
[stts.m.SvsDRAuto.P, stts.m.SvsDRAuto.H, stts.m.SvsDRAuto.stats]  = ranksum(sumdat.mSoloAuto.rsRaw, sumdat.mDuetAuto.rsRaw);
[stts.m.SvsDSAuto.P, stts.m.SvsDSAuto.H, stts.m.SvsDSAuto.stats]  = ranksum(sumdat.mSoloAuto.SPS, sumdat.mDuetAuto.SPS);

[stts.f.SvsDNAuto.P, stts.f.SvsDNAuto.H, stts.f.SvsDNAuto.stats]  = ranksum(sumdat.fSoloAuto.rsNorm, sumdat.fDuetAuto.rsNorm);
[stts.f.SvsDRAuto.P, stts.f.SvsDRAuto.H, stts.f.SvsDRAuto.stats]  = ranksum(sumdat.fSoloAuto.rsRaw, sumdat.fDuetAuto.rsRaw);
[stts.f.SvsDSAuto.P, stts.f.SvsDSAuto.H, stts.f.SvsDSAuto.stats]  = ranksum(sumdat.fSoloAuto.SPS, sumdat.fDuetAuto.SPS);

% Heterogenous duet RS significant from 0?

    stts.m.dNHetero.mean = mean(sumdat.mDuetHetero.rsNorm);
    [stts.m.dNHetero.H, stts.m.dNHetero.P, stts.m.dNHetero.CI, stts.m.dNHetero.stats]  = ttest(sumdat.mDuetHetero.rsNorm);
    stts.m.dRHetero.mean = mean(sumdat.mDuetHetero.rsRaw);
    [stts.m.dRHetero.H, stts.m.dRHetero.P, stts.m.dRHetero.CI, stts.m.dRHetero.stats]  = ttest(sumdat.mDuetHetero.rsRaw);
    
    stts.f.dNHetero.mean = mean(sumdat.fDuetHetero.rsNorm);
    [stts.f.dNHetero.H, stts.f.dNHetero.P, stts.f.dNHetero.CI, stts.f.dNHetero.stats]  = ttest(sumdat.fDuetHetero.rsNorm);
    stts.f.dRHetero.mean = mean(sumdat.fDuetHetero.rsRaw);
    [stts.f.dRHetero.H, stts.f.dRHetero.P, stts.f.dRHetero.CI, stts.f.dRHetero.stats]  = ttest(sumdat.fDuetHetero.rsRaw);
    
% Heterogenous Solo RS significant from 0?

    stts.m.sNHetero.mean = mean(sumdat.mSoloHetero.rsNorm);
    [stts.m.sNHetero.H, stts.m.sNHetero.P, stts.m.sNHetero.CI, stts.m.sNHetero.stats]  = ttest(sumdat.mSoloHetero.rsNorm);
    stts.m.sRHetero.mean = mean(sumdat.mSoloHetero.rsRaw);
    [stts.m.sRHetero.H, stts.m.sRHetero.P, stts.m.sRHetero.CI, stts.m.sRHetero.stats]  = ttest(sumdat.mSoloHetero.rsRaw);
    
    stts.f.sNHetero.mean = mean(sumdat.fSoloHetero.rsNorm);
    [stts.f.sNHetero.H, stts.f.sNHetero.P, stts.f.sNHetero.CI, stts.f.sNHetero.stats]  = ttest(sumdat.fSoloHetero.rsNorm);
    stts.f.sRHetero.mean = mean(sumdat.fSoloHetero.rsRaw);
    [stts.f.sRHetero.H, stts.f.sRHetero.P, stts.f.sRHetero.CI, stts.f.sRHetero.stats]  = ttest(sumdat.fSoloHetero.rsRaw);

% Difference between Heterogenous Duet and Solo RS motor?

[stts.m.SvsDNHetero.P, stts.m.SvsDNHetero.H, stts.m.SvsDNHetero.stats]  = ranksum(sumdat.mSoloHetero.rsNorm, sumdat.mDuetHetero.rsNorm);
[stts.m.SvsDRHetero.P, stts.m.SvsDRHetero.H, stts.m.SvsDRHetero.stats]  = ranksum(sumdat.mSoloHetero.rsRaw, sumdat.mDuetHetero.rsRaw);
[stts.m.SvsDSHetero.P, stts.m.SvsDSHetero.H, stts.m.SvsDSHetero.stats]  = ranksum(sumdat.mSoloHetero.SPS, sumdat.mDuetHetero.SPS);

[stts.f.SvsDNHetero.P, stts.f.SvsDNHetero.H, stts.f.SvsDNHetero.stats]  = ranksum(sumdat.fSoloHetero.rsNorm, sumdat.fDuetHetero.rsNorm);
[stts.f.SvsDRHetero.P, stts.f.SvsDRHetero.H, stts.f.SvsDRHetero.stats]  = ranksum(sumdat.fSoloHetero.rsRaw, sumdat.fDuetHetero.rsRaw);
[stts.f.SvsDSHetero.P, stts.f.SvsDSHetero.H, stts.f.SvsDSHetero.stats]  = ranksum(sumdat.fSoloHetero.SPS, sumdat.fDuetHetero.SPS);

[MvFHS.H, MvFHS.P, MvFHS.CI, MvFHS.stats]  = ttest2(sumdat.fSoloHetero.rsRaw, sumdat.mSoloHetero.rsRaw, 'vartype', 'unequal');

% Test for equal variance

if whichlist == 4
    MaleDuetSoloC = [sumdat.mDuetAuto.rsRaw, sumdat.mSoloAuto.rsRaw];
    MaleDuetSoloCIDX = [ones(1,length(sumdat.mDuetAuto.rsRaw)), 2*ones(1,length(sumdat.mSoloAuto.rsRaw))];
    [MaleVarAutoChronP, MaleVarAutoChronstats] = vartestn(MaleDuetSoloC', MaleDuetSoloCIDX','TestType','LeveneAbsolute');
    fprintf('Male Chronic Auto: Variance difference between Duet and Solo SPS? p = %1.5f \n', MaleVarAutoChronP);
end

if whichlist == 2
    FeMaleDuetSoloC = [sumdat.fDuetAuto.rsRaw, sumdat.fSoloAuto.rsRaw];
    FeMaleDuetSoloCIDX = [ones(1,length(sumdat.fDuetAuto.rsRaw)), 2*ones(1,length(sumdat.fSoloAuto.rsRaw))];
    [FemaleVarAutoChronP,FemaleVarAutoChronstats] = vartestn(FeMaleDuetSoloC', FeMaleDuetSoloCIDX','TestType','LeveneAbsolute');
    fprintf('Female Chronic Auto: Variance difference between Duet and Solo SPS? p = %1.5f \n', FemaleVarAutoChronP);
end

% FvsMAutoSoloC = [sumdat.fSoloAuto.rsRaw, sumdat.mSoloAuto.rsRaw];
% FvsMAutoSoloCIDX = [ones(1,length(sumdat.mSoloAuto.rsRaw)), 2*ones(1,length(sumdat.fSoloAuto.rsRaw))];
% [FvsMVarAutoChronP,FvsMVarAutoChronstats] = vartestn(FvsMAutoSoloC', FvsMAutoSoloCIDX','TestType','LeveneAbsolute');


% REPORT THE STATS FOR THE PAPER
%     fprintf('Male Auto Duet Raw RS vs Solo? p = %1.5f \n', stts.m.SvsDRAuto.P);
%     fprintf('Female Auto Duet Raw RS vs Solo? p = %1.5f \n', stts.f.SvsDRAuto.P);

if whichlist == 4
    fprintf('Male Auto SPS Duet vs Solo? p = %1.5f \n', stts.m.SvsDSAuto.P);
end
if whichlist == 2
    fprintf('Female Auto SPS Duet vs Solo? p = %1.5f \n', stts.f.SvsDSAuto.P);
end
    fprintf(' \n');
    
%     fprintf('Male Auto Duet Raw RS different from zero? p = %1.5f \n', stts.m.dRAuto.P);
%     fprintf('Male Auto Solo Raw RS different from zero? p = %1.5f \n', stts.m.sRAuto.P);
%     fprintf('Female Auto Duet Raw RS different from zero? p = %1.5f \n', stts.f.dRAuto.P);
%     fprintf('Female Auto Solo Raw RS different from zero? p = %1.5f \n', stts.f.sRAuto.P);
%     
%     fprintf(' \n');

%     fprintf('Male Hetero Duet Raw RS different from zero? p = %1.5f \n', stts.m.dRHetero.P);
%     fprintf('Male Hetero Solo Raw RS different from zero? p = %1.5f \n', stts.m.sRHetero.P);
%     fprintf('Female Hetero Duet Raw RS different from zero? p = %1.5f \n', stts.f.dRHetero.P);
%     fprintf('Female Hetero Solo Raw RS different from zero? p = %1.5f \n', stts.f.sRHetero.P);
%     
%     fprintf(' \n');

%     fprintf('Male Hetero Duet Raw RS vs Solo? p = %1.5f \n', stts.m.SvsDRHetero.P);
%     fprintf('Female Hetero Duet Raw RS vs Solo? p = %1.5f \n', stts.f.SvsDRHetero.P);
if whichlist == 4
    fprintf('Female Hetero Duet SPS vs Solo? p = %1.5f \n', stts.f.SvsDSHetero.P);
end
if whichlist == 2
    fprintf('Male Hetero Duet SPS vs Solo? p = %1.5f \n', stts.m.SvsDSHetero.P);
end
%
%     fprintf('Male versus Female Solo Hetero Raw RS different? p = %1.5f \n', MvFHS.P);
%     
     fprintf(' \n');

%     fprintf('FeMale Chronic Auto: Variance difference between Duet and Solo Raw RS? p = %1.5f \n', FemaleVarAutoChronP);
    
fprintf('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n');


%     fprintf('Male Auto Duet Norm RS vs Solo? p = %1.5f \n', stts.m.SvsDNAuto.P);
%     fprintf('Female Auto Duet Norm RS vs Solo? p = %1.5f \n', stts.f.SvsDNAuto.P);
%     
%     fprintf(' \n');

%     fprintf('Male Hetero Duet Norm RS different from zero? p = %1.5f \n', stts.m.dNHetero.P);
%     fprintf('Male Hetero Solo Norm RS different from zero? p = %1.5f \n', stts.m.sNHetero.P);
%     fprintf('Female Hetero Duet Norm RS different from zero? p = %1.5f \n', stts.f.dNHetero.P);
%     fprintf('Female Hetero Solo Norm RS different from zero? p = %1.5f \n', stts.f.sNHetero.P);
%     
%     fprintf(' \n');

%     fprintf('Male Hetero Duet Norm RS vs Solo? p = %1.5f \n', stts.m.SvsDNHetero.P);
%     fprintf('Female Hetero Duet Norm RS vs Solo? p = %1.5f \n', stts.f.SvsDNHetero.P);
%         
%     fprintf(' \n');

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%% Response Strength nested function

function qwe = rs(struc, syllabl, spontan, padme, ccuu)
 
   % Get spontaneous rate 
        sponSpikeCount = 0;  

if ccuu == 2 % CHRONIC
    spikeys = struc.Cspikes;
end
if ccuu == 1 % URETHANE
    spikeys = struc.Aspikes;
end
        for i=1:length(spikeys) % For each rep in acute, or each electrode in chronic
            sponSpikeCount = sponSpikeCount + length(find(spikeys{i} > spontan(1) & spikeys{i} < spontan(2)));
        end

        sponSPS = sponSpikeCount / (spontan(2) - spontan(1)); % This is spikes per second
        sponSPS = sponSPS/length(spikeys); % Divide by number of reps or number of electrodes
        
        
        % padme is positive for heterogenous and negative for autogenous
        
        % In CHRONIC, for PREMOTOR ACTIVITY, we want to shift the window
        % earlier (subtract time) and shift later (add time) for SENSORY.
        % When we shift the window later for heterogenous (SENSORY) in CHRONIC, the window will start
        % to overlap with the premotor from the next syllable, so we use
        % the variable skinny to solve that issue.
        
        % In ACUTE, we always want to shift the window later (add time).
        
 
            if ccuu == 1 % URETHANE
                padme = abs(padme); % ALWAYS add time to shift window later
                skinny = 0; % No special truncations
            end
            if ccuu == 2 % CHRONIC
                if padme > 0 % We have a heterogenous syllable and need to truncate the analysis. 
                    skinny = 2.5 * padme; % IMPORTANT, try 2 for a more conservative truncation.
                end
            end
        
        
    % Loop for each syllable
    
    for j = 1:length(syllabl)    
    
        % Get the spikes for that syllable
        stimSpikeCount = 0; 
    
        for i=1:length(spikeys) % For each electrode or repetition (padme shifts the window in seconds, negative earlier, positive later)
            stimSpikeCount = stimSpikeCount + length(find(spikeys{i} >= struc.syl(syllabl(j)).tim(1) + padme & spikeys{i} < struc.syl(syllabl(j)).tim(2) + padme - skinny));
        end
        
        stimSPS = stimSpikeCount / ((struc.syl(syllabl(j)).tim(2) - struc.syl(syllabl(j)).tim(1)) - skinny); % This is spikes per second
        stimSPS = stimSPS/length(spikeys); % Divide by number of reps or number of electrodes
        
        
        qwe(j).sylnum = syllabl(j);
        qwe(j).rsNorm = (stimSPS - sponSPS) / (sponSPS + 0.0000000000001); % Avoid divide by zero
        qwe(j).rsRaw = stimSPS - sponSPS;
        qwe(j).sponrate = sponSPS;
        qwe(j).spikerate = stimSPS;
        qwe(j).spontim = spontan;
        qwe(j).syltim = [struc.syl(syllabl(j)).tim(1), struc.syl(syllabl(j)).tim(2)];
        qwe(j).pad = padme;
        
    end
end

end







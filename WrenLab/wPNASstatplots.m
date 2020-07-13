function [out, sumdat, stts] = wPNASstatplots(in, padding)
% Usage: [out, sumdat, stts] = wRS_Chronic(w, padding)
% Calculates response strength to solo and duet syllables.
% Load the Chronic data structure first:
% load ChronicCompleat2019f.mat (Current as of 8-Jan-2020)
% w is the structure from that file with all of the data
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

analyzethis = [2, 4]; % The birdlist groups to analyze

for pp = analyzethis
    
pad = 0.050; 

% The user can specify the padding via an argin for convenience.
if nargin == 2; pad = padding; end


% Initializing variables

% Normalized and Raw response strengths for DUET syllables: Female and Male, Autogenous and Heterogenous 
    sumdat(max(analyzethis)).fDuetAuto.rsNorm = []; sumdat(max(analyzethis)).fDuetAuto.rsRaw = [];
    sumdat(max(analyzethis)).mDuetAuto.rsNorm = []; sumdat(max(analyzethis)).mDuetAuto.rsRaw = [];
    sumdat(max(analyzethis)).fDuetHetero.rsNorm = []; sumdat(max(analyzethis)).fDuetHetero.rsRaw = [];
    sumdat(max(analyzethis)).mDuetHetero.rsNorm = []; sumdat(max(analyzethis)).mDuetHetero.rsRaw = [];

% Normalized and Raw response strengths for SOLO syllables: Female and Male, Autogenous and Heterogenous 
    sumdat(max(analyzethis)).fSoloAuto.rsNorm = []; sumdat(max(analyzethis)).fSoloAuto.rsRaw = [];
    sumdat(max(analyzethis)).mSoloAuto.rsNorm = []; sumdat(max(analyzethis)).mSoloAuto.rsRaw = [];
    sumdat(max(analyzethis)).fSoloHetero.rsNorm = []; sumdat(max(analyzethis)).fSoloHetero.rsRaw = [];
    sumdat(max(analyzethis)).mSoloHetero.rsNorm = []; sumdat(max(analyzethis)).mSoloHetero.rsRaw = [];

% Spike rates (spikes per second) for solo and duet syllables, Female and Male, Autogenous and Heterogenous 
    sumdat(max(analyzethis)).mSoloHetero.SPS = []; sumdat(max(analyzethis)).fSoloHetero.SPS = [];
    sumdat(max(analyzethis)).mDuetHetero.SPS = []; sumdat(max(analyzethis)).fDuetHetero.SPS = [];
    sumdat(max(analyzethis)).mSoloAuto.SPS = []; sumdat(max(analyzethis)).fSoloAuto.SPS = [];
    sumdat(max(analyzethis)).mDuetAuto.SPS = []; sumdat(max(analyzethis)).fDuetAuto.SPS = [];

%% List of Chronic singing data with syllable indices and locations for spontaneous activity

[msolosyls, mduetsyls, fsolosyls, fduetsyls, Cspon, ~, birdlist] = wData;

%% Loop to calculate data values for each pair of wrens   

for curpair = birdlist{pp} 
    
    % Solo syllables MALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(msolosyls{curpair}) % Male sang solo syllables
        
        % Calculate RS values
        out(curpair).fSoloHetero = rs(in(curpair*2), msolosyls{curpair}, Cspon(:,curpair), pad);
        out(curpair).mSoloAuto = rs(in((curpair*2)-1), msolosyls{curpair}, Cspon(:,curpair), -pad);
        
        for kk = 1:length(msolosyls{curpair})
            sumdat(pp).fSoloHetero.rsNorm(end+1) = out(curpair).fSoloHetero(kk).rsNorm;
            sumdat(pp).fSoloHetero.rsRaw(end+1) = out(curpair).fSoloHetero(kk).rsRaw;
            sumdat(pp).fSoloHetero.SPS(end+1) = out(curpair).fSoloHetero(kk).spikerate;
            sumdat(pp).mSoloAuto.rsNorm(end+1) = out(curpair).mSoloAuto(kk).rsNorm;
            sumdat(pp).mSoloAuto.rsRaw(end+1) = out(curpair).mSoloAuto(kk).rsRaw;
            sumdat(pp).mSoloAuto.SPS(end+1) = out(curpair).mSoloAuto(kk).spikerate;
        end
    end
    
    % Solo syllables FEMALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(fsolosyls{curpair}) % Female sang solo syllables
        
        out(curpair).mSoloHetero = rs(in((curpair*2)-1), fsolosyls{curpair}, Cspon(:,curpair), pad);
        out(curpair).fSoloAuto = rs(in(curpair*2), fsolosyls{curpair}, Cspon(:,curpair), -pad);
        
        for kk = 1:length(fsolosyls{curpair})
            sumdat(pp).mSoloHetero.rsNorm(end+1) = out(curpair).mSoloHetero(kk).rsNorm;
            sumdat(pp).mSoloHetero.rsRaw(end+1) = out(curpair).mSoloHetero(kk).rsRaw;
            sumdat(pp).mSoloHetero.SPS(end+1) = out(curpair).mSoloHetero(kk).spikerate;
            sumdat(pp).fSoloAuto.rsNorm(end+1) = out(curpair).fSoloAuto(kk).rsNorm;
            sumdat(pp).fSoloAuto.rsRaw(end+1) = out(curpair).fSoloAuto(kk).rsRaw;
            sumdat(pp).fSoloAuto.SPS(end+1) = out(curpair).fSoloAuto(kk).spikerate;
        end
    end
    
    %% Duet syllables MALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(mduetsyls{curpair}) % Male sang duet syllables
        
        out(curpair).mDuetAuto = rs(in((curpair*2)-1), mduetsyls{curpair}, Cspon(:,curpair), -pad);
        out(curpair).fDuetHetero = rs(in(curpair*2), mduetsyls{curpair}, Cspon(:,curpair), pad);

        for kk = 1:length(mduetsyls{curpair})
            sumdat(pp).mDuetAuto.rsNorm(end+1) = out(curpair).mDuetAuto(kk).rsNorm; 
            sumdat(pp).mDuetAuto.rsRaw(end+1) = out(curpair).mDuetAuto(kk).rsRaw;
            sumdat(pp).mDuetAuto.SPS(end+1) = out(curpair).mDuetAuto(kk).spikerate;
            sumdat(pp).fDuetHetero.rsNorm(end+1) = out(curpair).fDuetHetero(kk).rsNorm;
            sumdat(pp).fDuetHetero.rsRaw(end+1) = out(curpair).fDuetHetero(kk).rsRaw;
            sumdat(pp).fDuetHetero.SPS(end+1) = out(curpair).fDuetHetero(kk).spikerate;
        end
    end 
    
    %% Duet syllables FEMALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(fduetsyls{curpair}) % Female sang solo syllables
        
        out(curpair).mDuetHetero = rs(in((curpair*2)-1), fduetsyls{curpair}, Cspon(:,curpair), pad);
        out(curpair).fDuetAuto = rs(in(curpair*2), fduetsyls{curpair}, Cspon(:,curpair), -pad);

        for kk = 1:length(fduetsyls{curpair})
            sumdat(pp).mDuetHetero.rsNorm(end+1) = out(curpair).mDuetHetero(kk).rsNorm;
            sumdat(pp).mDuetHetero.rsRaw(end+1) = out(curpair).mDuetHetero(kk).rsRaw;
            sumdat(pp).mDuetHetero.SPS(end+1) = out(curpair).mDuetHetero(kk).spikerate;
            sumdat(pp).fDuetAuto.rsNorm(end+1) = out(curpair).fDuetAuto(kk).rsNorm;
            sumdat(pp).fDuetAuto.rsRaw(end+1) = out(curpair).fDuetAuto(kk).rsRaw;
            sumdat(pp).fDuetAuto.SPS(end+1) = out(curpair).fDuetAuto(kk).spikerate;
        end
    end

end % End of calculations


%% Get means and std for MOTOR plots

% For Normalized RS data

    meanNorm(1) = mean(sumdat(pp).mDuetAuto.rsNorm); s(1) = std(sumdat(pp).mDuetAuto.rsNorm);
    meanNorm(2) = mean(sumdat(pp).mSoloAuto.rsNorm); s(2) = std(sumdat(pp).mSoloAuto.rsNorm);
    meanNorm(3) = mean(sumdat(pp).fDuetAuto.rsNorm); s(3) = std(sumdat(pp).fDuetAuto.rsNorm);
    meanNorm(4) = mean(sumdat(pp).fSoloAuto.rsNorm); s(4) = std(sumdat(pp).fSoloAuto.rsNorm);

% For raw RS data

    meanRaw(1) = mean(sumdat(pp).mDuetAuto.rsRaw); sraw(1) = std(sumdat(pp).mDuetAuto.rsRaw);
    meanRaw(2) = mean(sumdat(pp).mSoloAuto.rsRaw); sraw(2) = std(sumdat(pp).mSoloAuto.rsRaw);
    meanRaw(3) = mean(sumdat(pp).fDuetAuto.rsRaw); sraw(3) = std(sumdat(pp).fDuetAuto.rsRaw);
    meanRaw(4) = mean(sumdat(pp).fSoloAuto.rsRaw); sraw(4) = std(sumdat(pp).fSoloAuto.rsRaw);

% For SPS data
    
    meanSPS(1) = mean(sumdat(pp).mDuetAuto.SPS); sps(1) = std(sumdat(pp).mDuetAuto.SPS);
    meanSPS(2) = mean(sumdat(pp).mSoloAuto.SPS); sps(2) = std(sumdat(pp).mSoloAuto.SPS);
    meanSPS(3) = mean(sumdat(pp).fDuetAuto.SPS); sps(3) = std(sumdat(pp).fDuetAuto.SPS);
    meanSPS(4) = mean(sumdat(pp).fSoloAuto.SPS); sps(4) = std(sumdat(pp).fSoloAuto.SPS);
    
figure(pp); clf; % Spikes Per Second PLOTS    
hold on; title('Auto Spikes/Second'); 
    plot([2 1], meanSPS(1:2), 'b.', 'MarkerSize', 16); 
    errorbar([2 1], meanSPS(1:2), sps(1:2), 'b');
        for p=1:length(sumdat(pp).mDuetAuto.SPS); plot(2.1, sumdat(pp).mDuetAuto.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat(pp).mSoloAuto.SPS); plot(1.1, sumdat(pp).mSoloAuto.SPS(p), 'k.', 'MarkerSize', 8); end
    plot([4 3], meanSPS(3:4), 'm.', 'MarkerSize', 16); 
    errorbar([4 3], meanSPS(3:4), sps(3:4), 'm' );
        for p=1:length(sumdat(pp).fDuetAuto.SPS); plot(4.1, sumdat(pp).fDuetAuto.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat(pp).fSoloAuto.SPS); plot(3.1, sumdat(pp).fSoloAuto.SPS(p), 'k.', 'MarkerSize', 8); end
    %ylim([-5 40]); 
    xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');
    xticklabels({' ','S',' ','D',' ','S',' ','D',' '})

% figure(pp); clf; % RAW AUTOGENOUS PLOTS 
% hold on; title('Auto Raw RS'); 
%     plot([2 1], meanRaw(1:2), 'b.', 'MarkerSize', 16); 
%     errorbar([2 1], meanRaw(1:2), sraw(1:2), 'b');
%         for p=1:length(sumdat.mDuetAuto.rsRaw); plot(2.1, sumdat.mDuetAuto.rsRaw(p), 'k.', 'MarkerSize', 8); end
%         for p=1:length(sumdat.mSoloAuto.rsRaw); plot(1.1, sumdat.mSoloAuto.rsRaw(p), 'k.', 'MarkerSize', 8); end
%     plot([4 3], meanRaw(3:4), 'm.', 'MarkerSize', 16); 
%     errorbar([4 3], meanRaw(3:4), sraw(3:4), 'm' );
%         for p=1:length(sumdat.fDuetAuto.rsRaw); plot(4.1, sumdat.fDuetAuto.rsRaw(p), 'k.', 'MarkerSize', 8); end
%         for p=1:length(sumdat.fSoloAuto.rsRaw); plot(3.1, sumdat.fSoloAuto.rsRaw(p), 'k.', 'MarkerSize', 8); end
%     ylim([-10 65]); xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');
%     xticklabels({' ','S',' ','D',' ','S',' ','D',' '})
%     
% figure(pp); clf; % NORM AUTOGENOUS PLOTS    
% hold on; title('Auto Norm RS'); 
%     plot([2 1], meanNorm(1:2), 'b.', 'MarkerSize', 16); 
%     errorbar([2 1], meanNorm(1:2), s(1:2), 'b');
%         for p=1:length(sumdat.mDuetAuto.rsNorm); plot(2.1, sumdat.mDuetAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
%         for p=1:length(sumdat.mSoloAuto.rsNorm); plot(1.1, sumdat.mSoloAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
%     plot([4 3], meanNorm(3:4), 'm.', 'MarkerSize', 16); 
%     errorbar([4 3], meanNorm(3:4), s(3:4), 'm' );
%         for p=1:length(sumdat.fDuetAuto.rsNorm); plot(4.1, sumdat.fDuetAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
%         for p=1:length(sumdat.fSoloAuto.rsNorm); plot(3.1, sumdat.fSoloAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
%     ylim([-5 40]); xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');
%     xticklabels({' ','S',' ','D',' ','S',' ','D',' '})
    
%% Get means and std for SENSORY plots

% For Normalized RS data

    meanNorm(1) = mean(sumdat(pp).mDuetHetero.rsNorm); s(1) = std(sumdat(pp).mDuetHetero.rsNorm);
    meanNorm(2) = mean(sumdat(pp).mSoloHetero.rsNorm); s(2) = std(sumdat(pp).mSoloHetero.rsNorm);
    meanNorm(3) = mean(sumdat(pp).fDuetHetero.rsNorm); s(3) = std(sumdat(pp).fDuetHetero.rsNorm);
    meanNorm(4) = mean(sumdat(pp).fSoloHetero.rsNorm); s(4) = std(sumdat(pp).fSoloHetero.rsNorm);

% For raw RS data

    meanRaw(1) = mean(sumdat(pp).mDuetHetero.rsRaw); sraw(1) = std(sumdat(pp).mDuetHetero.rsRaw);
    meanRaw(2) = mean(sumdat(pp).mSoloHetero.rsRaw); sraw(2) = std(sumdat(pp).mSoloHetero.rsRaw);
    meanRaw(3) = mean(sumdat(pp).fDuetHetero.rsRaw); sraw(3) = std(sumdat(pp).fDuetHetero.rsRaw);
    meanRaw(4) = mean(sumdat(pp).fSoloHetero.rsRaw); sraw(4) = std(sumdat(pp).fSoloHetero.rsRaw);

% For SPS data
    
    meanSPS(1) = mean(sumdat(pp).mDuetHetero.SPS); sps(1) = std(sumdat(pp).mDuetHetero.SPS);
    meanSPS(2) = mean(sumdat(pp).mSoloHetero.SPS); sps(2) = std(sumdat(pp).mSoloHetero.SPS);
    meanSPS(3) = mean(sumdat(pp).fDuetHetero.SPS); sps(3) = std(sumdat(pp).fDuetHetero.SPS);
    meanSPS(4) = mean(sumdat(pp).fSoloHetero.SPS); sps(4) = std(sumdat(pp).fSoloHetero.SPS);
    
figure(pp+1); clf; % RAW HETEROGENOUS PLOTS
hold on; title('Hetero SPS');
    plot([2 1], meanSPS(1:2), 'b.', 'MarkerSize', 16); 
    errorbar([2 1], meanSPS(1:2), sps(1:2), 'b' );
        for p=1:length(sumdat(pp).mDuetHetero.SPS); plot(2.1, sumdat(pp).mDuetHetero.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat(pp).mSoloHetero.SPS); plot(1.1, sumdat(pp).mSoloHetero.SPS(p), 'k.', 'MarkerSize', 8); end
    plot([4 3], meanSPS(3:4), 'm.', 'MarkerSize', 16); 
    errorbar([4 3], meanSPS(3:4), sps(3:4), 'm' );
        for p=1:length(sumdat(pp).fDuetHetero.SPS); plot(4.1, sumdat(pp).fDuetHetero.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat(pp).fSoloHetero.SPS); plot(3.1, sumdat(pp).fSoloHetero.SPS(p), 'k.', 'MarkerSize', 8); end
    % ylim([-10 65]); 
    xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');
    xticklabels({' ','S',' ','D',' ','S',' ','D',' '})

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
%     
    
%% Compute stats     
    
% Autogenous duet RS significant from 0?
    stts.m(pp).dNAuto.mean = mean(sumdat(pp).mDuetAuto.rsNorm);
    [stts.m(pp).dNAuto.H, stts.m(pp).dNAuto.P, stts.m(pp).dNAuto.CI, stts.m(pp).dNAuto.stats]  = ttest(sumdat(pp).mDuetAuto.rsNorm);
    stts.m(pp).dRAuto.mean = mean(sumdat.mDuetAuto.rsRaw);
    [stts.m(pp).dRAuto.H, stts.m(pp).dRAuto.P, stts.m(pp).dRAuto.CI, stts.m(pp).dRAuto.stats]  = ttest(sumdat(pp).mDuetAuto.rsRaw);
%    stts.m.dSAuto.mean = mean(sumdat.mDuetAuto.SPS);
%    [stts.m.dSAuto.H, stts.m.dSAuto.P, stts.m.dSAuto.CI, stts.m.dSAuto.stats]  = ttest(sumdat.mDuetAuto.SPS);
    
    stts.f(pp).dNAuto.mean = mean(sumdat(pp).fDuetAuto.rsNorm);
    [stts.f(pp).dNAuto.H, stts.f(pp).dNAuto.P, stts.f(pp).dNAuto.CI, stts.f(pp).dNAuto.stats]  = ttest(sumdat(pp).fDuetAuto.rsNorm);
    stts.f(pp).dRAuto.mean = mean(sumdat(pp).fDuetAuto.rsRaw);
    [stts.f(pp).dRAuto.H, stts.f(pp).dRAuto.P, stts.f(pp).dRAuto.CI, stts.f(pp).dRAuto.stats]  = ttest(sumdat(pp).fDuetAuto.rsRaw);
%    stts.f.dSAuto.mean = mean(sumdat.fDuetAuto.SPS);
%    [stts.f.dSAuto.H, stts.f.dSAuto.P, stts.f.dSAuto.CI, stts.f.dSAuto.stats]  = ttest(sumdat.mDuetAuto.SPS);
    
% Autogenous Solo RS significant from 0?

    stts.m(pp).sNAuto.mean = mean(sumdat(pp).mSoloAuto.rsNorm);
    [stts.m(pp).sNAuto.H, stts.m(pp).sNAuto.P, stts.m(pp).sNAuto.CI, stts.m(pp).sNAuto.stats]  = ttest(sumdat(pp).mSoloAuto.rsNorm);
    stts.m(pp).sRAuto.mean = mean(sumdat(pp).mSoloAuto.rsRaw);
    [stts.m(pp).sRAuto.H, stts.m(pp).sRAuto.P, stts.m(pp).sRAuto.CI, stts.m(pp).sRAuto.stats]  = ttest(sumdat(pp).mSoloAuto.rsRaw);
    
    stts.f.sNAuto.mean = mean(sumdat.fSoloAuto.rsNorm);
    [stts.f(pp).sNAuto.H, stts.f(pp).sNAuto.P, stts.f(pp).sNAuto.CI, stts.f(pp).sNAuto.stats]  = ttest(sumdat(pp).fSoloAuto.rsNorm);
    stts.m(pp).sRAuto.mean = mean(sumdat(pp).fSoloAuto.rsRaw);
    [stts.f(pp).sRAuto.H, stts.f(pp).sRAuto.P, stts.f(pp).sRAuto.CI, stts.f(pp).sRAuto.stats]  = ttest(sumdat(pp).fSoloAuto.rsRaw);

% Difference between Autogenous Duet and Solo ?

% [stts.m.SvsDNAuto.H, stts.m.SvsDNAuto.P, stts.m.SvsDNAuto.CI, stts.m.SvsDNAuto.stats]  = ttest2(sumdat.mSoloAuto.rsNorm, sumdat.mDuetAuto.rsNorm, 'vartype', 'unequal');
% [stts.m.SvsDRAuto.H, stts.m.SvsDRAuto.P, stts.m.SvsDRAuto.CI, stts.m.SvsDRAuto.stats]  = ttest2(sumdat.mSoloAuto.rsRaw, sumdat.mDuetAuto.rsRaw, 'vartype', 'unequal');
% [stts.f.SvsDNAuto.H, stts.f.SvsDNAuto.P, stts.f.SvsDNAuto.CI, stts.f.SvsDNAuto.stats]  = ttest2(sumdat.fSoloAuto.rsNorm, sumdat.fDuetAuto.rsNorm, 'vartype', 'unequal');
% [stts.f.SvsDRAuto.H, stts.f.SvsDRAuto.P, stts.f.SvsDRAuto.CI, stts.f.SvsDRAuto.stats]  = ttest2(sumdat.fSoloAuto.rsRaw, sumdat.fDuetAuto.rsRaw, 'vartype', 'unequal');

[stts.m(pp).SvsDNAuto.P, stts.m(pp).SvsDNAuto.H, stts.m(pp).SvsDNAuto.stats]  = ranksum(sumdat(pp).mSoloAuto.rsNorm, sumdat(pp).mDuetAuto.rsNorm);
[stts.m(pp).SvsDRAuto.P, stts.m(pp).SvsDRAuto.H, stts.m(pp).SvsDRAuto.stats]  = ranksum(sumdat(pp).mSoloAuto.rsRaw, sumdat(pp).mDuetAuto.rsRaw);
[stts.m(pp).SvsDSAuto.P, stts.m(pp).SvsDSAuto.H, stts.m(pp).SvsDSAuto.stats]  = ranksum(sumdat(pp).mSoloAuto.SPS, sumdat(pp).mDuetAuto.SPS);

[stts.f(pp).SvsDNAuto.P, stts.f(pp).SvsDNAuto.H, stts.f(pp).SvsDNAuto.stats]  = ranksum(sumdat(pp).fSoloAuto.rsNorm, sumdat(pp).fDuetAuto.rsNorm);
[stts.f(pp).SvsDRAuto.P, stts.f(pp).SvsDRAuto.H, stts.f(pp).SvsDRAuto.stats]  = ranksum(sumdat(pp).fSoloAuto.rsRaw, sumdat(pp).fDuetAuto.rsRaw);
[stts.f(pp).SvsDSAuto.P, stts.f(pp).SvsDSAuto.H, stts.f(pp).SvsDSAuto.stats]  = ranksum(sumdat(pp).fSoloAuto.SPS, sumdat(pp).fDuetAuto.SPS);


% Heterogenous duet RS significant from 0?

    stts.m(pp).dNHetero.mean = mean(sumdat(pp).mDuetHetero.rsNorm);
    [stts.m(pp).dNHetero.H, stts.m(pp).dNHetero.P, stts.m(pp).dNHetero.CI, stts.m(pp).dNHetero.stats]  = ttest(sumdat(pp).mDuetHetero.rsNorm);
    stts.m(pp).dRHetero.mean = mean(sumdat.mDuetHetero.rsRaw);
    [stts.m(pp).dRHetero.H, stts.m(pp).dRHetero.P, stts.m(pp).dRHetero.CI, stts.m(pp).dRHetero.stats]  = ttest(sumdat(pp).mDuetHetero.rsRaw);
    
    stts.f(pp).dNHetero.mean = mean(sumdat(pp).fDuetHetero.rsNorm);
    [stts.f(pp).dNHetero.H, stts.f(pp).dNHetero.P, stts.f(pp).dNHetero.CI, stts.f(pp).dNHetero.stats]  = ttest(sumdat(pp).fDuetHetero.rsNorm);
    stts.f(pp).dRHetero.mean = mean(sumdat(pp).fDuetHetero.rsRaw);
    [stts.f(pp).dRHetero.H, stts.f(pp).dRHetero.P, stts.f(pp).dRHetero.CI, stts.f(pp).dRHetero.stats]  = ttest(sumdat(pp).fDuetHetero.rsRaw);
    
% Heterogenous Solo RS significant from 0?

    stts.m(pp).sNHetero.mean = mean(sumdat(pp).mSoloHetero.rsNorm);
    [stts.m(pp).sNHetero.H, stts.m(pp).sNHetero.P, stts.m(pp).sNHetero.CI, stts.m(pp).sNHetero.stats]  = ttest(sumdat(pp).mSoloHetero.rsNorm);
    stts.m(pp).sRHetero.mean = mean(sumdat(pp).mSoloHetero.rsRaw);
    [stts.m(pp).sRHetero.H, stts.m(pp).sRHetero.P, stts.m(pp).sRHetero.CI, stts.m(pp).sRHetero.stats]  = ttest(sumdat(pp).mSoloHetero.rsRaw);
    
    stts.f(pp).sNHetero.mean = mean(sumdat(pp).fSoloHetero.rsNorm);
    [stts.f(pp).sNHetero.H, stts.f(pp).sNHetero.P, stts.f(pp).sNHetero.CI, stts.f(pp).sNHetero.stats]  = ttest(sumdat(pp).fSoloHetero.rsNorm);
    stts.f(pp).sRHetero.mean = mean(sumdat(pp).fSoloHetero.rsRaw);
    [stts.f(pp).sRHetero.H, stts.f(pp).sRHetero.P, stts.f(pp).sRHetero.CI, stts.f(pp).sRHetero.stats]  = ttest(sumdat(pp).fSoloHetero.rsRaw);

% Difference between Heterogenous Duet and Solo RS motor?

% [stts.m.SvsDNHetero.H, stts.m.SvsDNHetero.P, stts.m.SvsDNHetero.CI, stts.m.SvsDNHetero.stats]  = ttest2(sumdat.mSoloHetero.rsNorm, sumdat.mDuetHetero.rsNorm, 'vartype', 'unequal');
% [stts.m.SvsDRHetero.H, stts.m.SvsDRHetero.P, stts.m.SvsDRHetero.CI, stts.m.SvsDRHetero.stats]  = ttest2(sumdat.mSoloHetero.rsRaw, sumdat.mDuetHetero.rsRaw, 'vartype', 'unequal');
% [stts.f.SvsDNHetero.H, stts.f.SvsDNHetero.P, stts.f.SvsDNHetero.CI, stts.f.SvsDNHetero.stats]  = ttest2(sumdat.fSoloHetero.rsNorm, sumdat.fDuetHetero.rsNorm, 'vartype', 'unequal');
% [stts.f.SvsDRHetero.H, stts.f.SvsDRHetero.P, stts.f.SvsDRHetero.CI, stts.f.SvsDRHetero.stats]  = ttest2(sumdat.fSoloHetero.rsRaw, sumdat.fDuetHetero.rsRaw, 'vartype', 'unequal');

[stts.m(pp).SvsDNHetero.P, stts.m(pp).SvsDNHetero.H, stts.m(pp).SvsDNHetero.stats]  = ranksum(sumdat(pp).mSoloHetero.rsNorm, sumdat(pp).mDuetHetero.rsNorm);
[stts.m(pp).SvsDRHetero.P, stts.m(pp).SvsDRHetero.H, stts.m(pp).SvsDRHetero.stats]  = ranksum(sumdat(pp).mSoloHetero.rsRaw, sumdat(pp).mDuetHetero.rsRaw);
[stts.m(pp).SvsDSHetero.P, stts.m(pp).SvsDSHetero.H, stts.m(pp).SvsDSHetero.stats]  = ranksum(sumdat(pp).mSoloHetero.SPS, sumdat(pp).mDuetHetero.SPS);

[stts.f(pp).SvsDNHetero.P, stts.f(pp).SvsDNHetero.H, stts.f(pp).SvsDNHetero.stats]  = ranksum(sumdat(pp).fSoloHetero.rsNorm, sumdat(pp).fDuetHetero.rsNorm);
[stts.f(pp).SvsDRHetero.P, stts.f(pp).SvsDRHetero.H, stts.f(pp).SvsDRHetero.stats]  = ranksum(sumdat(pp).fSoloHetero.rsRaw, sumdat(pp).fDuetHetero.rsRaw);
[stts.f(pp).SvsDSHetero.P, stts.f(pp).SvsDSHetero.H, stts.f(pp).SvsDSHetero.stats]  = ranksum(sumdat(pp).fSoloHetero.SPS, sumdat(pp).fDuetHetero.SPS);

[MvFHS.H, MvFHS.P, MvFHS.CI, MvFHS.stats]  = ttest2(sumdat(pp).fSoloHetero.rsRaw, sumdat(pp).mSoloHetero.rsRaw, 'vartype', 'unequal');

% Test for equal variance

MaleDuetSoloC = [sumdat(pp).mDuetAuto.rsRaw, sumdat(pp).mSoloAuto.rsRaw];
MaleDuetSoloCIDX = [ones(1,length(sumdat(pp).mDuetAuto.rsRaw)), 2*ones(1,length(sumdat(pp).mSoloAuto.rsRaw))];
[MaleVarAutoChronP,MaleVarAutoChronstats] = vartestn(MaleDuetSoloC', MaleDuetSoloCIDX','TestType','LeveneAbsolute');

FeMaleDuetSoloC = [sumdat(pp).fDuetAuto.rsRaw, sumdat(pp).fSoloAuto.rsRaw];
FeMaleDuetSoloCIDX = [ones(1,length(sumdat(pp).fDuetAuto.rsRaw)), 2*ones(1,length(sumdat(pp).fSoloAuto.rsRaw))];
[FemaleVarAutoChronP,FemaleVarAutoChronstats] = vartestn(FeMaleDuetSoloC', FeMaleDuetSoloCIDX','TestType','LeveneAbsolute');

FvsMAutoSoloC = [sumdat(pp).fSoloAuto.rsRaw, sumdat(pp).mSoloAuto.rsRaw];
FvsMAutoSoloCIDX = [ones(1,length(sumdat(pp).mSoloAuto.rsRaw)), 2*ones(1,length(sumdat(pp).fSoloAuto.rsRaw))];
[FvsMVarAutoChronP,FvsMVarAutoChronstats] = vartestn(FvsMAutoSoloC', FvsMAutoSoloCIDX','TestType','LeveneAbsolute');


    fprintf('Male Auto Duet Raw RS vs Solo? p = %1.5f \n', stts.m(pp).SvsDRAuto.P);
    fprintf('Female Auto Duet Raw RS vs Solo? p = %1.5f \n', stts.f(pp).SvsDRAuto.P);
    fprintf('Male Auto SPS Duet vs Solo? p = %1.5f \n', stts.m(pp).SvsDSAuto.P);
    fprintf('Female Auto SPS Duet vs Solo? p = %1.5f \n', stts.f(pp).SvsDSAuto.P);
    
    fprintf(' \n');
    
    fprintf('Male Auto Duet Raw RS different from zero? p = %1.5f \n', stts.m(pp).dRAuto.P);
    fprintf('Male Auto Solo Raw RS different from zero? p = %1.5f \n', stts.m(pp).sRAuto.P);
    fprintf('Female Auto Duet Raw RS different from zero? p = %1.5f \n', stts.f(pp).dRAuto.P);
    fprintf('Female Auto Solo Raw RS different from zero? p = %1.5f \n', stts.f(pp).sRAuto.P);
    
    fprintf(' \n');

    fprintf('Male Hetero Duet Raw RS different from zero? p = %1.5f \n', stts.m(pp).dRHetero.P);
    fprintf('Male Hetero Solo Raw RS different from zero? p = %1.5f \n', stts.m(pp).sRHetero.P);
    fprintf('Female Hetero Duet Raw RS different from zero? p = %1.5f \n', stts.f(pp).dRHetero.P);
    fprintf('Female Hetero Solo Raw RS different from zero? p = %1.5f \n', stts.f(pp).sRHetero.P);
    
    fprintf(' \n');

    fprintf('Male Hetero Duet Raw RS vs Solo? p = %1.5f \n', stts.m(pp).SvsDRHetero.P);
    fprintf('Female Hetero Duet Raw RS vs Solo? p = %1.5f \n', stts.f(pp).SvsDRHetero.P);
    fprintf('Male Hetero Duet SPS vs Solo? p = %1.5f \n', stts.m(pp).SvsDSHetero.P);
    fprintf('Female Hetero Duet SPS vs Solo? p = %1.5f \n', stts.f(pp).SvsDSHetero.P);
    
    fprintf('Male versus Female Solo Hetero Raw RS different? p = %1.5f \n', MvFHS.P);
    
    fprintf(' \n');

    fprintf('Male Chronic Auto: Variance difference between Duet and Solo Raw RS? p = %1.5f \n', MaleVarAutoChronP);
    fprintf('FeMale Chronic Auto: Variance difference between Duet and Solo Raw RS? p = %1.5f \n', FemaleVarAutoChronP);
    
fprintf('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n');


    fprintf('Male Auto Duet Norm RS vs Solo? p = %1.5f \n', stts.m(pp).SvsDNAuto.P);
    fprintf('Female Auto Duet Norm RS vs Solo? p = %1.5f \n', stts.f(pp).SvsDNAuto.P);
    
    fprintf(' \n');

    fprintf('Male Hetero Duet Norm RS different from zero? p = %1.5f \n', stts.m(pp).dNHetero.P);
    fprintf('Male Hetero Solo Norm RS different from zero? p = %1.5f \n', stts.m(pp).sNHetero.P);
    fprintf('Female Hetero Duet Norm RS different from zero? p = %1.5f \n', stts.f(pp).dNHetero.P);
    fprintf('Female Hetero Solo Norm RS different from zero? p = %1.5f \n', stts.f(pp).sNHetero.P);
    
    fprintf(' \n');

    fprintf('Male Hetero Duet Norm RS vs Solo? p = %1.5f \n', stts.m(pp).SvsDNHetero.P);
    fprintf('Female Hetero Duet Norm RS vs Solo? p = %1.5f \n', stts.f(pp).SvsDNHetero.P);
        
    fprintf(' \n');

end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%% Response Strength nested function

function qwe = rs(struc, syllabl, spontan, padme)
 
   % Get spontaneous rate 
        sponSpikeCount = 0;  
    
        for i=1:4 % 4 electrodes in a tetrode always
            sponSpikeCount = sponSpikeCount + length(find(struc.Cspikes{i} > spontan(1) & struc.Cspikes{i} < spontan(2)));
        end

        sponSPS = sponSpikeCount / (spontan(2) - spontan(1)); % This is spikes per second
        sponSPS = sponSPS/4; % Divide by 4 because we have 4 electrodes

        skinny = 0; 
        if padme > 0 % We have a heterogenous syllable. Because premotor is advanced and auditory feedback 
                     % is delayed, we need to truncate the window to avoid
                     % counting premotor spikes.  We are also truncating
                     % the solo heterogenous for a fair comparison.
            skinny = 2.5 * padme; % IMPORTANT, try 1.5 instead of 2 for a more conservative truncation.
        end
                     
    % Loop for each syllable
    
    for j = 1:length(syllabl)    
    
        % Get the spikes for that syllable
        stimSpikeCount = 0; 
    
        for i=1:4 % 4 electrodes in a tetrode always (padme shifts the window in seconds, negative earlier, positive later)
            stimSpikeCount = stimSpikeCount + length(find(struc.Cspikes{i} >= struc.syl(syllabl(j)).tim(1)+padme & struc.Cspikes{i} < struc.syl(syllabl(j)).tim(2)+padme-skinny));
        end
        
        stimSPS = stimSpikeCount / ((struc.syl(syllabl(j)).tim(2) - struc.syl(syllabl(j)).tim(1)) - skinny); % This is spikes per second
        stimSPS = stimSPS/4; % Divide by 4 because we have 4 electrodes
        
        
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





function [out, sumdat, stats] = wRS_Chronic(in, padding)
% Usage: Calculates response strength to solo and duet syllables.
% Relies on rs, a nested function below, to calculate Response Strength.
% Load the Chronic data structure first:
% load ChronicCompleat2017p.mat (OLD)
% load ChronicCompleat2018a.mat (Current as of 4-Aug-2018)

% 'pad' is a critical variable - it is the shift in the window around the
% clicked boundaries of the syllables for the calculation of RS. 
% The boundaries of each syllable are used for this analysis. This is
% inherently problematic as we expect pre-motor activity in awake animals to occur PRIOR to
% the sound and auditory activity in urethane-anesthetized animals to occur AFTER the sound.
% We are comfortable with using a value of '0' it is a rather unbiased. Change
% the value of padding in seconds (e.g. 0.020 or -0.030) to look at the effects on the results.
pad = 0.000; 

% The user can specify the padding via an argin for convenience.
if nargin == 2; pad = padding; end


% Initializing variables
sumdat.fDuetAuto.rsNorm = []; sumdat.fDuetAuto.rsRaw = [];
sumdat.mDuetAuto.rsNorm = []; sumdat.mDuetAuto.rsRaw = [];
sumdat.fDuetHetero.rsNorm = []; sumdat.fDuetHetero.rsRaw = [];
sumdat.mDuetHetero.rsNorm = []; sumdat.mDuetHetero.rsRaw = [];
sumdat.fSoloAuto.rsNorm = []; sumdat.fSoloAuto.rsRaw = [];
sumdat.mSoloAuto.rsNorm = []; sumdat.mSoloAuto.rsRaw = [];
sumdat.fSoloHetero.rsNorm = []; sumdat.fSoloHetero.rsRaw = [];
sumdat.mSoloHetero.rsNorm = []; sumdat.mSoloHetero.rsRaw = [];

sumdat.mSoloHetero.SPS = [];
sumdat.fSoloHetero.SPS = [];
sumdat.mDuetHetero.SPS = [];
sumdat.fDuetHetero.SPS = [];
sumdat.mSoloAuto.SPS = [];
sumdat.fSoloAuto.SPS = [];
sumdat.mDuetAuto.SPS = [];
sumdat.fDuetAuto.SPS = [];

%% List of Chronic singing data with syllable indices and locations for spontaneous activity

[msolosyls, mduetsyls, fsolosyls, fduetsyls, spon] = wData;

%% Loop to calculate RS values for each pair of wrens   
    
for curpair = 1:length(spon)
    
    % Solo syllables MALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(msolosyls{curpair}) % Male sang solo syllables
        
        % Calculate RS values
        out(curpair).fSoloHetero = rs(in(curpair*2), msolosyls{curpair}, spon(:,curpair), pad);
        out(curpair).mSoloAuto = rs(in((curpair*2)-1), msolosyls{curpair}, spon(:,curpair), pad);
        
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
        
        out(curpair).mSoloHetero = rs(in((curpair*2)-1), fsolosyls{curpair}, spon(:,curpair), pad);
        out(curpair).fSoloAuto = rs(in(curpair*2), fsolosyls{curpair}, spon(:,curpair), pad);
        
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
        
        out(curpair).mDuetAuto = rs(in((curpair*2)-1), mduetsyls{curpair}, spon(:,curpair), pad);
        out(curpair).fDuetHetero = rs(in(curpair*2), mduetsyls{curpair}, spon(:,curpair), pad);

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
        
        out(curpair).mDuetHetero = rs(in((curpair*2)-1), fduetsyls{curpair}, spon(:,curpair), pad);
        out(curpair).fDuetAuto = rs(in(curpair*2), fduetsyls{curpair}, spon(:,curpair), pad);

        for kk = 1:length(fduetsyls{curpair})
            sumdat.mDuetHetero.rsNorm(end+1) = out(curpair).mDuetHetero(kk).rsNorm;
            sumdat.mDuetHetero.rsRaw(end+1) = out(curpair).mDuetHetero(kk).rsRaw;
            sumdat.mDuetHetero.SPS(end+1) = out(curpair).mDuetHetero(kk).spikerate;
            sumdat.fDuetAuto.rsNorm(end+1) = out(curpair).fDuetAuto(kk).rsNorm;
            sumdat.fDuetAuto.rsRaw(end+1) = out(curpair).fDuetAuto(kk).rsRaw;
            sumdat.fDuetAuto.SPS(end+1) = out(curpair).fDuetAuto(kk).spikerate;
        end
    end

end


%% Plot MOTOR

    meanNorm(1) = mean(sumdat.mDuetAuto.rsNorm); s(1) = std(sumdat.mDuetAuto.rsNorm);
    meanNorm(2) = mean(sumdat.mSoloAuto.rsNorm); s(2) = std(sumdat.mSoloAuto.rsNorm);
    meanNorm(3) = mean(sumdat.fDuetAuto.rsNorm); s(3) = std(sumdat.fDuetAuto.rsNorm);
    meanNorm(4) = mean(sumdat.fSoloAuto.rsNorm); s(4) = std(sumdat.fSoloAuto.rsNorm);

% For raw RS data

    meanRaw(1) = mean(sumdat.mDuetAuto.rsRaw); sraw(1) = std(sumdat.mDuetAuto.rsRaw);
    meanRaw(2) = mean(sumdat.mSoloAuto.rsRaw); sraw(3) = std(sumdat.mSoloAuto.rsRaw);
    meanRaw(3) = mean(sumdat.fDuetAuto.rsRaw); sraw(2) = std(sumdat.fDuetAuto.rsRaw);
    meanRaw(4) = mean(sumdat.fSoloAuto.rsRaw); sraw(4) = std(sumdat.fSoloAuto.rsRaw);

% For raw SPS data
    
    meanSPS(1) = mean(sumdat.mDuetAuto.SPS); sps(1) = std(sumdat.mDuetAuto.SPS);
    meanSPS(2) = mean(sumdat.mSoloAuto.SPS); sps(2) = std(sumdat.mSoloAuto.SPS);
    meanSPS(3) = mean(sumdat.fDuetAuto.SPS); sps(3) = std(sumdat.fDuetAuto.SPS);
    meanSPS(4) = mean(sumdat.fSoloAuto.SPS); sps(4) = std(sumdat.fSoloAuto.SPS);
    
figure(21); clf; 

subplot(131); hold on; title('Auto Normalized RS');
    plot([1 2], meanNorm(1:2), 'bo'); 
    errorbar([1 2], meanNorm(1:2), s(1:2), 'b' );
        for p=1:length(sumdat.mDuetAuto.rsNorm); plot(1.1, sumdat.mDuetAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.mSoloAuto.rsNorm); plot(2.1, sumdat.mSoloAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
    plot([3 4], meanNorm(3:4), 'mo'); 
    errorbar([3 4], meanNorm(3:4), s(3:4), 'm' );
        for p=1:length(sumdat.fDuetAuto.rsNorm); plot(3.1, sumdat.fDuetAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.fSoloAuto.rsNorm); plot(4.1, sumdat.fSoloAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
    ylim([-5 40]); xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');

subplot(132); hold on; title('Auto Raw RS');
    plot([1 2], meanRaw(1:2), 'bo'); 
    errorbar([1 2], meanRaw(1:2), sraw(1:2), 'b' );
        for p=1:length(sumdat.mDuetAuto.rsRaw); plot(1.1, sumdat.mDuetAuto.rsRaw(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.mSoloAuto.rsRaw); plot(2.1, sumdat.mSoloAuto.rsRaw(p), 'k.', 'MarkerSize', 8); end
    plot([3 4], meanRaw(3:4), 'mo'); 
    errorbar([3 4], meanRaw(3:4), sraw(3:4), 'm' );
        for p=1:length(sumdat.fDuetAuto.rsRaw); plot(3.1, sumdat.fDuetAuto.rsRaw(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.fSoloAuto.rsRaw); plot(4.1, sumdat.fSoloAuto.rsRaw(p), 'k.', 'MarkerSize', 8); end
    ylim([-10 100]); xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');

subplot(133); hold on; title('Auto Spikes/Second');
    plot([1 2], meanSPS(1:2), 'bo'); 
    errorbar([1 2], meanSPS(1:2), sps(1:2), 'b' );
        for p=1:length(sumdat.mDuetAuto.SPS); plot(1.1, sumdat.mDuetAuto.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.mSoloAuto.SPS); plot(2.1, sumdat.mSoloAuto.SPS(p), 'k.', 'MarkerSize', 8); end
    plot([3 4], meanSPS(3:4), 'mo'); 
    errorbar([3 4], meanSPS(3:4), sps(3:4), 'm' );
        for p=1:length(sumdat.fDuetAuto.SPS); plot(3.1, sumdat.fDuetAuto.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.fSoloAuto.SPS); plot(4.1, sumdat.fSoloAuto.SPS(p), 'k.', 'MarkerSize', 8); end
    ylim([0 100]); xlim([0.5 4.5]); 

    
%% Plot SENSORY

    meanNorm(1) = mean(sumdat.mDuetHetero.rsNorm); s(1) = std(sumdat.mDuetHetero.rsNorm);
    meanNorm(2) = mean(sumdat.mSoloHetero.rsNorm); s(2) = std(sumdat.mSoloHetero.rsNorm);
    meanNorm(3) = mean(sumdat.fDuetHetero.rsNorm); s(3) = std(sumdat.fDuetHetero.rsNorm);
    meanNorm(4) = mean(sumdat.fSoloHetero.rsNorm); s(4) = std(sumdat.fSoloHetero.rsNorm);

% For raw RS data

    meanRaw(1) = mean(sumdat.mDuetHetero.rsRaw); sraw(1) = std(sumdat.mDuetHetero.rsRaw);
    meanRaw(2) = mean(sumdat.mSoloHetero.rsRaw); sraw(3) = std(sumdat.mSoloHetero.rsRaw);
    meanRaw(3) = mean(sumdat.fDuetHetero.rsRaw); sraw(2) = std(sumdat.fDuetHetero.rsRaw);
    meanRaw(4) = mean(sumdat.fSoloHetero.rsRaw); sraw(4) = std(sumdat.fSoloHetero.rsRaw);

% For raw SPS data
    
    meanSPS(1) = mean(sumdat.mDuetHetero.SPS); sps(1) = std(sumdat.mDuetHetero.SPS);
    meanSPS(2) = mean(sumdat.mSoloHetero.SPS); sps(2) = std(sumdat.mSoloHetero.SPS);
    meanSPS(3) = mean(sumdat.fDuetHetero.SPS); sps(3) = std(sumdat.fDuetHetero.SPS);
    meanSPS(4) = mean(sumdat.fSoloHetero.SPS); sps(4) = std(sumdat.fSoloHetero.SPS);
    
figure(22); clf; 

subplot(131); hold on; title('Hetero Normalized RS');
    plot([1 2], meanNorm(1:2), 'bo'); 
    errorbar([1 2], meanNorm(1:2), s(1:2), 'b' );
        for p=1:length(sumdat.mDuetHetero.rsNorm); plot(1.1, sumdat.mDuetHetero.rsNorm(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.mSoloHetero.rsNorm); plot(2.1, sumdat.mSoloHetero.rsNorm(p), 'k.', 'MarkerSize', 8); end
    plot([3 4], meanNorm(3:4), 'mo'); 
    errorbar([3 4], meanNorm(3:4), s(3:4), 'm' );
        for p=1:length(sumdat.fDuetHetero.rsNorm); plot(3.1, sumdat.fDuetHetero.rsNorm(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.fSoloHetero.rsNorm); plot(4.1, sumdat.fSoloHetero.rsNorm(p), 'k.', 'MarkerSize', 8); end
    ylim([-5 35]); xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');

subplot(132); hold on; title('Hetero Raw RS');
    plot([1 2], meanRaw(1:2), 'bo'); 
    errorbar([1 2], meanRaw(1:2), sraw(1:2), 'b' );
        for p=1:length(sumdat.mDuetHetero.rsRaw); plot(1.1, sumdat.mDuetHetero.rsRaw(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.mSoloHetero.rsRaw); plot(2.1, sumdat.mSoloHetero.rsRaw(p), 'k.', 'MarkerSize', 8); end
    plot([3 4], meanRaw(3:4), 'mo'); 
    errorbar([3 4], meanRaw(3:4), sraw(3:4), 'm' );
        for p=1:length(sumdat.fDuetHetero.rsRaw); plot(3.1, sumdat.fDuetHetero.rsRaw(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.fSoloHetero.rsRaw); plot(4.1, sumdat.fSoloHetero.rsRaw(p), 'k.', 'MarkerSize', 8); end
    ylim([-10 50]); xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');

subplot(133); hold on; title('Hetero Spikes/Second');
    plot([1 2], meanSPS(1:2), 'bo'); 
    errorbar([1 2], meanSPS(1:2), sps(1:2), 'b' );
        for p=1:length(sumdat.mDuetHetero.SPS); plot(1.1, sumdat.mDuetHetero.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.mSoloHetero.SPS); plot(2.1, sumdat.mSoloHetero.SPS(p), 'k.', 'MarkerSize', 8); end
    plot([3 4], meanSPS(3:4), 'mo'); 
    errorbar([3 4], meanSPS(3:4), sps(3:4), 'm' );
        for p=1:length(sumdat.fDuetHetero.SPS); plot(3.1, sumdat.fDuetHetero.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.fSoloHetero.SPS); plot(4.1, sumdat.fSoloHetero.SPS(p), 'k.', 'MarkerSize', 8); end
    ylim([-1 55]); xlim([0.5 4.5]); 


    stats= 0;
% subplot(132); hold on;
%     plot([1 2], mraw(1:2), 'b*'); 
%     errorbar([1 2], mraw(1:2), sraw(1:2), 'b' );% /sqrt(length(fChron)));
%     plot([3 4], mraw(3:4), 'm*'); hold on;
%     errorbar([3 4], mraw(3:4), sraw(3:4), 'm' );% /sqrt(length(mChron)));

%%%% SENSORY

% figure(2); clf; 
% 
% subplot(122); hold on;
%     plot(2, mean(sumdat.mAud.rsNorm), 'bo');
%         errorbar(2, mean(sumdat.mAud.rsNorm), std(sumdat.mAud.rsNorm), 'b');
%     plot(1, mean(sumdat.mduetHeterogenous.rsNorm), 'bo');
%         errorbar(1, mean(sumdat.mduetHeterogenous.rsNorm), std(sumdat.mduetHeterogenous.rsNorm), 'b');
%     plot(4, mean(sumdat.fAud.rsNorm), 'mo');
%         errorbar(4, mean(sumdat.fAud.rsNorm), std(sumdat.mAud.rsNorm), 'm');
%     plot(3, mean(sumdat.fduetHeterogenous.rsNorm), 'mo');
%         errorbar(3, mean(sumdat.fduetHeterogenous.rsNorm), std(sumdat.fduetHeterogenous.rsNorm), 'm');
%     ylim([-5 20]); xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');
% 
% %% Compute stats     
%     
% % Autogenous duet RS significant from 0?
%     stats.m.dNAuto.mean = mean(sumdat.mduetAutogenous.rsNorm);
%     [stats.m.dNAuto.H, stats.m.dNAuto.P, stats.m.dNAuto.CI, stats.m.dNAuto.stats]  = ttest(sumdat.mduetAutogenous.rsNorm);
%     stats.m.dRAuto.mean = mean(sumdat.mduetAutogenous.rsRaw);
%     [stats.m.dRAuto.H, stats.m.dRAuto.P, stats.m.dRAuto.CI, stats.m.dRAuto.stats]  = ttest(sumdat.mduetAutogenous.rsRaw);
%     stats.f.dNAuto.mean = mean(sumdat.fduetAutogenous.rsNorm);
%     [stats.f.dNAuto.H, stats.f.dNAuto.P, stats.f.dNAuto.CI, stats.f.dNAuto.stats]  = ttest(sumdat.fduetAutogenous.rsNorm);
%     stats.f.dRAuto.mean = mean(sumdat.fduetAutogenous.rsRaw);
%     [stats.f.dRAuto.H, stats.f.dRAuto.P, stats.f.dRAuto.CI, stats.f.dRAuto.stats]  = ttest(sumdat.fduetAutogenous.rsRaw);
%     
% % Autogenous Solo RS significant from 0?
% 
%     stats.m.sNAuto.mean = mean(sumdat.mSolo.rsNorm);
%     [stats.m.sNAuto.H, stats.m.sNAuto.P, stats.m.sNAuto.CI, stats.m.sNAuto.stats]  = ttest(sumdat.mSolo.rsNorm);
%     stats.m.sRAuto.mean = mean(sumdat.mSolo.rsRaw);
%     [stats.m.sRAuto.H, stats.m.sRAuto.P, stats.m.sRAuto.CI, stats.m.sRAuto.stats]  = ttest(sumdat.mSolo.rsRaw);
%     stats.f.sNAuto.mean = mean(sumdat.fSolo.rsNorm);
%     [stats.f.sNAuto.H, stats.f.sNAuto.P, stats.f.sNAuto.CI, stats.f.sNAuto.stats]  = ttest(sumdat.fSolo.rsNorm);
%     stats.m.sRAuto.mean = mean(sumdat.fSolo.rsRaw);
%     [stats.f.sRAuto.H, stats.f.sRAuto.P, stats.f.sRAuto.CI, stats.f.sRAuto.stats]  = ttest(sumdat.fSolo.rsRaw);
% 
% % Difference between Autogenous Duet and Solo RS motor?
% 
% [stats.m.SvsDNAuto.H, stats.m.SvsDNAuto.P, stats.m.SvsDNAuto.CI, stats.m.SvsDNAuto.stats]  = ttest2(sumdat.mSolo.rsNorm, sumdat.mduetAutogenous.rsNorm);
% [stats.m.SvsDRAuto.H, stats.m.SvsDRAuto.P, stats.m.SvsDRAuto.CI, stats.m.SvsDRAuto.stats]  = ttest2(sumdat.mSolo.rsRaw, sumdat.mduetAutogenous.rsRaw);
% [stats.f.SvsDNAuto.H, stats.f.SvsDNAuto.P, stats.f.SvsDNAuto.CI, stats.f.SvsDNAuto.stats]  = ttest2(sumdat.fSolo.rsNorm, sumdat.fduetAutogenous.rsNorm);
% [stats.f.SvsDRAuto.H, stats.f.SvsDRAuto.P, stats.f.SvsDRAuto.CI, stats.f.SvsDRAuto.stats]  = ttest2(sumdat.fSolo.rsRaw, sumdat.fduetAutogenous.rsRaw);
% 
% % Heterogenous duet RS significant from 0?
% 
%     stats.m.dNHetero.mean = mean(sumdat.mduetHeterogenous.rsNorm);
%     [stats.m.dNHetero.H, stats.m.dNHetero.P, stats.m.dNHetero.CI, stats.m.dNHetero.stats]  = ttest(sumdat.mduetHeterogenous.rsNorm);
%     stats.m.dRHetero.mean = mean(sumdat.mduetHeterogenous.rsRaw);
%     [stats.m.dRHetero.H, stats.m.dRHetero.P, stats.m.dRHetero.CI, stats.m.dRHetero.stats]  = ttest(sumdat.mduetHeterogenous.rsRaw);
%     stats.f.dNHetero.mean = mean(sumdat.fduetHeterogenous.rsNorm);
%     [stats.f.dNHetero.H, stats.f.dNHetero.P, stats.f.dNHetero.CI, stats.f.dNHetero.stats]  = ttest(sumdat.fduetHeterogenous.rsNorm);
%     stats.f.dRHetero.mean = mean(sumdat.fduetHeterogenous.rsRaw);
%     [stats.f.dRHetero.H, stats.f.dRHetero.P, stats.f.dRHetero.CI, stats.f.dRHetero.stats]  = ttest(sumdat.fduetHeterogenous.rsRaw);
%     
% % Heterogenous Solo RS significant from 0?
% 
%     stats.m.sNHetero.mean = mean(sumdat.mAud.rsNorm);
%     [stats.m.sNHetero.H, stats.m.sNHetero.P, stats.m.sNHetero.CI, stats.m.sNHetero.stats]  = ttest(sumdat.mAud.rsNorm);
%     stats.m.sRHetero.mean = mean(sumdat.mAud.rsRaw);
%     [stats.m.sRHetero.H, stats.m.sRHetero.P, stats.m.sRHetero.CI, stats.m.sRHetero.stats]  = ttest(sumdat.mAud.rsRaw);
%     stats.f.sNHetero.mean = mean(sumdat.fAud.rsNorm);
%     [stats.f.sNHetero.H, stats.f.sNHetero.P, stats.f.sNHetero.CI, stats.f.sNHetero.stats]  = ttest(sumdat.fAud.rsNorm);
%     stats.f.sRHetero.mean = mean(sumdat.fAud.rsRaw);
%     [stats.f.sRHetero.H, stats.f.sRHetero.P, stats.f.sRHetero.CI, stats.f.sRHetero.stats]  = ttest(sumdat.fAud.rsRaw);
% 
% % Difference between Heterogenous Duet and Solo RS motor?
% 
% [stats.m.SvsDNHetero.H, stats.m.SvsDNHetero.P, stats.m.SvsDNHetero.CI, stats.m.SvsDNHetero.stats]  = ttest2(sumdat.mAud.rsNorm, sumdat.mduetHeterogenous.rsNorm);
% [stats.m.SvsDRHetero.H, stats.m.SvsDRHetero.P, stats.m.SvsDRHetero.CI, stats.m.SvsDRHetero.stats]  = ttest2(sumdat.mAud.rsRaw, sumdat.mduetHeterogenous.rsRaw);
% [stats.f.SvsDNHetero.H, stats.f.SvsDNHetero.P, stats.f.SvsDNHetero.CI, stats.f.SvsDNHetero.stats]  = ttest2(sumdat.fAud.rsNorm, sumdat.fduetHeterogenous.rsNorm);
% [stats.f.SvsDRHetero.H, stats.f.SvsDRHetero.P, stats.f.SvsDRHetero.CI, stats.f.SvsDRHetero.stats]  = ttest2(sumdat.fAud.rsRaw, sumdat.fduetHeterogenous.rsRaw);


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
   
    % Loop for each syllable
    
        for j = 1:length(syllabl)    
    
        % Get the spikes for that syllable
        stimSpikeCount = 0; 
    
        for i=1:4 % 4 electrodes in a tetrode always
            stimSpikeCount = stimSpikeCount + length(find(struc.Cspikes{i} >= struc.syl(syllabl(j)).tim(1)-padme & struc.Cspikes{i} < struc.syl(syllabl(j)).tim(2)-padme));
        end
        
        stimSPS = stimSpikeCount / (struc.syl(syllabl(j)).tim(2) - struc.syl(syllabl(j)).tim(1)); % This is spikes per second
        stimSPS = stimSPS/4; % Divide by 4 because we have 4 electrodes
        
        
        qwe(j).sylnum = syllabl(j);
        qwe(j).rsNorm = (stimSPS - sponSPS) / sponSPS + 0.0000000000001; % Avoid divide by zero
        qwe(j).rsRaw = stimSPS - sponSPS;
        qwe(j).sponrate = sponSPS;
        qwe(j).spikerate = stimSPS;
        qwe(j).spontim = spontan;
        qwe(j).syltim = [struc.syl(syllabl(j)).tim(1), struc.syl(syllabl(j)).tim(2)];
        qwe(j).pad = padme;
        
    end
end

end




%% Award-winning (tm) code (not used!)
%     % Get firing rates during Auditory-only 
%         MaudSpikeCount = 0; Mauduration = 0; Msylcount = 0;
%         FaudSpikeCount = 0; Fauduration = 0; Fsylcount = 0;
% 
%     if ~isempty(fsolosyls{j})
%         for k = 1:length(fsolosyls{j})
%             for i=1:4
%                 MaudSpikeCount = MaudSpikeCount + length(find(w((j*2)-1).Cspikes{i} > w((j*2)-1).syl(fsolosyls{j}(k)).tim(1) & w((j*2)-1).Cspikes{i} < w((j*2)-1).syl(fsolosyls{j}(k)).tim(2)));
%             end    
%             Mauduration = Mauduration + (w((j*2)-1).syl(fsolosyls{j}(k)).tim(2) - w((j*2)-1).syl(fsolosyls{j}(k)).tim(1));
%             Msylcount = Msylcount + 1;
%         end
%         MaudSpikeRate(j) = MaudSpikeCount / Mauduration;
%         MnumSyls(j) = Msylcount;
%     end
%     if ~isempty(msolosyls{j})
%         for k = 1:length(msolosyls{j})
%             for i=1:4
%                 FaudSpikeCount = FaudSpikeCount + length(find(w((j*2)).Cspikes{i} > w((j*2)).syl(msolosyls{j}(k)).tim(1) & w((j*2)).Cspikes{i} < w((j*2)-1).syl(msolosyls{j}(k)).tim(2)));
%             end;
%             Fauduration = Fauduration + (w((j*2)).syl(msolosyls{j}(k)).tim(2) - w((j*2)).syl(msolosyls{j}(k)).tim(1));
%             Fsylcount = Fsylcount + 1;
%         end
%         FaudSpikeRate(j) = FaudSpikeCount / Fauduration;
%         FnumSyls(j) = Fsylcount;
%     end
% end

        
    
%end % End of embedded function RS
  


% figure(birdnum); clf; 

% xlimits = [1.0, 3.0]; % FOR 9-10
% xlimits = [1.2, 3.2]; % FOR 11-12
% xlimits = [2.6, 4.6]; % FOR 13-14 (no female urethane)
% 
% tt = find(w(birdnum).tim > xlimits(1) & w(birdnum).tim < xlimits(2));
% 
% mChist = bs_swPSTH(w(birdnum).Cspikes, [-5 10], 100, 0);
% fChist = bs_swPSTH(w(birdnum+1).Cspikes, [-5 10], 100, 0);
% 
% mAhist = bs_swPSTH(w(birdnum).Aspikes, [-5 10], 150, 0);
% fAhist = bs_swPSTH(w(birdnum+1).Aspikes, [-5 10], 150, 0);


% Show the first 5 syllables of a duet, color coded

% birdnum = 13;
% lastsyl = 8;
% 
%     figure(27); clf; 
%     subplot(212);
%     specgram(w(birdnum).duet, 1024, w(birdnum).Fs); xlim([0 (w(birdnum).syl(lastsyl).tim(1) - w(birdnum).tim(1))]);
%     subplot(211); hold on;
%     bs_raster(w(birdnum).Cspikes, 'b'); xlim([w(birdnum).tim(1) w(birdnum).syl(lastsyl).tim(1)]);
%     bs_raster(w(birdnum+1).Cspikes, 'm'); xlim([w(birdnum).tim(1) w(birdnum).syl(lastsyl).tim(1)]);


% subplot(411); specgram(w(i).duet(tt), 512, w(i).Fs); ylim([500 4500]); caxis([-50 20]);
% subplot(412); bs_raster(w(i).Cspikes, 'b'); xlim(xlimits);
% subplot(413); bs_raster(w(i+1).Cspikes, 'm'); xlim(xlimits);
% subplot(414); hold on; 
%     plot(mChist.tim, 1-(mChist.spers/max(mChist.spers)), 'b');
%     plot(fChist.tim, fChist.spers/max(fChist.spers), 'm');
%     xlim(xlimits);

    
% subplot(311); specgram(w(birdnum).duet(tt), 512, w(birdnum).Fs, [], 500); ylim([500 4500]); caxis([-25 15]);
%     colormap(flipud(gray));
% subplot(312); hold on; 
%     plot(mChist.tim, 1-(mChist.spers/max(mChist.spers)), 'b');
%     plot(fChist.tim, fChist.spers/max(fChist.spers), 'm');
%     xlim(xlimits);
% subplot(313); hold on; 
%     plot(mAhist.tim, 1-(mAhist.spers/max(mAhist.spers)), 'b');
%     plot(fAhist.tim, fAhist.spers/max(fAhist.spers), 'm');
%     xlim(xlimits);
    
    
%end;
    

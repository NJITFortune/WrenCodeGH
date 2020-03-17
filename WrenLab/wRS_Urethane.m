function [out, sumdat, stts] = wRS_Urethane(in, padding)
% Usage: Calculates response strength to solo and duet syllables.
% Relies on rs, a nested function below, to calculate Response Strength.
% Load the Chronic data structure first:
% load ChronicCompleat2017p.mat (OLD)
% load ChronicCompleat2018d.mat (Current as of 4-Jan-2019)

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

[msolosyls, mduetsyls, fsolosyls, fduetsyls, ~, Aspon] = wData;

%% Loop to calculate RS values for each pair of wrens   
    
for curpair = [1 2 3 4 5 6] % 1:length(Aspon) % length(spon)
    
    % Solo syllables MALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(msolosyls{curpair}) % Male sang solo syllables
        
        % Calculate RS values
        out(curpair).fSoloHetero = rs(in(curpair*2), msolosyls{curpair}, Aspon(:,curpair), pad);
        out(curpair).mSoloAuto = rs(in((curpair*2)-1), msolosyls{curpair}, Aspon(:,curpair), pad);
        
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
        
        out(curpair).mSoloHetero = rs(in((curpair*2)-1), fsolosyls{curpair}, Aspon(:,curpair), pad);
        out(curpair).fSoloAuto = rs(in(curpair*2), fsolosyls{curpair}, Aspon(:,curpair), pad);
        
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
        
        out(curpair).mDuetAuto = rs(in((curpair*2)-1), mduetsyls{curpair}, Aspon(:,curpair), pad);
        out(curpair).fDuetHetero = rs(in(curpair*2), mduetsyls{curpair}, Aspon(:,curpair), pad);

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
        
        out(curpair).mDuetHetero = rs(in((curpair*2)-1), fduetsyls{curpair}, Aspon(:,curpair), pad);
        out(curpair).fDuetAuto = rs(in(curpair*2), fduetsyls{curpair}, Aspon(:,curpair), pad);

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


%% Plot MALE

    meanNorm(1) = mean(sumdat.mDuetAuto.rsNorm); s(1) = std(sumdat.mDuetAuto.rsNorm);
    meanNorm(2) = mean(sumdat.mDuetHetero.rsNorm); s(2) = std(sumdat.mDuetHetero.rsNorm);
    meanNorm(3) = mean(sumdat.mSoloAuto.rsNorm); s(3) = std(sumdat.mSoloAuto.rsNorm);
    meanNorm(4) = mean(sumdat.mSoloHetero.rsNorm); s(4) = std(sumdat.mSoloHetero.rsNorm);

% For raw RS data

    meanRaw(1) = mean(sumdat.mDuetAuto.rsRaw); sraw(1) = std(sumdat.mDuetAuto.rsRaw);
    meanRaw(2) = mean(sumdat.mDuetHetero.rsRaw); sraw(2) = std(sumdat.mDuetHetero.rsRaw);
    meanRaw(3) = mean(sumdat.mSoloAuto.rsRaw); sraw(3) = std(sumdat.mSoloAuto.rsRaw);
    meanRaw(4) = mean(sumdat.mSoloHetero.rsRaw); sraw(4) = std(sumdat.mSoloHetero.rsRaw);

% For raw SPS data
    
    meanSPS(1) = mean(sumdat.mDuetAuto.SPS); sps(1) = std(sumdat.mDuetAuto.SPS);
    meanSPS(2) = mean(sumdat.mDuetHetero.SPS); sps(2) = std(sumdat.mDuetHetero.SPS);
    meanSPS(3) = mean(sumdat.mSoloAuto.SPS); sps(3) = std(sumdat.mSoloAuto.SPS);
    meanSPS(4) = mean(sumdat.mSoloHetero.SPS); sps(4) = std(sumdat.mSoloHetero.SPS);
    
figure(23); clf; 

subplot(131); hold on; title('Male Normalized RS');
    plot([1 2], meanNorm(1:2), 'bo'); 
    errorbar([1 2], meanNorm(1:2), s(1:2), 'b', 'LineWidth', 2);
        for p=1:length(sumdat.mDuetAuto.rsNorm); plot(1.1, sumdat.mDuetAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.mDuetHetero.rsNorm); plot(2.1, sumdat.mDuetHetero.rsNorm(p), 'k.', 'MarkerSize', 8); end
    plot([3 4], meanNorm(3:4), 'mo'); 
    errorbar([3 4], meanNorm(3:4), s(3:4), 'b', 'LineWidth', 1);
        for p=1:length(sumdat.mSoloAuto.rsNorm); plot(3.1, sumdat.mSoloAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.mSoloHetero.rsNorm); plot(4.1, sumdat.mSoloHetero.rsNorm(p), 'k.', 'MarkerSize', 8); end
    % ylim([-5 40]); 
    xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');

subplot(132); hold on; title('Male Raw RS');
    plot([1 2], meanRaw(1:2), 'bo'); 
    errorbar([1 2], meanRaw(1:2), sraw(1:2), 'b', 'LineWidth', 2);
        for p=1:length(sumdat.mDuetAuto.rsRaw); plot(1.1, sumdat.mDuetAuto.rsRaw(p), 'w.', 'MarkerSize', 8); end
        for p=1:length(sumdat.mDuetHetero.rsRaw); plot(2.1, sumdat.mDuetHetero.rsRaw(p), 'w.', 'MarkerSize', 8); end
        for p=1:length(sumdat.mDuetAuto.rsRaw); text(1.1, sumdat.mDuetAuto.rsRaw(p), num2str(p)); end
        for p=1:length(sumdat.mDuetHetero.rsRaw); text(2.1, sumdat.mDuetHetero.rsRaw(p), num2str(p)); end
    plot([3 4], meanRaw(3:4), 'mo'); 
    errorbar([3 4], meanRaw(3:4), sraw(3:4), 'b', 'LineWidth', 1);
%         for p=1:length(sumdat.mSoloAuto.rsRaw); plot(3.1, sumdat.mSoloAuto.rsRaw(p), 'k.', 'MarkerSize', 8); end
%         for p=1:length(sumdat.mSoloHetero.rsRaw); plot(4.1, sumdat.mSoloHetero.rsRaw(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.mSoloAuto.rsRaw); text(3.1, sumdat.mSoloAuto.rsRaw(p), num2str(p)); end
        for p=1:length(sumdat.mSoloHetero.rsRaw); text(4.1, sumdat.mSoloHetero.rsRaw(p), num2str(p)); end
    % ylim([-10 100]); 
    xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');
    xticklabels({'DA','DH','SA','SH'})

subplot(133); hold on; title('Male Spikes/Second');
    plot([1 2], meanSPS(1:2), 'bo'); 
    errorbar([1 2], meanSPS(1:2), sps(1:2), 'b', 'LineWidth', 2);
        for p=1:length(sumdat.mDuetAuto.SPS); plot(1.1, sumdat.mDuetAuto.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.mDuetHetero.SPS); plot(2.1, sumdat.mDuetHetero.SPS(p), 'k.', 'MarkerSize', 8); end
    plot([3 4], meanSPS(3:4), 'mo'); 
    errorbar([3 4], meanSPS(3:4), sps(3:4), 'b', 'LineWidth', 1);
        for p=1:length(sumdat.mSoloAuto.SPS); plot(3.1, sumdat.mSoloAuto.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.mSoloHetero.SPS); plot(4.1, sumdat.mSoloHetero.SPS(p), 'k.', 'MarkerSize', 8); end
    % ylim([0 100]); 
    xlim([0.5 4.5]); 

    
%% Plot FEMALE

    meanNorm(1) = mean(sumdat.fDuetAuto.rsNorm); s(1) = std(sumdat.fDuetAuto.rsNorm);
    meanNorm(2) = mean(sumdat.fDuetHetero.rsNorm); s(2) = std(sumdat.fDuetHetero.rsNorm);
    meanNorm(3) = mean(sumdat.fSoloAuto.rsNorm); s(3) = std(sumdat.fSoloAuto.rsNorm);
    meanNorm(4) = mean(sumdat.fSoloHetero.rsNorm); s(4) = std(sumdat.fSoloHetero.rsNorm);

% For raw RS data

    meanRaw(1) = mean(sumdat.fDuetAuto.rsRaw); sraw(1) = std(sumdat.fDuetAuto.rsRaw);
    meanRaw(2) = mean(sumdat.fDuetHetero.rsRaw); sraw(2) = std(sumdat.fDuetHetero.rsRaw);
    meanRaw(3) = mean(sumdat.fSoloAuto.rsRaw); sraw(3) = std(sumdat.fSoloAuto.rsRaw);
    meanRaw(4) = mean(sumdat.fSoloHetero.rsRaw); sraw(4) = std(sumdat.fSoloHetero.rsRaw);

% For raw SPS data
    
    meanSPS(1) = mean(sumdat.fDuetAuto.SPS); sps(1) = std(sumdat.fDuetAuto.SPS);
    meanSPS(2) = mean(sumdat.fDuetHetero.SPS); sps(2) = std(sumdat.fDuetHetero.SPS);
    meanSPS(3) = mean(sumdat.fSoloAuto.SPS); sps(3) = std(sumdat.fSoloAuto.SPS);
    meanSPS(4) = mean(sumdat.fSoloHetero.SPS); sps(4) = std(sumdat.fSoloHetero.SPS);

    
figure(24); clf; 

subplot(131); hold on; title('Female Normalized RS');
    plot([1 2], meanNorm(1:2), 'bo'); 
    errorbar([1 2], meanNorm(1:2), s(1:2), 'm', 'LineWidth', 2);
        for p=1:length(sumdat.fDuetAuto.rsNorm); plot(1.1, sumdat.fDuetAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.fDuetHetero.rsNorm); plot(2.1, sumdat.fDuetHetero.rsNorm(p), 'k.', 'MarkerSize', 8); end
    plot([3 4], meanNorm(3:4), 'mo'); 
    errorbar([3 4], meanNorm(3:4), s(3:4), 'm', 'LineWidth', 1);
        for p=1:length(sumdat.fSoloAuto.rsNorm); plot(3.1, sumdat.fSoloAuto.rsNorm(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.fSoloHetero.rsNorm); plot(4.1, sumdat.fSoloHetero.rsNorm(p), 'k.', 'MarkerSize', 8); end
    % ylim([-5 35]); 
    xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');

subplot(132); hold on; title('Female Raw RS');
    plot([1 2], meanRaw(1:2), 'bo'); 
    errorbar([1 2], meanRaw(1:2), sraw(1:2), 'm', 'LineWidth', 2);
        for p=1:length(sumdat.fDuetAuto.rsRaw); plot(1.1, sumdat.fDuetAuto.rsRaw(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.fDuetHetero.rsRaw); plot(2.1, sumdat.fDuetHetero.rsRaw(p), 'k.', 'MarkerSize', 8); end
    plot([3 4], meanRaw(3:4), 'mo'); 
    errorbar([3 4], meanRaw(3:4), sraw(3:4), 'm', 'LineWidth', 1);
        for p=1:length(sumdat.fSoloAuto.rsRaw); plot(3.1, sumdat.fSoloAuto.rsRaw(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.fSoloHetero.rsRaw); plot(4.1, sumdat.fSoloHetero.rsRaw(p), 'k.', 'MarkerSize', 8); end
    % ylim([-10 50]); 
    xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');
    xticklabels({'DA','DH','SA','SH'})

subplot(133); hold on; title('Female Spikes/Second');
    plot([1 2], meanSPS(1:2), 'bo'); 
    errorbar([1 2], meanSPS(1:2), sps(1:2), 'm', 'LineWidth', 2);
        for p=1:length(sumdat.fDuetAuto.SPS); plot(1.1, sumdat.fDuetAuto.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.fDuetHetero.SPS); plot(2.1, sumdat.fDuetHetero.SPS(p), 'k.', 'MarkerSize', 8); end
    plot([3 4], meanSPS(3:4), 'mo'); 
    errorbar([3 4], meanSPS(3:4), sps(3:4), 'm', 'LineWidth', 1);
        for p=1:length(sumdat.fSoloAuto.SPS); plot(3.1, sumdat.fSoloAuto.SPS(p), 'k.', 'MarkerSize', 8); end
        for p=1:length(sumdat.fSoloHetero.SPS); plot(4.1, sumdat.fSoloHetero.SPS(p), 'k.', 'MarkerSize', 8); end
    % ylim([-1 55]); 
    xlim([0.5 4.5]); 
    
% subplot(132); hold on;
%     plot([1 2], mraw(1:2), 'b*'); 
%     errorbar([1 2], mraw(1:2), sraw(1:2), 'b' );% /sqrt(length(fChron)));
%     plot([3 4], mraw(3:4), 'm*'); hold on;
%     errorbar([3 4], mraw(3:4), sraw(3:4), 'm' );% /sqrt(length(mChron)));

%%%% SENSORY

% figure(2); clf; 
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
    
% Autogenous duet RS significant from 0?
    stts.m.dNAuto.mean = mean(sumdat.mDuetAuto.rsNorm);
    [stts.m.dNAuto.H, stts.m.dNAuto.P, stts.m.dNAuto.CI, stts.m.dNAuto.stats]  = ttest(sumdat.mDuetAuto.rsNorm);
    stts.m.dRAuto.mean = mean(sumdat.mDuetAuto.rsRaw);
    [stts.m.dRAuto.H, stts.m.dRAuto.P, stts.m.dRAuto.CI, stts.m.dRAuto.stats]  = ttest(sumdat.mDuetAuto.rsRaw);
    
    stts.f.dNAuto.mean = mean(sumdat.fDuetAuto.rsNorm);
    [stts.f.dNAuto.H, stts.f.dNAuto.P, stts.f.dNAuto.CI, stts.f.dNAuto.stats]  = ttest(sumdat.fDuetAuto.rsNorm);
    stts.f.dRAuto.mean = mean(sumdat.fDuetAuto.rsRaw);
    [stts.f.dRAuto.H, stts.f.dRAuto.P, stts.f.dRAuto.CI, stts.f.dRAuto.stats]  = ttest(sumdat.fDuetAuto.rsRaw);
    
% Autogenous Solo RS significant from 0?

    stts.m.sNAuto.mean = mean(sumdat.mSoloAuto.rsNorm);
    [stts.m.sNAuto.H, stts.m.sNAuto.P, stts.m.sNAuto.CI, stts.m.sNAuto.stats]  = ttest(sumdat.mSoloAuto.rsNorm);
    stts.m.sRAuto.mean = mean(sumdat.mSoloAuto.rsRaw);
    [stts.m.sRAuto.H, stts.m.sRAuto.P, stts.m.sRAuto.CI, stts.m.sRAuto.stats]  = ttest(sumdat.mSoloAuto.rsRaw);
    
    stts.f.sNAuto.mean = mean(sumdat.fSoloAuto.rsNorm);
    [stts.f.sNAuto.H, stts.f.sNAuto.P, stts.f.sNAuto.CI, stts.f.sNAuto.stats]  = ttest(sumdat.fSoloAuto.rsNorm);
    stts.f.sRAuto.mean = mean(sumdat.fSoloAuto.rsRaw);
    [stts.f.sRAuto.H, stts.f.sRAuto.P, stts.f.sRAuto.CI, stts.f.sRAuto.stats]  = ttest(sumdat.fSoloAuto.rsRaw);

% Difference between Autogenous Duet and Solo RS motor?

[stts.m.SvsDNAuto.H, stts.m.SvsDNAuto.P, stts.m.SvsDNAuto.CI, stts.m.SvsDNAuto.stats]  = ttest2(sumdat.mSoloAuto.rsNorm, sumdat.mDuetAuto.rsNorm);
[stts.m.SvsDRAuto.H, stts.m.SvsDRAuto.P, stts.m.SvsDRAuto.CI, stts.m.SvsDRAuto.stats]  = ttest2(sumdat.mSoloAuto.rsRaw, sumdat.mDuetAuto.rsRaw);
[stts.f.SvsDNAuto.H, stts.f.SvsDNAuto.P, stts.f.SvsDNAuto.CI, stts.f.SvsDNAuto.stats]  = ttest2(sumdat.fSoloAuto.rsNorm, sumdat.fDuetAuto.rsNorm);
[stts.f.SvsDRAuto.H, stts.f.SvsDRAuto.P, stts.f.SvsDRAuto.CI, stts.f.SvsDRAuto.stats]  = ttest2(sumdat.fSoloAuto.rsRaw, sumdat.fDuetAuto.rsRaw);

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

[stts.m.SvsDNHetero.H, stts.m.SvsDNHetero.P, stts.m.SvsDNHetero.CI, stts.m.SvsDNHetero.stats]  = ttest2(sumdat.mSoloHetero.rsNorm, sumdat.mDuetHetero.rsNorm);
[stts.m.SvsDRHetero.H, stts.m.SvsDRHetero.P, stts.m.SvsDRHetero.CI, stts.m.SvsDRHetero.stats]  = ttest2(sumdat.mSoloHetero.rsRaw, sumdat.mDuetHetero.rsRaw);
[stts.f.SvsDNHetero.H, stts.f.SvsDNHetero.P, stts.f.SvsDNHetero.CI, stts.f.SvsDNHetero.stats]  = ttest2(sumdat.fSoloHetero.rsNorm, sumdat.fDuetHetero.rsNorm);
[stts.f.SvsDRHetero.H, stts.f.SvsDRHetero.P, stts.f.SvsDRHetero.CI, stts.f.SvsDRHetero.stats]  = ttest2(sumdat.fSoloHetero.rsRaw, sumdat.fDuetHetero.rsRaw);


% Difference between Auto and Hetero in the duet?

[stts.m.AvsHNDuet.H, stts.m.AvsHNDuet.P, stts.m.AvsHNDuet.CI, stts.m.AvsHNDuet.stats]  = ttest2(sumdat.mDuetHetero.rsNorm, sumdat.mDuetAuto.rsNorm);
[stts.m.AvsHRDuet.H, stts.m.AvsHRDuet.P, stts.m.AvsHRDuet.CI, stts.m.AvsHRDuet.stats]  = ttest2(sumdat.mDuetHetero.rsRaw, sumdat.mDuetAuto.rsRaw);
[stts.f.AvsHNDuet.H, stts.f.AvsHNDuet.P, stts.f.AvsHNDuet.CI, stts.f.AvsHNDuet.stats]  = ttest2(sumdat.fDuetHetero.rsNorm, sumdat.fDuetAuto.rsNorm);
[stts.f.AvsHRDuet.H, stts.f.AvsHRDuet.P, stts.f.AvsHRDuet.CI, stts.f.AvsHRDuet.stats]  = ttest2(sumdat.fDuetHetero.rsRaw, sumdat.fDuetAuto.rsRaw);

% Difference between Auto and Hetero in the Solo?

[stts.m.AvsHNSolo.H, stts.m.AvsHNSolo.P, stts.m.AvsHNSolo.CI, stts.m.AvsHNSolo.stats]  = ttest2(sumdat.mSoloHetero.rsNorm, sumdat.mSoloAuto.rsNorm);
[stts.m.AvsHRSolo.H, stts.m.AvsHRSolo.P, stts.m.AvsHRSolo.CI, stts.m.AvsHRSolo.stats]  = ttest2(sumdat.mSoloHetero.rsRaw, sumdat.mSoloAuto.rsRaw);
[stts.f.AvsHNSolo.H, stts.f.AvsHNSolo.P, stts.f.AvsHNSolo.CI, stts.f.AvsHNSolo.stats]  = ttest2(sumdat.fSoloHetero.rsNorm, sumdat.fSoloAuto.rsNorm);
[stts.f.AvsHRSolo.H, stts.f.AvsHRSolo.P, stts.f.AvsHRSolo.CI, stts.f.AvsHRSolo.stats]  = ttest2(sumdat.fSoloHetero.rsRaw, sumdat.fSoloAuto.rsRaw);


    fprintf('Male Auto Duet Raw RS different from zero? p = %1.2e \n', stts.m.dRAuto.P);
    fprintf('Male Hetero Duet Raw RS different from zero? p = %1.2e \n', stts.m.dRHetero.P);
    fprintf('Male Auto Solo Raw RS different from zero? p = %1.2e \n', stts.m.sRAuto.P);
    fprintf('Male Hetero Solo Raw RS different from zero? p = %1.2e \n', stts.m.sRHetero.P);
    
    fprintf(' \n');
    
    fprintf('Female Auto Duet Raw RS different from zero? p = %1.2e \n', stts.f.dRAuto.P);
    fprintf('Female Hetero Duet Raw RS different from zero? p = %1.2e \n', stts.f.dRHetero.P);
    fprintf('Female Auto Solo Raw RS different from zero? p = %1.2e \n', stts.f.sRAuto.P);
    fprintf('Female Hetero Solo Raw RS different from zero? p = %1.2e \n', stts.f.sRHetero.P);

    fprintf(' \n');

    fprintf('Male Auto Duet vs Hetero? p = %1.2e \n', stts.m.AvsHRDuet.P);
    fprintf('Female Auto Duet vs Hetero? p = %1.2e \n', stts.f.AvsHRDuet.P);

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%% Response Strength nested function

function qwe = rs(struc, syllabl, spontan, padme)
 
   % Get spontaneous rate 
        sponSpikeCount = 0;  
    
        for i=1:length(struc.Aspikes) % For each repetition of the stimulus
            sponSpikeCount = sponSpikeCount + length(find(struc.Aspikes{i} > spontan(1) & struc.Aspikes{i} < spontan(2)));
        end

        sponSPS = sponSpikeCount / (spontan(2) - spontan(1)); % This is spikes per second
        sponSPS = sponSPS/length(struc.Aspikes); % Divide by number of repetitions
   
    % Loop for each syllable
    
    for j = 1:length(syllabl)    
    
        % Get the spikes for that syllable
        stimSpikeCount = 0; 
    
        for i=1:length(struc.Aspikes) % For each repetition of the stimulus
            % (padme shifts the window in seconds, negative earlier, positive later)
%             if j < length(syllabl)
%                 sylender = struc.syl(syllabl(j+1)).tim(1)+padme; % Beginning of next syllable              
%             end
%             if j == length(syllabl)
%                 sylender = struc.syl(syllabl(j)).tim(2)+padme; % End of the current syllable
%             end
            
            stimSpikeCount = stimSpikeCount + length(find(struc.Aspikes{i} >= struc.syl(syllabl(j)).tim(1)+padme & struc.Aspikes{i} < struc.syl(syllabl(j)).tim(2)+padme));
%            stimSpikeCount = stimSpikeCount + length(find(struc.Aspikes{i} >= struc.syl(syllabl(j)).tim(1)+padme & struc.Aspikes{i} < sylender));

        end
        
        stimSPS = stimSpikeCount / (struc.syl(syllabl(j)).tim(2) - struc.syl(syllabl(j)).tim(1)); % This is spikes per second
        stimSPS = stimSPS/length(struc.Aspikes); % Divide by number of repetitions
        
        
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




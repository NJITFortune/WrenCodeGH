function [rsData, dsiData, sts] = wRS(in, padding)
% Usage: function [RSchron, RSacute, dsiAcute, dsiChron, sts] = wRS(in, padding)
% Calculates response strength to solo and duet syllables.
%
sts = 0;

% load ChronicCompleat2017p.mat

pad = 0; 

if nargin == 2; pad = padding; end

% The boundaries of each syllable are used for this analysis. This is
% inherently problematic as we expect pre-motor activity in awake animals to occur PRIOR to
% the sound and auditory activity in urethane-anesthetized animals to occur AFTER the sound.
% We are comfortable with this approach because it is a rather unbiased. Change
% the value of padding in seconds (e.g. 0.005 or -0.003) to look at the effects.
    
%% List of duets with syllable indices and locations for spontaneous activity

% m17: March 2017, male index 1, female index 2

    msolosyls{1} = [1 2 3 4]; 
    mduetsyls{1} = [5 7 9 11 13];
    fsolosyls{1} = []; 
    fduetsyls{1} = [6 8 10 12 14];
    spon(:,1) = [-5.5, -0.5]; % This is a mess. 
    
% j160806: January 2016, 08:06am, male index 3, female index 4

    msolosyls{2} = []; 
    mduetsyls{2} = [5 8 10 13];
    fsolosyls{2} = [1 2 3 4]; 
    fduetsyls{2} = [6 7 9 11 12 14];    
    spon(:,2) = [-5.0, 0.0]; 
    
% j160807: January 2016, 08:07am, male index 5, female index 6
    
    msolosyls{3} = []; 
    mduetsyls{3} = [3 6 8 11 13 16 18 21 23];
    fsolosyls{3} = [1 2]; 
    fduetsyls{3} = [4 5 7 9 10 13 14 15 17 19 20 22 24];
    spon(:,3) = [-5.0, 0.0];
    
% j160815: January 2016, 08:15am, male index 7, female index 8
    
    msolosyls{4} = []; 
    mduetsyls{4} = [3 5 7];
    fsolosyls{4} = [1 2]; 
    fduetsyls{4} = [4 6 8];
    spon(:,4) = [-2.5, -1.0];

% j161009: January 2016, 10:09am, male index 9, female index 10
    
    msolosyls{5} = []; 
    mduetsyls{5} = [3 5 7 9 11 14 16 19 21];
    fsolosyls{5} = [1 2]; 
    fduetsyls{5} = [4 6 8 10 12 13 15 17 18 20];
    spon(:,5) = [-5.0, 0.0];

% j161022: January 2016, 10:22am, male index 11, female index 12
    
    msolosyls{6} = []; mduetsyls{6} = [4 6 8 10 12 14];
    fsolosyls{6} = [1 2 3]; fduetsyls{6} = [5 7 9 11 13];
    spon(:,6) = [-3.5, 0.0];

% j17060848: 06 January 2017, 08:48am, male index 13, female index 14
    
    msolosyls{7} = [1 2 3 4 5]; mduetsyls{7} = [8 10 13 15 18 20];
    fsolosyls{7} = []; fduetsyls{7} = [6 7 9 11 12 14 17 19];
    spon(:,7) = [-3.0, 0.0];

% j17081733: 08 January 2017, 5:33pm, male index 15, female index 16

    msolosyls{8} = []; mduetsyls{8} = [3 8];
    fsolosyls{8} = 1; fduetsyls{8} = [2 7];
    spon(:,8) = [-5, 0];
    
    
%% Loop to calculate RS values for each pair of wrens  
% 
        tmp.sylnum = [];
        tmp.rsNorm = [];
        tmp.rsRaw = [];
        tmp.sponSPS = [];
        tmp.stimSPS = [];
        tmp.spontim = [];
        tmp.syltim = [];
        tmp.pad = [];

fHsoloC = tmp; mAsoloC = tmp;
mHsoloC = tmp; fAsoloC = tmp;
fHduetC = tmp; mAduetC = tmp;
mHduetC = tmp; fAduetC = tmp;

fHsoloA = tmp; mAsoloA = tmp;
mHsoloA = tmp; fAsoloA = tmp;
fHduetA = tmp; mAduetA = tmp;
mHduetA = tmp; fAduetA = tmp;

for curpair = 1:length(spon)

    % Solo syllables MALE
    if ~isempty(msolosyls{curpair}) % Male sang solo syllables        
        for kk = 1:length(msolosyls{curpair}) % For each syllable
            [fHsoloC(end+1), fHsoloA(end+1)] = rs(in(curpair*2), msolosyls{curpair}(kk), spon(:,curpair), pad);
            [mAsoloC(end+1), mAsoloA(end+1)] = rs(in((curpair*2)-1), msolosyls{curpair}(kk), spon(:,curpair), pad);            
        end
    end
    
    % Solo syllables FEMALE
    if ~isempty(fsolosyls{curpair}) % Female sang solo syllables
        for kk = 1:length(fsolosyls{curpair})
            [mHsoloC(end+1), mHsoloA(end+1)] = rs(in((curpair*2)-1), fsolosyls{curpair}(kk), spon(:,curpair), pad);
            [fAsoloC(end+1), fAsoloA(end+1)] = rs(in(curpair*2), fsolosyls{curpair}(kk), spon(:,curpair), pad);
        end
    end
    
    % Duet syllables MALE
    if ~isempty(mduetsyls{curpair}) % Male sang solo syllables        
        for kk = 1:length(mduetsyls{curpair}) % For each syllable
            [fHduetC(end+1), fHduetA(end+1)] = rs(in(curpair*2), mduetsyls{curpair}(kk), spon(:,curpair), pad);
            [mAduetC(end+1), mAduetA(end+1)] = rs(in((curpair*2)-1), mduetsyls{curpair}(kk), spon(:,curpair), pad);            
        end
    end
    
    % Duet syllables FEMALE
    if ~isempty(fduetsyls{curpair}) % Female sang solo syllables
        for kk = 1:length(fsolosyls{curpair})
            [mHduetC(end+1), mHduetA(end+1)] = rs(in((curpair*2)-1), fduetsyls{curpair}(kk), spon(:,curpair), pad);
            [fAduetC(end+1), fAduetA(end+1)] = rs(in(curpair*2), fduetsyls{curpair}(kk), spon(:,curpair), pad);
        end
    end 
end

% Assemble into the output structure (so lazy!)

rsData.fc.Hsolo = fHsoloC(2:end); 
rsData.mc.Asolo = mAsoloC(2:end);
rsData.mc.Hsolo = mHsoloC(2:end); 
rsData.fc.Asolo = fAsoloC(2:end);
rsData.fc.Hduet = fHduetC(2:end); 
rsData.mc.Aduet = mAduetC(2:end);
rsData.mc.Hduet = mHduetC(2:end); 
rsData.fc.Aduet = fAduetC(2:end);

rsData.ma.Aduet = mAduetA(2:end);
rsData.ma.Hduet = mHduetA(2:end); 
rsData.ma.Asolo = mAsoloA(2:end);
rsData.ma.Hsolo = mHsoloA(2:end); 
rsData.fa.Hsolo = fHsoloA(2:end-1); 
rsData.fa.Asolo = fAsoloA(2:end-1);
rsData.fa.Hduet = fHduetA(2:end-1); 
rsData.fa.Aduet = fAduetA(2:end-1);

clear fAduetA mHduetA mAduetA fHduetA fAsoloA mHsoloA mAsoloA fHsoloA
clear fAduetC mHduetC mAduetC fHduetC fAsoloC mHsoloC mAsoloC fHsoloC

figure; clf;

axx(1) = subplot(121); hold on;
        plot(0.9, mean([rsData.mc.Aduet.rsNorm]), 'b*');
            errorbar(0.9, mean([rsData.mc.Aduet.rsNorm]), std([rsData.mc.Aduet.rsNorm]), 'b' );% /sqrt(length(mChron)));
        plot(1.1, mean([rsData.mc.Hduet.rsNorm]), 'b*');
            errorbar(1.1, mean([rsData.mc.Hduet.rsNorm]), std([rsData.mc.Hduet.rsNorm]), 'b' );% /sqrt(length(mChron)));
        text(0.9, 5, 'A'); text(1.1, 5, 'H');

        plot(1.9, mean([rsData.mc.Asolo.rsNorm]), 'b*');
            errorbar(1.9, mean([rsData.mc.Asolo.rsNorm]), std([rsData.mc.Asolo.rsNorm]), 'b' );% /sqrt(length(mChron)));
        plot(2.1, mean([rsData.mc.Hsolo.rsNorm]), 'b*');
            errorbar(2.1, mean([rsData.mc.Hsolo.rsNorm]), std([rsData.mc.Hsolo.rsNorm]), 'b' );% /sqrt(length(mChron)));
        text(1.9, 5, 'A'); text(2.1, 5, 'H');

        plot(2.9, mean([rsData.fc.Aduet.rsNorm]), 'm*');
            errorbar(2.9, mean([rsData.fc.Aduet.rsNorm]), std([rsData.fc.Aduet.rsNorm]), 'm' );% /sqrt(length(mChron)));
        plot(3.1, mean([rsData.fc.Hduet.rsNorm]), 'm*');
            errorbar(3.1, mean([rsData.fc.Hduet.rsNorm]), std([rsData.fc.Hduet.rsNorm]), 'm' );% /sqrt(length(mChron)));
        text(2.9, 5, 'A'); text(3.1, 5, 'H');

        plot(3.9, mean([rsData.fc.Asolo.rsNorm]), 'm*');
            errorbar(3.9, mean([rsData.fc.Asolo.rsNorm]), std([rsData.fc.Asolo.rsNorm]), 'm' );% /sqrt(length(mChron)));
        plot(4.1, mean([rsData.fc.Hsolo.rsNorm]), 'm*');
            errorbar(4.1, mean([rsData.fc.Hsolo.rsNorm]), std([rsData.fc.Hsolo.rsNorm]), 'm' );% /sqrt(length(mChron)));
        text(3.9, 5, 'A'); text(4.1, 5, 'H');

        text(2.5, -1.2, 'Chronic', 'HorizontalAlignment', 'center');     
        plot([0.5 4.5], [0 0], 'k-');
        ylabel('Normalized Response Strength');

axx(2) = subplot(122); hold on;
        plot(0.9, mean([rsData.ma.Aduet.rsNorm]), 'b*');
            errorbar(0.9, mean([rsData.ma.Aduet.rsNorm]), std([rsData.ma.Aduet.rsNorm]), 'b' );% /sqrt(length(mChron)));
        plot(1.1, mean([rsData.ma.Hduet.rsNorm]), 'b*');
            errorbar(1.1, mean([rsData.ma.Hduet.rsNorm]), std([rsData.ma.Hduet.rsNorm]), 'b' );% /sqrt(length(mChron)));
        text(0.9, 5, 'A'); text(1.1, 5, 'H');

        plot(1.9, mean([rsData.ma.Asolo.rsNorm]), 'b*');
            errorbar(1.9, mean([rsData.ma.Asolo.rsNorm]), std([rsData.ma.Asolo.rsNorm]), 'b' );% /sqrt(length(mChron)));
        plot(2.1, mean([rsData.ma.Hsolo.rsNorm]), 'b*');
            errorbar(2.1, mean([rsData.ma.Hsolo.rsNorm]), std([rsData.ma.Hsolo.rsNorm]), 'b' );% /sqrt(length(mChron)));
        text(1.9, 5, 'A'); text(2.1, 5, 'H');
            
        plot(2.9, mean([rsData.fa.Aduet.rsNorm]), 'm*');
            errorbar(2.9, mean([rsData.fa.Aduet.rsNorm]), std([rsData.fa.Aduet.rsNorm]), 'm' );% /sqrt(length(mChron)));
        plot(3.1, mean([rsData.fa.Hduet.rsNorm]), 'm*');
            errorbar(3.1, mean([rsData.fa.Hduet.rsNorm]), std([rsData.fa.Hduet.rsNorm]), 'm' );% /sqrt(length(mChron)));
        text(2.9, 5, 'A'); text(3.1, 5, 'H');
            
        plot(3.9, mean([rsData.fa.Asolo.rsNorm]), 'm*');
            errorbar(3.9, mean([rsData.fa.Asolo.rsNorm]), std([rsData.fa.Asolo.rsNorm]), 'm' );% /sqrt(length(mChron)));
        plot(4.1, mean([rsData.fa.Hsolo.rsNorm]), 'm*');
            errorbar(4.1, mean([rsData.fa.Hsolo.rsNorm]), std([rsData.fa.Hsolo.rsNorm]), 'm' );% /sqrt(length(mChron)));
        text(3.9, 5, 'A'); text(4.1, 5, 'H');

        text(2.5, -1.2, 'Acute', 'HorizontalAlignment', 'center');
        plot([0.5 4.5], [0 0], 'k-');

        linkaxes(axx, 'xy');
        xlim([0.5 4.5]);         




%% Loop to calculate DSI values (sanity check for the wDSI more restrictive analysis)

    dsiData.MC.solo = dsi(rsData.mc.Hsolo, rsData.mc.Asolo);
    dsiData.MC.duet = dsi(rsData.mc.Hduet, rsData.mc.Aduet);
    dsiData.MA.solo = dsi(rsData.ma.Hsolo, rsData.ma.Asolo);
    dsiData.MA.duet = dsi(rsData.ma.Hduet, rsData.ma.Aduet);
    dsiData.FC.solo = dsi(rsData.fc.Hsolo, rsData.fc.Asolo);
    dsiData.FC.duet = dsi(rsData.fc.Hduet, rsData.fc.Aduet);
    dsiData.FA.solo = dsi(rsData.fa.Hsolo, rsData.fa.Asolo);
    dsiData.FA.duet = dsi(rsData.fa.Hduet, rsData.fa.Aduet);
    
    
% Plot dsi figure
    figure; clf; 
        ax(1) = subplot(121); hold on;
        ylabel('Duet Selectivity Index');
        plot(1, dsiData.MC.duet, 'b*');
        plot(2, dsiData.MC.solo, 'b*');
        plot(3, dsiData.FC.duet, 'm*');
        plot(4, dsiData.FC.solo, 'm*');
        plot([0.5 4.5], [0 0], 'k-', 'LineWidth', 0.1);
        text(2.5, -0.8, 'Chronic', 'HorizontalAlignment', 'center');

        ax(2) = subplot(122); hold on;
        plot(1, dsiData.MA.duet, 'b*');
        plot(2, dsiData.MA.solo, 'b*');
        plot(3, dsiData.FA.duet, 'm*');
        plot(4, dsiData.FA.solo, 'm*');
        plot([0.5 4.5], [0 0], 'k-', 'LineWidth', 0.1);
        text(2.5, -0.8, 'Acute', 'HorizontalAlignment', 'center');

        linkaxes(ax, 'xy');
        xlim([0.5 4.5]); ylim([-1 1]);

%% Plot 

%%%% MOTOR

% For the Normalized data. Not doing STD because of assymetry of
% distribution.
% m(1) = mean(sumdat.mduetAutogenous.rsNorm); s(1) = std(sumdat.mduetAutogenous.rsNorm);
% m(2) = mean(sumdat.mSolo.rsNorm); s(2) = std(sumdat.mSolo.rsNorm);
% m(3) = mean(sumdat.fduetAutogenous.rsNorm); s(3) = std(sumdat.fduetAutogenous.rsNorm);
% m(4) = mean(sumdat.fSolo.rsNorm); s(4) = std(sumdat.fSolo.rsNorm);
% 
% % For raw data
% mraw(1) = mean(sumdat.mduetAutogenous.rsNorm); sraw(1) = std(sumdat.mduetAutogenous.rsRaw);
% mraw(2) = mean(sumdat.mSolo.rsNorm); sraw(3) = std(sumdat.mSolo.rsRaw);
% mraw(3) = mean(sumdat.fduetAutogenous.rsNorm); sraw(2) = std(sumdat.fduetAutogenous.rsRaw);
% mraw(4) = mean(sumdat.fSolo.rsNorm); sraw(4) = std(sumdat.fSolo.rsRaw);
% 
% figure(1); clf; 
% 
% subplot(121); hold on;
%     plot([1 2], m(1:2), 'bo'); 
%     errorbar([1 2], m(1:2), s(1:2), 'b' );
%     plot([3 4], m(3:4), 'mo'); 
%     errorbar([3 4], m(3:4), s(3:4), 'm' );
%     ylim([-5 20]);
% 
% xlim([0.5 4.5]); 
% 
% % subplot(132); hold on;
% %     plot([1 2], mraw(1:2), 'b*'); 
% %     errorbar([1 2], mraw(1:2), sraw(1:2), 'b' );% /sqrt(length(fChron)));
% %     plot([3 4], mraw(3:4), 'm*'); hold on;
% %     errorbar([3 4], mraw(3:4), sraw(3:4), 'm' );% /sqrt(length(mChron)));
% 
% %%%% SENSORY
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
%     ylim([-5 20]);
% 
% %% Compute stats     
%     
% % Autogenous duet RS significant from 0?
% stats.m.dNAuto.mean = mean(sumdat.mduetAutogenous.rsNorm);
% [stats.m.dNAuto.H, stats.m.dNAuto.P, stats.m.dNAuto.CI, stats.m.dNAuto.stats]  = ttest(sumdat.mduetAutogenous.rsNorm);
% stats.m.dRAuto.mean = mean(sumdat.mduetAutogenous.rsRaw);
% [stats.m.dRAuto.H, stats.m.dRAuto.P, stats.m.dRAuto.CI, stats.m.dRAuto.stats]  = ttest(sumdat.mduetAutogenous.rsRaw);
% stats.f.dNAuto.mean = mean(sumdat.fduetAutogenous.rsNorm);
% [stats.f.dNAuto.H, stats.f.dNAuto.P, stats.f.dNAuto.CI, stats.f.dNAuto.stats]  = ttest(sumdat.fduetAutogenous.rsNorm);
% stats.f.dRAuto.mean = mean(sumdat.fduetAutogenous.rsRaw);
% [stats.f.dRAuto.H, stats.f.dRAuto.P, stats.f.dRAuto.CI, stats.f.dRAuto.stats]  = ttest(sumdat.fduetAutogenous.rsRaw);
%     
% % Autogenous Solo RS significant from 0?
% 
% stats.m.sNAuto.mean = mean(sumdat.mSolo.rsNorm);
% [stats.m.sNAuto.H, stats.m.sNAuto.P, stats.m.sNAuto.CI, stats.m.sNAuto.stats]  = ttest(sumdat.mSolo.rsNorm);
% stats.m.sRAuto.mean = mean(sumdat.mSolo.rsRaw);
% [stats.m.sRAuto.H, stats.m.sRAuto.P, stats.m.sRAuto.CI, stats.m.sRAuto.stats]  = ttest(sumdat.mSolo.rsRaw);
% stats.f.sNAuto.mean = mean(sumdat.fSolo.rsNorm);
% [stats.f.sNAuto.H, stats.f.sNAuto.P, stats.f.sNAuto.CI, stats.f.sNAuto.stats]  = ttest(sumdat.fSolo.rsNorm);
% stats.m.sRAuto.mean = mean(sumdat.fSolo.rsRaw);
% [stats.f.sRAuto.H, stats.f.sRAuto.P, stats.f.sRAuto.CI, stats.f.sRAuto.stats]  = ttest(sumdat.fSolo.rsRaw);
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
% stats.m.dNHetero.mean = mean(sumdat.mduetHeterogenous.rsNorm);
% [stats.m.dNHetero.H, stats.m.dNHetero.P, stats.m.dNHetero.CI, stats.m.dNHetero.stats]  = ttest(sumdat.mduetHeterogenous.rsNorm);
% stats.m.dRHetero.mean = mean(sumdat.mduetHeterogenous.rsRaw);
% [stats.m.dRHetero.H, stats.m.dRHetero.P, stats.m.dRHetero.CI, stats.m.dRHetero.stats]  = ttest(sumdat.mduetHeterogenous.rsRaw);
% stats.f.dNHetero.mean = mean(sumdat.fduetHeterogenous.rsNorm);
% [stats.f.dNHetero.H, stats.f.dNHetero.P, stats.f.dNHetero.CI, stats.f.dNHetero.stats]  = ttest(sumdat.fduetHeterogenous.rsNorm);
% stats.f.dRHetero.mean = mean(sumdat.fduetHeterogenous.rsRaw);
% [stats.f.dRHetero.H, stats.f.dRHetero.P, stats.f.dRHetero.CI, stats.f.dRHetero.stats]  = ttest(sumdat.fduetHeterogenous.rsRaw);
%     
% % Heterogenous Solo RS significant from 0?
% 
% stats.m.sNHetero.mean = mean(sumdat.mAud.rsNorm);
% [stats.m.sNHetero.H, stats.m.sNHetero.P, stats.m.sNHetero.CI, stats.m.sNHetero.stats]  = ttest(sumdat.mAud.rsNorm);
% stats.m.sRHetero.mean = mean(sumdat.mAud.rsRaw);
% [stats.m.sRHetero.H, stats.m.sRHetero.P, stats.m.sRHetero.CI, stats.m.sRHetero.stats]  = ttest(sumdat.mAud.rsRaw);
% stats.f.sNHetero.mean = mean(sumdat.fAud.rsNorm);
% [stats.f.sNHetero.H, stats.f.sNHetero.P, stats.f.sNHetero.CI, stats.f.sNHetero.stats]  = ttest(sumdat.fAud.rsNorm);
% stats.f.sRHetero.mean = mean(sumdat.fAud.rsRaw);
% [stats.f.sRHetero.H, stats.f.sRHetero.P, stats.f.sRHetero.CI, stats.f.sRHetero.stats]  = ttest(sumdat.fAud.rsRaw);
% 
% % Difference between Heterogenous Duet and Solo RS motor?
% 
% [stats.m.SvsDNHetero.H, stats.m.SvsDNHetero.P, stats.m.SvsDNHetero.CI, stats.m.SvsDNHetero.stats]  = ttest2(sumdat.mAud.rsNorm, sumdat.mduetHeterogenous.rsNorm);
% [stats.m.SvsDRHetero.H, stats.m.SvsDRHetero.P, stats.m.SvsDRHetero.CI, stats.m.SvsDRHetero.stats]  = ttest2(sumdat.mAud.rsRaw, sumdat.mduetHeterogenous.rsRaw);
% [stats.f.SvsDNHetero.H, stats.f.SvsDNHetero.P, stats.f.SvsDNHetero.CI, stats.f.SvsDNHetero.stats]  = ttest2(sumdat.fAud.rsNorm, sumdat.fduetHeterogenous.rsNorm);
% [stats.f.SvsDRHetero.H, stats.f.SvsDRHetero.P, stats.f.SvsDRHetero.CI, stats.f.SvsDRHetero.stats]  = ttest2(sumdat.fAud.rsRaw, sumdat.fduetHeterogenous.rsRaw);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Response Strength nested function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function [cRs, aRs] = rs(struc, syllabl, spontan, padme)

   % Get spontaneous rate 
    CsponSpikeCount = 0;  
    AsponSpikeCount = 0;  
    CstimSpikeCount = 0; 
    AstimSpikeCount = 0; 
    
        for i=1:4 % 4 electrodes in a tetrode always
            CsponSpikeCount = CsponSpikeCount + length(find(struc.Cspikes{i} > spontan(1) & struc.Cspikes{i} < spontan(2)));
            CstimSpikeCount = CstimSpikeCount + length(find(struc.Cspikes{i} > struc.syl(syllabl).tim(1)-padme & struc.Cspikes{i} < struc.syl(syllabl).tim(2)-padme));
        end
            Csponrate = CsponSpikeCount / (spontan(2) - spontan(1));
            Cstimrate = CstimSpikeCount / (struc.syl(syllabl).tim(2) - struc.syl(syllabl).tim(1));        

            
        for i=1:length(struc.Aspikes) % Number of replays under urethane
            AsponSpikeCount = AsponSpikeCount + length(find(struc.Aspikes{i} > spontan(1) & struc.Aspikes{i} < spontan(2)));
            AstimSpikeCount = AstimSpikeCount + length(find(struc.Aspikes{i} > struc.syl(syllabl).tim(1)-padme & struc.Aspikes{i} < struc.syl(syllabl).tim(2)-padme));
        end
            Asponrate = AsponSpikeCount / (spontan(2) - spontan(1));
            Astimrate = AstimSpikeCount / (struc.syl(syllabl).tim(2) - struc.syl(syllabl).tim(1));        
        
        cRs.sylnum = syllabl;
        if Csponrate > 0.01
            cRs.rsNorm = (Cstimrate - Csponrate) / Csponrate;
        else
            cRs.rsNorm = [];
        end
        cRs.rsRaw = Cstimrate - Csponrate;
        cRs.sponSPS = Csponrate;
        cRs.stimSPS = Cstimrate;
        cRs.spontim = spontan;
        cRs.syltim = [struc.syl(syllabl).tim(1), struc.syl(syllabl).tim(2)];
        cRs.pad = padme;

        aRs.sylnum = syllabl;
        if Asponrate > 0.01
            aRs.rsNorm = (Astimrate - Asponrate) / Asponrate;
        else
            aRs.rsNorm = [];
        end
        aRs.rsRaw = Astimrate - Asponrate;
        aRs.sponSPS = Asponrate;
        aRs.stimSPS = Astimrate;
        aRs.spontim = spontan;
        aRs.syltim = [struc.syl(syllabl).tim(1), struc.syl(syllabl).tim(2)];
        aRs.pad = padme;
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Duet Selectivity Index nested function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

function dsiValue = dsi(Het, Aut)
        
    for j=2:length(Het); H(j) = Het(j).stimSPS; end
    for j=2:length(Aut); A(j) = Aut(j).stimSPS; end

    dsiValue = (mean(A) - mean(H)) / (mean(A) + mean(H));            

end



end

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
    

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
sumdat.fduetAutogenous.rsNorm = []; sumdat.fduetAutogenous.rsRaw = [];
sumdat.mduetAutogenous.rsNorm = []; sumdat.mduetAutogenous.rsRaw = [];
sumdat.fduetHeterogenous.rsNorm = []; sumdat.fduetHeterogenous.rsRaw = [];
sumdat.mduetHeterogenous.rsNorm = []; sumdat.mduetHeterogenous.rsRaw = [];
sumdat.fSolo.rsNorm = []; sumdat.fSolo.rsRaw = [];
sumdat.mSolo.rsNorm = []; sumdat.mSolo.rsRaw = [];
sumdat.fAud.rsNorm = []; sumdat.fAud.rsRaw = [];
sumdat.mAud.rsNorm = []; sumdat.mAud.rsRaw = [];

    
%% List of Chronic singing data with syllable indices and locations for spontaneous activity

[msolosyls, mduetsyls, fsolosyls, fduetsyls, spon] = wData;
% % 1-2: m17
% 
%     msolosyls{1} = [1 2 3 4]; 
%     mduetsyls{1} = [5 7 9 11 13];
%     fsolosyls{1} = []; 
%     fduetsyls{1} = [6 8 10 12 14];
%     spon(:,1) = [-5.5, -0.5]; % This is a mess. 
%     
% % 3-4: j160806
% 
%     msolosyls{2} = [2]; 
%     mduetsyls{2} = [4 6 8 10 12];
%     fsolosyls{2} = [1]; 
%     fduetsyls{2} = [3 5 7 9 11 13];    
%     spon(:,2) = [-5.0, 0.0]; 
%     
% % 5-6: j160807
%     
%     msolosyls{3} = []; 
%     mduetsyls{3} = [3 6 8 11 13 16 18 21 23];
%     fsolosyls{3} = [1 2]; 
%     fduetsyls{3} = [4 5 7 9 10 13 14 15 17 19 20 22 24];
%     spon(:,3) = [-5.0, 0.0];
%     
% % 7-8: j160815
%     
%     msolosyls{4} = []; 
%     mduetsyls{4} = [3 5 7];
%     fsolosyls{4} = [1 2]; 
%     fduetsyls{4} = [4 6 8];
%     spon(:,4) = [-2.5, -1.0];
% 
% % 9-10: j161009
%     
%     msolosyls{5} = []; 
%     mduetsyls{5} = [3 5 7 9 11 14 16 19 21];
%     fsolosyls{5} = [1 2]; 
%     fduetsyls{5} = [4 6 8 10 12 13 15 17 18 20];
%     spon(:,5) = [-5.0, 0.0];
% 
% % 11-12: j161022
%     
%     msolosyls{6} = []; 
%     mduetsyls{6} = [4 6 8 10 12 14];
%     fsolosyls{6} = [1 2 3]; 
%     fduetsyls{6} = [5 7 9 11 13];
%     spon(:,6) = [-3.5, 0.0];
% 
% % 13-14: j17060848
%     
%     msolosyls{7} = [1 2 3 4 5]; 
%     mduetsyls{7} = [8 10 13 15 18 20];
%     fsolosyls{7} = []; 
%     fduetsyls{7} = [6 7 9 11 12 14 17 19];
%     spon(:,7) = [-3.0, 0.0];
% 
% % 15-16: j170081733 
% 
% %     msolosyls{8} = []; 
% %     mduetsyls{8} = [3 8];
% %     fsolosyls{8} = 1; 
% %     fduetsyls{8} = [2 4 5 7 9 10];
% %     spon(:,8) = [-5, 0];

%% Loop to calculate RS values for each pair of wrens   
    
for curpair = 1:length(spon)
    curpair
    % Solo syllables MALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(msolosyls{curpair}) % Male sang solo syllables
        out(curpair).fAud = rs(in(curpair*2), msolosyls{curpair}, spon(:,curpair), pad);
        out(curpair).mSolo = rs(in((curpair*2)-1), msolosyls{curpair}, spon(:,curpair), pad);
        for kk = 1:length(msolosyls{curpair})
            sumdat.fAud.rsNorm(end+1) = out(curpair).fAud(kk).rsNorm;
            sumdat.fAud.rsRaw(end+1) = out(curpair).fAud(kk).rsRaw;
            sumdat.mSolo.rsNorm(end+1) = out(curpair).mSolo(kk).rsNorm;
            sumdat.mSolo.rsRaw(end+1) = out(curpair).mSolo(kk).rsRaw;
        end
    end
    
    % Solo syllables FEMALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(fsolosyls{curpair}) % Female sang solo syllables
        out(curpair).mAud = rs(in((curpair*2)-1), fsolosyls{curpair}, spon(:,curpair), pad);
        out(curpair).fSolo = rs(in(curpair*2), fsolosyls{curpair}, spon(:,curpair), pad);
        
        for kk = 1:length(fsolosyls{curpair})
            sumdat.mAud.rsNorm(end+1) = out(curpair).mAud(kk).rsNorm;
            sumdat.mAud.rsRaw(end+1) = out(curpair).mAud(kk).rsRaw;
            sumdat.fSolo.rsNorm(end+1) = out(curpair).fSolo(kk).rsNorm;
            sumdat.fSolo.rsRaw(end+1) = out(curpair).fSolo(kk).rsRaw;
        end
    end
    
    % Motor activity (Autogenous) during duet (bird's own syllables) %%%%%%
    out(curpair).mduetA = rs(in((curpair*2)-1), mduetsyls{curpair}, spon(:,curpair), pad);
    out(curpair).fduetA = rs(in(curpair*2), fduetsyls{curpair}, spon(:,curpair), pad);
    
    for kk = 1:length(out(curpair).mduetA)
        sumdat.mduetAutogenous.rsNorm(end+1) = out(curpair).mduetA(kk).rsNorm;
        sumdat.mduetAutogenous.rsRaw(end+1) = out(curpair).mduetA(kk).rsRaw;
    end
    for kk = 1:length(out(curpair).fduetA)
        sumdat.fduetAutogenous.rsNorm(end+1) = out(curpair).fduetA(kk).rsNorm;
        sumdat.fduetAutogenous.rsRaw(end+1) = out(curpair).fduetA(kk).rsRaw;
    end
    
    % Auditory activity (Heterogenous) during duet (other bird's syllables)
    out(curpair).mduetH = rs(in((curpair*2)-1), fduetsyls{curpair}, spon(:,curpair), pad);
    out(curpair).fduetH = rs(in(curpair*2), mduetsyls{curpair}, spon(:,curpair), pad);
    
    for kk = 1:length(out(curpair).mduetH)
        sumdat.mduetHeterogenous.rsNorm(end+1) = out(curpair).mduetH(kk).rsNorm;
        sumdat.mduetHeterogenous.rsRaw(end+1) = out(curpair).mduetH(kk).rsRaw;
    end
    for kk = 1:length(out(curpair).fduetH)
        sumdat.fduetHeterogenous.rsNorm(end+1) = out(curpair).fduetH(kk).rsNorm;
        sumdat.fduetHeterogenous.rsRaw(end+1) = out(curpair).fduetH(kk).rsRaw;
    end

end


%% Plot 

%%%% MOTOR

% For the Normalized data. Not doing STD because of asymetry of
% distribution.
    m(1) = mean(sumdat.mduetAutogenous.rsNorm); s(1) = std(sumdat.mduetAutogenous.rsNorm);
    m(2) = mean(sumdat.mSolo.rsNorm); s(2) = std(sumdat.mSolo.rsNorm);
%     m(3) = mean(sumdat.fduetAutogenous.rsNorm); s(3) = std(sumdat.fduetAutogenous.rsNorm);
%     m(4) = mean(sumdat.fSolo.rsNorm); s(4) = std(sumdat.fSolo.rsNorm);
    m(3) = mean(sumdat.fduetAutogenous.rsNorm(6:end-8)); s(3) = std(sumdat.fduetAutogenous.rsNorm(6:end-8));
    m(4) = mean(sumdat.fSolo.rsNorm); s(4) = std(sumdat.fSolo.rsNorm);



% For raw data
    mraw(1) = mean(sumdat.mduetAutogenous.rsNorm); sraw(1) = std(sumdat.mduetAutogenous.rsRaw);
    mraw(2) = mean(sumdat.mSolo.rsNorm); sraw(3) = std(sumdat.mSolo.rsRaw);
    mraw(3) = mean(sumdat.fduetAutogenous.rsNorm); sraw(2) = std(sumdat.fduetAutogenous.rsRaw);
    mraw(4) = mean(sumdat.fSolo.rsNorm); sraw(4) = std(sumdat.fSolo.rsRaw);

figure(1); clf; 

subplot(121); hold on;
    plot([1 2], m(1:2), 'bo'); 
    errorbar([1 2], m(1:2), s(1:2), 'b' );
    plot([3 4], m(3:4), 'mo'); 
    errorbar([3 4], m(3:4), s(3:4), 'm' );
    ylim([-5 20]); xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');

% subplot(132); hold on;
%     plot([1 2], mraw(1:2), 'b*'); 
%     errorbar([1 2], mraw(1:2), sraw(1:2), 'b' );% /sqrt(length(fChron)));
%     plot([3 4], mraw(3:4), 'm*'); hold on;
%     errorbar([3 4], mraw(3:4), sraw(3:4), 'm' );% /sqrt(length(mChron)));

%%%% SENSORY

subplot(122); hold on;
    plot(2, mean(sumdat.mAud.rsNorm), 'bo');
        errorbar(2, mean(sumdat.mAud.rsNorm), std(sumdat.mAud.rsNorm), 'b');
    plot(1, mean(sumdat.mduetHeterogenous.rsNorm), 'bo');
        errorbar(1, mean(sumdat.mduetHeterogenous.rsNorm), std(sumdat.mduetHeterogenous.rsNorm), 'b');
    plot(4, mean(sumdat.fAud.rsNorm), 'mo');
        errorbar(4, mean(sumdat.fAud.rsNorm), std(sumdat.mAud.rsNorm), 'm');
    plot(3, mean(sumdat.fduetHeterogenous.rsNorm), 'mo');
        errorbar(3, mean(sumdat.fduetHeterogenous.rsNorm), std(sumdat.fduetHeterogenous.rsNorm), 'm');
    ylim([-5 20]); xlim([0.5 4.5]); plot([1,4], [0,0], 'k-');

%% Compute stats     
    
% Autogenous duet RS significant from 0?
    stats.m.dNAuto.mean = mean(sumdat.mduetAutogenous.rsNorm);
    [stats.m.dNAuto.H, stats.m.dNAuto.P, stats.m.dNAuto.CI, stats.m.dNAuto.stats]  = ttest(sumdat.mduetAutogenous.rsNorm);
    stats.m.dRAuto.mean = mean(sumdat.mduetAutogenous.rsRaw);
    [stats.m.dRAuto.H, stats.m.dRAuto.P, stats.m.dRAuto.CI, stats.m.dRAuto.stats]  = ttest(sumdat.mduetAutogenous.rsRaw);
    stats.f.dNAuto.mean = mean(sumdat.fduetAutogenous.rsNorm);
    [stats.f.dNAuto.H, stats.f.dNAuto.P, stats.f.dNAuto.CI, stats.f.dNAuto.stats]  = ttest(sumdat.fduetAutogenous.rsNorm);
    stats.f.dRAuto.mean = mean(sumdat.fduetAutogenous.rsRaw);
    [stats.f.dRAuto.H, stats.f.dRAuto.P, stats.f.dRAuto.CI, stats.f.dRAuto.stats]  = ttest(sumdat.fduetAutogenous.rsRaw);
    
% Autogenous Solo RS significant from 0?

    stats.m.sNAuto.mean = mean(sumdat.mSolo.rsNorm);
    [stats.m.sNAuto.H, stats.m.sNAuto.P, stats.m.sNAuto.CI, stats.m.sNAuto.stats]  = ttest(sumdat.mSolo.rsNorm);
    stats.m.sRAuto.mean = mean(sumdat.mSolo.rsRaw);
    [stats.m.sRAuto.H, stats.m.sRAuto.P, stats.m.sRAuto.CI, stats.m.sRAuto.stats]  = ttest(sumdat.mSolo.rsRaw);
    stats.f.sNAuto.mean = mean(sumdat.fSolo.rsNorm);
    [stats.f.sNAuto.H, stats.f.sNAuto.P, stats.f.sNAuto.CI, stats.f.sNAuto.stats]  = ttest(sumdat.fSolo.rsNorm);
    stats.m.sRAuto.mean = mean(sumdat.fSolo.rsRaw);
    [stats.f.sRAuto.H, stats.f.sRAuto.P, stats.f.sRAuto.CI, stats.f.sRAuto.stats]  = ttest(sumdat.fSolo.rsRaw);

% Difference between Autogenous Duet and Solo RS motor?

[stats.m.SvsDNAuto.H, stats.m.SvsDNAuto.P, stats.m.SvsDNAuto.CI, stats.m.SvsDNAuto.stats]  = ttest2(sumdat.mSolo.rsNorm, sumdat.mduetAutogenous.rsNorm);
[stats.m.SvsDRAuto.H, stats.m.SvsDRAuto.P, stats.m.SvsDRAuto.CI, stats.m.SvsDRAuto.stats]  = ttest2(sumdat.mSolo.rsRaw, sumdat.mduetAutogenous.rsRaw);
[stats.f.SvsDNAuto.H, stats.f.SvsDNAuto.P, stats.f.SvsDNAuto.CI, stats.f.SvsDNAuto.stats]  = ttest2(sumdat.fSolo.rsNorm, sumdat.fduetAutogenous.rsNorm);
[stats.f.SvsDRAuto.H, stats.f.SvsDRAuto.P, stats.f.SvsDRAuto.CI, stats.f.SvsDRAuto.stats]  = ttest2(sumdat.fSolo.rsRaw, sumdat.fduetAutogenous.rsRaw);

% Heterogenous duet RS significant from 0?

    stats.m.dNHetero.mean = mean(sumdat.mduetHeterogenous.rsNorm);
    [stats.m.dNHetero.H, stats.m.dNHetero.P, stats.m.dNHetero.CI, stats.m.dNHetero.stats]  = ttest(sumdat.mduetHeterogenous.rsNorm);
    stats.m.dRHetero.mean = mean(sumdat.mduetHeterogenous.rsRaw);
    [stats.m.dRHetero.H, stats.m.dRHetero.P, stats.m.dRHetero.CI, stats.m.dRHetero.stats]  = ttest(sumdat.mduetHeterogenous.rsRaw);
    stats.f.dNHetero.mean = mean(sumdat.fduetHeterogenous.rsNorm);
    [stats.f.dNHetero.H, stats.f.dNHetero.P, stats.f.dNHetero.CI, stats.f.dNHetero.stats]  = ttest(sumdat.fduetHeterogenous.rsNorm);
    stats.f.dRHetero.mean = mean(sumdat.fduetHeterogenous.rsRaw);
    [stats.f.dRHetero.H, stats.f.dRHetero.P, stats.f.dRHetero.CI, stats.f.dRHetero.stats]  = ttest(sumdat.fduetHeterogenous.rsRaw);
    
% Heterogenous Solo RS significant from 0?

    stats.m.sNHetero.mean = mean(sumdat.mAud.rsNorm);
    [stats.m.sNHetero.H, stats.m.sNHetero.P, stats.m.sNHetero.CI, stats.m.sNHetero.stats]  = ttest(sumdat.mAud.rsNorm);
    stats.m.sRHetero.mean = mean(sumdat.mAud.rsRaw);
    [stats.m.sRHetero.H, stats.m.sRHetero.P, stats.m.sRHetero.CI, stats.m.sRHetero.stats]  = ttest(sumdat.mAud.rsRaw);
    stats.f.sNHetero.mean = mean(sumdat.fAud.rsNorm);
    [stats.f.sNHetero.H, stats.f.sNHetero.P, stats.f.sNHetero.CI, stats.f.sNHetero.stats]  = ttest(sumdat.fAud.rsNorm);
    stats.f.sRHetero.mean = mean(sumdat.fAud.rsRaw);
    [stats.f.sRHetero.H, stats.f.sRHetero.P, stats.f.sRHetero.CI, stats.f.sRHetero.stats]  = ttest(sumdat.fAud.rsRaw);

% Difference between Heterogenous Duet and Solo RS motor?

[stats.m.SvsDNHetero.H, stats.m.SvsDNHetero.P, stats.m.SvsDNHetero.CI, stats.m.SvsDNHetero.stats]  = ttest2(sumdat.mAud.rsNorm, sumdat.mduetHeterogenous.rsNorm);
[stats.m.SvsDRHetero.H, stats.m.SvsDRHetero.P, stats.m.SvsDRHetero.CI, stats.m.SvsDRHetero.stats]  = ttest2(sumdat.mAud.rsRaw, sumdat.mduetHeterogenous.rsRaw);
[stats.f.SvsDNHetero.H, stats.f.SvsDNHetero.P, stats.f.SvsDNHetero.CI, stats.f.SvsDNHetero.stats]  = ttest2(sumdat.fAud.rsNorm, sumdat.fduetHeterogenous.rsNorm);
[stats.f.SvsDRHetero.H, stats.f.SvsDRHetero.P, stats.f.SvsDRHetero.CI, stats.f.SvsDRHetero.stats]  = ttest2(sumdat.fAud.rsRaw, sumdat.fduetHeterogenous.rsRaw);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%% Response Strength nested function

function qwe = rs(struc, syllabl, spontan, padme)
 
   % Get spontaneous rate 
    sponSpikeCount = 0;  
    
        for i=1:4 % 4 electrodes in a tetrode always
            sponSpikeCount = sponSpikeCount + length(find(struc.Cspikes{i} > spontan(1) & struc.Cspikes{i} < spontan(2)));
        end

        sponrate = sponSpikeCount / (spontan(2) - spontan(1));
        sponrate = sponrate/4;
   
    % Loop for each syllable
    
    for j = 1:length(syllabl)    
    
    % Get the spikes for that syllable
    stimSpikeCount = 0; 
    
        for i=1:4 % 4 electrodes in a tetrode always
            stimSpikeCount = stimSpikeCount + length(find(struc.Cspikes{i} > struc.syl(syllabl(j)).tim(1)-padme & struc.Cspikes{i} < struc.syl(syllabl(j)).tim(2)-padme));
        end
        
        stimrate = stimSpikeCount / (struc.syl(syllabl(j)).tim(2) - struc.syl(syllabl(j)).tim(1));
        stimrate = stimrate/4;
        
        
        qwe(j).sylnum = syllabl(j);
        qwe(j).rsNorm = (stimrate - sponrate) / sponrate + 0.0000000000001;
        qwe(j).rsRaw = stimrate - sponrate;
        qwe(j).sponrate = sponrate;
        qwe(j).spikerate = stimrate;
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
    

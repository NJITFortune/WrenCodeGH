function [out, sumdat, stts] = wRS_Acute(in, padding)
% Usage: Calculates response strength to solo and duet syllables.
% in the Urethane data
% load ChronicCompleat2018d.mat (Current as of 4-Jan-2019)


sumdat.fDuetAuto.rsNorm = []; sumdat.fDuetAuto.rsRaw = [];
sumdat.mDuetAuto.rsNorm = []; sumdat.mDuetAuto.rsRaw = [];
sumdat.fDuetHetero.rsNorm = []; sumdat.fDuetHetero.rsRaw = [];
sumdat.mDuetHetero.rsNorm = []; sumdat.mDuetHetero.rsRaw = [];

sumdat.fSoloAuto.rsNorm = []; sumdat.fSoloAuto.rsRaw = [];
sumdat.mSoloAuto.rsNorm = []; sumdat.mSoloAuto.rsRaw = [];
sumdat.fSoloHetero.rsNorm = []; sumdat.fSoloHetero.rsRaw = [];
sumdat.mSoloHetero.rsNorm = []; sumdat.mSoloHetero.rsRaw = [];

% load ChronicCompleat2017p.mat

pad = 0; 

if nargin == 2; pad = padding; end
% The boundaries of each syllable are used for this analysis. This is
% inherently problematic as we expect pre-motor activity in awake animals to occur PRIOR to
% the sound and auditory activity in urethane-anesthetized animals to occur AFTER the sound.
% We are comfortable with this approach because it is a rather unbiased. Change
% the value of padding in seconds (e.g. 0.005 or -0.003) to look at the effects.
    

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
    
%     % Solo syllables MALE
%     if ~isempty(msolosyls{curpair}) % Male sang solo syllables
%         if ~isempty(in(curpair*2).Aspikes) % Female physiology exists
%             out(curpair).fSoloHetero = rs(in(curpair*2), msolosyls{curpair}, spon(:,curpair), pad);
%             for kk = 1:length(msolosyls{curpair})
%                 sumdat.fSoloHetero.rsNorm(end+1) = out(curpair).fSoloHetero(kk).rsNorm;
%                 sumdat.fSoloHetero.rsRaw(end+1) = out(curpair).fSoloHetero(kk).rsRaw;
%             end
%         end        
%         if ~isempty(in((curpair*2)-1).Aspikes) % Male physiology exists
%             out(curpair).mSoloAuto = rs(in((curpair*2)-1), msolosyls{curpair}, spon(:,curpair), pad);
%             for kk = 1:length(msolosyls{curpair})
%                 sumdat.mSoloAuto.rsNorm(end+1) = out(curpair).mSoloAuto(kk).rsNorm;
%                 sumdat.mSoloAuto.rsRaw(end+1) = out(curpair).mSoloAuto(kk).rsRaw;
%             end
%         end
%     end
    
    % Solo syllables FEMALE
    if ~isempty(fsolosyls{curpair}) % Female sang solo syllables
        if ~isempty(in(curpair*2).Aspikes) % Female physiology exists
            out(curpair).fSoloAuto = rs(in(curpair*2), fsolosyls{curpair}, spon(:,curpair), pad);
            for kk = 1:length(fsolosyls{curpair})
                sumdat.fSoloAuto.rsNorm(end+1) = out(curpair).fSoloAuto(kk).rsNorm;
                sumdat.fSoloAuto.rsRaw(end+1) = out(curpair).fSoloAuto(kk).rsRaw;
            end
        end        
        if ~isempty(in((curpair*2)-1).Aspikes) % Male physiology exists
            out(curpair).mSoloHetero = rs(in((curpair*2)-1), fsolosyls{curpair}, spon(:,curpair), pad);
            for kk = 1:length(fsolosyls{curpair})
                sumdat.mSoloHetero.rsNorm(end+1) = out(curpair).mSoloHetero(kk).rsNorm;
                sumdat.mSoloHetero.rsRaw(end+1) = out(curpair).mSoloHetero(kk).rsRaw;
            end
        end
    end
    
    % Duet syllables FEMALE
    if ~isempty(fduetsyls{curpair}) % Female sang duet syllables
        if ~isempty(in(curpair*2).Aspikes) % Female physiology exists
            out(curpair).fDuetAuto = rs(in(curpair*2), fduetsyls{curpair}, spon(:,curpair), pad);
            for kk = 1:length(fduetsyls{curpair})
                sumdat.fDuetAuto.rsNorm(end+1) = out(curpair).fDuetAuto(kk).rsNorm;
                sumdat.fDuetAuto.rsRaw(end+1) = out(curpair).fDuetAuto(kk).rsRaw;
            end
        end        
        if ~isempty(in((curpair*2)-1).Aspikes) % Male physiology exists
            out(curpair).mDuetHetero = rs(in((curpair*2)-1), fduetsyls{curpair}, spon(:,curpair), pad);
            for kk = 1:length(fduetsyls{curpair})
                sumdat.mDuetHetero.rsNorm(end+1) = out(curpair).mDuetHetero(kk).rsNorm;
                sumdat.mDuetHetero.rsRaw(end+1) = out(curpair).mDuetHetero(kk).rsRaw;
            end
        end
    end
    
    % Duet syllables MALE
    if ~isempty(mduetsyls{curpair}) % Male sang duet syllables
        if ~isempty(in(curpair*2).Aspikes) % Female physiology exists
            out(curpair).fDuetHetero = rs(in(curpair*2), mduetsyls{curpair}, spon(:,curpair), pad);
            for kk = 1:length(mduetsyls{curpair})
                sumdat.fDuetHetero.rsNorm(end+1) = out(curpair).fDuetHetero(kk).rsNorm;
                sumdat.fDuetHetero.rsRaw(end+1) = out(curpair).fDuetHetero(kk).rsRaw;
            end
        end        
        if ~isempty(in((curpair*2)-1).Aspikes) % Male physiology exists
            out(curpair).mDuetAuto = rs(in((curpair*2)-1), mduetsyls{curpair}, spon(:,curpair), pad);
            for kk = 1:length(mduetsyls{curpair})
                sumdat.mDuetAuto.rsNorm(end+1) = out(curpair).mDuetAuto(kk).rsNorm;
                sumdat.mDuetAuto.rsRaw(end+1) = out(curpair).mDuetAuto(kk).rsRaw;
            end
        end
    end
    
%     % Autogenous activity during duet (bird's own syllables)
%     if ~isempty(in((curpair*2)-1).Aspikes)
%         out(curpair).mduetA = rs(in((curpair*2)-1), mduetsyls{curpair}, spon(:,curpair), pad);
%     end
%     if ~isempty(in(curpair*2).Aspikes)
%         out(curpair).fduetA = rs(in(curpair*2), fduetsyls{curpair}, spon(:,curpair), pad);
%     end
%     %FIX
%     for kk = 1:length(out(curpair).mduetA);
%         sumdat.mduetAutogenous.rsNorm(end+1) = out(curpair).mduetA(kk).rsNorm;
%         sumdat.mduetAutogenous.rsRaw(end+1) = out(curpair).mduetA(kk).rsRaw;
%     end
%     for kk = 1:length(out(curpair).fduetA)
%         sumdat.fduetAutogenous.rsNorm(end+1) = out(curpair).fduetA(kk).rsNorm;
%         sumdat.fduetAutogenous.rsRaw(end+1) = out(curpair).fduetA(kk).rsRaw;
%     end
%     
%     % Heterogenous activity during duet (other bird's syllables)
%     if ~isempty(in((curpair*2)-1).Aspikes)
%         out(curpair).mduetH = rs(in((curpair*2)-1), fduetsyls{curpair}, spon(:,curpair), pad);
%     end
%     if ~isempty(in(curpair*2).Aspikes)
%         out(curpair).fduetH = rs(in(curpair*2), mduetsyls{curpair}, spon(:,curpair), pad);
%     end
%     %FIX
%     for kk = 1:length(out(curpair).mduetH)
%         sumdat.mduetHeterogenous.rsNorm(end+1) = out(curpair).mduetH(kk).rsNorm;
%         sumdat.mduetHeterogenous.rsRaw(end+1) = out(curpair).mduetH(kk).rsRaw;
%     end
%     for kk = 1:length(out(curpair).fduetH)
%         sumdat.fduetHeterogenous.rsNorm(end+1) = out(curpair).fduetH(kk).rsNorm;
%         sumdat.fduetHeterogenous.rsRaw(end+1) = out(curpair).fduetH(kk).rsRaw;
%     end

end


%% Plot 

% For the Normalized data. 

mNorm(1) = mean(sumdat.mDuetAuto.rsNorm); sNorm(1) = std(sumdat.mDuetAuto.rsNorm);
mNorm(2) = mean(sumdat.mSoloAuto.rsNorm); sNorm(2) = std(sumdat.mSoloAuto.rsNorm);
mNorm(3) = mean(sumdat.fDuetAuto.rsNorm); sNorm(3) = std(sumdat.fDuetAuto.rsNorm);
mNorm(4) = mean(sumdat.fSoloAuto.rsNorm); sNorm(4) = std(sumdat.fSoloAuto.rsNorm);

% For raw data

mraw(1) = mean(sumdat.mDuetAuto.rsRaw); sraw(1) = std(sumdat.mDuetAuto.rsRaw);
mraw(2) = mean(sumdat.mSoloAuto.rsRaw); sraw(3) = std(sumdat.mSoloAuto.rsRaw);
mraw(3) = mean(sumdat.fDuetAuto.rsRaw); sraw(2) = std(sumdat.fDuetAuto.rsRaw);
mraw(4) = mean(sumdat.fSoloAuto.rsRaw); sraw(4) = std(sumdat.fSoloAuto.rsRaw);

% AUTOGENOUS
figure(1); clf; 
subplot(121); hold on; % NORMALIZED
    plot([1 2], mNorm(1:2), 'bo'); 
        errorbar([1 2], mNorm(1:2), sNorm(1:2), 'b' );
    plot([3 4], mNorm(3:4), 'mo'); 
        errorbar([3 4], mNorm(3:4), sNorm(3:4), 'm' );
    xlim([0.5 4.5]); 
    
subplot(122); hold on; % RAW
    plot([1 2], mraw(1:2), 'bo'); 
        errorbar([1 2], mraw(1:2), sraw(1:2), 'b' );
    plot([3 4], mraw(3:4), 'mo'); 
        errorbar([3 4], mraw(3:4), sraw(3:4), 'm' );
    xlim([0.5 4.5]); 

% HETEROGENOUS
figure(2); clf; 
subplot(121); hold on; % NORMALIZED
    plot(1, mean(sumdat.mDuetHetero.rsNorm), 'bo');
        errorbar(1, mean(sumdat.mDuetHetero.rsNorm), std(sumdat.mDuetHetero.rsNorm), 'b');
    plot(2, mean(sumdat.mSoloHetero.rsNorm), 'bo');
        errorbar(2, mean(sumdat.mSoloHetero.rsNorm), std(sumdat.mSoloHetero.rsNorm), 'b');
    plot(3, mean(sumdat.fDuetHetero.rsNorm), 'mo');
        errorbar(3, mean(sumdat.fDuetHetero.rsNorm), std(sumdat.fDuetHetero.rsNorm), 'm');
    plot(4, mean(sumdat.fSoloHetero.rsNorm), 'mo');
        errorbar(4, mean(sumdat.fSoloHetero.rsNorm), std(sumdat.fSoloHetero.rsNorm), 'm');
    xlim([0.5 4.5]);

subplot(122); hold on; % RAW
    plot(1, mean(sumdat.mDuetHetero.rsRaw), 'bo');
        errorbar(1, mean(sumdat.mDuetHetero.rsRaw), std(sumdat.mDuetHetero.rsRaw), 'b');
    plot(2, mean(sumdat.mSoloHetero.rsRaw), 'bo');
        errorbar(2, mean(sumdat.mSoloHetero.rsRaw), std(sumdat.mSoloHetero.rsRaw), 'b');
    plot(3, mean(sumdat.fDuetHetero.rsRaw), 'mo');
        errorbar(3, mean(sumdat.fDuetHetero.rsRaw), std(sumdat.fDuetHetero.rsRaw), 'm');
    plot(4, mean(sumdat.fSoloHetero.rsRaw), 'mo');
        errorbar(4, mean(sumdat.fSoloHetero.rsRaw), std(sumdat.fSoloHetero.rsRaw), 'm');
    xlim([0.5 4.5]);
    
%% Compute stats     
    
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
stts.m.sRAuto.mean = mean(sumdat.fSoloAuto.rsRaw);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%% Response Strength nested function

function qwe = rs(struc, syllabl, spontan, padme)
 
   % Get spontaneous rate 
    sponSpikeCount = 0;  
    
        tt = length(struc.Aspikes);
        for i=1:tt % each rep 
            sponSpikeCount = sponSpikeCount + length(find(struc.Aspikes{i} > spontan(1) & struc.Aspikes{i} < spontan(2)));
        end

        sponrate = sponSpikeCount / (spontan(2) - spontan(1)); 
        sponrate = sponrate/tt; % Should be spikes per second.
   
    % Loop for each syllable
    
    for j = 1:length(syllabl)    
    
    % Get the spikes for that syllable
    stimSpikeCount = 0; 
    
        for i=1:tt % for each rep
            stimSpikeCount = stimSpikeCount + length(find(struc.Aspikes{i} > struc.syl(syllabl(j)).tim(1)-padme & struc.Aspikes{i} < struc.syl(syllabl(j)).tim(2)-padme));
        end
        
        stimrate = stimSpikeCount / (struc.syl(syllabl(j)).tim(2) - struc.syl(syllabl(j)).tim(1));
        stimrate = stimrate/tt;
        
        
        qwe(j).sylnum = syllabl(j);
        qwe(j).rsNorm = (stimrate - sponrate) / sponrate;
        qwe(j).rsRaw = stimrate - sponrate;
        qwe(j).sponrate = sponrate;
        qwe(j).spikerate = stimrate;
        qwe(j).spontim = spontan;
        qwe(j).syltim = [struc.syl(syllabl(j)).tim(1), struc.syl(syllabl(j)).tim(2)];
        qwe(j).pad = padme;
        
    end
end % End of embedded function rs

end % End of the function


%     % Get firing rates during Auditory-only 
%         MaudSpikeCount = 0; Mauduration = 0; Msylcount = 0;
%         FaudSpikeCount = 0; Fauduration = 0; Fsylcount = 0;
% 
%     if ~isempty(fsolosyls{j})
%         for k = 1:length(fsolosyls{j})
%             for i=1:4
%                 MaudSpikeCount = MaudSpikeCount + length(find(w((j*2)-1).Aspikes{i} > w((j*2)-1).syl(fsolosyls{j}(k)).tim(1) & w((j*2)-1).Aspikes{i} < w((j*2)-1).syl(fsolosyls{j}(k)).tim(2)));
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
%                 FaudSpikeCount = FaudSpikeCount + length(find(w((j*2)).Aspikes{i} > w((j*2)).syl(msolosyls{j}(k)).tim(1) & w((j*2)).Aspikes{i} < w((j*2)-1).syl(msolosyls{j}(k)).tim(2)));
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
% mChist = bs_swPSTH(w(birdnum).Aspikes, [-5 10], 100, 0);
% fChist = bs_swPSTH(w(birdnum+1).Aspikes, [-5 10], 100, 0);
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
%     bs_raster(w(birdnum).Aspikes, 'b'); xlim([w(birdnum).tim(1) w(birdnum).syl(lastsyl).tim(1)]);
%     bs_raster(w(birdnum+1).Aspikes, 'm'); xlim([w(birdnum).tim(1) w(birdnum).syl(lastsyl).tim(1)]);


% subplot(411); specgram(w(i).duet(tt), 512, w(i).Fs); ylim([500 4500]); caxis([-50 20]);
% subplot(412); bs_raster(w(i).Aspikes, 'b'); xlim(xlimits);
% subplot(413); bs_raster(w(i+1).Aspikes, 'm'); xlim(xlimits);
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
    

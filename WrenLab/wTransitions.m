function [M, F] = wTransitions(in, wwindow)
% Usage out = wTransitions(w, window)
% This generates plots that focus on the transitions between female and male syllables. 
% w is the data structure 'w' from ChronicCompleat2019f.mat (Current as of 8-Jan-2020)
% window is the time before and after the start of the 2nd syllable in 
% pair for plotting. Default is 0.250 (250msec).  300 looks good too.

%% Preparations
% Default window width for histogram if user didn't specify window

    widow = 0.250; % Time before and after transition. 300 msec looks pretty good. 250 for publication
    numbins = 5; % How many bins before and after the onset of our focal syllable?
    
    % If user specified a different window length, use that.
    if nargin > 1; widow = wwindow(1); numbins = wwindow(2); end   
    
    windur = widow / numbins;

%% Load the list of Chronic singing data with syllable indices and locations for spontaneous activity

[msolosyls, mduetsyls, fsolosyls, fduetsyls, Cspon, Aspon] = wData;

% Some variables because I am not a talented coder.
    Fsyldur = []; Msyldur = [];  % A side calculation - durations of syllables
    M2FISI = []; F2MISI = [];    % Another side calculation - durations of InterSyllable Intervals
    
    % Keeping track of how many solo syllables that have produced data
    Mwhichmalesolosyl = 0; Mwhichfemalesolosyl = 0; 
    Fwhichmalesolosyl = 0; Fwhichfemalesolosyl = 0;

    % More variables that we need
    MAHU(1).SPS = []; MAHU(1).RSraw = []; MAHU(1).RSnorm = [];

    MAHC = MAHU; MHAU = MAHU; MHAC = MAHU; 
    MSAU = MAHU; MSAC = MAHU; MSHU = MAHU; MSHC = MAHU;

    FAHU = MAHU; FAHC = MAHU; FHAU = MAHU; FHAC = MAHU; 
    FSAU = MAHU; FSAC = MAHU; FSHU = MAHU; FSHC = MAHU;

%% Cycle for each bird in our list

% Choose which data to analyze
%    birdlist = 16:-1:1; % All compleat data with both ACUTE (urethane) AND CHRONIC (awake)
birdlist = 1:length(in);

for ff = birdlist
        
    sylstrdx = ceil(ff/2); % Apologies. The syllable indices from wData.m 
                           % each refer to two entries in w, one for each male datum (odd entries)
                           % and one for each female datum (even entries). This
                           % just resolves that indexing issue. 
    
% Variables to hold sylable start times for the current entry

    currM2Fsyltim = []; % Male to Female transition
    currF2Msyltim = []; % Female to Male transition
    currMsolosyltims = []; % Male Solo syllable
    currFsolosyltims = []; % Female Solo syllable

% Calculate spontaneous rates for current entry

    ChronSpon = 0; % CHRONIC DATA
    
    for z = 1:length(in(ff).Cspikes)
        ChronSpon = ChronSpon + length(find(in(ff).Cspikes{z} > Cspon(1,sylstrdx) & in(ff).Cspikes{z} < Cspon(2,sylstrdx)));
    end
        ChronSpon = ChronSpon / (Cspon(2,sylstrdx) - Cspon(1,sylstrdx)); % Divide by duration, SPIKES PER SECOND 
        ChronSpon = ChronSpon / length(in(ff).Cspikes); % Divide by number of reps (always 4 for Chronic)
            
    AcuteSpon = 0; % ACUTE DATA
    
    if sylstrdx < 7 % We only have complete Urethane data for the first 3 pairs of birds
        for z = 1:length(in(ff).Aspikes)
            AcuteSpon = AcuteSpon + length(find(in(ff).Aspikes{z} > Aspon(1,sylstrdx) & in(ff).Aspikes{z} < Aspon(2,sylstrdx)));
        end
            AcuteSpon = AcuteSpon / length(in(ff).Aspikes); % Divide by number of reps 
            AcuteSpon = AcuteSpon / (Aspon(2,sylstrdx) - Aspon(1,sylstrdx)); % Divide by duration, SPIKES PER SECOND
    end

% Find all Male-female and female-male duet syllable transitions
        
    % Find male to female syllable transitions M2F
    for t = 1:length(mduetsyls{sylstrdx}) % For every male syllable in the entry...
       if ~isempty(find(fduetsyls{sylstrdx} == mduetsyls{sylstrdx}(t)+1, 1)) % If the next syllable is female
           
           currentFemaleIndex = mduetsyls{sylstrdx}(t)+1; % Not a necessary step
           currM2Fsyltim(end+1) = in(ff).syl(currentFemaleIndex).tim(1); % Time of the start of the female syllable
           
           Fsyldur(end+1) = in(ff).syl(currentFemaleIndex).tim(2) - in(ff).syl(currentFemaleIndex).tim(1); % Female syllable duration
           M2FISI(end+1) = in(ff).syl(currentFemaleIndex).tim(1) - in(ff).syl(currentFemaleIndex-1).tim(2); % Duration of ISI
           
       end
    end
    
    % Find female to male syllable transitions F2M
    for t = 1:length(fduetsyls{sylstrdx}) % For every female syllable in the entry...
       if ~isempty(find(mduetsyls{sylstrdx} == fduetsyls{sylstrdx}(t)+1, 1)) % If the next syllable is male
           
           currentMaleIndex = fduetsyls{ceil(ff/2)}(t)+1; % Not a necessary stp
           currF2Msyltim(end+1) = in(ff).syl(currentMaleIndex).tim(1); % Time of the start of the male syllable

           Msyldur(end+1) = in(ff).syl(currentMaleIndex).tim(2) -in(ff).syl(currentMaleIndex).tim(1); % Male syllable duration
           F2MISI(end+1) = in(ff).syl(currentMaleIndex).tim(1) - in(ff).syl(currentMaleIndex-1).tim(2); % Duration of ISI
       end
    end

    % Find Solo Male syllables
    for t = 1:length(msolosyls{sylstrdx}) % For every male syllable
         currMsolosyltims(end+1) = in(ff).syl(msolosyls{sylstrdx}(t)).tim(1); % Time of the start of the male syllable
    end
    
    % Find Solo Female syllables
    for t = 1:length(fsolosyls{sylstrdx}) % For every female syllable
         currFsolosyltims(end+1) = in(ff).syl(fsolosyls{sylstrdx}(t)).tim(1); % Time of the start of the female syllable
    end

%% Generate the transition histograms

if in(ff).sexy == 1 % This is a male
        
    if ~isempty(mduetsyls{sylstrdx}) % For songs with male duet syllables...   
        %  Acute
        if ~isempty(in(ff).Aspikes)
            [tmp, M.bintims] = wPhaseHist(in(ff).Aspikes, currM2Fsyltim, widow, numbins, AcuteSpon);
                for kk = 1:length(tmp); MAHU(end+1) = tmp(kk); end; clear tmp;
            [tmp, ~] = wPhaseHist(in(ff).Aspikes, currF2Msyltim, widow, numbins, AcuteSpon);        
                for kk = 1:length(tmp); MHAU(end+1) = tmp(kk); end; clear tmp;
        end
        % Chronic    
            [tmp, ~] = wPhaseHist(in(ff).Cspikes, currM2Fsyltim, widow, numbins, ChronSpon);        
                for kk = 1:length(tmp); MAHC(end+1) = tmp(kk); end; clear tmp;   
            [tmp, ~] = wPhaseHist(in(ff).Cspikes, currF2Msyltim, widow, numbins, ChronSpon);
                for kk = 1:length(tmp); MHAC(end+1) = tmp(kk); end; clear tmp;
    end
    
    if ~isempty(msolosyls{sylstrdx}) % For songs with male solo syllable
        % Acute
        if ~isempty(in(ff).Aspikes)
        [tmp, ~] = wPhaseHist(in(ff).Aspikes, currMsolosyltims, widow, numbins, AcuteSpon);
            for kk = 1:length(tmp); MSAU(end+1) = tmp(kk); end; clear tmp; 
        end
        % Chronic
        [tmp, ~] = wPhaseHist(in(ff).Cspikes, currMsolosyltims, widow, numbins, ChronSpon);
            for kk = 1:length(tmp); MSAC(end+1) = tmp(kk); end; clear tmp;
    end
    
    if ~isempty(fsolosyls{sylstrdx}) % % For songs with female solo syllable
        % Acute
        if ~isempty(in(ff).Aspikes)        
        [tmp, ~] = wPhaseHist(in(ff).Aspikes, currFsolosyltims, widow, numbins, AcuteSpon);
            for kk = 1:length(tmp); MSHU(end+1) = tmp(kk); end; clear tmp; 
        end
        % Chronic
        [tmp, ~] = wPhaseHist(in(ff).Cspikes, currFsolosyltims, widow, numbins, ChronSpon);
            for kk = 1:length(tmp); MSHC(end+1) = tmp(kk); end; clear tmp;           
    end

end % End of male

if in(ff).sexy == 2 % This is a female
    
    if ~isempty(mduetsyls{sylstrdx}) % For songs with male duet syllables...   
    
    % Acute
    if ~isempty(in(ff).Aspikes)
        [tmp, F.bintims] = wPhaseHist(in(ff).Aspikes, currF2Msyltim, widow, numbins, AcuteSpon);
            for kk = 1:length(tmp); FAHU(end+1) = tmp(kk); end; clear tmp;
        [tmp, ~] = wPhaseHist(in(ff).Aspikes, currM2Fsyltim, widow, numbins, AcuteSpon);
            for kk = 1:length(tmp); FHAU(end+1) = tmp(kk); end; clear tmp;
    end
    
    % Chronic
        if ~isempty(currF2Msyltim)
            [tmp, F.bintims] = wPhaseHist(in(ff).Cspikes, currF2Msyltim, widow, numbins, ChronSpon);
             for kk = 1:length(tmp)
                 FAHC(end+1).SPS = tmp(kk).SPS; 
                 FAHC(end).RSraw = tmp(kk).RSraw; 
                 FAHC(end).RSnorm = tmp(kk).RSnorm; 
             end
             clear tmp; 
        end
        if ~isempty(currM2Fsyltim)
            [tmp, ~] = wPhaseHist(in(ff).Cspikes, currM2Fsyltim, widow, numbins, ChronSpon);
             for kk = 1:length(tmp) 
                 FHAC(end+1).SPS = tmp(kk).SPS; 
                 FHAC(end).RSraw = tmp(kk).RSraw; 
                 FHAC(end).RSnorm = tmp(kk).RSnorm; 
                    figure(27); axr(1) = subplot(121); title('Female Autogenous Duet'); hold on; plot(tmp(kk).SPS);
             end
             clear tmp;
        end
    
    bins4plot = (F(1).bintims(2:end) + F(1).bintims(1:end-1)) /2; % Time bins adjusted for proper plotting
    
    end
    
    if ~isempty(msolosyls{sylstrdx}) % For songs with male solo syllable
        if ~isempty(in(ff).Aspikes)
            [tmp, ~] = wPhaseHist(in(ff).Aspikes, currMsolosyltims, widow, numbins, AcuteSpon);
            for kk = 1:length(tmp); FSHU(end+1) = tmp(kk); end; clear tmp;
        end
            [tmp, ~] = wPhaseHist(in(ff).Cspikes, currMsolosyltims, widow, numbins, ChronSpon);
             for kk = 1:length(tmp); FSHC(end+1) = tmp(kk); end; clear tmp;
    end
    
    if ~isempty(fsolosyls{sylstrdx}) % For songs with female solo syllable
        if ~isempty(in(ff).Aspikes)
        [tmp, ~] = wPhaseHist(in(ff).Aspikes, currFsolosyltims, widow, numbins, AcuteSpon);
            for kk = 1:length(tmp); FSAU(end+1) = tmp(kk); end; clear tmp;
        end
         [tmp, ~] = wPhaseHist(in(ff).Cspikes, currFsolosyltims, widow, numbins, ChronSpon);
            for kk = 1:length(tmp) 
                FSAC(end+1) = tmp(kk); 
                figure(27); axr(2) = subplot(122); title('Female Autogenous Solo');hold on; plot(tmp(kk).SPS);
                linkaxes(axr, 'y');
            end; clear tmp;
    end
    
end % End of female
                
end % End of cycling for every bird


%% Build the Cool fill plots for the DUET data

    msFHAC = concatHist(FHAC);                        
    msFHAU = concatHist(FHAU);
    msFAHC = concatHist(FAHC);
    msFAHU = concatHist(FAHU);

    msMHAC = concatHist(MHAC);                        
    msMHAU = concatHist(MHAU);
    msMAHC = concatHist(MAHC);
    msMAHU = concatHist(MAHU);


figure(1); clf; set(gcf, 'Color', [1,1,1]);

% PLOT M2F DATA 

axc(1) = subplot(221); hold on; title('M2F Chronic'); 
    plot([0 0], [-10 55], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
% Male
%fill([bins4plot bins4plot(end:-1:1)], [msMAHC.meanRSraw - msMAHC.steRSraw, msMAHC.meanRSraw(end:-1:1) + msMAHC.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
%plot(bins4plot, msMAHC.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
fill([bins4plot bins4plot(end:-1:1)], [msMAHC.meanSPS - msMAHC.steSPS, msMAHC.meanSPS(end:-1:1) + msMAHC.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMAHC.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
%fill([bins4plot bins4plot(end:-1:1)], [msFHAC.meanRSraw - msFHAC.steRSraw, msFHAC.meanRSraw(end:-1:1) + msFHAC.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
%plot(bins4plot, msFHAC.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
fill([bins4plot bins4plot(end:-1:1)], [msFHAC.meanSPS - msFHAC.steSPS, msFHAC.meanSPS(end:-1:1) + msFHAC.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFHAC.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

axu(1) = subplot(223); hold on; title('M2F Urethane'); 
    plot([0 0], [-10 55], 'k-', 'LineWidth', 2); %plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    set(axu(1),'Color', [0.9, 0.9, 0.9]);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMAHU.meanSPS - msMAHU.steSPS, msMAHU.meanSPS(end:-1:1) + msMAHU.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMAHU.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFHAU.meanSPS - msFHAU.steSPS, msFHAU.meanSPS(end:-1:1) + msFHAU.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFHAU.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    
% PLOT F2M DATA

axc(2) = subplot(222); hold on; title('F2M Chronic'); 
    plot([0 0], [-10 55], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
% Female
%fill([bins4plot bins4plot(end:-1:1)], [msFAHC.meanRSraw - msFAHC.steRSraw, msFAHC.meanRSraw(end:-1:1) + msFAHC.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
%plot(bins4plot, msFAHC.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
fill([bins4plot bins4plot(end:-1:1)], [msFAHC.meanSPS - msFAHC.steSPS, msFAHC.meanSPS(end:-1:1) + msFAHC.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFAHC.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Male
%fill([bins4plot bins4plot(end:-1:1)], [msMHAC.meanRSraw - msMHAC.steRSraw, msMHAC.meanRSraw(end:-1:1) + msMHAC.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
%plot(bins4plot, msMHAC.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
fill([bins4plot bins4plot(end:-1:1)], [msMHAC.meanSPS - msMHAC.steSPS, msMHAC.meanSPS(end:-1:1) + msMHAC.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMHAC.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

axu(2) = subplot(224); hold on; title('F2M Urethane'); 
    plot([0 0], [-10 55], 'k-', 'LineWidth', 2); %plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    set(axu(2),'Color', [0.9, 0.9, 0.9]);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFAHU.meanSPS - msFAHU.steSPS, msFAHU.meanSPS(end:-1:1) + msFAHU.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFAHU.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMHAU.meanSPS - msMHAU.steSPS, msMHAU.meanSPS(end:-1:1) + msMHAU.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMHAU.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

linkaxes(axc, 'xy'); figure(1); subplot(221);  ylim([-8 60]);  xlim([-widow-0.0001, widow+0.0001]);
linkaxes(axu, 'xy'); figure(1); subplot(223);  ylim([0 55]);  xlim([-widow-0.0001, widow+0.0001]);

% figure(3); clf; set(gcf, 'Color', [1,1,1]);
% xax(1) = subplot(221); hold on; title('M2F Chronic'); 
%     plot([0 0], [-1 5], 'k-', 'LineWidth', 1); plot([min(bins4plot), max(bins4plot)], [0, 0], 'k-', 'LineWidth', 1);
% plot(bins4plot, msMAHC.RSN, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% plot(bins4plot, msFHAC.RSN, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% 
% xax(2) = subplot(223); hold on; title('M2F Urethane'); 
% set(xax(2),'Color', [0.9, 0.9, 0.9]);
%     plot([0 0], [-1 50], 'k-', 'LineWidth', 1); ylim([-1 50]);
% plot(bins4plot, msMAHU.RSN, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% plot(bins4plot, msFHAU.RSN, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% 
% xax(3) = subplot(222); hold on; title('F2M Chronic'); 
%     plot([0 0], [-1 5], 'k-', 'LineWidth', 1); plot([min(bins4plot), max(bins4plot)], [0, 0], 'k-', 'LineWidth', 1);
% plot(bins4plot, msFAHC.RSN, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% plot(bins4plot, msMHAC.RSN, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% 
% xax(4) = subplot(224); hold on; title('F2M Urethane'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
% set(xax(4),'Color', [0.9, 0.9, 0.9]);
%     plot([0 0], [-1 50], 'k-', 'LineWidth', 1); ylim([-1 50]);
% plot(bins4plot, msFAHU.RSN, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% plot(bins4plot, msMHAU.RSN, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% 
% linkaxes(xax, 'x');


%% Do the solo syllable calculations

    msFSAC = concatHist(FSAC);                        
    msFSAU = concatHist(FSAU);                        
    msFSHC = concatHist(FSHC);                        
    msFSHU = concatHist(FSHU);                        

    msMSAC = concatHist(MSAC);                        
    msMSAU = concatHist(MSAU);                        
    msMSHC = concatHist(MSHC);                        
    msMSHU = concatHist(MSHU);                        

figure(2); clf; set(gcf, 'Color', [1,1,1]);

% PLOT M Solo DATA
axxc(1) = subplot(221); hold on; title('M Solo Chronic'); 
    plot([0 0], [-10 55], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
% Male
%fill([bins4plot bins4plot(end:-1:1)], [msMSAC.meanRSraw - msMSAC.steRSraw, msMSAC.meanRSraw(end:-1:1) + msMSAC.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
%plot(bins4plot, msMSAC.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
fill([bins4plot bins4plot(end:-1:1)], [msMSAC.meanSPS - msMSAC.steSPS, msMSAC.meanSPS(end:-1:1) + msMSAC.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSAC.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
%fill([bins4plot bins4plot(end:-1:1)], [msFSHC.meanRSraw - msFSHC.steRSraw, msFSHC.meanRSraw(end:-1:1) + msFSHC.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
%plot(bins4plot, msFSHC.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
fill([bins4plot bins4plot(end:-1:1)], [msFSHC.meanSPS - msFSHC.steSPS, msFSHC.meanSPS(end:-1:1) + msFSHC.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSHC.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

axxu(1) = subplot(223); hold on; title('M Solo Urethane'); 
    plot([0 0], [-10 55], 'k-', 'LineWidth', 2); %plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    set(axxu(1),'Color', [0.9, 0.9, 0.9]);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMSAU.meanSPS - msMSAU.steSPS, msMSAU.meanSPS(end:-1:1) + msMSAU.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSAU.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSHU.meanSPS - msFSHU.steSPS, msFSHU.meanSPS(end:-1:1) + msFSHU.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSHU.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

% PLOT F Solo DATA
axxc(2) = subplot(222); hold on; title('F Solo Chronic'); 
    plot([0 0], [-10 55], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
% Female
%fill([bins4plot bins4plot(end:-1:1)], [msFSAC.meanRSraw - msFSAC.steRSraw, msFSAC.meanRSraw(end:-1:1) + msFSAC.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
%plot(bins4plot, msFSAC.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
fill([bins4plot bins4plot(end:-1:1)], [msFSAC.meanSPS - msFSAC.steSPS, msFSAC.meanSPS(end:-1:1) + msFSAC.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSAC.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Male
%fill([bins4plot bins4plot(end:-1:1)], [msMSHC.meanRSraw - msMSHC.steRSraw, msMSHC.meanRSraw(end:-1:1) + msMSHC.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
%plot(bins4plot, msMSHC.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
fill([bins4plot bins4plot(end:-1:1)], [msMSHC.meanSPS - msMSHC.steSPS, msMSHC.meanSPS(end:-1:1) + msMSHC.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSHC.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

axxu(2) = subplot(224); hold on; title('F Solo Urethane'); 
    plot([0 0], [-10 55], 'k-', 'LineWidth', 2); %plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    set(axxu(2),'Color', [0.9, 0.9, 0.9]);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSAU.meanSPS - msFSAU.steSPS, msFSAU.meanSPS(end:-1:1) + msFSAU.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSAU.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMSHU.meanSPS - msMSHU.steSPS, msMSHU.meanSPS(end:-1:1) + msMSHU.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSHU.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

linkaxes(axxc, 'xy'); figure(2); subplot(221); ylim([-8 60]);
linkaxes(axxu, 'xy'); figure(2); subplot(223);  ylim([0 55]);

%% PLOTS FOR MANUSCRIPTS

figure(7); % For luck
clf;
% SOLO AUTOGENOUS CHRONIC PANEL
fmc(1) = subplot(221); hold on; title('Auto Solo Chronic'); 
    plot([0 0], [-10 70], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    
    % Male
%    fill([bins4plot bins4plot(end:-1:1)], [msMSAC.meanRSraw - msMSAC.steRSraw, msMSAC.meanRSraw(end:-1:1) + msMSAC.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
%    plot(bins4plot, msMSAC.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    fill([bins4plot bins4plot(end:-1:1)], [msMSAC.meanSPS - msMSAC.steSPS, msMSAC.meanSPS(end:-1:1) + msMSAC.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
    plot(bins4plot, msMSAC.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    % Female
%    fill([bins4plot bins4plot(end:-1:1)], [msFSAC.meanRSraw - msFSAC.steRSraw, msFSAC.meanRSraw(end:-1:1) + msFSAC.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
%    plot(bins4plot, msFSAC.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    fill([bins4plot bins4plot(end:-1:1)], [msFSAC.meanSPS - msFSAC.steSPS, msFSAC.meanSPS(end:-1:1) + msFSAC.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
    plot(bins4plot, msFSAC.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

% DUET AUTOGENOUS CHRONIC PANEL
fmc(2) = subplot(222); hold on; title('Auto Duet Chronic'); 
    plot([0 0], [-10 70], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    
    % Male
%    fill([bins4plot bins4plot(end:-1:1)], [msMHAC.meanRSraw - msMHAC.steRSraw, msMHAC.meanRSraw(end:-1:1) + msMHAC.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
%    plot(bins4plot, msMHAC.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    fill([bins4plot bins4plot(end:-1:1)], [msMHAC.meanSPS - msMHAC.steSPS, msMHAC.meanSPS(end:-1:1) + msMHAC.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
    plot(bins4plot, msMHAC.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    % Female
%    fill([bins4plot bins4plot(end:-1:1)], [msFHAC.meanRSraw - msFHAC.steRSraw, msFHAC.meanRSraw(end:-1:1) + msFHAC.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
%    plot(bins4plot, msFHAC.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    fill([bins4plot bins4plot(end:-1:1)], [msFHAC.meanSPS - msFHAC.steSPS, msFHAC.meanSPS(end:-1:1) + msFHAC.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
    plot(bins4plot, msFHAC.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    
% SOLO HETEROGENOUS CHRONIC PANEL
fmc(3) = subplot(223); hold on; title('Hetero Solo Chronic'); 
    plot([0 0], [-10 70], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    
    % Male
%    fill([bins4plot bins4plot(end:-1:1)], [msMSHC.meanRSraw - msMSHC.steRSraw, msMSHC.meanRSraw(end:-1:1) + msMSHC.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
%    plot(bins4plot, msMSHC.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    fill([bins4plot bins4plot(end:-1:1)], [msMSHC.meanSPS - msMSHC.steSPS, msMSHC.meanSPS(end:-1:1) + msMSHC.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
    plot(bins4plot, msMSHC.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    % Female
%    fill([bins4plot bins4plot(end:-1:1)], [msFSHC.meanRSraw - msFSHC.steRSraw, msFSHC.meanRSraw(end:-1:1) + msFSHC.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
%    plot(bins4plot, msFSHC.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    fill([bins4plot bins4plot(end:-1:1)], [msFSHC.meanSPS - msFSHC.steSPS, msFSHC.meanSPS(end:-1:1) + msFSHC.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
    plot(bins4plot, msFSHC.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

% DUET HETEROGENOUS CHRONIC PANEL
fmc(4) = subplot(224); hold on; title('Hetero Duet Chronic'); 
    plot([0 0], [-10 70], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    
    % Male
%    fill([bins4plot bins4plot(end:-1:1)], [msMAHC.meanRSraw - msMAHC.steRSraw, msMAHC.meanRSraw(end:-1:1) + msMAHC.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
%    plot(bins4plot, msMAHC.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    fill([bins4plot bins4plot(end:-1:1)], [msMAHC.meanSPS - msMAHC.steSPS, msMAHC.meanSPS(end:-1:1) + msMAHC.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
    plot(bins4plot, msMAHC.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    % Female
%    fill([bins4plot bins4plot(end:-1:1)], [msFAHC.meanRSraw - msFAHC.steRSraw, msFAHC.meanRSraw(end:-1:1) + msFAHC.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
%    plot(bins4plot, msFAHC.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    fill([bins4plot bins4plot(end:-1:1)], [msFAHC.meanSPS - msFAHC.steSPS, msFAHC.meanSPS(end:-1:1) + msFAHC.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
    plot(bins4plot, msFAHC.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

linkaxes(fmc, 'xy'); figure(7); subplot(221); ylim([-8 63]); xlim([-widow-0.0001, widow+0.0001]);

figure(8); clf;

fmu(1) = subplot(221); hold on; title('HA Duet Acute');
    plot([0 0], [0 55], 'k-', 'LineWidth', 2); %plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    set(fmu(1),'Color', [0.9, 0.9, 0.9]);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMHAU.meanSPS - msMHAU.steSPS, msMHAU.meanSPS(end:-1:1) + msMHAU.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMHAU.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFHAU.meanSPS - msFHAU.steSPS, msFHAU.meanSPS(end:-1:1) + msFHAU.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFHAU.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

fmu(2) = subplot(223); hold on; title('Solo Female Syllable Acute');
    plot([0 0], [0 55], 'k-', 'LineWidth', 2); %plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    set(fmu(2),'Color', [0.9, 0.9, 0.9]);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMSHU.meanSPS - msMSHU.steSPS, msMSHU.meanSPS(end:-1:1) + msMSHU.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSHU.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSAU.meanSPS - msFSAU.steSPS, msFSAU.meanSPS(end:-1:1) + msFSAU.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSAU.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

fmu(3) = subplot(224); hold on; title('Solo Male Syllable Acute');
    plot([0 0], [0 55], 'k-', 'LineWidth', 2); %plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    set(fmu(3),'Color', [0.9, 0.9, 0.9]);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMSAU.meanSPS - msMSAU.steSPS, msMSAU.meanSPS(end:-1:1) + msMSAU.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSAU.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSHU.meanSPS - msFSHU.steSPS, msFSHU.meanSPS(end:-1:1) + msFSHU.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSHU.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);


fmu(4) = subplot(222); hold on; title('AH Duet Acute');
    plot([0 0], [0 55], 'k-', 'LineWidth', 2); %plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    set(fmu(4),'Color', [0.9, 0.9, 0.9]);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMAHU.meanSPS - msMAHU.steSPS, msMAHU.meanSPS(end:-1:1) + msMAHU.steSPS(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMAHU.meanSPS, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFAHU.meanSPS - msFAHU.steSPS, msFAHU.meanSPS(end:-1:1) + msFAHU.steSPS(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFAHU.meanSPS, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);


linkaxes(fmu, 'xy'); figure(8); subplot(223);  ylim([0 50]); xlim([-widow-0.0001, widow+0.0001]);

%% Put everything into the output structures

M(1).MAHC = MAHC; M(1).MAHU = MAHU;
M(1).MHAC = MHAC; M(1).MHAU = MHAU;
M(1).MSAC = MSAC; M(1).MSAU = MSAU;
M(1).MSHC = MSHC; M(1).MSHU = MSHU;
M(1).msMAHC = msMAHC; M(1).msMAHU = msMAHU;
M(1).msMHAC = msMHAC; M(1).msMHAU = msMHAU;
M(1).msMSAC = msMSAC; M(1).msMSAU = msMSAU;
M(1).msMSHC = msMSHC; M(1).msMSHU = msMSHU;

F(1).FAHC = FAHC; F(1).FAHU = FAHU;
F(1).FHAC = FHAC; F(1).FHAU = FHAU;
F(1).FSAC = FSAC; F(1).FSAU = FSAU;
F(1).FSHC = FSHC; F(1).FSHU = FSHU;
F(1).msFAHC = msFAHC; F(1).msFAHU = msFAHU;
F(1).msFHAC = msFHAC; F(1).msFHAU = msFHAU;
F(1).msFSAC = msFSAC; F(1).msFSAU = msFSAU;
F(1).msFSHC = msFSHC; F(1).msFSHU = msFSHU;


%% Statistics, if you please.
    
fprintf('The mean and std for Male syllable duration is  %1.3f %1.3f \n', mean(Msyldur), std(Msyldur));
fprintf('The mean and std for Female syllable duration is  %1.3f %1.3f \n', mean(Fsyldur), std(Fsyldur));
fprintf('The mean and std for M2F ISI is  %1.3f %1.3f \n', mean(M2FISI), std(M2FISI));
fprintf('The mean and std for F2M ISI is  %1.3f %1.3f \n', mean(F2MISI), std(F2MISI));

%% Embedded Concatonation function
function tuo = concatHist(xin)
        
    for qq = length(xin):-1:1
        
        save xin.mat xin
        
                if ~exist('cSPS', 'var')
                    if ~isempty(xin(qq).SPS)
                    cSPS(:,1) = xin(qq).SPS;
                    cRSrw(:,1) = xin(qq).RSraw;
                    cRSnrm(:,1) = xin(qq).RSnorm;
                    end
                else
                    if ~isempty(xin(qq).SPS)
                    cSPS(:,end+1) = xin(qq).SPS;
                    cRSnrm(:,end+1) = xin(qq).RSnorm;
                    cRSrw(:,end+1) = xin(qq).RSraw;
                    end
                end
    end
        
    for j = 1:length(cSPS(:,1))
        meanSPS(j) = mean(cSPS(j,:)); 
        stdSPS(j) = std(cSPS(j,:)); 
            steSPS(j) = stdSPS(j) / sqrt(length(xin));
        meanRSnorm(j) = mean(cRSnrm(j,:)); 
        stdRSnorm(j) = std(cRSnrm(j,:));
            steRSnorm(j) = stdRSnorm(j) / sqrt(length(xin));
        meanRSraw(j) = mean(cRSrw(j,:)); 
        stdRSraw(j) = std(cRSrw(j,:));
            steRSraw(j) = stdRSraw(j) / sqrt(length(xin));
    end
        
    tuo.meanSPS = meanSPS;
    tuo.meanRSnorm = meanRSnorm;
    tuo.meanRSraw = meanRSraw;
    tuo.stdSPS = stdSPS;
    tuo.stdRSnorm = stdRSnorm;
    tuo.stdRSraw = stdRSraw;
    tuo.steSPS = steSPS;
    tuo.steRSnorm = steRSnorm;
    tuo.steRSraw = steRSraw;
    
end

%% Embedded histogram function
function [out, bintims] = wPhaseHist(spiketimes, tims, wid, numbin, sponSPS)

    out = [];
    
        % wid is window in msec before and after start of the syllable at the end of the focal ISI.
        binwid = wid / numbin; % Width of each bin
        
        % Specify the OVERLAP percentage here
        overlap = 50; % Overlap is XX% of previous window
        
        overlap = 1-(overlap/100); % Converts to step size for advancing the window
        
        bintims = -wid:binwid*overlap:wid; % List of bin times centered on zero

        

    for j = length(tims):-1:1 % For each syllable in the list
        
        bintimstarts = tims(j)-wid:binwid*overlap:tims(j)+wid-(binwid*overlap); % Start and end times for current syllable
        
        spkcnts = zeros(1, length(-wid:binwid*overlap:wid-(binwid*overlap)));
                
        for k = 1:length(spiketimes) % For each row of spikes
            for m = 1:length(bintimstarts) % For each bin of our PSTH
                spkcnts(m) = spkcnts(m) + length(find(spiketimes{k} > bintimstarts(m) & spiketimes{k} < bintimstarts(m)+binwid));
            end
        end
         
        spikearray(:,j) = spkcnts / length(spiketimes); % mean spikes per rep in each bin
 
    end
    
    % Convert raw spike counts useful measures
        
%         for k = 1:length(spikearray(:,1))
%             SPShist(k) = sum(spikearray(k,:)); 
%             SPShist(k) = SPShist(k) / totalreps; % Divide by number of 'reps'
%             SPShist(k) = SPShist(k) / binwid; % Divide by length of bin
%             SPShist(k) = SPShist(k) / length(spikearray(1,:)); %Divide by number of transitions
% 
%             RSrawhist(k) = SPShist(k) - sponSPS; % Subtract Spontaneous rate
%             RSnorm(k) = RSrawhist(k) / (sponSPS + 0.0001); % Divide by Spontaneous rate
%         end

        goodlist = [];
        for j=1:length(spikearray(1,:))
            if ~isnan(spikearray(1,j))
                goodlist = [goodlist j];
            end
        end

        for j=1:length(goodlist)
            out(j).SPS = spikearray(:, goodlist(j)) / binwid; 
            out(j).RSraw = (spikearray(:, goodlist(j)) / binwid) - sponSPS;
            out(j).RSnorm = ((spikearray(:, goodlist(j)) / binwid) - sponSPS) / sponSPS;
        end
        
        
end

end



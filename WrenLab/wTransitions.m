function [M, F] = wTransitions(in)
% Usage out = wTransitions(w, window)
% This generates plots that focus on the transitions between female and male syllables. 
% w is the data structure 'w' from ChronicCompleat2019f.mat (Current as of 8-Jan-2020)
% window is the time before and after the start of the 2nd syllable in 
% pair for plotting. Default is 0.250 (250msec).  300 looks good too.

%% Preparations
% Default window width for histogram if user didn't specify window
    widow = 0.500; % 300 msec looks pretty good with numbins 10 and overlap 50
    numbins = 5; % How many bins before and after the onset of our focal syllable?
    
    windur = widow / numbins;

%% Load the list of Chronic singing data with syllable indices and locations for spontaneous activity

[msolosyls, mduetsyls, fsolosyls, fduetsyls, Cspon, Aspon] = wData;

% Some variables because I am not a talented coder.
    Fsyldur = []; Msyldur = [];  % A side calculation - durations of syllables
    M2FISI = []; F2MISI = [];    % Another side calculation - durations of InterSyllable Intervals
    
    % Keeping track of how many solo syllables that have produced data
    Mwhichmalesolosyl = 0; Mwhichfemalesolosyl = 0; 
    Fwhichmalesolosyl = 0; Fwhichfemalesolosyl = 0;
    Mwhichduet = 0; Fwhichduet = 0;

% Choose which data to analyze
%    birdlist = 16:-1:1; % All compleat data with both ACUTE (urethane) AND CHRONIC (awake)
birdlist = 1:length(in);
    
%% For each bird in our list

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

MAHU(1).SPS = []; MAHU(1).RSraw = []; MAHU(1).RSnorm = [];
MAHC = MAHU; MHAU = MAHU; MHAC = MAHU; 
MSAU = MAHU; MSAC = MAHU; MSHU = MAHU; MSHC = MAHU;

FAHU = MAHU; FAHC = MAHU; FHAU = MAHU; FHAC = MAHU; 
FSAU = MAHU; FSAC = MAHU; FSHU = MAHU; FSHC = MAHU;


if in(ff).sexy == 1 % This is a male
    
    % Use the PhaseCut embedded function to create the histogram for Duet
    % data.  A=Autogenous, H=Heterogenous, U=Urethane, C=Chronic
    
    if ~isempty(mduetsyls{sylstrdx}) % For every male duet syllable...   
    Mwhichduet = Mwhichduet + 1;
    
    [tmp, M(Mwhichduet).bintims] = wPhaseHist(in(ff).Aspikes, currM2Fsyltim, widow, numbins, AcuteSpon);
        for kk = length(tmp) 
            MAHU(end+1).SPS = tmp(kk).SPS; MAHU(end+1).RSnorm = tmp(kk).RSnorm; MAHU(end+1).RSraw = tmp(kk).RSraw;
        end
    [tmp, ~] = wPhaseHist(in(ff).Cspikes, currM2Fsyltim, widow, numbins, ChronSpon);        
        for kk = length(tmp); MAHC(end+1) = tmp(kk); end
    [tmp, ~] = wPhaseHist(in(ff).Aspikes, currF2Msyltim, widow, numbins, AcuteSpon);        
        for kk = length(tmp); MHAU(end+1) = tmp(kk); end
    [tmp, ~] = wPhaseHist(in(ff).Cspikes, currF2Msyltim, widow, numbins, ChronSpon);
        for kk = length(tmp); MHAC(end+1) = tmp(kk); end
    
    bins4plot = (M(Mwhichduet).bintims(2:end) + M(Mwhichduet).bintims(1:end-1))/2; % Time bins adjusted for proper plotting
    end
    
    % Solo data S=Solo, F=Female, M=Male, A=Autogenous, H=Heterogenous,
    % U=Urethane, C=Chronic
    if ~isempty(msolosyls{sylstrdx})
        Mwhichmalesolosyl = Mwhichmalesolosyl+1; % We are using a different indexing here
        [tmp, ~] = wPhaseHist(in(ff).Aspikes, currMsolosyltims, widow, numbins, AcuteSpon);
         for kk = length(tmp); MSAU(end+1) = tmp(kk); end
       [tmp, ~] = wPhaseHist(in(ff).Cspikes, currMsolosyltims, widow, numbins, ChronSpon);
        for kk = length(tmp); MSAC(end+1) = tmp(kk); end
    end
    if ~isempty(fsolosyls{sylstrdx})
        Mwhichfemalesolosyl = Mwhichfemalesolosyl +1;
        [tmp, ~] = wPhaseHist(in(ff).Aspikes, currFsolosyltims, widow, numbins, AcuteSpon);
         for kk = length(tmp); MSHU(end+1) = tmp(kk); end
        [tmp, ~] = wPhaseHist(in(ff).Cspikes, currFsolosyltims, widow, numbins, ChronSpon);
         for kk = length(tmp); MSHC(end+1) = tmp(kk); end
    end

end % End of male

if in(ff).sexy == 2 % This is a female
    
    % Use the PhaseCut embedded function to create the histogram for Duet
    % data.  A=Autogenous, H=Heterogenous, U=Urethane, C=Chronic
    if ~isempty(fduetsyls{sylstrdx})    
    Fwhichduet = Fwhichduet + 1;
    
    if ~isempty(in(ff).Aspikes)
        [tmp, F(Fwhichduet).bintims] = wPhaseHist(in(ff).Aspikes, currF2Msyltim, widow, numbins, AcuteSpon);
         for kk = length(tmp); FAHU(end+1) = tmp(kk); end
        [tmp, ~] = wPhaseHist(in(ff).Aspikes, currM2Fsyltim, widow, numbins, AcuteSpon);
         for kk = length(tmp); FHAU(end+1) = tmp(kk); end
    end
    
    if ~isempty(in(ff).Cspikes)
    if ~isempty(currF2Msyltim)
        [tmp, F(Fwhichduet).bintims] = wPhaseHist(in(ff).Cspikes, currF2Msyltim, widow, numbins, ChronSpon);
         for kk = length(tmp); FAHC(end+1) = tmp(kk); end
    end
    if ~isempty(currM2Fsyltim)
        [tmp, ~] = wPhaseHist(in(ff).Cspikes, currM2Fsyltim, widow, numbins, ChronSpon);
         for kk = length(tmp); FHAC(end+1) = tmp(kk); end
    end
    end
    
    bins4plot = (F(1).bintims(2:end) + F(1).bintims(1:end-1)) /2; % Time bins adjusted for proper plotting
    
    end
    
    % Solo data S=Solo, F=Female, M=Male, A=Autogenous, H=Heterogenous,
    % U=Urethane, C=Chronic
    if ~isempty(msolosyls{sylstrdx})
        Fwhichmalesolosyl = Fwhichmalesolosyl +1; % Again, we are using a different indexing for solo data
        [tmp, ~] = wPhaseHist(in(ff).Aspikes, currMsolosyltims, widow, numbins, AcuteSpon);
         for kk = length(tmp); FSHU(end+1) = tmp(kk); end
        [tmp, ~] = wPhaseHist(in(ff).Cspikes, currMsolosyltims, widow, numbins, ChronSpon);
         for kk = length(tmp); FSHC(end+1) = tmp(kk); end
    end
    if ~isempty(fsolosyls{sylstrdx})
        Fwhichfemalesolosyl = Fwhichfemalesolosyl +1;
        [tmp, ~] = wPhaseHist(in(ff).Aspikes, currFsolosyltims, widow, numbins, AcuteSpon);
         for kk = length(tmp); FSAU(end+1) = tmp(kk); end
        [tmp, ~] = wPhaseHist(in(ff).Cspikes, currFsolosyltims, widow, numbins, ChronSpon);
         for kk = length(tmp); FSAC(end+1) = tmp(kk); end
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
ax(1) = subplot(221); hold on; title('M2F Chronic'); 
    plot([0 0], [-10 80], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMAHC.meanRSraw - msMAHC.steRSraw, msMAHC.meanRSraw(end:-1:1) + msMAHC.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMAHC.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFHAC.meanRSraw - msFHAC.steRSraw, msFHAC.meanRSraw(end:-1:1) + msFHAC.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFHAC.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

ax(2) = subplot(223); hold on; title('M2F Urethane'); 
    plot([0 0], [-10 80], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
set(ax(2),'Color', [0.9, 0.9, 0.9]);

% Male
fill([bins4plot bins4plot(end:-1:1)], [msMAHU.meanRSraw - msMAHU.steRSraw, msMAHU.meanRSraw(end:-1:1) + msMAHU.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMAHU.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFHAU.meanRSraw - msFHAU.steRSraw, msFHAU.meanRSraw(end:-1:1) + msFHAU.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFHAU.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    
% PLOT F2M DATA
ax(3) = subplot(222); hold on; title('F2M Chronic'); 
    plot([0 0], [-10 80], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFAHC.meanRSraw - msFAHC.steRSraw, msFAHC.meanRSraw(end:-1:1) + msFAHC.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFAHC.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMHAC.meanRSraw - msMHAC.steRSraw, msMHAC.meanRSraw(end:-1:1) + msMHAC.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMHAC.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

ax(4) = subplot(224); hold on; title('F2M Urethane'); 
    plot([0 0], [-10 80], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
set(ax(4),'Color', [0.9, 0.9, 0.9]);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFAHU.meanRSraw - msFAHU.steRSraw, msFAHU.meanRSraw(end:-1:1) + msFAHU.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFAHU.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMHAU.meanRSraw - msMHAU.steRSraw, msMHAU.meanRSraw(end:-1:1) + msMHAU.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMHAU.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

linkaxes(ax, 'xy'); figure(1); subplot(222); % ylim([0 1]);

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
axx(1) = subplot(221); hold on; title('M Solo Chronic'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMSAC.meanRSraw - msMSAC.steRSraw, msMSAC.meanRSraw(end:-1:1) + msMSAC.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSAC.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSHC.meanRSraw - msFSHC.steRSraw, msFSHC.meanRSraw(end:-1:1) + msFSHC.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSHC.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

axx(2) = subplot(223); hold on; title('M Solo Urethane'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
set(axx(2),'Color', [0.9, 0.9, 0.9]);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMSAU.meanRSraw - msMSAU.steRSraw, msMSAU.meanRSraw(end:-1:1) + msMSAU.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSAU.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSHU.meanRSraw - msFSHU.steRSraw, msFSHU.meanRSraw(end:-1:1) + msFSHU.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSHU.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

% PLOT F Solo DATA
axx(3) = subplot(222); hold on; title('F Solo Chronic'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSAC.meanRSraw - msFSAC.steRSraw, msFSAC.meanRSraw(end:-1:1) + msFSAC.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSAC.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMSHC.meanRSraw - msMSHC.steRSraw, msMSHC.meanRSraw(end:-1:1) + msMSHC.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSHC.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

axx(4) = subplot(224); hold on; title('F Solo Urethane'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
set(axx(4),'Color', [0.9, 0.9, 0.9]);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSAU.meanRSraw - msFSAU.steRSraw, msFSAU.meanRSraw(end:-1:1) + msFSAU.steRSraw(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSAU.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMSHU.meanRSraw - msMSHU.steRSraw, msMSHU.meanRSraw(end:-1:1) + msMSHU.steRSraw(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSHU.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

linkaxes(axx, 'xy'); 

figure(2); subplot(222); % ylim([0 1]);


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
    
    xin
    
    for qq = length(xin):-1:1
            SPS(:,qq) = xin(qq).SPS;
            RSnrm(:,qq) = xin(qq).RSnorm;
            RSrw(:,qq) = xin(qq).RSraw;
    end
        
    for j = length(xin(1).SPS):-1:1
        meanSPS(j) = mean(SPS(j,:)); stdSPS(j) = std(SPS(j,:)); 
            steSPS(j) = stdSPS(j) / sqrt(length(xin));
        meanRSnorm(j) = mean(RSnrm(j,:)); stdRSnorm(j) = std(RSnrm(j,:));
            steRSnorm(j) = stdRSnorm(j) / sqrt(length(xin));
        meanRSraw(j) = mean(SPS(j,:)); stdRSraw(j) = std(RSrw(j,:));
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

        % wid is window in msec before and after start of the syllable at the end of the focal ISI.
        binwid = wid / numbin; % Width of each bin
        
        % Specify the OVERLAP percentage here
        overlap = 75; % Overlap is XX% of previous window
        
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
        
        out.SPS = spikearray / binwid;
        out.RSraw = out.SPS - sponSPS;
        out.RSnorm = out.RSraw / sponSPS;
        
        
end

end



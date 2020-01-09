function [M, F] = wTransitions(in)
% Usage out = wTransitions(w, window)
% This generates plots that focus on the transitions between female and male syllables. 
% w is the data structure 'w' from ChronicCompleat2019f.mat (Current as of 8-Jan-2020)
% window is the time before and after the start of the 2nd syllable in 
% pair for plotting. Default is 0.250 (250msec).  300 looks good too.

%% Preparations
% Default window width for histogram if user didn't specify window
    widow = 0.250; % 300 msec looks pretty good with numbins 10 and overlap 50
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

    ChronSpon(ff) = 0; % CHRONIC DATA
    
    for z = 1:length(in(ff).Cspikes)
        ChronSpon(ff) = ChronSpon(ff) + length(in(ff).Cspikes{z} > Cspon(1,sylstrdx) & in(ff).Cspikes{z} < Cspon(2,sylstrdx));
    end
        ChronSpon(ff) = ChronSpon(ff) / length(in(ff).Cspikes); % Divide by number of reps (always 4 for Chronic)
        ChronSpon(ff) = ChronSpon(ff) / (Cspon(2,sylstrdx) - Cspon(1,sylstrdx)); % Divide by duration, SPIKES PER SECOND 
        
    
    AcuteSpon(ff) = 0; % ACUTE DATA
    
    if sylstrdx < 7 % We only have complete Urethane data for the first 3 pairs of birds
        for z = 1:length(in(ff).Aspikes)
            AcuteSpon(ff) = AcuteSpon(ff) + length(in(ff).Aspikes{z} > Aspon(1,sylstrdx) & in(ff).Aspikes{z} < Aspon(2,sylstrdx));
        end
            AcuteSpon(ff) = AcuteSpon(ff) / length(in(ff).Aspikes); % Divide by number of reps 
            AcuteSpon(ff) = AcuteSpon(ff) / (Aspon(2,sylstrdx) - Aspon(1,sylstrdx)); % Divide by duration, SPIKES PER SECOND
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
    
    % Use the PhaseCut embedded function to create the histogram for Duet
    % data.  A=Autogenous, H=Heterogenous, U=Urethane, C=Chronic
    
    if ~isempty(mduetsyls{sylstrdx}) % For every male duet syllable...   
    Mwhichduet = Mwhichduet + 1;
    
    [MAHU(Mwhichduet), M(Mwhichduet).bintims] = wPhaseHist(in(ff).Aspikes, currM2Fsyltim, widow, numbins, AcuteSpon(ff));        
    [MAHC(Mwhichduet), ~] = wPhaseHist(in(ff).Cspikes, currM2Fsyltim, widow, numbins, ChronSpon(ff));        
    [MHAU(Mwhichduet), ~] = wPhaseHist(in(ff).Aspikes, currF2Msyltim, widow, numbins, AcuteSpon(ff));        
    [MHAC(Mwhichduet), ~] = wPhaseHist(in(ff).Cspikes, currF2Msyltim, widow, numbins, ChronSpon(ff));
    
    bins4plot = (M(Mwhichduet).bintims(2:end) + M(Mwhichduet).bintims(1:end-1))/2; % Time bins adjusted for proper plotting
    end
    
    % Solo data S=Solo, F=Female, M=Male, A=Autogenous, H=Heterogenous,
    % U=Urethane, C=Chronic
    if ~isempty(msolosyls{sylstrdx})
        Mwhichmalesolosyl = Mwhichmalesolosyl+1; % We are using a different indexing here
        [MSAU(Mwhichmalesolosyl), ~] = wPhaseHist(in(ff).Aspikes, currMsolosyltims, widow, numbins, AcuteSpon(ff));
        [MSAC(Mwhichmalesolosyl), ~] = wPhaseHist(in(ff).Cspikes, currMsolosyltims, widow, numbins, ChronSpon(ff));
    end
    if ~isempty(fsolosyls{sylstrdx})
        Mwhichfemalesolosyl = Mwhichfemalesolosyl +1;
        [MSHU(Mwhichfemalesolosyl), ~] = wPhaseHist(in(ff).Aspikes, currFsolosyltims, widow, numbins, AcuteSpon(ff));
        [MSHC(Mwhichfemalesolosyl), ~] = wPhaseHist(in(ff).Cspikes, currFsolosyltims, widow, numbins, ChronSpon(ff));
    end

end % End of male

if in(ff).sexy == 2 % This is a female
    
    % Use the PhaseCut embedded function to create the histogram for Duet
    % data.  A=Autogenous, H=Heterogenous, U=Urethane, C=Chronic
    if ~isempty(fduetsyls{sylstrdx})    
    Fwhichduet = Fwhichduet + 1;
    
    if ~isempty(in(ff).Aspikes)
        [FAHU(Fwhichduet), F(Fwhichduet).bintims] = wPhaseHist(in(ff).Aspikes, currF2Msyltim, widow, numbins, AcuteSpon(ff));
        [FHAU(Fwhichduet), ~] = wPhaseHist(in(ff).Aspikes, currM2Fsyltim, widow, numbins, AcuteSpon(ff));
    end
    
    if ~isempty(in(ff).Cspikes)
    if ~isempty(currF2Msyltim)
        [FAHC(Fwhichduet), F(Fwhichduet).bintims] = wPhaseHist(in(ff).Cspikes, currF2Msyltim, widow, numbins, ChronSpon(ff));
    end
    if ~isempty(currM2Fsyltim)
        [FHAC(Fwhichduet), ~] = wPhaseHist(in(ff).Cspikes, currM2Fsyltim, widow, numbins, ChronSpon(ff));
    end
    end
    
    bins4plot = (F(1).bintims(2:end) + F(1).bintims(1:end-1)) /2; % Time bins adjusted for proper plotting
    
    end
    
    % Solo data S=Solo, F=Female, M=Male, A=Autogenous, H=Heterogenous,
    % U=Urethane, C=Chronic
    if ~isempty(msolosyls{sylstrdx})
        Fwhichmalesolosyl = Fwhichmalesolosyl +1; % Again, we are using a different indexing for solo data
        [FSHU(Fwhichmalesolosyl), ~] = wPhaseHist(in(ff).Aspikes, currMsolosyltims, widow, numbins, AcuteSpon(ff));
        [FSHC(Fwhichmalesolosyl), ~] = wPhaseHist(in(ff).Cspikes, currMsolosyltims, widow, numbins, ChronSpon(ff));
    end
    if ~isempty(fsolosyls{sylstrdx})
        Fwhichfemalesolosyl = Fwhichfemalesolosyl +1;
        [FSAU(Fwhichfemalesolosyl), ~] = wPhaseHist(in(ff).Aspikes, currFsolosyltims, widow, numbins, AcuteSpon(ff));
        [FSAC(Fwhichfemalesolosyl), ~] = wPhaseHist(in(ff).Cspikes, currFsolosyltims, widow, numbins, ChronSpon(ff));
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
ax(1) = subplot(221); hold on; title('M2F Chronic'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMAHC.meanRSraw - msMAHC.stdRSraw/2, msMAHC.meanRSraw(end:-1:1) + msMAHC.stdRSraw(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMAHC.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFHAC.meanRSraw - msFHAC.stdRSraw/2, msFHAC.meanRSraw(end:-1:1) + msFHAC.stdRSraw(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFHAC.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

ax(2) = subplot(223); hold on; title('M2F Urethane'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
set(ax(2),'Color', [0.9, 0.9, 0.9]);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMAHU.mean - msMAHU.std/2, msMAHU.mean(end:-1:1) + msMAHU.std(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMAHU.mean, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFHAU.mean - msFHAU.std/2, msFHAU.mean(end:-1:1) + msFHAU.std(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFHAU.mean, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
    
% PLOT F2M DATA
ax(3) = subplot(222); hold on; title('F2M Chronic'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFAHC.mean - msFAHC.std/2, msFAHC.mean(end:-1:1) + msFAHC.std(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFAHC.mean, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMHAC.mean - msMHAC.std/2, msMHAC.mean(end:-1:1) + msMHAC.std(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMHAC.mean, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

ax(4) = subplot(224); hold on; title('F2M Urethane'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
set(ax(4),'Color', [0.9, 0.9, 0.9]);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFAHU.mean - msFAHU.std/2, msFAHU.mean(end:-1:1) + msFAHU.std(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFAHU.mean, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMHAU.mean - msMHAU.std/2, msMHAU.mean(end:-1:1) + msMHAU.std(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMHAU.mean, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

linkaxes(ax, 'xy'); figure(1); subplot(222); ylim([0 1]);


figure(3); clf; set(gcf, 'Color', [1,1,1]);
xax(1) = subplot(221); hold on; title('M2F Chronic'); 
    plot([0 0], [-1 5], 'k-', 'LineWidth', 1); plot([min(bins4plot), max(bins4plot)], [0, 0], 'k-', 'LineWidth', 1);
plot(bins4plot, msMAHC.RSN, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
plot(bins4plot, msFHAC.RSN, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

xax(2) = subplot(223); hold on; title('M2F Urethane'); 
set(xax(2),'Color', [0.9, 0.9, 0.9]);
    plot([0 0], [-1 50], 'k-', 'LineWidth', 1); ylim([-1 50]);
plot(bins4plot, msMAHU.RSN, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
plot(bins4plot, msFHAU.RSN, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

xax(3) = subplot(222); hold on; title('F2M Chronic'); 
    plot([0 0], [-1 5], 'k-', 'LineWidth', 1); plot([min(bins4plot), max(bins4plot)], [0, 0], 'k-', 'LineWidth', 1);
plot(bins4plot, msFAHC.RSN, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
plot(bins4plot, msMHAC.RSN, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

xax(4) = subplot(224); hold on; title('F2M Urethane'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
set(xax(4),'Color', [0.9, 0.9, 0.9]);
    plot([0 0], [-1 50], 'k-', 'LineWidth', 1); ylim([-1 50]);
plot(bins4plot, msFAHU.RSN, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
plot(bins4plot, msMHAU.RSN, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

linkaxes(xax, 'x');


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
fill([bins4plot bins4plot(end:-1:1)], [msMSAC.meanRSraw - msMSAC.stdRSraw/2, msMSAC.meanRSraw(end:-1:1) + msMSAC.stdRSraw(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSAC.meanRSraw, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSHC.meanRSraw - msFSHC.stdRSraw/2, msFSHC.meanRSraw(end:-1:1) + msFSHC.stdRSraw(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSHC.meanRSraw, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

axx(2) = subplot(223); hold on; title('M Solo Urethane'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
set(axx(2),'Color', [0.9, 0.9, 0.9]);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMSAU.mean - msMSAU.std/2, msMSAU.mean(end:-1:1) + msMSAU.std(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSAU.mean, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSHU.mean - msFSHU.std/2, msFSHU.mean(end:-1:1) + msFSHU.std(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSHU.mean, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

% PLOT F Solo DATA
axx(3) = subplot(222); hold on; title('F Solo Chronic'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSAC.mean - msFSAC.std/2, msFSAC.mean(end:-1:1) + msFSAC.std(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSAC.mean, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMSHC.mean - msMSHC.std/2, msMSHC.mean(end:-1:1) + msMSHC.std(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSHC.mean, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

axx(4) = subplot(224); hold on; title('F Solo Urethane'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
set(axx(4),'Color', [0.9, 0.9, 0.9]);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSAU.mean - msFSAU.std/2, msFSAU.mean(end:-1:1) + msFSAU.std(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSAU.mean, 'm-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMSHU.mean - msMSHU.std/2, msMSHU.mean(end:-1:1) + msMSHU.std(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSHU.mean, 'b-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);

linkaxes(axx, 'xy'); figure(2); subplot(222); ylim([0 1]);


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
    
    for qq = 1:length(xin)
            SPS(:,qq) = xin(qq).SPS;
            RSnrm(:,qq) = xin(qq).RSnorm;
            RSrw(:,qq) = xin(qq).RSraw;
    end
    
    for j = length(SPS(:,1))
        meanSPS(j) = mean(SPS(j,:)); stdSPS(j) = std(SPS(j,:));
        meanRSnorm(j) = mean(RSnrm(j,:)); stdRSnorm(j) = std(RSnrm(j,:));
        meanRSraw(j) = mean(SPS(j,:)); stdRSraw(j) = std(RSrw(j,:));
    end
    
    tuo.meanSPS = meanSPS;
    tuo.meanRSnorm = meanRSnorm;
    tuo.meanRSraw = meanRSraw;
    tuo.stdSPS = stdSPS;
    tuo.stdRSnorm = stdRSnorm;
    tuo.stdRSraw = stdRSraw;
    
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
         
        spikearray(:,j) = spkcnts;
 
    end
    
    % Convert raw spike counts useful measures
        
        for k = 1:length(spikearray(:,1))
            SPShist(k) = sum(spikearray(k,:)) / length(spikearray(k,:)); % Divide by number of 'reps'
            SPShist(k) = SPShist(k) / binwid; % Divide by length of bin

            RSrawhist(k) = SPShist(k) - sponSPS; % Subtract Spontaneous rate
            RSnorm(k) = RSrawhist(k) / sponSPS; % Divide by Spontaneous rate
        end

        out.SPS = SPShist;
        out.RSraw = RSrawhist;
        out.RSnorm = RSnorm;

end


end



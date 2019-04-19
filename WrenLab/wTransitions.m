function [M, F] = wTransitions(in, widow)
% Usage out = wTransitions(in, window)
% This generates plots that focus on the transitions between
% female and male syllables. 
% in is the data structure w
% window is the time before and after the start of the 2nd syllable in 
% pair for plotting.

%% Preparations
% Default window width for histogram if user didn't specify window
if nargin < 2; widow = 0.500; end % 500 msec looks pretty good with numbins 10 and overlap 50



%% Load the list of Chronic singing data with syllable indices and locations for spontaneous activity

[msolosyls, mduetsyls, fsolosyls, fduetsyls, ~, ~] = wData;

% Some variables because I am not a talented coder.
    Fsyldur = []; Msyldur = [];  % A side calculation - durations of syllables
    M2FISI = []; F2MISI = []; % Another side calculation - durations of InterSyllable Intervals
    
    % Keeping track of how many solo syllables that have produced data
    Mwhichmalesolosyl = 0; Mwhichfemalesolosyl = 0; 
    Fwhichmalesolosyl = 0; Fwhichfemalesolosyl = 0;
    Mwhichduet = 0; Fwhichduet = 0;

% Choose which data to analyze
    birdlist = 16:-1:1; % All compleat data with both ACUTE (urethane) AND CHRONIC (awake)

    
%% For each bird in our list

for ff = birdlist
        
% Get all male-female and female-male duet syllable transitions

    % These will have the list of sylable start times for the current song
    currM2Fsyltim = []; currF2Msyltim = []; currMsolosyltims = []; currFsolosyltims = [];
    
    sylstrdx = ceil(ff/2); % Apologies, but the syllable indices from wData.m 
                               % each refer to two entries in w. This
                               % just resolves that indexing issue.
    
    % Find male to female syllable transitions
    for t = 1:length(mduetsyls{sylstrdx}) % For every male syllable
       if ~isempty(find(fduetsyls{sylstrdx} == mduetsyls{sylstrdx}(t)+1, 1)) % If the next syllable is female
           
           currentFemaleIndex = mduetsyls{sylstrdx}(t)+1; % Not a necessary step
           currM2Fsyltim(end+1) = in(ff).syl(currentFemaleIndex).tim(1); % Time of the start of the female syllable
           
           Fsyldur(end+1) = in(ff).syl(currentFemaleIndex).tim(2) - in(ff).syl(currentFemaleIndex).tim(1); % Female syllable duration
           M2FISI(end+1) = in(ff).syl(currentFemaleIndex).tim(1) - in(ff).syl(currentFemaleIndex-1).tim(2); % Duration of ISI
       end
    end
    
    % Find female to male syllable transitions
    for t = 1:length(fduetsyls{sylstrdx}) % For every female syllable
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
    
%% Fetch the histograms

if in(ff).sexy == 1 % This is a male
    
    % Use the PhaseCut embedded function to create the histogram for Duet
    % data.  A=Autogenous, H=Heterogenous, U=Urethane, C=Chronic
    Mwhichduet = Mwhichduet + 1;
    [MAHU(Mwhichduet).spkcnt, M(Mwhichduet).bintims] = wPhaseCut(in(ff).Aspikes, currM2Fsyltim, widow);
    [MAHC(Mwhichduet).spkcnt, ~] = wPhaseCut(in(ff).Cspikes, currM2Fsyltim, widow);
    [MHAU(Mwhichduet).spkcnt, ~] = wPhaseCut(in(ff).Aspikes, currF2Msyltim, widow);
    [MHAC(Mwhichduet).spkcnt, ~] = wPhaseCut(in(ff).Cspikes, currF2Msyltim, widow);
    
    bins4plot = (M(Mwhichduet).bintims(2:end) + M(Mwhichduet).bintims(1:end-1))/2; % Time bins adjusted for proper plotting
    
    % Solo data S=Solo, F=Female, M=Male, A=Autogenous, H=Heterogenous,
    % U=Urethane, C=Chronic
    if ~isempty(msolosyls{sylstrdx})
        Mwhichmalesolosyl = Mwhichmalesolosyl+1; % We are using a different indexing here
        [MSAU(Mwhichmalesolosyl).spkcnt, ~] = wPhaseCut(in(ff).Aspikes, currMsolosyltims, widow);
        [MSAC(Mwhichmalesolosyl).spkcnt, ~] = wPhaseCut(in(ff).Cspikes, currMsolosyltims, widow);
    end
    if ~isempty(fsolosyls{sylstrdx})
        Mwhichfemalesolosyl = Mwhichfemalesolosyl +1;
        [MSHU(Mwhichfemalesolosyl).spkcnt, ~] = wPhaseCut(in(ff).Aspikes, currFsolosyltims, widow);
        [MSHC(Mwhichfemalesolosyl).spkcnt, ~] = wPhaseCut(in(ff).Cspikes, currFsolosyltims, widow);
    end

end % End of male

if in(ff).sexy == 2 % This is a female
    
    % Use the PhaseCut embedded function to create the histogram for Duet
    % data.  A=Autogenous, H=Heterogenous, U=Urethane, C=Chronic
    Fwhichduet = Fwhichduet + 1;
    [FAHU(Fwhichduet).spkcnt, F(Fwhichduet).bintims] = wPhaseCut(in(ff).Aspikes, currF2Msyltim, widow);
    [FAHC(Fwhichduet).spkcnt, ~] = wPhaseCut(in(ff).Cspikes, currF2Msyltim, widow);
    [FHAU(Fwhichduet).spkcnt, ~] = wPhaseCut(in(ff).Aspikes, currM2Fsyltim, widow);
    [FHAC(Fwhichduet).spkcnt, ~] = wPhaseCut(in(ff).Cspikes, currM2Fsyltim, widow);

    bins4plot = (F(Fwhichduet).bintims(2:end) + F(Fwhichduet).bintims(1:end-1)) /2; % Time bins adjusted for proper plotting

    % Solo data S=Solo, F=Female, M=Male, A=Autogenous, H=Heterogenous,
    % U=Urethane, C=Chronic
    if ~isempty(msolosyls{sylstrdx})
        Fwhichmalesolosyl = Fwhichmalesolosyl +1; % Again, we are using a different indexing for solo data
        [FSHU(Fwhichmalesolosyl).spkcnt, ~] = wPhaseCut(in(ff).Aspikes, currMsolosyltims, widow);
        [FSHC(Fwhichmalesolosyl).spkcnt, ~] = wPhaseCut(in(ff).Cspikes, currMsolosyltims, widow);
    end
    if ~isempty(fsolosyls{sylstrdx})
        Fwhichfemalesolosyl = Fwhichfemalesolosyl +1;
        [FSAU(Fwhichfemalesolosyl).spkcnt, ~] = wPhaseCut(in(ff).Aspikes, currFsolosyltims, widow);
        [FSAC(Fwhichfemalesolosyl).spkcnt, ~] = wPhaseCut(in(ff).Cspikes, currFsolosyltims, widow);
    end
    
end % End of female

end % End of cycling for every bird


%% Build the Cool fill plots for the DUET data

    msFHAC = concatHist(FHAC, length(FHAC(1).spkcnt(:,1)));                        
    msFHAU = concatHist(FHAU, length(FHAU(1).spkcnt(:,1)));
    msFAHC = concatHist(FAHC, length(FAHC(1).spkcnt(:,1)));
    msFAHU = concatHist(FAHU, length(FAHU(1).spkcnt(:,1)));

    msMHAC = concatHist(MHAC, length(MHAC(1).spkcnt(:,1)));                        
    msMHAU = concatHist(MHAU, length(MHAU(1).spkcnt(:,1)));
    msMAHC = concatHist(MAHC, length(MAHC(1).spkcnt(:,1)));
    msMAHU = concatHist(MAHU, length(MAHU(1).spkcnt(:,1)));


figure(3); clf; % PLOT M2F DATA
ax(1) = subplot(211); hold on; title('M2F Chronic'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFHAC.mean - msFHAC.std/2, msFHAC.mean(end:-1:1) + msFHAC.std(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFHAC.mean, 'm-o', 'LineWidth', 2);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMAHC.mean - msMAHC.std/2, msMAHC.mean(end:-1:1) + msMAHC.std(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMAHC.mean, 'b-o', 'LineWidth', 2);

ax(2) = subplot(212); hold on; title('M2F Urethane'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFHAU.mean - msFHAU.std/2, msFHAU.mean(end:-1:1) + msFHAU.std(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFHAU.mean, 'm-o', 'LineWidth', 2);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMAHU.mean - msMAHU.std/2, msMAHU.mean(end:-1:1) + msMAHU.std(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMAHU.mean, 'b-o', 'LineWidth', 2);

linkaxes(ax, 'xy'); ylim([0 1]);
    
figure(4); clf; % PLOT F2M DATA
axx(1) = subplot(211); hold on; title('F2M Chronic'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFAHC.mean - msFAHC.std/2, msFAHC.mean(end:-1:1) + msFAHC.std(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFAHC.mean, 'm-o', 'LineWidth', 2);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMHAC.mean - msMHAC.std/2, msMHAC.mean(end:-1:1) + msMHAC.std(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMHAC.mean, 'b-o', 'LineWidth', 2);

axx(2) = subplot(212); hold on; title('F2M Urethane'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFAHU.mean - msFAHU.std/2, msFAHU.mean(end:-1:1) + msFAHU.std(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFAHU.mean, 'm-o', 'LineWidth', 2);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMHAU.mean - msMHAU.std/2, msMHAU.mean(end:-1:1) + msMHAU.std(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMHAU.mean, 'b-o', 'LineWidth', 2);

linkaxes(axx, 'xy'); ylim([0 1]);

%% Do the solo syllable stuff

    msFSAC = concatHist(FSAC, length(FSAC(1).spkcnt(:,1)));                        
    msFSAU = concatHist(FSAU, length(FSAC(1).spkcnt(:,1)));                        
    msFSHC = concatHist(FSHC, length(FSHC(1).spkcnt(:,1)));                        
    msFSHU = concatHist(FSHU, length(FSHC(1).spkcnt(:,1)));                        

    msMSAC = concatHist(MSAC, length(MSAC(1).spkcnt(:,1)));                        
    msMSAU = concatHist(MSAU, length(MSAC(1).spkcnt(:,1)));                        
    msMSHC = concatHist(MSHC, length(MSHC(1).spkcnt(:,1)));                        
    msMSHU = concatHist(MSHU, length(MSHC(1).spkcnt(:,1)));                        

figure(5); clf; % PLOT M Solo DATA
hxx(1) = subplot(211); hold on; title('M Solo Chronic'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSHC.mean - msFSHC.std/2, msFSHC.mean(end:-1:1) + msFSHC.std(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSHC.mean, 'm-o', 'LineWidth', 2);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMSAC.mean - msMSAC.std/2, msMSAC.mean(end:-1:1) + msMSAC.std(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSAC.mean, 'b-o', 'LineWidth', 2);

hxx(2) = subplot(212); hold on; title('M Solo Urethane'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSHU.mean - msFSHU.std/2, msFSHU.mean(end:-1:1) + msFSHU.std(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSHU.mean, 'm-o', 'LineWidth', 2);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMSAU.mean - msMSAU.std/2, msMSAU.mean(end:-1:1) + msMSAU.std(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSAU.mean, 'b-o', 'LineWidth', 2);

linkaxes(hxx, 'xy'); ylim([0 1]);
    
figure(6); clf; % PLOT F Solo DATA
hyx(1) = subplot(211); hold on; title('F Solo Chronic'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSAC.mean - msFSAC.std/2, msFSAC.mean(end:-1:1) + msFSAC.std(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSAC.mean, 'm-o', 'LineWidth', 2);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMSHC.mean - msMSHC.std/2, msMSHC.mean(end:-1:1) + msMSHC.std(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSHC.mean, 'b-o', 'LineWidth', 2);

hyx(2) = subplot(212); hold on; title('F Solo Urethane'); plot([0 0], [0 1], 'k-', 'LineWidth', 2);
% Female
fill([bins4plot bins4plot(end:-1:1)], [msFSAU.mean - msFSAU.std/2, msFSAU.mean(end:-1:1) + msFSAU.std(end:-1:1)/2], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, msFSAU.mean, 'm-o', 'LineWidth', 2);
% Male
fill([bins4plot bins4plot(end:-1:1)], [msMSHU.mean - msMSHU.std/2, msMSHU.mean(end:-1:1) + msMSHU.std(end:-1:1)/2], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, msMSHU.mean, 'b-o', 'LineWidth', 2);

linkaxes(hyx, 'xy'); ylim([0 1]);
    
    

fprintf('The mean and std for Male syllable duration is  %1.3f %1.3f \n', mean(Msyldur), std(Msyldur));
fprintf('The mean and std for Female syllable duration is  %1.3f %1.3f \n', mean(Fsyldur), std(Fsyldur));
fprintf('The mean and std for M2F ISI is  %1.3f %1.3f \n', mean(M2FISI), std(M2FISI));
fprintf('The mean and std for F2M ISI is  %1.3f %1.3f \n', mean(F2MISI), std(F2MISI));

%% Embedded Concatonation function
function out = concatHist(in, len)

    dat(1,:) = zeros(1,len);
    
    for qq = 1:length(in)
        
        for ww = 1:length(in(qq).spkcnt(1,:))
            if sum(in(qq).spkcnt(:,ww)) > 1
            dat(end+1,:) = in(qq).spkcnt(:,ww) / max(in(qq).spkcnt(:,ww));
            end
        end
    end
    
    out.mean = sum(dat) / (length(dat(:,1))-1); % Mean normalized data
    out.std = std(dat);
    out.N = length(dat(:,1));
    clear dat;
    
end


%% Embedded histogram function
function [spikearray, bintims] = wPhaseCut(spiketimes, tims, wid)

        % wid is window in msec before and after start of the syllable at the end of the focal ISI.
        numbins = 10; % How many bins before and after the onset of our focal syllable?
        binwid = wid / numbins; % Width of each bin
        
        % Specify the OVERLAP percentage here
        overlap = 50; % Overlap is 80% of previous window
        
        overlap = 1-(overlap/100); % Converts to step size for advancing the window
        
        bintims = -wid:binwid*overlap:wid; % List of bin times centered on zero
    
    for j = length(tims):-1:1 % For each syllable
        
        bintimstarts = tims(j)-wid:binwid*overlap:tims(j)+wid-(binwid*overlap); % Start and end times for current syllable
        spkcnts = zeros(1, length(-wid:binwid*overlap:wid-(binwid*overlap)));
                
        for k = 1:length(spiketimes) % For each row of spikes
            
            for m = 1:length(bintimstarts) % For each bin of our PSTH
                spkcnts(m) = spkcnts(m) + length(find(spiketimes{k} > bintimstarts(m) & spiketimes{k} < bintimstarts(m)+binwid));
            end
        end
         
        spikearray(:,j) = spkcnts;
 
    end

end


end



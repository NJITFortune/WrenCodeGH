function [M, F] = wTransitionsPlot(in, wwindow)
% Usage out = wTransitionsPlot(w, window)
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

[msolosyls, mduetsyls, fsolosyls, fduetsyls, Cspon, Aspon, birdlist] = wData;

analyzethis = [2, 4]; % The birdlist groups to analyze (2 is duets with female solo syllables and 4 with male solo syllables)

%% Cycle for each bird group for transition plots
for jj=1:length(analyzethis)

    % Some variables because I am not a talented coder.
    Fsyldur = []; Msyldur = [];  % A side calculation - durations of syllables
    M2FISI = []; F2MISI = [];    % Another side calculation - durations of InterSyllable Intervals
    
    clear MAHU;
    
    % More variables that we need
    MAHU(1).SPS = []; MAHU(1).RSraw = []; MAHU(1).RSnorm = [];

    MAHC = MAHU; MHAU = MAHU; MHAC = MAHU; 
    MSAU = MAHU; MSAC = MAHU; MSHU = MAHU; MSHC = MAHU;

    FAHU = MAHU; FAHC = MAHU; FHAU = MAHU; FHAC = MAHU; 
    FSAU = MAHU; FSAC = MAHU; FSHU = MAHU; FSHC = MAHU;
    
    
    % Cycle for each bird
    for ff = birdlist{analyzethis(jj)}
        
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
    
% Generate the transition histograms

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
            for kk = 1:length(tmp); MSAC(end+1) = tmp(kk); 
                % figure(27); axr(2) = subplot(122); title('Autogenous Solo'); hold on; plot(tmp(kk).SPS, 'b');
            end; clear tmp;
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
                    %figure(27); axr(1) = subplot(121); title('Female Autogenous Duet'); hold on; plot(tmp(kk).SPS);
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
                %figure(27); axr(2) = subplot(122); title('Autogenous Solo'); hold on; plot(tmp(kk).SPS, 'm');
                %linkaxes(axr, 'y');
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

    msFSAC = concatHist(FSAC);                        
    msFSAU = concatHist(FSAU);                        
    msFSHC = concatHist(FSHC);                        
    msFSHU = concatHist(FSHU);                        

    msMSAC = concatHist(MSAC);                        
    msMSAU = concatHist(MSAU);                        
    msMSHC = concatHist(MSHC);                        
    msMSHU = concatHist(MSHU);                        


typeofplot = 2; % 1 is SPS, 2 is rawRS, 3 is normRS

% figure(1+((jj-1)*10)); clf; set(gcf, 'Color', [1,1,1]); 
% 
% % PLOT Male Autogenous DATA 
% 
% axc(1) = subplot(221); hold on; title('Solo Male Chronic'); 
%     plot([0 0], [-10 100], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
%     % Female
%         plotmasteryoda(msFSHC, bins4plot, typeofplot, 2);
%     % Male
%         plotmasteryoda(msMSAC, bins4plot, typeofplot, 1);
% 
% axc(2) = subplot(222); hold on; title('F2M Chronic'); 
%     plot([0 0], [-10 100], 'k-', 'LineWidth', 2); %plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
%     % Female
%         plotmasteryoda(msFAHC, bins4plot, typeofplot, 2);
%     % Male
%         plotmasteryoda(msMHAC, bins4plot, typeofplot, 1);
%     
% % PLOT Female Autogenous DATA
% 
% axc(3) = subplot(223); hold on; title('Solo Female Chronic'); 
%     plot([0 0], [-10 100], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
%     % Male
%         plotmasteryoda(msMSHC, bins4plot, typeofplot, 1);
%     % Female
%         plotmasteryoda(msFSAC, bins4plot, typeofplot, 2);
% 
% axc(4) = subplot(224); hold on; title('M2F Chronic'); 
%     plot([0 0], [-10 100], 'k-', 'LineWidth', 2); %plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
%     % Male
%         plotmasteryoda(msMAHC, bins4plot, typeofplot, 1);
%     % Female
%         plotmasteryoda(msFHAC, bins4plot, typeofplot, 2);

        
% Custom figure
if analyzethis(jj) == 2 % Female data
    
%% AWAKE PLOTS 
figure(1); clf; 
    
axf(1) = subplot(222); hold on; title('Solo Female Chronic'); 
    plot([0 0], [-10 100], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    % Male
        plotmasteryoda(msMSHC, bins4plot, typeofplot, 1);
    % Female
        plotmasteryoda(msFSAC, bins4plot, typeofplot, 2);

axf(2) = subplot(224); hold on; title('M2F Chronic'); 
    plot([0 0], [-10 100], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    % Male
        plotmasteryoda(msMAHC, bins4plot, typeofplot, 1);
    % Female
        plotmasteryoda(msFHAC, bins4plot, typeofplot, 2);

linkaxes(axf, 'xy'); figure(1); subplot(222); ylim([-5 38]);  xlim([-widow, widow]);


%% URETHANE PLOTS (not used in manuscript)
figure(2); clf;
    
axxf(1) = subplot(222); hold on; title('Solo Female Acute'); 
    plot([0 0], [-10 100], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    % Male
        plotmasteryoda(msMSHU, bins4plot, typeofplot, 1);
    % Female
        plotmasteryoda(msFSAU, bins4plot, typeofplot, 2);

axxf(2) = subplot(224); hold on; title('M2F Acute'); 
    plot([0 0], [-10 100], 'k-', 'LineWidth', 2); %plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    % Male
        plotmasteryoda(msMAHU, bins4plot, typeofplot, 1);
    % Female
        plotmasteryoda(msFHAU, bins4plot, typeofplot, 2);

linkaxes(axxf, 'xy'); figure(2); subplot(222); ylim([-10 50]); xlim([-widow, widow]);

end

if analyzethis(jj) == 4 % Male data
    
figure(1);

axm(1) = subplot(221); hold on; title('Solo Male Chronic'); 
    plot([0 0], [-10 100], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    % Female
        plotmasteryoda(msFSHC, bins4plot, typeofplot, 2);
    % Male
        plotmasteryoda(msMSAC, bins4plot, typeofplot, 1);
        
axm(2) = subplot(223); hold on; title('F2M Chronic'); 
    plot([0 0], [-10 100], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    % Female
        plotmasteryoda(msFAHC, bins4plot, typeofplot, 2);
    % Male
        plotmasteryoda(msMHAC, bins4plot, typeofplot, 1);
        
linkaxes(axm, 'xy'); figure(1); subplot(223); ylim([-10 78]);  xlim([-widow, widow]);

figure(2);

axxm(1) = subplot(221); hold on; title('Solo Male Acute'); 
    plot([0 0], [-10 100], 'k-', 'LineWidth', 2); plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    % Female
        plotmasteryoda(msFSHU, bins4plot, typeofplot, 2);
    % Male
        plotmasteryoda(msMSAU, bins4plot, typeofplot, 1);
        
axxm(2) = subplot(223); hold on; title('F2M Acute'); 
    plot([0 0], [-10 100], 'k-', 'LineWidth', 2); %plot([-widow widow], [0 0], 'k-', 'LineWidth', 2);
    % Female
        plotmasteryoda(msFAHU, bins4plot, typeofplot, 2);
    % Male
        plotmasteryoda(msMHAU, bins4plot, typeofplot, 1);
        
linkaxes(axxm, 'xy'); figure(2); subplot(223); ylim([0 50]);  xlim([-widow, widow]);

end

% linkaxes(axc, 'xy'); figure(1+((jj-1)*10)); subplot(221); xlim([-widow-0.0001, widow+0.0001]);
%     
%     if analyzethis(jj)==2; ylim([0 40]); end
%     if analyzethis(jj)==3; ylim([0 40]); end
%     if analyzethis(jj)==4; ylim([0 75]); end
    


end % Cycling for each group














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

        if sponSPS == 0; sponSPS = 1; end
        
        for j=1:length(goodlist)
            out(j).SPS = spikearray(:, goodlist(j)) / binwid; 
            out(j).RSraw = (spikearray(:, goodlist(j)) / binwid) - sponSPS;
            out(j).RSnorm = ((spikearray(:, goodlist(j)) / binwid) - sponSPS) / sponSPS;
        end
        
        
end

%% Embedded plotting function
function plotmasteryoda(data, bins, typ, sx)

        if sx == 2 % This is a female
            col = [0.9, 0.7, 0.9];
            mrkr = 'm-';
        elseif sx == 1 % This is a male
            col = [0.6, 0.9, 0.9];
            mrkr = 'b-';
        end
        
        if typ == 1 % spikes per second
fill([bins bins(end:-1:1)], [data.meanSPS - data.steSPS, data.meanSPS(end:-1:1) + data.steSPS(end:-1:1)], col, 'LineStyle', 'none');
plot(bins, data.meanSPS, mrkr, 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 16);
        end
        if typ == 2 % Raw Response Strength
fill([bins bins(end:-1:1)], [data.meanRSraw - data.steRSraw, data.meanRSraw(end:-1:1) + data.steRSraw(end:-1:1)], col, 'LineStyle', 'none');
plot(bins, data.meanRSraw, mrkr, 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 16);
        end
        if typ == 3 % Normalized Response Strength
fill([bins bins(end:-1:1)], [data.meanRSnorm - data.steRSnorm, data.meanRSnorm(end:-1:1) + data.steRSnorm(end:-1:1)], col, 'LineStyle', 'none');
plot(bins, data.meanRSnorm, mrkr, 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 16);
        end

    end


end

function [M, F] = wTransitions(in, widow)
% Usage out = wTransitions(in, window)
% This generates plots that focus on the transitions between
% female and male syllables. 
% in is the data structure w
% window is the time before and after the start of the 2nd syllable in 
% pair for plotting.

% Default window width for histogram if user didn't specify window
if nargin < 2; widow = 0.250; end

figure(1); clf; subplot(121); hold on; title('M2F Chronic'); subplot(122); hold on; title('M2F Urethane'); 
figure(2); clf; subplot(121); hold on;  title('F2M Chronic'); subplot(122); hold on;  title('F2M Urethane'); 

figure(3); clf; subplot(121); hold on; title('M2F Chronic'); subplot(122); hold on; title('M2F Urethane'); 
figure(4); clf; subplot(121); hold on;  title('F2M Chronic'); subplot(122); hold on;  title('F2M Urethane'); 

%% List of Chronic singing data with syllable indices and locations for spontaneous activity

[msolosyls, mduetsyls, fsolosyls, fduetsyls, ~, ~] = wData;

% Choose which data to analyze
    birdlist = 14:-1:1; % All compleat data with both ACUTE (urethane) AND CHRONIC (awake)

    Fsyldur = []; Msyldur = [];
    M2FISI = []; F2MISI = [];
    Mwhichmalesolosyl = 0; Mwhichfemalesolosyl = 0;
    Fwhichmalesolosyl = 0; Fwhichfemalesolosyl = 0;
%% For each bird
for ff = birdlist
        
%% Get all male-female and female-male duet syllable transitions
    M2Fsyltim = []; F2Msyltim = []; Msolosyltims = []; Fsolosyltims = [];
    
    % Male to female transitions
    for t = 1:length(mduetsyls{ceil(ff/2)}) % For every male syllable
       if ~isempty(find(fduetsyls{ceil(ff/2)} == mduetsyls{ceil(ff/2)}(t)+1, 1)) % If the next syllable is female
           currentFemaleIndex = mduetsyls{ceil(ff/2)}(t)+1;
           M2Fsyltim(end+1) = in(ff).syl(currentFemaleIndex).tim(1); % Time of the start of the female syllable
           Fsyldur(end+1) = in(ff).syl(currentFemaleIndex).tim(2) - in(ff).syl(currentFemaleIndex).tim(1);
           M2FISI(end+1) = in(ff).syl(currentFemaleIndex).tim(1) - in(ff).syl(currentFemaleIndex-1).tim(2);
       end
    end
    % Female to male transitions
    for t = 1:length(fduetsyls{ceil(ff/2)}) % For every female syllable
       if ~isempty(find(mduetsyls{ceil(ff/2)} == fduetsyls{ceil(ff/2)}(t)+1, 1)) % If the next syllable is male
           currentMaleIndex = fduetsyls{ceil(ff/2)}(t)+1;
           F2Msyltim(end+1) = in(ff).syl(currentMaleIndex).tim(1); % Time of the start of the male syllable
           Msyldur(end+1) = in(ff).syl(currentMaleIndex).tim(2) -in(ff).syl(currentMaleIndex).tim(1);
           F2MISI(end+1) = in(ff).syl(currentMaleIndex).tim(1) - in(ff).syl(currentMaleIndex-1).tim(2);           
       end
    end

    % All Solo Male syllables
    for t = 1:length(msolosyls{ceil(ff/2)}) % For every male syllable
         Msolosyltims(end+1) = in(ff).syl(msolosyls{ceil(ff/2)}(t)).tim(1); % Time of the start of the male syllable
    end
    
    % All Solo Female syllables
    for t = 1:length(fsolosyls{ceil(ff/2)}) % For every female syllable
         Fsolosyltims(end+1) = in(ff).syl(fsolosyls{ceil(ff/2)}(t)).tim(1); % Time of the start of the female syllable
    end
    
%% Fetch the histograms

if in(ff).sexy == 1 % This is a male
    [M(ff).AHU.spkcnt, M(ff).bintims] = wPhaseCut(in(ff).Aspikes, M2Fsyltim, widow);
    [M(ff).AHC.spkcnt, ~] = wPhaseCut(in(ff).Cspikes, M2Fsyltim, widow);
    [M(ff).HAU.spkcnt, ~] = wPhaseCut(in(ff).Aspikes, F2Msyltim, widow);
    [M(ff).HAC.spkcnt, ~] = wPhaseCut(in(ff).Cspikes, F2Msyltim, widow);
    
    bins4plot = (M(ff).bintims(2:end) + M(ff).bintims(1:end-1)) /2;
    
    if ~isempty(msolosyls{ceil(ff/2)})
        Mwhichmalesolosyl = Mwhichmalesolosyl +1;
        [MSAU(:,Mwhichmalesolosyl), ~] = wPhaseCut(in(ff).Aspikes, Msolosyltims, widow);
        [MSAC(:,Mwhichmalesolosyl), ~] = wPhaseCut(in(ff).Cspikes, Msolosyltims, widow);
    end
    if ~isempty(fsolosyls{ceil(ff/2)})
        Mwhichfemalesolosyl = Mwhichfemalesolosyl +1;
        [MSHU(:,Mwhichfemalesolosyl), ~] = wPhaseCut(in(ff).Aspikes, Fsolosyltims, widow);
        [MSHC(:,Mwhichfemalesolosyl), ~] = wPhaseCut(in(ff).Cspikes, Fsolosyltims, widow);
    end

% RAW    
%     figure(1); subplot(121); plot(M(ff).bintims(1:end-1), M(ff).AHC.spkcnt, 'b');
%     figure(1); subplot(122); plot(M(ff).bintims(1:end-1), M(ff).AHU.spkcnt, 'b');
%     figure(2); subplot(121); plot(M(ff).bintims(1:end-1), M(ff).HAC.spkcnt, 'b');
%     figure(2); subplot(122); plot(M(ff).bintims(1:end-1), M(ff).HAU.spkcnt, 'b');
% NORMALIZED    
    figure(1); subplot(121); plot(bins4plot, M(ff).AHC.spkcnt/max(M(ff).AHC.spkcnt), 'b*-');
    figure(1); subplot(122); plot(bins4plot, M(ff).AHU.spkcnt/max(M(ff).AHU.spkcnt), 'b*-');
    figure(2); subplot(121); plot(bins4plot, M(ff).HAC.spkcnt/max(M(ff).HAC.spkcnt), 'b*-');
    figure(2); subplot(122); plot(bins4plot, M(ff).HAU.spkcnt/max(M(ff).HAU.spkcnt), 'b*-');
end

if in(ff).sexy == 2 % This is a female
    [F(ff).AHU.spkcnt, F(ff).bintims] = wPhaseCut(in(ff).Aspikes, F2Msyltim, widow);
    [F(ff).AHC.spkcnt, ~] = wPhaseCut(in(ff).Cspikes, F2Msyltim, widow);
    [F(ff).HAU.spkcnt, ~] = wPhaseCut(in(ff).Aspikes, M2Fsyltim, widow);
    [F(ff).HAC.spkcnt, ~] = wPhaseCut(in(ff).Cspikes, M2Fsyltim, widow);

    bins4plot = (F(ff).bintims(2:end) + F(ff).bintims(1:end-1)) /2;

    if ~isempty(msolosyls{ceil(ff/2)})
        Fwhichmalesolosyl = Fwhichmalesolosyl +1;
        [FSHU(:,Fwhichmalesolosyl), ~] = wPhaseCut(in(ff).Aspikes, Msolosyltims, widow);
        [FSHC(:,Fwhichmalesolosyl), ~] = wPhaseCut(in(ff).Cspikes, Msolosyltims, widow);
    end
    if ~isempty(fsolosyls{ceil(ff/2)})
        Fwhichfemalesolosyl = Fwhichfemalesolosyl +1;
        [FSAU(:,Fwhichfemalesolosyl), ~] = wPhaseCut(in(ff).Aspikes, Fsolosyltims, widow);
        [FSAC(:,Fwhichfemalesolosyl), ~] = wPhaseCut(in(ff).Cspikes, Fsolosyltims, widow);
    end
    
% RAW
%     figure(1); subplot(121); plot(F(ff).bintims(1:end-1), F(ff).HAC.spkcnt, 'm');
%     figure(1); subplot(122); plot(F(ff).bintims(1:end-1), F(ff).HAU.spkcnt, 'm');
%     figure(2); subplot(121); plot(F(ff).bintims(1:end-1), F(ff).AHC.spkcnt, 'm');
%     figure(2); subplot(122); plot(F(ff).bintims(1:end-1), F(ff).AHU.spkcnt, 'm');
% NORMALIZED    
    figure(1); subplot(121); plot(bins4plot, F(ff).HAC.spkcnt/max(F(ff).HAC.spkcnt), 'm*-');
    figure(1); subplot(122); plot(bins4plot, F(ff).HAU.spkcnt/max(F(ff).HAU.spkcnt), 'm*-');
    figure(2); subplot(121); plot(bins4plot, F(ff).AHC.spkcnt/max(F(ff).AHC.spkcnt), 'm*-');
    figure(2); subplot(122); plot(bins4plot, F(ff).AHU.spkcnt/max(F(ff).AHU.spkcnt), 'm*-');
end

end % Cycle for every bird

%% Build the Cool fill plots

for kk = length(F):-2:2
    FHAC(:,kk) = F(kk).HAC.spkcnt/max(F(kk).HAC.spkcnt);
    FHAU(:,kk) = F(kk).HAU.spkcnt/max(F(kk).HAU.spkcnt);
    FAHC(:,kk) = F(kk).AHC.spkcnt/max(F(kk).AHC.spkcnt);
    FAHU(:,kk) = F(kk).AHU.spkcnt/max(F(kk).AHU.spkcnt);
end
for kk = length(M):-2:1
    MHAC(:,kk) = M(kk).HAC.spkcnt/max(M(kk).HAC.spkcnt);
    MHAU(:,kk) = M(kk).HAU.spkcnt/max(M(kk).HAU.spkcnt);
    MAHC(:,kk) = M(kk).AHC.spkcnt/max(M(kk).AHC.spkcnt);
    MAHU(:,kk) = M(kk).AHU.spkcnt/max(M(kk).AHU.spkcnt);
end

for jj = length(FHAC(:,1)):-1:1
    mFHAC(jj) = mean(FHAC(jj,:)); vFHAC(jj) = std(FHAC(jj,:));
    mFHAU(jj) = mean(FHAU(jj,:)); vFHAU(jj) = std(FHAU(jj,:));
    mFAHC(jj) = mean(FAHC(jj,:)); vFAHC(jj) = std(FAHC(jj,:));
    mFAHU(jj) = mean(FAHU(jj,:)); vFAHU(jj) = std(FAHU(jj,:));
end    
for jj = length(MHAC(:,1)):-1:1
    mMHAC(jj) = mean(MHAC(jj,:)); vMHAC(jj) = std(MHAC(jj,:));
    mMHAU(jj) = mean(MHAU(jj,:)); vMHAU(jj) = std(MHAU(jj,:));
    mMAHC(jj) = mean(MAHC(jj,:)); vMAHC(jj) = std(MAHC(jj,:));
    mMAHU(jj) = mean(MAHU(jj,:)); vMAHU(jj) = std(MAHU(jj,:));
end

%% Do the solo syllable stuff

% NORMALIZE
for jj = 1:length(FSAU(1,:))
    FSAU(:,jj) = FSAU(:,jj) / max(FSAU(:,jj));
end
for jj = 1:length(FSAC(1,:))
    FSAC(:,jj) = FSAC(:,jj) / max(FSAC(:,jj));
end
for jj = 1:length(FSHU(1,:))
    FSHU(:,jj) = FSHU(:,jj) / max(FSHU(:,jj));
end
for jj = 1:length(FSHC(1,:))
    FSHC(:,jj) = FSHC(:,jj) / max(FSHC(:,jj));
end
for jj = 1:length(FSAU(1,:))
    FSAU(:,jj) = FSAU(:,jj) / max(FSAU(:,jj));
end

for jj = 1:length(MSAU(1,:))
    MSAU(:,jj) = MSAU(:,jj) / max(MSAU(:,jj));
end
for jj = 1:length(MSAC(1,:))
    MSAC(:,jj) = MSAC(:,jj) / max(MSAC(:,jj));
end
for jj = 1:length(MSHU(1,:))
    MSHU(:,jj) = MSHU(:,jj) / max(MSHU(:,jj));
end
for jj = 1:length(MSHC(1,:))
    MSHC(:,jj) = MSHC(:,jj) / max(MSHC(:,jj));
end

% AVERAGE
for jj = 1:length(FSAU(:,1))
    mFSAU(jj) = mean(FSAU(jj,:)); vFSAU(jj) = std(FSAU(jj,:));
end
for jj = 1:length(FSHU(:,1))
    mFSHU(jj) = mean(FSHU(jj,:)); vFSHU(jj) = std(FSHU(jj,:));
end
for jj = 1:length(FSAC(:,1))
    mFSAC(jj) = mean(FSAC(jj,:)); vFSAC(jj) = std(FSAC(jj,:));
end
for jj = 1:length(FSHC(:,1))
    mFSHC(jj) = mean(FSHC(jj,:)); vFSHC(jj) = std(FSHC(jj,:));
end
for jj = 1:length(MSAU(:,1))
    mMSAU(jj) = mean(MSAU(jj,:)); vMSAU(jj) = std(MSAU(jj,:));
end
for jj = 1:length(MSHU(:,1))
    mMSHU(jj) = mean(MSHU(jj,:)); vMSHU(jj) = std(MSHU(jj,:));
end
for jj = 1:length(MSAC(:,1))
    mMSAC(jj) = mean(MSAC(jj,:)); vMSAC(jj) = std(MSAC(jj,:));
end
for jj = 1:length(MSHC(:,1))
    mMSHC(jj) = mean(MSHC(jj,:)); vMSHC(jj) = std(MSHC(jj,:));
end

figure(127); clf; 
    subplot(121); hold on; title('FemaleSolo Chronic'); plot(mFSAC, 'm'); plot(mMSHC, 'b') 
    subplot(122); hold on; title('FemaleSolo Urethane'); plot(mFSAU, 'm'); plot(mMSHU, 'b') 
figure(128); clf; 
    subplot(121); hold on; title('MaleSolo Chronic'); plot(mFSHC, 'm'); plot(mMSAC, 'b') 
    subplot(122); hold on; title('MaleSolo Urethane'); plot(mFSHU, 'm'); plot(mMSAU, 'b') 

vMAHC = vMAHC/2; vMAHU = vMAHU/2; vMHAC = vMHAC/2; vMHAU = vMHAU/2;
vFAHC = vFAHC/2; vFAHU = vFAHU/2; vFHAC = vFHAC/2; vFHAU = vFHAU/2;

figure(3); subplot(121); plot([0 0], [0 0.8], 'k-', 'LineWidth', 2);
fill([bins4plot bins4plot(end:-1:1)], [mMAHC - vMAHC, mMAHC(end:-1:1) + vMAHC(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, mMAHC, 'b-o', 'LineWidth', 2);
figure(3); subplot(121);
fill([bins4plot bins4plot(end:-1:1)], [mFHAC - vFHAC, mFHAC(end:-1:1) + vFHAC(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, mFHAC, 'm-o', 'LineWidth', 2);

figure(3); subplot(122); plot([0 0], [0 0.8], 'k-', 'LineWidth', 2);
fill([bins4plot bins4plot(end:-1:1)], [mMAHU - vMAHU, mMAHU(end:-1:1) + vMAHU(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, mMAHU, 'b-o', 'LineWidth', 2);
figure(3); subplot(122);
fill([bins4plot bins4plot(end:-1:1)], [mFHAU - vFHAU, mFHAU(end:-1:1) + vFHAU(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, mFHAU, 'm-o', 'LineWidth', 2);


figure(4); subplot(121); plot([0 0], [0 0.8], 'k-', 'LineWidth', 2);
fill([bins4plot bins4plot(end:-1:1)], [mMHAC - vMHAC, mMHAC(end:-1:1) + vMHAC(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, mMHAC, 'b-o', 'LineWidth', 2);
figure(4); subplot(121);
fill([bins4plot bins4plot(end:-1:1)], [mFAHC - vFAHC, mFAHC(end:-1:1) + vFAHC(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, mFAHC, 'm-o', 'LineWidth', 2);

figure(4); subplot(122);
fill([bins4plot bins4plot(end:-1:1)], [mFAHU - vFAHU, mFAHU(end:-1:1) + vFAHU(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, mFAHU, 'm-o', 'LineWidth', 2);
figure(4); subplot(122); plot([0 0], [0 0.8], 'k-', 'LineWidth', 2);
fill([bins4plot bins4plot(end:-1:1)], [mMHAU - vMHAU, mMHAU(end:-1:1) + vMHAU(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, mMHAU, 'b-o', 'LineWidth', 2);

fprintf('The mean and std for Male syllable duration is  %1.3f %1.3f \n', mean(Msyldur), std(Msyldur));
fprintf('The mean and std for Female syllable duration is  %1.3f %1.3f \n', mean(Fsyldur), std(Fsyldur));
fprintf('The mean and std for M2F ISI is  %1.3f %1.3f \n', mean(M2FISI), std(M2FISI));
fprintf('The mean and std for F2M ISI is  %1.3f %1.3f \n', mean(F2MISI), std(F2MISI));




%% Embedded histogram function
function [spkcnts, bintims] = wPhaseCut(spiketimes, tims, wid)

        % wid is window in msec before and after start of the syllable at the end of the focal ISI.
        numbins = 10; % How many bins before and after the onset of our focal syllable?
        binwid = wid / numbins; % Width of each bin
        
        % Specify the OVERLAP percentage here
        overlap = 75; % Overlap is 80% of previous window
        
        overlap = 1-(overlap/100); % Converts to step size for advancing the window
        
        bintims = -wid:binwid*overlap:wid; % List of bin times centered on zero

        spkcnts = zeros(1, length(-wid:binwid*overlap:wid-(binwid*overlap)));

    
    for j = length(tims):-1:1 % For each syllable
        
        bintimstarts = tims(j)-wid:binwid*overlap:tims(j)+wid-(binwid*overlap); % Start and end times for current syllable
                
        for k = 1:length(spiketimes) % For each row of spikes
            
            for m = 1:length(bintimstarts) % For each bin of our PSTH
                spkcnts(m) = spkcnts(m) + length(find(spiketimes{k} > bintimstarts(m) & spiketimes{k} < bintimstarts(m)+binwid));
            end
         end
 
    end

end


end



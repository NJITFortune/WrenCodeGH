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

[~, mduetsyls, ~, fduetsyls, ~, ~] = wData;

% Choose which data to analyze
    birdlist = 12:-1:1; % All compleat data

    Fsyldur = []; Msyldur = [];
    M2FISI = []; F2MISI = [];

%% For each bird
for ff = birdlist
        
%% Get all male-female and female-male duet syllable transitions
    M2Fsyltim = []; F2Msyltim = [];
    
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

%% Fetch the histograms

if in(ff).sexy == 1 % This is a male
    [M(ff).AHU.spkcnt, M(ff).bintims] = wPhaseCut(in(ff).Aspikes, M2Fsyltim, widow);
    [M(ff).AHC.spkcnt, ~] = wPhaseCut(in(ff).Cspikes, M2Fsyltim, widow);
    [M(ff).HAU.spkcnt, ~] = wPhaseCut(in(ff).Aspikes, F2Msyltim, widow);
    [M(ff).HAC.spkcnt, ~] = wPhaseCut(in(ff).Cspikes, F2Msyltim, widow);
    
    bins4plot = (M(ff).bintims(2:end) + M(ff).bintims(1:end-1)) /2;

% RAW    
    figure(1); subplot(121); plot(M(ff).bintims(1:end-1), M(ff).AHC.spkcnt, 'b');
    figure(1); subplot(122); plot(M(ff).bintims(1:end-1), M(ff).AHU.spkcnt, 'b');
    figure(2); subplot(121); plot(M(ff).bintims(1:end-1), M(ff).HAC.spkcnt, 'b');
    figure(2); subplot(122); plot(M(ff).bintims(1:end-1), M(ff).HAU.spkcnt, 'b');
% NORMALIZED    
%     figure(1); subplot(121); plot(bins4plot, M(ff).AHC.spkcnt/max(M(ff).AHC.spkcnt), 'b*-');
%     figure(1); subplot(122); plot(bins4plot, M(ff).AHU.spkcnt/max(M(ff).AHU.spkcnt), 'b*-');
%     figure(2); subplot(121); plot(bins4plot, M(ff).HAC.spkcnt/max(M(ff).HAC.spkcnt), 'b*-');
%     figure(2); subplot(122); plot(bins4plot, M(ff).HAU.spkcnt/max(M(ff).HAU.spkcnt), 'b*-');
end

if in(ff).sexy == 2 % This is a female
    [F(ff).AHU.spkcnt, F(ff).bintims] = wPhaseCut(in(ff).Aspikes, F2Msyltim, widow);
    [F(ff).AHC.spkcnt, ~] = wPhaseCut(in(ff).Cspikes, F2Msyltim, widow);
    [F(ff).HAU.spkcnt, ~] = wPhaseCut(in(ff).Aspikes, M2Fsyltim, widow);
    [F(ff).HAC.spkcnt, ~] = wPhaseCut(in(ff).Cspikes, M2Fsyltim, widow);

    bins4plot = (F(ff).bintims(2:end) + F(ff).bintims(1:end-1)) /2;

% RAW    
    figure(1); subplot(121); plot(F(ff).bintims(1:end-1), F(ff).HAC.spkcnt, 'm');
    figure(1); subplot(122); plot(F(ff).bintims(1:end-1), F(ff).HAU.spkcnt, 'm');
    figure(2); subplot(121); plot(F(ff).bintims(1:end-1), F(ff).AHC.spkcnt, 'm');
    figure(2); subplot(122); plot(F(ff).bintims(1:end-1), F(ff).AHU.spkcnt, 'm');
% NORMALIZED    
%     figure(1); subplot(121); plot(bins4plot, F(ff).HAC.spkcnt/max(F(ff).HAC.spkcnt), 'm*-');
%     figure(1); subplot(122); plot(bins4plot, F(ff).HAU.spkcnt/max(F(ff).HAU.spkcnt), 'm*-');
%     figure(2); subplot(121); plot(bins4plot, F(ff).AHC.spkcnt/max(F(ff).AHC.spkcnt), 'm*-');
%     figure(2); subplot(122); plot(bins4plot, F(ff).AHU.spkcnt/max(F(ff).AHU.spkcnt), 'm*-');
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

vMAHC = vMAHC/2; vMAHU = vMAHU/2; vMHAC = vMHAC/2; vMHAU = vMHAU/2;
vFAHC = vFAHC/2; vFAHU = vFAHU/2; vFHAC = vFHAC/2; vFHAU = vFHAU/2;

figure(3); subplot(121); plot([0 0], [0 0.8], 'k-', 'LineWidth', 2);
fill([bins4plot bins4plot(end:-1:1)], [mMAHC - vMAHC, mMAHC(end:-1:1) + vMAHC(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, mMAHC, 'b', 'LineWidth', 2);
figure(3); subplot(121);
fill([bins4plot bins4plot(end:-1:1)], [mFHAC - vFHAC, mFHAC(end:-1:1) + vFHAC(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, mFHAC, 'm', 'LineWidth', 2);

figure(3); subplot(122); plot([0 0], [0 0.8], 'k-', 'LineWidth', 2);
fill([bins4plot bins4plot(end:-1:1)], [mMAHU - vMAHU, mMAHU(end:-1:1) + vMAHU(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, mMAHU, 'b', 'LineWidth', 2);
figure(3); subplot(122);
fill([bins4plot bins4plot(end:-1:1)], [mFHAU - vFHAU, mFHAU(end:-1:1) + vFHAU(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, mFHAU, 'm', 'LineWidth', 2);


figure(4); subplot(121); plot([0 0], [0 0.8], 'k-', 'LineWidth', 2);
fill([bins4plot bins4plot(end:-1:1)], [mMHAC - vMHAC, mMHAC(end:-1:1) + vMHAC(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, mMHAC, 'b', 'LineWidth', 2);
figure(4); subplot(121);
fill([bins4plot bins4plot(end:-1:1)], [mFAHC - vFAHC, mFAHC(end:-1:1) + vFAHC(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, mFAHC, 'm', 'LineWidth', 2);

figure(4); subplot(122);
fill([bins4plot bins4plot(end:-1:1)], [mFAHU - vFAHU, mFAHU(end:-1:1) + vFAHU(end:-1:1)], [0.9, 0.7, 0.9], 'LineStyle', 'none');
plot(bins4plot, mFAHU, 'm', 'LineWidth', 2);
figure(4); subplot(122); plot([0 0], [0 0.8], 'k-', 'LineWidth', 2);
fill([bins4plot bins4plot(end:-1:1)], [mMHAU - vMHAU, mMHAU(end:-1:1) + vMHAU(end:-1:1)], [0.6, 0.9, 0.9], 'LineStyle', 'none');
plot(bins4plot, mMHAU, 'b', 'LineWidth', 2);

fprintf('The mean and std for Male syllable duration is  %1.3f %1.3f \n', mean(Msyldur), std(Msyldur));
fprintf('The mean and std for Female syllable duration is  %1.3f %1.3f \n', mean(Fsyldur), std(Fsyldur));
fprintf('The mean and std for M2F ISI is  %1.3f %1.3f \n', mean(M2FISI), std(M2FISI));
fprintf('The mean and std for F2M ISI is  %1.3f %1.3f \n', mean(F2MISI), std(F2MISI));




%% Embedded histogram function
function [spkcnts, bintims] = wPhaseCut(spiketimes, tims, wid)

        % wid is window in msec before and after the middle of the ISI.
        numbins = 5; % How many bins before and after the onset of our focal syllable?
        binwid = wid / numbins; % Width of each bin
        
        bintims = -wid:binwid:wid; % List of bin times centered on zero
        
                spkcnts = zeros(1,numbins*2);

    
    for j = length(tims):-1:1 % For each syllable
        
        realbintims = tims(j)-wid:binwid:tims(j)+wid; % Start and end times for current syllable
                
        for k = 1:length(spiketimes) % For each row of spikes
            
            for m = 1:length(realbintims)-1 % For each bin of our PSTH
                spkcnts(m) = spkcnts(m) + length(find(spiketimes{k} > realbintims(m) & spiketimes{k} < realbintims(m+1)));
            end
         end
 
    end

end


end



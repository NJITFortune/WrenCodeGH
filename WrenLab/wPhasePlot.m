function [out, stts] = wPhasePlot(in, pad)
% Usage out = wPhasePlot(in)
% This makes a phase plot and calculates vector strength for AH and HA
% intervals. 

if nargin < 2; pad = 0; end

%% List of Chronic singing data with syllable indices and locations for spontaneous activity

[msolosyls, mduetsyls, fsolosyls, fduetsyls, Cspon, Aspon] = wData;

% User can choose which data to analyze
% birdlist = [3 4 9 10 11 12 13 14];
birdlist = 1:12; % All compleat data


%% For each bird
for ff = birdlist
        
% List of all duet syllables for this bird
    allduetsyls = sort([mduetsyls{ceil(ff/2)}, fduetsyls{ceil(ff/2)}]);

%% Get all male-female and female-male duet syllable transitions
    goodmalesyls = []; goodfemsyls = [];

    % Male to female transitions
    for t = 1:length(mduetsyls{ceil(ff/2)}) 
       if ~isempty(find(allduetsyls == mduetsyls{ceil(ff/2)}(t)+1, 1))
           goodmalesyls(end+1) = mduetsyls{ceil(ff/2)}(t);
       end
    end
    % Female to male transitions
    for t = 1:length(fduetsyls{ceil(ff/2)})
       if ~isempty(find(allduetsyls == fduetsyls{ceil(ff/2)}(t)+1, 1))
           goodfemsyls(end+1) = fduetsyls{ceil(ff/2)}(t);
       end
    end

%% Get all of the syllable times (previous start, middle of interval, next end)
for zz = length(goodmalesyls):-1:1   
    Msyltims(zz,1) = in(ff).syl(goodmalesyls(zz)).tim(1)-pad;
    Msyltims(zz,2) = mean([in(ff).syl(goodmalesyls(zz)+1).tim(1), in(ff).syl(goodmalesyls(zz)).tim(2) ]);
    Msyltims(zz,3) = in(ff).syl(goodmalesyls(zz)+1).tim(2)+pad;
end
for zz = length(goodfemsyls):-1:1   
    Fsyltims(zz,1) = in(ff).syl(goodfemsyls(zz)).tim(1)-pad;
    Fsyltims(zz,2) = mean([in(ff).syl(goodfemsyls(zz)+1).tim(1), in(ff).syl(goodfemsyls(zz)).tim(2) ]);
    Fsyltims(zz,3) = in(ff).syl(goodfemsyls(zz)+1).tim(2)+pad;
end

%% Fetch the vector strengths and plot

Rplotwid = 0.500; Rplotwind = 10; % For raw data
plotwid = pi*2; plotwind = 200; % For Normalized data

if in(ff).sexy == 1 % This is a male
    [out(ff).AH.Uvs, out(ff).AH.uPhaseSpikes, out(ff).AH.uSpikes] = wPhasor(in(ff).Aspikes, Msyltims);
    [out(ff).AH.Cvs, out(ff).AH.cPhaseSpikes, out(ff).AH.cSpikes] = wPhasor(in(ff).Cspikes, Msyltims);
    [out(ff).HA.Uvs, out(ff).HA.uPhaseSpikes, out(ff).HA.uSpikes] = wPhasor(in(ff).Aspikes, Fsyltims);
    [out(ff).HA.Cvs, out(ff).HA.cPhaseSpikes, out(ff).HA.cSpikes] = wPhasor(in(ff).Cspikes, Fsyltims);
end

if in(ff).sexy == 2 % This is a female
    [out(ff).AH.Uvs, out(ff).AH.uPhaseSpikes, out(ff).AH.uSpikes] = wPhasor(in(ff).Aspikes, Fsyltims);
    [out(ff).AH.Cvs, out(ff).AH.cPhaseSpikes, out(ff).AH.cSpikes] = wPhasor(in(ff).Cspikes, Fsyltims);
    [out(ff).HA.Uvs, out(ff).HA.uPhaseSpikes, out(ff).HA.uSpikes] = wPhasor(in(ff).Aspikes, Msyltims);
    [out(ff).HA.Cvs, out(ff).HA.cPhaseSpikes, out(ff).HA.cSpikes] = wPhasor(in(ff).Cspikes, Msyltims);    
end

end % Cycle for every bird

%% Now that we have all birds, compare vector strength values for males and females and do stats

for ff = birdlist
    % Which urethane vector strength is better?    
    out(ff).maxUV = max([out(ff).HA.Uvs out(ff).AH.Uvs]);     
end

[stts.H, stts.P, stts.CI] = ttest2([out(1:2:end).maxUV], [out(2:2:end).maxUV]);

pltplt = [0 0 0]; % No plotting of sliding window PSTHs
% pltplt = [1 6 2]; % Blue for male, black for urethane, magenta for female

%% Use this to plot each bird and calculate the mean
for ff = birdlist
    if in(ff).sexy == 1 % This is a male
%        figure(ff); clf;
%        subplot(411); 
            mAHcP(ff) = bs_swPSTH(out(ff).AH.cPhaseSpikes,[0, plotwid], plotwind, pltplt(1)); 
%        subplot(412); 
            mAHuP(ff) = bs_swPSTH(out(ff).AH.uPhaseSpikes,[0, plotwid], plotwind, pltplt(2));
%        subplot(413); 
            mHAcP(ff) = bs_swPSTH(out(ff).HA.cPhaseSpikes,[0, plotwid], plotwind, pltplt(1));
%        subplot(414); 
            mHAuP(ff) = bs_swPSTH(out(ff).HA.uPhaseSpikes,[0, plotwid], plotwind, pltplt(2));
        
%        figure(ff+100); clf;
%        subplot(411); 
            mAHcR(ff) = bs_swPSTH(out(ff).AH.cSpikes,[0, Rplotwid], Rplotwind, pltplt(1)); 
%        subplot(412); 
            mAHuR(ff) = bs_swPSTH(out(ff).AH.uSpikes,[0, Rplotwid], Rplotwind, pltplt(2));
%        subplot(413); 
            mHAcR(ff) = bs_swPSTH(out(ff).HA.cSpikes,[0, Rplotwid], Rplotwind, pltplt(1));
%        subplot(414); 
            mHAuR(ff) = bs_swPSTH(out(ff).HA.uSpikes,[0, Rplotwid], Rplotwind, pltplt(2));
    end

    if in(ff).sexy == 2 % This is a female
%    figure(ff); clf;
%        subplot(411); 
            fAHcP(ff) = bs_swPSTH(out(ff).AH.cPhaseSpikes,[0, plotwid], plotwind, pltplt(3)); 
%        subplot(412); 
            fAHuP(ff) = bs_swPSTH(out(ff).AH.uPhaseSpikes,[0, plotwid], plotwind, pltplt(2));
%        subplot(413); 
            fHAcP(ff) = bs_swPSTH(out(ff).HA.cPhaseSpikes,[0, plotwid], plotwind, pltplt(3));
%        subplot(414); 
            fHAuP(ff) = bs_swPSTH(out(ff).HA.uPhaseSpikes,[0, plotwid], plotwind, pltplt(2));
        
%        figure(ff+100); clf;
%        subplot(411); 
            fAHcR(ff) = bs_swPSTH(out(ff).AH.cSpikes,[0, Rplotwid], Rplotwind, pltplt(3)); 
%        subplot(412); 
            fAHuR(ff) = bs_swPSTH(out(ff).AH.uSpikes,[0, Rplotwid], Rplotwind, pltplt(2));
%        subplot(413); 
            fHAcR(ff) = bs_swPSTH(out(ff).HA.cSpikes,[0, Rplotwid], Rplotwind, pltplt(3));
%        subplot(414); 
            fHAuR(ff) = bs_swPSTH(out(ff).HA.uSpikes,[0, Rplotwid], Rplotwind, pltplt(2));
    end
end

% Initialize males
ChronicMaleRealHA = mHAcR(1).spers / max(mHAcR(1).spers);
ChronicMalePhaseHA = mHAcP(1).spers / max(mHAcP(1).spers);
ChronicMaleRealAH = mAHcR(1).spers / max(mAHcR(1).spers);
ChronicMalePhaseAH = mAHcP(1).spers / max(mAHcP(1).spers);

UrethaneMaleRealHA = mHAuR(1).spers / max(mHAuR(1).spers);
UrethaneMalePhaseHA = mHAuP(1).spers / max(mHAuP(1).spers);
UrethaneMaleRealAH = mAHuR(1).spers / max(mAHuR(1).spers);
UrethaneMalePhaseAH = mAHuP(1).spers / max(mAHuP(1).spers);

ChronicFemaleRealHA = fHAcR(2).spers / max(fHAcR(2).spers);
ChronicFemalePhaseHA = fHAcP(2).spers / max(fHAcP(2).spers);
ChronicFemaleRealAH = fAHcR(2).spers / max(fAHcR(2).spers);
ChronicFemalePhaseAH = fAHcP(2).spers / max(fAHcP(2).spers);

UrethaneFemaleRealHA = fHAuR(2).spers / max(fHAuR(2).spers);
UrethaneFemalePhaseHA = fHAuP(2).spers / max(fHAuP(2).spers);
UrethaneFemaleRealAH = fAHuR(2).spers / max(fAHuR(2).spers);
UrethaneFemalePhaseAH = fAHuP(2).spers / max(fAHuP(2).spers);

for tt = 3:2:length(birdlist) % The remaining male data
    ChronicMalePhaseAH = ChronicMalePhaseAH + mAHcP(tt).spers / max(mAHcP(tt).spers);
    ChronicMaleRealAH = ChronicMaleRealAH + mAHcR(tt).spers / max(mAHcR(tt).spers);
    ChronicMalePhaseHA = ChronicMalePhaseHA + mHAcP(tt).spers / max(mHAcP(tt).spers);
    ChronicMaleRealHA = ChronicMaleRealHA + mHAcR(tt).spers / max(mHAcR(tt).spers);
    
    UrethaneMalePhaseAH = UrethaneMalePhaseAH + mAHuP(tt).spers / max(mAHuP(tt).spers);
    UrethaneMaleRealAH = UrethaneMaleRealAH + mAHuR(tt).spers / max(mAHuR(tt).spers);
    UrethaneMalePhaseHA = UrethaneMalePhaseHA + mHAuP(tt).spers / max(mHAuP(tt).spers);
    UrethaneMaleRealHA = UrethaneMaleRealHA + mHAuR(tt).spers / max(mHAuR(tt).spers);
end

for tt = 4:2:length(birdlist) % The remaining female data
    ChronicFemalePhaseAH = ChronicFemalePhaseAH + fAHcP(tt).spers / max(fAHcP(tt).spers);
    ChronicFemaleRealAH = ChronicFemaleRealAH + fAHcR(tt).spers / max(fAHcR(tt).spers);   
    ChronicFemalePhaseHA = ChronicFemalePhaseHA + fHAcP(tt).spers / max(fHAcP(tt).spers);
    ChronicFemaleRealHA = ChronicFemaleRealHA + fHAcR(tt).spers / max(fHAcR(tt).spers);    

    UrethaneFemalePhaseAH = UrethaneFemalePhaseAH + fAHuP(tt).spers / max(fAHuP(tt).spers);
    UrethaneFemaleRealAH = UrethaneFemaleRealAH + fAHuR(tt).spers / max(fAHuR(tt).spers);   
    UrethaneFemalePhaseHA = UrethaneFemalePhaseHA + fHAuP(tt).spers / max(fHAuP(tt).spers);
    UrethaneFemaleRealHA = UrethaneFemaleRealHA + fHAuR(tt).spers / max(fHAuR(tt).spers);    
end

figure(27); clf; % Female followed by Male PHASE
subplot(2,1,1); plot(ChronicMalePhaseHA, 'b'); hold on;
subplot(2,1,1); plot(ChronicFemalePhaseAH, 'm');
xlim([0 100*pi]);
subplot(2,1,2); plot(UrethaneMalePhaseHA, 'b'); hold on;
subplot(2,1,2); plot(UrethaneFemalePhaseAH, 'm');
xlim([0 100*pi]);

figure(28); clf; % Male followed by Female PHASE
subplot(2,1,1); plot(ChronicMalePhaseAH, 'b'); hold on;
subplot(2,1,1); plot(ChronicFemalePhaseHA, 'm');
xlim([0 200*pi]);
subplot(2,1,2); plot(UrethaneMalePhaseAH, 'b'); hold on;
subplot(2,1,2); plot(UrethaneFemalePhaseHA, 'm');
xlim([0 200*pi]);


% xlim([0 200*pi]);
% subplot(2,1,2); plot(UrethaneMaleRealAH, 'b'); hold on;
% subplot(2,1,2); plot(UrethaneFemaleRealHA, 'm');
% xlim([0 500]);



%% Embedded vector strength function
function [vector_strength, phasespikes, regularspikes] = wPhasor(spiketimes, tims)

    wid = 0.250; paddington = 0.050;
    
    for j = length(tims):-1:1     
        
        prespikes = []; postspikes = [];
        PHASEprespikes = []; PHASEpostspikes = [];
        
        for k = 1:length(spiketimes)

        prespikes =  [prespikes,  ((spiketimes{k}(spiketimes{k} > tims(j,2)-wid & spiketimes{k} < tims(j,2))) - (tims(j,2)-wid) )'];
        postspikes = [postspikes, ((spiketimes{k}(spiketimes{k} > tims(j,2) & spiketimes{k} < tims(j,2)+wid)) - (tims(j,2)-wid) )'];
        
        PHASEprespikes = [PHASEprespikes, (pi * ((spiketimes{k}(spiketimes{k} > tims(j,1)-paddington & spiketimes{k} < tims(j,2)))-tims(j,1)) / (tims(j,2) - tims(j,1)))'];
        PHASEpostspikes = [PHASEpostspikes, (pi + (pi * ((spiketimes{k}(spiketimes{k} > tims(j,2) & spiketimes{k} < tims(j,3)+paddington ))-tims(j,2)) / (tims(j,3) - tims(j,2))))' ];
        
        end
    
        regularspikes{j} = sort([prespikes postspikes]);
        
        phasespikes{j} = sort([PHASEprespikes PHASEpostspikes]);
        phasespikes{j} = phasespikes{j}(~isnan(phasespikes{j}));
        
        vector_strength(j) = sqrt(mean(cos(phasespikes{j})).^2 + mean(sin(phasespikes{j})).^2);

    end

% figure(27); clf; bs_raster(phasespikes);
    
    
end


end
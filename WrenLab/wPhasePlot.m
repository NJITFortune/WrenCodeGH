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
plotwid = pi*2; plotwind = 100; % For Normalized data

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

%% Use this to plot each bird and calculate the mean
for ff = birdlist
    if in(ff).sexy == 1 % This is a male
        figure(ff); clf;
        subplot(411); bs_swPSTH(out(ff).AH.cPhaseSpikes,[0, plotwid], plotwind, 1); 
        subplot(412); mAHuP(ff) = bs_swPSTH(out(ff).AH.uPhaseSpikes,[0, plotwid], plotwind, 6);
        subplot(413); mHAcP(ff) = bs_swPSTH(out(ff).HA.cPhaseSpikes,[0, plotwid], plotwind, 1);
        subplot(414); mHAuP(ff) = bs_swPSTH(out(ff).HA.uPhaseSpikes,[0, plotwid], plotwind, 6);
        
        figure(ff+100); clf;
        subplot(411); bs_swPSTH(out(ff).AH.cSpikes,[0, Rplotwid], Rplotwind, 1); 
        subplot(412); mAHuR(ff) = bs_swPSTH(out(ff).AH.uSpikes,[0, Rplotwid], Rplotwind, 6);
        subplot(413); mHAcR(ff) = bs_swPSTH(out(ff).HA.cSpikes,[0, Rplotwid], Rplotwind, 1);
        subplot(414); mHAuR(ff) = bs_swPSTH(out(ff).HA.uSpikes,[0, Rplotwid], Rplotwind, 6);
    end

    if in(ff).sexy == 2 % This is a female
    figure(ff); clf;
        subplot(411); fAHcP(ff) = bs_swPSTH(out(ff).AH.cPhaseSpikes,[0, plotwid], plotwind, 2); 
        subplot(412); fAHuP(ff) = bs_swPSTH(out(ff).AH.uPhaseSpikes,[0, plotwid], plotwind, 6);
        subplot(413); bs_swPSTH(out(ff).HA.cPhaseSpikes,[0, plotwid], plotwind, 2);
        subplot(414); fHAuP(ff) = bs_swPSTH(out(ff).HA.uPhaseSpikes,[0, plotwid], plotwind, 6);
        
        figure(ff+100); clf;
        subplot(411); fAHcR(ff) = bs_swPSTH(out(ff).AH.cSpikes,[0, Rplotwid], Rplotwind, 1); 
        subplot(412); fAHuR(ff) = bs_swPSTH(out(ff).AH.uSpikes,[0, Rplotwid], Rplotwind, 6);
        subplot(413); bs_swPSTH(out(ff).HA.cSpikes,[0, Rplotwid], Rplotwind, 1);
        subplot(414); fHAuR(ff) = bs_swPSTH(out(ff).HA.uSpikes,[0, Rplotwid], Rplotwind, 6);
    end
end

UrethaneMalePhase = mHAuP(1).spers / max(mHAuP(1).spers);
UrethaneFemalePhase = fAHuP(2).spers / max(fAHuP(2).spers);

UrethaneMaleReal = mHAuR(1).spers / max(mHAuR(1).spers);
UrethaneFemaleReal = fAHuR(2).spers / max(fAHuR(2).spers);

for tt = 3:2:length(birdlist)

    UrethaneMalePhase = UrethaneMalePhase + mHAuP(tt).spers / max(mHAuP(tt).spers);
    UrethaneMaleReal = UrethaneMaleReal + mHAuR(tt).spers / max(mHAuR(tt).spers);
end

for tt = 4:2:length(birdlist)

    UrethaneFemalePhase = UrethaneFemalePhase + fAHuP(tt).spers / max(fAHuP(tt).spers);
    UrethaneFemaleReal = UrethaneFemaleReal + fAHuR(tt).spers / max(fAHuR(tt).spers);   
    
end

figure(28); clf;
subplot(2,1,1); plot(UrethaneMalePhase, 'b'); hold on;
subplot(2,1,1); plot(UrethaneFemalePhase, 'm');
xlim([0 200*pi]);
subplot(2,1,2); plot(UrethaneMaleReal, 'b'); hold on;
subplot(2,1,2); plot(UrethaneFemaleReal, 'm');
xlim([0 500]);




%% Embedded vector strength function
function [vector_strength, phasespikes, regularspikes] = wPhasor(spiketimes, tims)

    wid = 0.250;
    
    for j = length(tims):-1:1     
        
        prespikes = []; postspikes = [];
        PHASEprespikes = []; PHASEpostspikes = [];
        
        for k = 1:length(spiketimes)

        prespikes =  [prespikes,  ((spiketimes{k}(spiketimes{k} > tims(j,2)-wid & spiketimes{k} < tims(j,2))) - (tims(j,2)-wid) )'];
        postspikes = [postspikes, ((spiketimes{k}(spiketimes{k} > tims(j,2) & spiketimes{k} < tims(j,2)+wid)) - (tims(j,2)-wid) )'];
        
        PHASEprespikes = [PHASEprespikes, (pi * ((spiketimes{k}(spiketimes{k} > tims(j,1) & spiketimes{k} < tims(j,2)))-tims(j,1)) / (tims(j,2) - tims(j,1)))'];
        PHASEpostspikes = [PHASEpostspikes, (pi + (pi * ((spiketimes{k}(spiketimes{k} > tims(j,2) & spiketimes{k} < tims(j,3)))-tims(j,2)) / (tims(j,3) - tims(j,2))))' ];
        
        end
    
        regularspikes{j} = sort([prespikes postspikes]);
        
        phasespikes{j} = sort([PHASEprespikes PHASEpostspikes]);
        phasespikes{j} = phasespikes{j}(~isnan(phasespikes{j}));
        
        vector_strength(j) = sqrt(mean(cos(phasespikes{j})).^2 + mean(sin(phasespikes{j})).^2);

    end

% figure(27); clf; bs_raster(phasespikes);
    
    
end


end
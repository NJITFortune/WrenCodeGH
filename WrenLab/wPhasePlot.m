function [out, stts] = wPhasePlot(in, pad)
% Usage out = wPhasePlot(in)
% This makes a phase plot and calculates vector strength for AH and HA
% intervals. 

if nargin < 2; pad = 0; end

%% List of Chronic singing data with syllable indices and locations for spontaneous activity

[msolosyls, mduetsyls, fsolosyls, fduetsyls, Cspon, Aspon] = wData;

% birdlist = [3 4 9 10 11 12 13 14];

birdlist = 1:12;

%% For each bird
for ff = birdlist
    
    
% List of all duet syllables
allduetsyls = sort([mduetsyls{ceil(ff/2)}, fduetsyls{ceil(ff/2)}]);

%% Get all male-female and female-male duet syllable transitions
goodmalesyls = []; goodfemsyls = [];

for t = 1:length(mduetsyls{ceil(ff/2)})
   if ~isempty(find(allduetsyls == mduetsyls{ceil(ff/2)}(t)+1, 1))
       goodmalesyls(end+1) = mduetsyls{ceil(ff/2)}(t);
   end
end
for t = 1:length(fduetsyls{ceil(ff/2)})
   if ~isempty(find(allduetsyls == fduetsyls{ceil(ff/2)}(t)+1, 1))
       goodfemsyls(end+1) = fduetsyls{ceil(ff/2)}(t);
   end
end

%% Get all of the syllable times
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

plotwid = 0.300; plotwind = 100; % For raw data
% plotwid = pi*2; plotwind = 100; % For Normalized data

if in(ff).sexy == 1 % This is a male
    [out(ff).AH.Uvs, out(ff).AH.uPhaseSpikes, out(ff).AH.uSpikes] = wPhasor(in(ff).Aspikes, Msyltims);
    [out(ff).AH.Cvs, out(ff).AH.cPhaseSpikes, out(ff).AH.cSpikes] = wPhasor(in(ff).Cspikes, Msyltims);
    [out(ff).HA.Uvs, out(ff).HA.uPhaseSpikes, out(ff).HA.uSpikes] = wPhasor(in(ff).Aspikes, Fsyltims);
    [out(ff).HA.Cvs, out(ff).HA.cPhaseSpikes, out(ff).HA.cSpikes] = wPhasor(in(ff).Cspikes, Fsyltims);   
    
figure(ff); clf;
subplot(411); bs_swPSTH(out(ff).AH.cSpikes,[0, plotwid], plotwind, 1); 
subplot(412); bs_swPSTH(out(ff).AH.uSpikes,[0, plotwid], plotwind, 6);
subplot(413); bs_swPSTH(out(ff).HA.cSpikes,[0, plotwid], plotwind, 1);
subplot(414); bs_swPSTH(out(ff).HA.uSpikes,[0, plotwid], plotwind, 6);
end

if in(ff).sexy == 2 % This is a female
    [out(ff).AH.Uvs, out(ff).AH.uPhaseSpikes, out(ff).AH.uSpikes] = wPhasor(in(ff).Aspikes, Fsyltims);
    [out(ff).AH.Cvs, out(ff).AH.cPhaseSpikes, out(ff).AH.cSpikes] = wPhasor(in(ff).Cspikes, Fsyltims);
    [out(ff).HA.Uvs, out(ff).HA.uPhaseSpikes, out(ff).HA.uSpikes] = wPhasor(in(ff).Aspikes, Msyltims);
    [out(ff).HA.Cvs, out(ff).HA.cPhaseSpikes, out(ff).HA.cSpikes] = wPhasor(in(ff).Cspikes, Msyltims);    

figure(ff); clf;
subplot(411); bs_swPSTH(out(ff).AH.cSpikes,[0, plotwid], plotwind, 2); 
subplot(412); bs_swPSTH(out(ff).AH.uSpikes,[0, plotwid], plotwind, 6);
subplot(413); bs_swPSTH(out(ff).HA.cSpikes,[0, plotwid], plotwind, 2);
subplot(414); bs_swPSTH(out(ff).HA.uSpikes,[0, plotwid], plotwind, 6);
end



%     subplot(411); bs_raster(out(ff).AH.cp); 
%     subplot(412); bs_raster(out(ff).AH.up); 
%     subplot(413); bs_raster(out(ff).HA.cp); 
%     subplot(414); bs_raster(out(ff).HA.up); 



end % Cycle for every bird

%% Now that we have all birds, compare vector strength values for males and females and do stats

for ff = birdlist
    % Which urethane vector strength is better?    
    out(ff).maxUV = max([out(ff).HA.Uvs out(ff).AH.Uvs]);     
end


[stts.H, stts.P, stts.CI] = ttest2([out(1:2:end).maxUV], [out(2:2:end).maxUV]) 



function [vector_strength, phasespikes, regularspikes] = wPhasor(spiketimes, tims)

    wid = 0.150;
    
    for j = length(tims):-1:1     
        
        prespikes = []; postspikes = [];
        PHASEprespikes = []; PHASEpostspikes = [];
        
        for k = 1:length(spiketimes)

        prespikes =  [prespikes,  ((spiketimes{k}(spiketimes{k} > tims(j,2)-wid & spiketimes{k} < tims(j,2))) - tims(j,2)-wid)'];
        postspikes = [postspikes, ((spiketimes{k}(spiketimes{k} > tims(j,2) & spiketimes{k} < tims(j,2)+wid)) - tims(j,2))'];
        
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
function [out, stts] = wPhasePlot(in, pad)
% Usage out = wPhasePlot(in)
% This makes a phase plot and calculates vector strength for AH and HA
% intervals. 

if nargin < 2; pad = 0; end

%% List of Chronic singing data with syllable indices and locations for spontaneous activity

[msolosyls, mduetsyls, fsolosyls, fduetsyls, Cspon, Aspon] = wData;

% birdlist = [3 4 9 10 11 12 13 14];

birdlist = 1:12;

% For each bird
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

if in(ff).sexy == 1 % This is a male
    [out(ff).AH.uv, out(ff).AH.up] = wPhasor(in(ff).Aspikes, Msyltims);
    [out(ff).AH.cv, out(ff).AH.cp] = wPhasor(in(ff).Cspikes, Msyltims);
    [out(ff).HA.uv, out(ff).HA.up] = wPhasor(in(ff).Aspikes, Fsyltims);
    [out(ff).HA.cv, out(ff).HA.cp] = wPhasor(in(ff).Cspikes, Fsyltims);   
    
figure(ff); clf;
subplot(411); bs_swPSTH(out(ff).AH.cp,[0, pi*2], 100, 1); 
subplot(412); bs_swPSTH(out(ff).AH.up,[0, pi*2], 100, 6);
subplot(413); bs_swPSTH(out(ff).HA.cp,[0, pi*2], 100, 1);
subplot(414); bs_swPSTH(out(ff).HA.up,[0, pi*2], 100, 6);
end

if in(ff).sexy == 2 % This is a female
    [out(ff).AH.uv, out(ff).AH.up] = wPhasor(in(ff).Aspikes, Fsyltims);
    [out(ff).AH.cv, out(ff).AH.cp] = wPhasor(in(ff).Cspikes, Fsyltims);
    [out(ff).HA.uv, out(ff).HA.up] = wPhasor(in(ff).Aspikes, Msyltims);
    [out(ff).HA.cv, out(ff).HA.cp] = wPhasor(in(ff).Cspikes, Msyltims);    

figure(ff); clf;
subplot(411); bs_swPSTH(out(ff).AH.cp,[0, pi*2], 100, 2); 
subplot(412); bs_swPSTH(out(ff).AH.up,[0, pi*2], 100, 6);
subplot(413); bs_swPSTH(out(ff).HA.cp,[0, pi*2], 100, 2);
subplot(414); bs_swPSTH(out(ff).HA.up,[0, pi*2], 100, 6);
end



%     subplot(411); bs_raster(out(ff).AH.cp); 
%     subplot(412); bs_raster(out(ff).AH.up); 
%     subplot(413); bs_raster(out(ff).HA.cp); 
%     subplot(414); bs_raster(out(ff).HA.up); 



end % Cycle for every bird

%% Now that we have all birds, compare vector strength values for males and females and do stats

for ff = birdlist
    % Which urethane vector strength is better?    
    out(ff).maxUV = max([out(ff).HA.uv out(ff).AH.uv]);     
end


[stts.H, stts.P, stts.CI] = ttest2([out(1:2:end).maxUV], [out(2:2:end).maxUV]) 



function [vector_strength, phasespikes] = wPhasor(spiketimes, tims)
    
    
    for j = length(tims):-1:1     
        prespikes = []; postspikes = [];
        
        for k = 1:length(spiketimes)

        prespikes = [prespikes, (pi * ((spiketimes{k}(spiketimes{k} > tims(j,1) & spiketimes{k} < tims(j,2)))-tims(j,1)) / (tims(j,2) - tims(j,1)))'];
        postspikes = [postspikes, (pi + (pi * ((spiketimes{k}(spiketimes{k} > tims(j,2) & spiketimes{k} < tims(j,3)))-tims(j,2)) / (tims(j,3) - tims(j,2))))' ];
        
        end
    
        phasespikes{j} = sort([prespikes postspikes]);
        phasespikes{j} = phasespikes{j}(~isnan(phasespikes{j}));
        vector_strength(j) = sqrt(mean(cos(phasespikes{j})).^2 + mean(sin(phasespikes{j})).^2);

    end

% figure(27); clf; bs_raster(phasespikes);
    
    
end


end
function out = wPhasor(in)
% Usage out = wPhasePlot(in)
% This makes a phase plot and calculates vector strength for AH and HA
% intervals during duets. 

%% List of Chronic singing data with syllable indices and locations for spontaneous activity

[~, mduetsyls, ~, fduetsyls, ~, ~] = wData;

% Choose which data to analyze
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
    Msyltims(zz,1) = in(ff).syl(goodmalesyls(zz)).tim(1);
    Msyltims(zz,2) = mean([in(ff).syl(goodmalesyls(zz)+1).tim(1), in(ff).syl(goodmalesyls(zz)).tim(2) ]);
    Msyltims(zz,3) = in(ff).syl(goodmalesyls(zz)+1).tim(2);
end
for zz = length(goodfemsyls):-1:1   
    Fsyltims(zz,1) = in(ff).syl(goodfemsyls(zz)).tim(1);
    Fsyltims(zz,2) = mean([in(ff).syl(goodfemsyls(zz)+1).tim(1), in(ff).syl(goodfemsyls(zz)).tim(2) ]);
    Fsyltims(zz,3) = in(ff).syl(goodfemsyls(zz)+1).tim(2);
end

%% Fetch the vector strengths and plot

if in(ff).sexy == 1 % This is a male
    [out(ff).AH.Uvs, out(ff).AH.uPhaseSpikes, out(ff).AH.uSpikes] = wPhasor(in(ff).Aspikes, Msyltims);
    [out(ff).AH.Cvs, out(ff).AH.cPhaseSpikes, out(ff).AH.cSpikes] = wPhasor(in(ff).Cspikes, Msyltims);
    [out(ff).HA.Uvs, out(ff).HA.uPhaseSpikes, out(ff).HA.uSpikes] = wPhasor(in(ff).Aspikes, Fsyltims);
    [out(ff).HA.Cvs, out(ff).HA.cPhaseSpikes, out(ff).HA.cSpikes] = wPhasor(in(ff).Cspikes, Fsyltims);
    out(ff).sexy = 1;
end

if in(ff).sexy == 2 % This is a female
    [out(ff).AH.Uvs, out(ff).AH.uPhaseSpikes, out(ff).AH.uSpikes] = wPhasor(in(ff).Aspikes, Fsyltims);
    [out(ff).AH.Cvs, out(ff).AH.cPhaseSpikes, out(ff).AH.cSpikes] = wPhasor(in(ff).Cspikes, Fsyltims);
    [out(ff).HA.Uvs, out(ff).HA.uPhaseSpikes, out(ff).HA.uSpikes] = wPhasor(in(ff).Aspikes, Msyltims);
    [out(ff).HA.Cvs, out(ff).HA.cPhaseSpikes, out(ff).HA.cSpikes] = wPhasor(in(ff).Cspikes, Msyltims);    
    out(ff).sexy = 2;
end

end % Cycle for every bird

%% Embedded vector strength function
function [vector_strength, phasespikes, regularspikes] = wPhasor(spiketimes, tims)

    wid = 0.250; % This is window in msec before and after the middle of the ISI.
    paddington = 0.050; % Pad the window so that we minimize edge effects during analysis.
    
    for j = length(tims):-1:1     
        
        prespikes = []; postspikes = [];
        PHASEprespikes = []; PHASEpostspikes = [];
        
        for k = 1:length(spiketimes)

        prespikes =  [prespikes,  ((spiketimes{k}(spiketimes{k} > tims(j,2)-(wid+paddington) & spiketimes{k} < tims(j,2))) - (tims(j,2)-wid) )'];
        postspikes = [postspikes, ((spiketimes{k}(spiketimes{k} > tims(j,2) & spiketimes{k} < tims(j,2)+(wid+paddington))) - (tims(j,2)-wid) )'];
        
        PHASEprespikes = [PHASEprespikes, (pi * ((spiketimes{k}(spiketimes{k} > tims(j,1)-paddington & spiketimes{k} < tims(j,2)))-tims(j,1)) / (tims(j,2) - tims(j,1)))'];
        PHASEpostspikes = [PHASEpostspikes, (pi + (pi * ((spiketimes{k}(spiketimes{k} > tims(j,2) & spiketimes{k} < tims(j,3)+paddington ))-tims(j,2)) / (tims(j,3) - tims(j,2))))' ];
        
        end
    
        regularspikes{j} = sort([prespikes postspikes]);
        
        phasespikes{j} = sort([PHASEprespikes PHASEpostspikes]);
        phasespikes{j} = phasespikes{j}(~isnan(phasespikes{j}));
        
        vector_strength(j) = sqrt(mean(cos(phasespikes{j})).^2 + mean(sin(phasespikes{j})).^2);

    end

end


end



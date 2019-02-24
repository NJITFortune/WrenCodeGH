function foo = wPhasorPlot(data)

% There are two ways to handle this... 
% 1) Concantenate all of the spikes from all of the examples and calculate
% the histogram from that. The problem with this approach is that different
% recordings have different firing rates, which biases the outcome towards
% neurons with the highest firing rates.
% 2) Normalize all of the histograms and add them together. This treats
% each recording equally. The problem with this approach we will likely be
% averaging across potentially different categories of neurons.
% With luck, both approaches will yield similar results...

    ldat = length(data); 
    Pwin = [0 2*pi];
    %Rwin = [0 0.500];
    Pbinwid = 0.2;

%% Concatenate method
mcSPIKES = {}; muSPIKES = {};
fcSPIKES = {}; fuSPIKES = {};

for mdat = 1:2:ldat
    for kk = 1:length(data(mdat).HA.cPhaseSpikes)
        mcSPIKES{end+1} = data(mdat).HA.cPhaseSpikes{kk};
    end
    for kk = 1:length(data(mdat).HA.uPhaseSpikes)
        muSPIKES{end+1} = data(mdat).HA.uPhaseSpikes{kk};
    end
end
for fdat = 2:2:ldat
    for kk = 1:length(data(fdat).AH.cPhaseSpikes)
        fcSPIKES{end+1} = data(fdat).AH.cPhaseSpikes{kk};
    end
    for kk = 1:length(data(fdat).AH.uPhaseSpikes)
        fuSPIKES{end+1} = data(fdat).AH.uPhaseSpikes{kk};
    end
end

% mcSPIKES = mcSPIKES{2:end}; muSPIKES = muSPIKES{2:end};
% fcSPIKES = fcSPIKES{2:end}; fuSPIKES = fuSPIKES{2:end};

[foo.SmcSPER, foo.Smctim, foo.SmcFs] = swPSTH(mcSPIKES, Pwin, Pbinwid);
[foo.SmuSPER, foo.Smutim, foo.SmuFs] = swPSTH(muSPIKES, Pwin, Pbinwid);
[foo.SfcSPER, foo.Sfctim, foo.SfcFs] = swPSTH(fcSPIKES, Pwin, Pbinwid);
[foo.SfuSPER, foo.Sfutim, foo.SfuFs] = swPSTH(fuSPIKES, Pwin, Pbinwid);

% Normalize
foo.SmcSPER = foo.SmcSPER / max(foo.SmcSPER);
foo.SmuSPER = foo.SmuSPER / max(foo.SmuSPER);
foo.SfcSPER = foo.SfcSPER / max(foo.SfcSPER);
foo.SfuSPER = foo.SfuSPER / max(foo.SfuSPER);


figure(28); clf;
    xa(1) = subplot(211); hold on;
    plot(foo.Smctim, foo.SmcSPER, 'b', 'LineWidth', 2);
    plot(foo.Sfctim, foo.SfcSPER, 'm', 'LineWidth', 2);
    plot([pi, pi], [0 1.1], 'k');
    xa(2) = subplot(212); hold on;
    plot(foo.Smutim, foo.SmuSPER, 'b', 'LineWidth', 2);
    plot(foo.Sfutim, foo.SfuSPER, 'm', 'LineWidth', 2);
    plot([pi, pi], [0 1.1], 'k');
    
linkaxes(xa, 'xy'); xlim([foo.Smctim(1), foo.Smctim(end)]); ylim([0, 1.05]);




%% Average PSTHs

% Concatenate the Male data
for mdat = 1:2:ldat
    % Get the data in Spikes per Second
    [mcSPER(:,ceil(mdat/2)), foo.mctim, foo.mcFs] = swPSTH(data(mdat).HA.cPhaseSpikes, Pwin, Pbinwid);
    [muSPER(:,ceil(mdat/2)), foo.mutim, foo.muFs] = swPSTH(data(mdat).HA.uPhaseSpikes, Pwin, Pbinwid);
    % Normalize the data
    mcSPER(:,ceil(mdat/2)) = mcSPER(:,ceil(mdat/2)) / max(mcSPER(:,ceil(mdat/2)));
    muSPER(:,ceil(mdat/2)) = muSPER(:,ceil(mdat/2)) / max(muSPER(:,ceil(mdat/2)));
end
% Concatenate the Female data
for fdat = 2:2:ldat
    % Get the data in Spikes per Second
    [fcSPER(:,ceil(fdat/2)), foo.fctim, foo.fcFs] = swPSTH(data(fdat).AH.cPhaseSpikes, Pwin, Pbinwid);
    [fuSPER(:,ceil(fdat/2)), foo.futim, foo.fuFs] = swPSTH(data(fdat).AH.uPhaseSpikes, Pwin, Pbinwid);
    % Normalize the data
    fcSPER(:,ceil(fdat/2)) = fcSPER(:,ceil(fdat/2)) / max(fcSPER(:,ceil(fdat/2)));
    fuSPER(:,ceil(fdat/2)) = fuSPER(:,ceil(fdat/2)) / max(fuSPER(:,ceil(fdat/2)));
end

for jj = length(fcSPER(:,1)):-1:1
    
    foo.mcSPER(jj) = sum(mcSPER(jj,:));
        foo.mcSPERstd(jj) = std(mcSPER(jj,:));
    foo.fcSPER(jj) = sum(fcSPER(jj,:));
        foo.fcSPERstd(jj) = std(fcSPER(jj,:));
    foo.muSPER(jj) = sum(muSPER(jj,:));
        foo.muSPERstd(jj) = std(muSPER(jj,:));
    foo.fuSPER(jj) = sum(fuSPER(jj,:));
        foo.fuSPERstd(jj) = std(fuSPER(jj,:));
    
end

% And normalize again!
    foo.mcSPER = foo.mcSPER / max(foo.mcSPER);
        foo.mcSPERstd = foo.mcSPERstd / max(foo.mcSPER);
    foo.muSPER = foo.muSPER / max(foo.muSPER);
        foo.fcSPERstd = foo.fcSPERstd / max(foo.fcSPER);
    foo.fcSPER = foo.fcSPER / max(foo.fcSPER);
        foo.muSPERstd = foo.muSPERstd / max(foo.muSPER);
    foo.fuSPER = foo.fuSPER / max(foo.fuSPER);
        foo.fuSPERstd = foo.fuSPERstd / max(foo.fuSPER);

figure(27); clf;
    ax(1) = subplot(211); hold on;
    plot(foo.mctim, foo.mcSPER, 'b', 'LineWidth', 2);
    plot(foo.fctim, foo.fcSPER, 'm', 'LineWidth', 2);
    plot([pi, pi], [0 1.1], 'k');
    ax(2) = subplot(212); hold on;
    plot(foo.mutim, foo.muSPER, 'b', 'LineWidth', 2);
    plot(foo.futim, foo.fuSPER, 'm', 'LineWidth', 2);
    plot([pi, pi], [0 1.1], 'k');
    
linkaxes(ax, 'xy'); xlim([foo.futim(1), foo.futim(end)]); ylim([0, 1.05]);



function [spers, tims, Fs] = swPSTH(spikes, win, binwidth)

%% Preparations

% This is the amount of overlap from bin to bin
overlap = 0.1;

% len is window over which we will do the analysis.
% We add one bin to either side to reduce edge issues.

    len = (win(2) + binwidth) - (win(1) + binwidth);

% Calculate how many bins we will need - divide the length (len) by the
% stepsize (binwidth*overlap)

    binnum = round( len / (binwidth*overlap) );

% spb is spikes per bin and tims will be the midpoint of each bin, which
% will be set below

    spers = zeros(binnum,1);
    tims = spers;

%% Fill the bins, cycle by bin

for i = 1:binnum
    
    % Get the start of the bin, end of the bin, and the midpoint for each bin
    binstart = (i-1)*(binwidth*overlap) + win(1);
    binend   =     binstart+binwidth;   
    tims(i) = binstart + binwidth/2;
    
    % For each bin we now cycle through each rep of the stimulus
    
    for j = 1:length(spikes)  % length(spikes) gives us the number of reps     
        % get number of spikes (sum) over duration of the current bin
        repbinspikes = sum(spikes{j} > binstart & spikes{j} <= binend);
        % For the current bin (i) add the number of spikes (repbinisi)
        spers(i) = spers(i) + repbinspikes;
    end
end

% Normalize to spikes per second. 
    spers = spers / (binwidth * length(spikes));
    Fs = binnum/len;

end


end

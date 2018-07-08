function spikes = stimphasplot(spiketimes, stimtimes, binum, prt)
% stimphasplot(spiketimes, stimtimes, binum, prt)
% spiketimes - when the spikes happened in time (seconds)
% stimtimes - when the cycles start (must be periodic in seconds)
% binnum is optional: Number of bins you want for 0-360.
% prt is optional: 1 generates a PDF file, 0 is no PDF
% This depends on bs_raster

%% Setup

if nargin < 3; binum = 90; end; % 90 bins is default
if nargin < 4; prt = 0; end; % Default no print

set(gcf, 'Color', [1,1,1]); % This sets the background to white

%% Extract the spikes for each cycle and convert to degrees

    for i = length(stimtimes):-1:2;
        spikes{i-1} = (360 * (spiketimes(spiketimes > stimtimes(i-1) & spiketimes <= stimtimes(i)) - stimtimes(i-1))) / (stimtimes(i) - stimtimes(i-1));
    end;

%% The raster plot

subplot(2,1,1);

bs_raster(spikes);
    xlim([0 360]);
    h = gca; set(h, 'XTickLabel', [], 'Ydir', 'reverse', 'YAxisLocation', 'left');
    set(h, 'Position', [0.1300, 0.5100, 0.7750, 0.4400]);

%% The Histogram

subplot(2,1,2);

% How many bins?

    a = zeros(binum,1);
    stps = 0:360/binum:360;

% Fill the bins, cycle by bin

for i = 1:binum-1;
    binstart = stps(i);
    binend   = stps(i+1);   
    
    % and cycle here by reps of the stimulus
    
    for j = 1:length(spikes);       
        repbinisi = sum (spikes{j} > binstart & spikes{j} <= binend);
        a(i) = a(i) + repbinisi;
    end;
end;

% And plot!

    bar(stps(2:end),a,'k');
    xlim([0 360]);
    h = gca; set(h, 'XTickLabel', [], 'Box', 'off');
    
    
%% Printing to PDF
    
    if prt == 1;
        fn = [struct(stimnum).birdname struct(stimnum).date struct(stimnum).sex struct(stimnum).site '_' struct(stimnum).unit '-' num2str(stimnum)];
        print(gcf,'-dpng', fn);
    end;    

    
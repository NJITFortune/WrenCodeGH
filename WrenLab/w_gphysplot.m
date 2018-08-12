function w_gphysplot(in, win, binwidth, plt_type)
% w_gphysplot(spikes, stim, Fs, win, binwidth, plt_type, prt)
% win: [start stop] in Seconds (pre-onset times are negative!!!!)
% binwidth: The histogram bindwidth in milliseconds
% plt_type is optional: 0 for osc and 1 for specgram
% prt is optional: 1 generates a PDF file, 0 is no PDF
% This depends on bs_raster, takes data from bs_converter.

%% Setup

if nargin < 4; plt_type = 0; end % Default osc plot

signal = in.duet;
Fs = in.Fs;
binwidth = binwidth/1000;
tim = 1/Fs:1/Fs:length(signal)/Fs;

set(gcf, 'Color', [1,1,1]); % This sets the background to white

%% The raster plot

% The axs bit is so that we can link axes so the user can zoom in and out
% and all of the subplots will stay aligned in time
axs(1)=subplot(2,1,1);

% This calls our plotting function "bs_raster"
    bs_raster(in.Cspikes);
    
    % set the window width to what the user asked for in "win"
    xlim([win(1) win(2)]);
    % Position the plot appropriately
    h = gca; set(h, 'XTickLabel', [], 'Ydir', 'reverse', 'YAxisLocation', 'left');
    set(h, 'Position', [0.1300, 0.5100, 0.7750, 0.4400]);

%% The Histogram

axs(2)=subplot(4,1,3);

% How many bins?

    len = win(2) - win(1);
    binnum = round( len / binwidth);
    a = zeros(binnum,1);
    tims = a;

% Fill the bins, cycle by bin

for i = 1:binnum
    binstart = (i-1)*binwidth + win(1);
    binend   =     i*binwidth + win(1);   
    tims(i) = binstart;
    
    % and cycle here by reps of the stimulus
    
    for j = 1:length(in.Cspikes)       
        repbinisi = sum (in.Cspikes{j} > binstart & in.Cspikes{j} <= binend);
        a(i) = a(i) + repbinisi;
    end
end

% And now add the overlapping histogram (because we can)

    swpsthdata = bs_swPSTH(in.Cspikes, win, binwidth*1000, 0);

% And plot!

    bar(tims + binwidth/2,a,'k');
    hold on; 
    plot(swpsthdata.tim, swpsthdata.spers, 'r');
    xlim([win(1) win(2)]);
    h = gca; set(h, 'XTickLabel', [], 'Box', 'off');
    
%% Plot the stimulus

axs(3)=subplot(4,1,4);

    tt = find(tim > win(1) & tim < win(2));

% Option to plot an oscillogram
    if plt_type == 0
        plot(tim(tt),signal(tt));
        h = gca; set(h, 'Ycolor', [1,1,1], 'Box', 'off'); 
        f = axis;
        axis([win(1) win(2) f(3) f(4)]);
    end

% Option to plot an spectrogram
    if plt_type == 1
        specgram(signal(tt),2048,Fs,[],2000);
        xx = [abs(win(1)), abs(win(1))]; yy = [600, 4900];
        hold on; plot(xx, yy, 'b', 'LineWidth', 1); hold off;
        ylim([500 5000]); colormap(flipud(hot)); caxis([0 50]);
%         h = gca; ts = str2num(get(h,'XTickLabel')); set(h, 'Box', 'off');
%         ts = ts + win(1);
%         set(h,'XTickLabel', ts);
    end

    linkaxes(axs,'x');
        
    
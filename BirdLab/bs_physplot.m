function bs_physplot(struct, stimnum, win, binwidth, plt_type, prt)
% bs_physplot(struct, stimnum, win, binwidth, plt_type, prt)
% struct: The output from bs_converter
% stimnum: Which stimulus we want to process
% win: [start stop] in Seconds (pre-onset times are negative!!!!)
% binwidth: The histogram bindwidth in milliseconds
% plt_type is optional: 0 for osc and 1 for specgram
% prt is optional: 1 generates a PDF file, 0 is no PDF
% This depends on bs_raster, takes data from bs_converter.

%% Setup

if nargin < 5; plt_type = 0; end; % Default osc plot
if nargin < 6; prt = 0; end; % Default no print

binwidth = binwidth/1000;

%spikes = struct.stim(stimnum).spikes;
%signal = struct.stim(stimnum).stim; 
%tim = struct.stim(stimnum).tim;
% stimname = struct.stim(stimnum).StimName;
spikes = struct(stimnum).spikes;
signal = struct(stimnum).stim;
tim = struct(stimnum).tim;
stimname = struct(stimnum).StimName;

set(gcf, 'Color', [1,1,1]); % This sets the background to white

%% The raster plot

% The axs bit is so that we can link axes so the user can zoom in and out
% and all of the subplots will stay aligned in time
axs(1)=subplot(2,1,1);

% This calls our plotting function "bs_raster"
    bs_raster(spikes);
    
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

for i = 1:binnum;
    binstart = (i-1)*binwidth + win(1);
    binend   =     i*binwidth + win(1);   
    tims(i) = binstart;
    
    % and cycle here by reps of the stimulus
    
    for j = 1:length(spikes);       
        repbinisi = sum (spikes{j} > binstart & spikes{j} <= binend);
        a(i) = a(i) + repbinisi;
    end;
end;

% And now add the overlapping histogram (because we can)

    swpsthdata = bs_swPSTH(spikes, win, binwidth*1000, 0);

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
    if plt_type == 0;
        plot(tim(tt),signal(tt));
        h = gca; set(h, 'Ycolor', [1,1,1], 'Box', 'off'); 
        f = axis;
        axis([win(1) win(2) f(3) f(4)]);
        text(win(1)+1, -2, stimname, 'BackgroundColor',[.7 .9 .7]);
    end;

% Option to plot an spectrogram
    if plt_type == 1;
        specgram(signal(tt),2048,struct.stim(stimnum).Fs,[],2000);
        xx = [abs(win(1)), abs(win(1))]; yy = [600, 4900];
        text(1, 4000, struct(stimnum).StimName, 'BackgroundColor',[.7 .9 .7]);
        hold on; plot(xx, yy, 'b', 'LineWidth', 1); hold off;
        ylim([500 5000]); colormap(flipud(hot)); caxis([0 50]);
        h = gca; ts = str2num(get(h,'XTickLabel')); set(h, 'Box', 'off');
        ts = ts + win(1);
        set(h,'XTickLabel', ts);
    end;

    linkaxes(axs,'x');
    
% Label the plot with relevant information
    lbl = [struct(stimnum).birdname ', ' struct(stimnum).date ', ' struct(stimnum).sex ', site:' struct(stimnum).site ', unit:' struct(stimnum).unit];
    if struct(stimnum).sex == 'M';
    xlabel(lbl, 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0 0 1]);
    end;
    if struct(stimnum).sex == 'F';
    xlabel(lbl, 'FontSize', 16, 'FontWeight', 'bold', 'Color', [1 0 1]);
    end;
    
%% Printing to PDF
    
    if prt == 1;
        fn = [struct(stimnum).birdname struct(stimnum).date struct(stimnum).sex struct(stimnum).site '_' struct(stimnum).unit '-' num2str(stimnum)];
        print(gcf,'-dpng', fn);
    end;    

    
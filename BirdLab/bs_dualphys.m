function bs_dualphys(fspikes, mspikes, signal, tim, win, binwidth, plt_type, prt)
% bs_physplot(fspikes, mspikes, signal, tim, win, binwidth, plt_type, prt)
% struct: The output from bs_converter
% stimnum: Which stimulus we want to process
% win: [start stop] in Seconds (pre-onset times are negative!!!!)
% binwidth: The histogram bindwidth in milliseconds
% plt_type is optional: 0 for osc and 1 for specgram
% prt is optional: 1 generates a PDF file, 0 is no PDF
% This depends on bs_raster, takes data from bs_converter.

%% Setup

binwidth = binwidth/1000; % convert to seconds
Fs = 1/(tim(2) - tim(1));

if nargin < 7; plt_type = 0; end; % Default osc plot
if nargin < 8; prt = 0; end; % Default no print

% spikes = struct(stimnum).spikes;
% signal = struct(stimnum).stim; 
% tim = struct(stimnum).tim;

set(gcf, 'Color', [1,1,1]); % This sets the background to white

%% The raster plot

ax(1) = subplot(2,1,1);

bs_raster(mspikes, 'b');
bs_raster(fspikes, 'm');
    xlim([win(1) win(2)]);
    h = gca; set(h, 'XTickLabel', [], 'Ydir', 'reverse', 'YAxisLocation', 'left');
    set(h, 'Position', [0.1300, 0.5100, 0.7750, 0.4400]);

%% The Histogram

ax(2) = subplot(4,1,3);

% How many bins?

    len = win(2) - win(1);
    stepsize = binwidth/2;
    binnum = round( len / stepsize);
    fhist = zeros(binnum,1);
    mhist = fhist;
    tims = fhist;

% Fill the bins, cycle by bin

for i = 1:binnum;
    binstart = (i-1)*stepsize + win(1);
    binend   =     binstart + binwidth;   
    tims(i) = binstart;
    
    % and cycle here by reps of the stimulus
    
    for j = 1:length(fspikes);       
        frepbincnt = sum (fspikes{j} > binstart & fspikes{j} <= binend);
        fhist(i) = fhist(i) + frepbincnt;
    end;
    for j = 1:length(mspikes);       
        mrepbincnt = sum (mspikes{j} > binstart & mspikes{j} <= binend);
        mhist(i) = mhist(i) + mrepbincnt;
    end;    
end;

% And plot!

    plot(tims + binwidth/2, fhist,'m-', 'LineWidth', 2); 
    hold on;
    plot(tims + binwidth/2, mhist,'b-', 'LineWidth', 2);
    hold off;
    xlim([win(1) win(2)]);
    h = gca; set(h, 'XTickLabel', [], 'Box', 'off');
    
%% Plot the stimulus

    tt = find(tim > win(1) & tim < win(2));

% Option to plot an oscillogram
    if plt_type == 0;
        ax(3) = subplot(4,1,4);
        plot(tim(tt),signal(tt), 'k');
        h = gca; set(h, 'Ycolor', [1,1,1], 'Box', 'off'); 
        f = axis;
        axis([win(1) win(2) f(3) f(4)]);
        %text(win(1)+1, -2, struct(stimnum).StimName, 'BackgroundColor',[.7 .9 .7]);
    end;

% Option to plot an spectrogram
    if plt_type == 1;
        subplot(4,1,4);
        specgram(signal(tt),2048,Fs,[],2000);
        xx = [abs(win(1)), abs(win(1))]; yy = [600, 4900];
        % text(1, 4000, struct(stimnum).StimName, 'BackgroundColor',[.7 .9 .7]);
        hold on; plot(xx, yy, 'b', 'LineWidth', 1); hold off;
        ylim([500 5000]); colormap(flipud(hot)); caxis([0 50]);
        h = gca; ts = str2num(get(h,'XTickLabel')); set(h, 'Box', 'off');
        ts = ts + win(1);
        set(h,'XTickLabel', ts);
    end;

linkaxes(ax,'x');
    
% Label the plot with relevant information
%     lbl = [struct(stimnum).birdname ', ' struct(stimnum).date ', ' struct(stimnum).sex ', site:' struct(stimnum).site ', unit:' struct(stimnum).unit];
%     if struct(stimnum).sex == 'M';
%     xlabel(lbl, 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0 0 1]);
%     end;
%     if struct(stimnum).sex == 'F';
%     xlabel(lbl, 'FontSize', 16, 'FontWeight', 'bold', 'Color', [1 0 1]);
%     end;
    
%% Printing to PDF
    
    if prt == 1;
%        fn = [struct(stimnum).birdname struct(stimnum).date struct(stimnum).sex struct(stimnum).site '_' struct(stimnum).unit '-' num2str(stimnum)];
        fn = 'dualplot'
        print(gcf,'-dpng', fn);
    end;    

    
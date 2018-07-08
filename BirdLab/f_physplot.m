function f_physplot(struct, stimnum, win, binwidth, plt_type, prt)
% f_physplot(struct, stimnum, win, binwidth, plt_type, prt)
% struct: The output from bs_converter
% stimnum: Which stimulus we want to process
% win: [start stop] in Seconds (pre-ons\et times are negative!!!!)
% binwidth: The histogram bindwidth in milliseconds
% plt_type is optional: 0 for osc and 1 for specgram
% prt is optional: 1 generates a PDF file, 0 is no PDF
% This depends on bs_raster, takes data from bs_converter.

%% Setup

if nargin < 5; plt_type = 0; end; % Default pos plot
if nargin < 6; prt = 0; end; % Default no print

%spikes = struct.stim(stimnum).spikes;
%signal = struct.stim(stimnum).stim; 
%tim = struct.stim(stimnum).tim;
% stimname = struct.stim(stimnum).StimName;
spikes = struct(stimnum).spikes;
pos = struct(stimnum).stim;
tim = struct(stimnum).tim;
vel = struct(stimnum).vel;
acc = struct(stimnum).acc;
vatim = struct(stimnum).vatim;

stimname = struct(stimnum).StimName;

set(gcf, 'Color', [1,1,1]); % This sets the background to white

%% The raster plot

ax(1) = subplot(2,1,1);

bs_raster(spikes);
    xlim([win(1) win(2)]);
    h = gca; set(h, 'XTickLabel', [], 'Ydir', 'reverse', 'YAxisLocation', 'left');
    set(h, 'Position', [0.1300, 0.5100, 0.7750, 0.4400]);

%% The Histogram

ax(2) = subplot(4,1,3);

% How many bins?

    len = win(2) - win(1);
    binnum = round( len*1000 / binwidth);
    a = zeros(binnum,1);
    tims = a;

% Fill the bins, cycle by bin

for i = 1:binnum;
    binstart = (i-1)*binwidth + win(1)*1000;
    binend   =     i*binwidth + win(1)*1000;   
    tims(i) = binstart;
    
    % and cycle here by reps of the stimulus
    
    for j = 1:length(spikes);       
        repbinisi = sum (spikes{j}*1000 > binstart & spikes{j}*1000 <= binend);
        a(i) = a(i) + repbinisi;
    end;
end;

% And plot!

    bar(tims + binwidth/2,a,'k');
    xlim([win(1)*1000 win(2)*1000]);
    h = gca; set(h, 'XTickLabel', [], 'Box', 'off');
    
%% Plot the stimulus

    ax(3) = subplot(4,1,4);

    tt = find(tim > win(1) & tim < win(2));

% Option to plot position
    if plt_type == 0;
        plot(tim(tt),pos(tt));
        h = gca; set(h, 'Ycolor', [1,1,1], 'Box', 'off'); 
        f = axis;
        axis([win(1) win(2) f(3) f(4)]);
        text(win(1)+1, -2, stimname, 'BackgroundColor',[.7 .9 .7]);
    end;

% Option to plot velocity and acceleration
    if plt_type == 1;
        vatt = find(vatim > win(1) & vatim < win(2));
        plot(vatim(vatt),[vel(vatt) acc(vatt)*10]);
        h = gca; set(h, 'Ycolor', [1,1,1], 'Box', 'off'); 
        f = axis;
        axis([win(1) win(2) f(3) f(4)]);
        text(win(1)+1, -2, stimname, 'BackgroundColor',[.7 .9 .7]);
    end;

% Label the plot with relevant information
    lbl = [struct(stimnum).date ', site:' struct(stimnum).site ', unit:' struct(stimnum).unit];
    xlabel(lbl, 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0 0 1]);
    
    
%% Printing to PDF
    
    if prt == 1;
        fn = [struct(stimnum).birdname struct(stimnum).date struct(stimnum).sex struct(stimnum).site '_' struct(stimnum).unit '-' num2str(stimnum)];
        print(gcf,'-dpng', fn);
    end;    

    %linkaxes(ax, 'x');
    
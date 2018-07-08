function wPhysPlot(data)
% wPhysPlot(data)
% This script is hard coded to plot all the chronic duet data.  Please edit 
% the section "User parameters" to change what is plotted.

% h = msgbox('This script is hard coded to plot the chronic duet data. Please edit the section "User parameters" in the code to change what is plotted', 'wPhysPlot');
% uiwait(h); % Forces Matlab to wait until the user has clicked the OK button

%% User parameters

% Alter the parameters below to customize the plotting.

% How to plot the duet audio recording.  Plot type: 0 for oscillogram, 1 for spectrogram
plt_type = 1; 

% Do you want to plot the raster plot for the acute data?  Usually this
% should be no, which is Arasplot = 0;  This can take A LONG TIME.

Arasplot = 0;

% The duration, in seconds, of the bins used for the PSTH. binwidth = 0.050 is a good value.
binwidth = 0.10; % Usually 0.050

% Which data to plot? There are 16 entries (8 duets, Male and Female physiology). 
% Use 'Widx = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];' for all duets in the structure
% Widx = [11 12 13 14 15 16];
Widx = [3 4 5 6 7 8 9 10 11 12 13 14 15 16];
% These are the time window around the start of the duet that will be
% plotted. [-2 8] plots 2 seconds prior and 8 seconds after the onset of
% the duet song. These can be adjusted to zoom in and zoom out.

    % March 2017. 2 motifs with male solo prior to the duet
        win(:,1) = [-1 5.6]; % Male [-2 8] [0.1 6.6]
        win(:,2) = [-1 5.6]; % Female [-2 8] [0.1 6.6]
        
    % January 2016, 08:06. 2.5 motifs with female solo prior to the duet
        win(:,3) = [0 7]; % Male [-2 11]
        win(:,4) = [0 7]; % Female [-2 11]
        
    % January 2016, 08:07. 4 motifs with female solo prior to duet
        win(:,5) = [-2 12]; % Male [-2 10]
        win(:,6) = [-2 12]; % Female [-2 10]
        
    % January 2016, 08:15. 2 motifs with female solo prior to duet
        win(:,7) = [-0.3 7]; % Male [-2 6]
        win(:,8) = [-0.3 7]; % Female [-2 6]
        
    % January 2016, 10:09. 4 motifs with female change (motifs 1,2 and 3,4 are similar). Female solo prior to duet.
        win(:,9) = [-1 11]; % Male  [-1 11]
        win(:,10) = [-1 11]; % Female [-1 11]
        
    % January 2016, 10:22. 3 motifs with female solo prior to duet
        win(:,11) = [0 7]; % Male [-2 8]
        win(:,12) = [0 7]; % Female [-2 8]
        
    % January 2017, 08:48. 3 motifs with male solo prior to duet
        win(:,13) = [-2 14]; % Male [-2 14]
        win(:,14) = [-2 14]; % Female [-2 14]
        
    % January 2017, 17:33. 2 motifs with overlap syllable and gap
        win(:,15) = [-2 10]; % Male [-2 10]
        win(:,16) = [-2 10]; % Female [-2 10]
 
prt = 0; % Print to file? 0 for no, 1 for PNG, 2 for PDF. 


%% Loop through selected birds for plotting

for kk = 1:length(Widx)

% Copy the data into local variables 
    Cspikes = data(Widx(kk)).Cspikes;
    Aspikes = data(Widx(kk)).Aspikes;
    duet = data(Widx(kk)).duet;
    tim = data(Widx(kk)).tim;
    Fs = data(Widx(kk)).Fs;
    id = data(Widx(kk)).id;
    %sx = wrensex;

figure(kk); clf; set(gcf, 'Color', [1,1,1]); % This sets the background to white

% Chronic raster plot
    axs(kk).axs(1) = subplot(4,1,1);
    w_raster(Cspikes);
    
    % set the window width to what the user asked for in "win"
    xlim([win(1,Widx(kk)) win(2,Widx(kk))]);
    % Position the plot appropriately    
    h = gca; set(h, 'XTickLabel', [], 'Ydir', 'reverse', 'YAxisLocation', 'left');
%    set(h, 'Position', [0.1300, 0.5100, 0.7750, 0.4400]);

%% The sliding window Histogram

    Cpsthdata = w_swPSTH(Cspikes, win(:,Widx(kk)), binwidth);
    axs(kk).axs(2) = subplot(4,1,2); hold on;
    plot(Cpsthdata.tim, Cpsthdata.spers/max(Cpsthdata.spers), 'k-', 'LineWidth', 1.5);
    xlim([win(1,Widx(kk)) win(2,Widx(kk))]);
    h = gca; set(h, 'XTickLabel', [], 'Box', 'off');
    text(win(1,Widx(kk))+0.5, 0.8, 'Awake');

    axs(kk).axs(3) = subplot(4,1,3); hold on;
    Apsthdata = w_swPSTH(Aspikes, win(:,Widx(kk)), binwidth);
    plot(Apsthdata.tim, Apsthdata.spers/max(Apsthdata.spers), 'k-.', 'LineWidth', 1.5);
    xlim([win(1,Widx(kk)) win(2,Widx(kk))]);
    h = gca; set(h, 'XTickLabel', [], 'Box', 'off');
    text(win(1,Widx(kk))+0.5, 0.8, 'Urethane');

% Add the syllable times lines

for slbls = 1:length(data(Widx(kk)).sylsex)
    subplot(411); 
        plot([data(Widx(kk)).syl(slbls).tim(1) data(Widx(kk)).syl(slbls).tim(1)], [0 5], 'g-');
        plot([data(Widx(kk)).syl(slbls).tim(2) data(Widx(kk)).syl(slbls).tim(2)], [0 5], 'r-');
    subplot(412); 
        plot([data(Widx(kk)).syl(slbls).tim(1) data(Widx(kk)).syl(slbls).tim(1)], [0 1], 'g-');
        plot([data(Widx(kk)).syl(slbls).tim(2) data(Widx(kk)).syl(slbls).tim(2)], [0 1], 'r-');
    subplot(413); 
        plot([data(Widx(kk)).syl(slbls).tim(1) data(Widx(kk)).syl(slbls).tim(1)], [0 1], 'g-');
        plot([data(Widx(kk)).syl(slbls).tim(2) data(Widx(kk)).syl(slbls).tim(2)], [0 1], 'r-');    
end
    
    
%% Plot the stimulus

    tt = find(tim > win(1,Widx(kk)) & tim < win(2,Widx(kk)));

% Option to plot an oscillogram
    if plt_type == 0
        axs(kk).axs(4) = subplot(4,1,4); hold on;
        plot(tim(tt),duet(tt));
        h = gca; set(h, 'Ycolor', [1,1,1], 'Box', 'off'); 
        f = axis;
        axis([win(1,Widx(kk)) win(2,Widx(kk)) f(3) f(4)]);
        
        % Add the syllable lines to the oscillogram and add the colors
        % indicated sex. This code is duplicated below for the sonogram option.
        for lblss = 1:length(data(Widx(kk)).sylsex) % add syllable lines
            plot([data(Widx(kk)).syl(lblss).tim(1) data(Widx(kk)).syl(lblss).tim(1)], [-0.9 0.9], 'g-');
            plot([data(Widx(kk)).syl(lblss).tim(2) data(Widx(kk)).syl(lblss).tim(2)], [-0.9 0.9], 'r-');    

            if data(Widx(kk)).sylsex(lblss) == 1 
                subplot(411); plot([data(Widx(kk)).syl(lblss).tim(1) data(Widx(kk)).syl(lblss).tim(2)], [0.1 0.1], 'b-', 'LineWidth', 2); 
                subplot(414); plot([data(Widx(kk)).syl(lblss).tim(1) data(Widx(kk)).syl(lblss).tim(2)], [-0.9 -0.9], 'b-', 'LineWidth', 2); 
            end 
            if data(Widx(kk)).sylsex(lblss) == 2 
                subplot(411); plot([data(Widx(kk)).syl(lblss).tim(1) data(Widx(kk)).syl(lblss).tim(2)], [0.1 0.1], 'm-', 'LineWidth', 2); 
                subplot(414); plot([data(Widx(kk)).syl(lblss).tim(1) data(Widx(kk)).syl(lblss).tim(2)], [-0.9 -0.9], 'm-', 'LineWidth', 2); 
            end
            
        end
        text(win(1,Widx(kk))+0.5, 0.5, id, 'BackgroundColor',[.7 .9 .7]);
    end;

% Option to plot an spectrogram
    if plt_type == 1
        subplot(4,1,4); hold on;
        specgram(duet(tt), 1024, Fs, [], 1000);
        offsetim = abs(win(1,Widx(kk)));
        %xx = [offsetim, offsetim]; 
        %yy = [600, 4900];
        %hold on; plot(xx, yy, 'k', 'LineWidth', 1); 
        text(0.5, 4000, id, 'BackgroundColor',[.7 .9 .7]);
        %ylim([500 5000]); colormap(flipud('HOT')); caxis([0 40]);
        ylim([500 5000]); colormap(flipud(gray)); caxis([-10 20]);
        
        % Add the syllable lines to the oscillogram and add the colors
        % indicated sex. This code is duplicated above for the oscillogram option.
        for lblss = 1:length(data(Widx(kk)).sylsex) % add syllable lines
            plot([data(Widx(kk)).syl(lblss).tim(1)+offsetim data(Widx(kk)).syl(lblss).tim(1)+offsetim], [600 4500], 'g-');
            plot([data(Widx(kk)).syl(lblss).tim(2)+offsetim data(Widx(kk)).syl(lblss).tim(2)+offsetim], [600 4500], 'r-');    

            if data(Widx(kk)).sylsex(lblss) == 1 
                subplot(411); plot([data(Widx(kk)).syl(lblss).tim(1) data(Widx(kk)).syl(lblss).tim(2)], [0.1 0.1], 'b-', 'LineWidth', 2); 
                subplot(414); plot([data(Widx(kk)).syl(lblss).tim(1)+offsetim data(Widx(kk)).syl(lblss).tim(2)+offsetim], [600 600], 'b-', 'LineWidth', 2); 
            end 
            if data(Widx(kk)).sylsex(lblss) == 2 
                subplot(411); plot([data(Widx(kk)).syl(lblss).tim(1) data(Widx(kk)).syl(lblss).tim(2)], [0.1 0.1], 'm-', 'LineWidth', 2); 
                subplot(414); plot([data(Widx(kk)).syl(lblss).tim(1)+offsetim data(Widx(kk)).syl(lblss).tim(2)+offsetim], [600 600], 'm-', 'LineWidth', 2); 
            end
            
        end
        
        xlim([0 (abs(win(1,Widx(kk)))+win(2,Widx(kk)))]);
%         h = gca; ts = str2num(get(h,'XTickLabel')); set(h, 'Box', 'off');
%         ts = ts + win(1);
%         set(h,'XTickLabel', ts);
    end

    linkaxes(axs(kk).axs,'x');
    
% Label the plot with relevant information
    lbl = [id ', ' data(Widx(kk)).wrensex ];
    if data(Widx(kk)).wrensex == 'M'
    xlabel(lbl, 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0 0 1]);
    end;
    if data(Widx(kk)).wrensex == 'F'
    xlabel(lbl, 'FontSize', 16, 'FontWeight', 'bold', 'Color', [1 0 1]);
    end;
    
% Printing to PNG or PDF
    
    if prt ~= 0
        fn = [data(Widx(kk)).birdname data(Widx(kk)).date data(Widx(kk)).sex data(Widx(kk)).site '_' data(Widx(kk)).unit '-' num2str(Widx(kk))];
        if prt == 1; print(gcf,'-dpng', fn); end;
        if prt == 2; print(gcf,'-dpdf', fn); end;
    end    

% Do the Acute Raster plot 
if Arasplot ~= 0
   
% Acute raster plot in a separate window

    figure(kk+20); clf; set(gcf, 'Color', [1,1,1]); % This sets the background to white

    w_raster(Aspikes);
    
    % set the window width to what the user asked for in "win"
    xlim([win(1,Widx(kk)) win(2,Widx(kk))]);
    
end
    
    
end % End of main plot loop 
    
%% Function w_raster    
function w_raster(spiketimes, col)
% w_raster(spiketimes, col)
% spiketimes is a cell array of reps with spike times
% col is the color, e.g. 'b' or 'm'
% Version 18 September 2017. Original J.A.Bender 7/25/07

if nargin < 2; col = 'k'; end; % Default to black tick marks

hold on;

for ttt = 1:length( spiketimes )
   for rr = 1:length( spiketimes{ttt} )
      plot( ones(1,2).*spiketimes{ttt}(rr), [ttt-0.3 ttt+0.3], col, 'LineWidth', 0.5);
      % plot( ones(1,2).*spiketimes{tt}(ss), [tt-0.3 tt+0.3], 'k', 'LineWidth', 1.5);
   end
end

ax = axis;
axis( [ax(1:2) 0 length( spiketimes )+1] )

end % End of w_raster

%% Function w_swPSTH
function out = w_swPSTH(spikes, wind, binwidth)
% out = bs_srast(spikes, window, binwidth, plot)
% This produces a smoothed histogram.
% spikes is the spike array from bs_converter
% window is the time relative to the start of the stimulus (0) e.g. [-10 30]
% binwidth is the duration of the sum in milliseconds, e.g. 50
% plot is optional - use plot == 0 to suppress plotting.
% 1 (default) is blue, 2 magenta, 3 red, 4 green, 5 cyan, 6 black

% Preparations

% This is the amount of overlap from bin to bin
    overlap = 0.1;

% This is the length of the window over which we will do the analysis.
% win values are given by the user.
    len = wind(2) - wind(1);

% Calculate how many bins we will need - divide the length (len) by the
% stepsize (binwidth*overlap)
    binnum = round( len / (binwidth * overlap) );

% For speed, we make an array of zeros that is number of bins we will use 
% spb is spikes per bin and tims will be the midpoint of each bin, which
% will be set below
    spb = zeros(binnum,1);
    tims = spb;

% Fill the bins, cycle by bin

for i = 1:binnum
    % Get the start of the bin, end of the bin, and the midpoint for each bin
    binstart = (i-1)*(binwidth*overlap) + wind(1);
    binend   =     binstart+binwidth;   
    tims(i) = binstart + binwidth/2;
    
    % For each bin we now cycle through each rep of the stimulus
    
    for j = 1:length(spikes)  % length(spikes) gives us the number of reps     
        % get number of spikes (sum) over duration of the current bin
        % repbinisi is a temporary variable for the current rep (j)
        repbinisi = sum (spikes{j} > binstart & spikes{j} <= binend);
        % For the current bin (i) add the number of spikes (repbinisi)
        spb(i) = spb(i) + repbinisi;
    end;
end;

% Normalize to spikes per second. 
spb = spb / (binwidth * length(spikes));

% And put the data into our output structure
out.tim = tims;
out.spers = spb;
out.Fs = binnum/len;

end % End of w_swPSTH 




end % End of wPhysPlot
    
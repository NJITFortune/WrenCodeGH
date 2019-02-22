function out = bs_swPSTH(spikes, win, binwidth, plt)
% out = bs_srast(spikes, window, binwidth, plot)
% This produces a smoothed histogram.
% spikes is the spike array from bs_converter
% window is the time relative to the start of the stimulus (0) e.g. [-10 30]
% binwidth is the duration of the sum in milliseconds, e.g. 50
% plot is optional - use plot == 0 to suppress plotting.
% 1 (default) is blue, 2 magenta, 3 red, 4 green, 5 cyan, 6 black

%% Preparations

if nargin < 4; plt = 1; end
clrs(1,:)='b-'; clrs(2,:)='m-'; clrs(3,:)='r-'; 
clrs(4,:)='g-'; clrs(5,:)='c-'; clrs(6,:)='k-';

% This is the amount of overlap from bin to bin
overlap = 0.1;

% Convert from milliseconds to seconds (our spike times are recorded in seconds)
binwidth = binwidth/1000;

% This is the length of the window over which we will do the analysis.
% win values are given by the user.
    len = win(2) - win(1);

% Calculate how many bins we will need - divide the length (len) by the
% stepsize (binwidth*overlap)
    binnum = round( len / (binwidth*overlap) );

% For speed, we make an array of zeros that is number of bins we will use 
% spb is spikes per bin and tims will be the midpoint of each bin, which
% will be set below
    spb = zeros(binnum,1);
    tims = spb;

%% Fill the bins, cycle by bin

for i = 1:binnum
    % Get the start of the bin, end of the bin, and the midpoint for each bin
    binstart = (i-1)*(binwidth*overlap) + win(1);
    binend   =     binstart+binwidth;   
    tims(i) = binstart + binwidth/2;
    
    % For each bin we now cycle through each rep of the stimulus
    
    for j = 1:length(spikes)  % length(spikes) gives us the number of reps     
        % get number of spikes (sum) over duration of the current bin
        % repbinisi is a temporary variable for the current rep (j)
        repbinisi = sum (spikes{j} > binstart & spikes{j} <= binend);
        % For the current bin (i) add the number of spikes (repbinisi)
        spb(i) = spb(i) + repbinisi;
    end
end

% Normalize to spikes per second. 
spb = spb / (binwidth * length(spikes));

%% And plot!
if plt ~= 0 
    plot(tims,spb, clrs(plt,:), 'LineWidth', 2); xlim([win(1) win(2)]);
end

%% And put the data into our output structure
out.tim = tims;
out.spers = spb;
out.Fs = binnum/len;


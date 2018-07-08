function out = getspiketimes(signal, tim, ud)
% out = getspiketimes(spikes, time, updown);
% spikes is the raw data for threshold
% time is the time base
% updown is 1 for positive-going and 0 for negative-going spikes
% You should have filtered and prepared the spike data before using this
% script.

tmp = zeros(1,length(signal)); % pre-allocate 

% Invert the signal on request
	if ud == 0; signal = -signal; end;

% Plot the spikes
	figure(1); clf; plot(tim, signal);

% Wait until the user has zoomed or done whatever needs to be done...
	a = input('Hit return once the plot is prepared for your click');

% Get the y value of the single click
	[~,thresh] = ginput(1);

% Set all of the data above the threshold to 1 (formerly 0).
	tmp(signal > thresh) = 1;

% Get the spike times as the positive-going 'starts' of each above threshold epoch
	out = tim(diff(tmp) == 1);


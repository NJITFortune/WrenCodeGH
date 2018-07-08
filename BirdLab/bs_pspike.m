function [spikes,psp,t] = bs_pspike(icell)

mpd = 3e-3; % minimum peak distance in seconds 
mph  = 6;  % minimum peak height (to be multiplied by standard deviation)

v = icell.values;
Ts = icell.interval;
Fs =  1/Ts;



% % Low pass filter for PSPs 
% order = 5;
% Fcutoff = 10;
% Wnlow = Fcutoff * 2 / Fs;
% [b,a] = butter(order, Wnlow, 'low');
% psp = filtfilt(b, a, v);


Tmedian = 15e-3; % 5 ms

% median filter to eliminate spikes
window = round(Tmedian/Ts)
if mod(window,2) == 0;
    window = window + 1;
end
psp = RankOrderFilter(v,window,30);


[~,spike_idx] = findpeaks(v-psp,'minpeakheight',mph*std(psp), ... 
                   'minpeakdistance',round(mpd/Ts));

t = (0:length(v)-1)*Ts;

spikes = t(spike_idx);

% Get the spikes


function oscson(sig,Fs, opt)
% Usage: oscon(signal, Fs, [color thresholds]);
% Plots a specgram with parameters for wren songs with an oscillogram
% imposed on it. 
% [low hi] for color threshold levels on spectrogram.

%% Housekeeping
% Set up a time series
tim = 1/Fs:1/Fs:length(sig)/Fs;
oscplt = sig - min(sig); % Add the minimum to the signal - min is now zero.
cntr = 500; % This is the frequency around which the oscillo will be centered
oscplt = oscplt * cntr/max(oscplt); %% Center the signal at desired f.

%% Plotting
% specgram is depricated. spectrogram is the replacement with new syntax.
% Use appropriate for your Matlab version.
specgram(sig,1024,Fs,[],1000); ylim([0 4500]); 
% spectrogram(sig,1024,1000,1024,Fs,'yaxis'); ylim([0 5000]);

% Lower spectral resolution
% specgram(sig,512,Fs,[],500); ylim([0 5000]); 

%colormap("hot"); % Use for Octave
colormap('HOT'); % Use for Matlab

%% Set the colormap
% Get the original colormap range
[lo, hi] = caxis; 

% If user doesn't tell us what to do, we adjust automagically to ensure
% that you can see the signal
if nargin == 2;
    lo = hi - (hi-lo)*0.25; caxis([lo hi]);
end;

if nargin ==3;
    lo = opt(1); hi = opt(2); caxis([lo hi]);
end;

%% Now that we've adjusted the colors, plot the oscillogram on top.
hold on;
plot(tim, oscplt, 'w');
hold off;

end

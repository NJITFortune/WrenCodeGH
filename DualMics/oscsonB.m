function oscsonB(sig1,sig2, Fs, opt, topper)
% Usage: oscon(signal, Fs, [color thresholds]);
% Plots a specgram with parameters for wren songs with an oscillogram
% imposed on it. 
% [low hi] for color threshold levels on spectrogram.

%% Housekeeping
% Set up a time series
tim1 = 1/Fs:1/Fs:length(sig1)/Fs;
tim2 = 1/Fs:1/Fs:length(sig2)/Fs;

oscplt1 = sig1 - min(sig1); % Add the minimum to the signal - min is now zero.
oscplt2= sig2 - min(sig2);

cntr = 1000; % This is the frequency on the specgram around which the oscillo will be centered

oscplt1 = oscplt1 * cntr/max(oscplt1); %% Center the signal at desired f.
oscplt2 = oscplt2 * cntr/max(oscplt2); %% Center the signal at desired f.

    if topper == 'F'
        sex{1} = 'Female'; sex{2} = 'Male';
    else 
        sex{1} = 'Male'; sex{2} = 'Female';
    end


%% Plotting
% specgram is depricated. spectrogram is the replacement with new syntax.
% Use appropriate for your Matlab version.
figure(1)
ax1=subplot(2,1,1);
specgram(sig1,1024,Fs,[],1000); ylim([0 4500]);
hold on;

ax2=subplot(2,1,2);
specgram(sig2,1024,Fs,[],1000); ylim([0 4500]);
hold on;
% spectrogram(sig,1024,1000,1024,Fs,'yaxis'); ylim([0 5000]);

% Lower spectral resolution
% specgram(sig,512,Fs,[],500); ylim([0 5000]); 

%colormap("hot"); % Use for Octave
%colormap('parula'); % Use for Matlab

%% Set the colormap
% Get the original colormap range
[lo, hi] = caxis; 

% If user doesn't tell us what to do, we adjust automagically to ensure
% that you can see the signal
if nargin == 2
    lo = hi - (hi-lo)*0.25; caxis([lo hi]);
end

if nargin ==3
    lo = opt(1); hi = opt(2); caxis([lo hi]);
end

%% Now that we've adjusted the colors, plot the oscillogram on top.

colormap(ax2,'gray');
subplot(2,1,2)
hold on;
    plot(ax2, tim2,oscplt2,'k');
    text(0.1, 3000, sex{2});
hold off;

subplot(2,1,1)
hold on;
    plot(ax1, tim1, oscplt1, 'k');
    text(0.1, 3000, sex{1});
hold off


end

%% Load data
[ww, Fs] = audioread('/Users/eric/Downloads/Gordon2017tmp/wrenagade.wav');

%% Setup
overlap = 0.95;

clrmap = flipud(gray);
cmap = [-24 33];

%% Figure 1
figure(1); clf; set(gcf, 'Renderer', 'painters');
% 1024
nfft = 1024; 
subplot(411);
    specgram(ww, nfft, Fs, [], ceil(nfft*overlap)); 
    colormap(clrmap);
    ylim([500 5000]);
    caxis(cmap);
    text(0.2, 4400, 'Fs=44100, nFFT=1024, coloraxis=[-24 33]', 'FontSize', 18, 'Color', 'r');
    
% 512
nfft = 512; 
subplot(412);
    specgram(ww, nfft, Fs, [], ceil(nfft*overlap)); 
    colormap(clrmap);
    ylim([500 5000]);
    caxis(cmap);
    text(0.2, 4400, 'Fs=44100, nFFT=512, coloraxis=[-24 33]', 'FontSize', 18, 'Color', 'r');

maxtim = length(ww)/Fs;
    newFs = 10000;    
    newtim = 1/newFs:1/newFs:maxtim;
    tim = 1/Fs:1/Fs:maxtim;
 neww = resample(ww, tim, newFs); 
 
 % 1024
nfft = 1024; 
subplot(413);
    specgram(neww, nfft, newFs, [], ceil(nfft*overlap)); 
    colormap(clrmap);
    ylim([500 5000]);
    caxis(cmap);
    text(0.2, 4400, 'Fs=10000, nFFT=1024, coloraxis=[-24 33]', 'FontSize', 18, 'Color', 'b');
    
% 512
nfft = 512; 
subplot(414);
    specgram(neww, nfft, newFs, [], ceil(nfft*overlap)); 
    colormap(clrmap);
    ylim([500 5000]);
    caxis(cmap);
    text(0.2, 4400, 'Fs=10000, nFFT=512, coloraxis=[-24 33]', 'FontSize', 18, 'Color', 'b');


    
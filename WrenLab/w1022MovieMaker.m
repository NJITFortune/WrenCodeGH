%% Load the data and plot an initial confirmation

load /Users/eric/Sync/Wren/ChronicCompleat2019f.mat

[Svideo, Fs] = audioread('~/Sync/Wren/cVideo/ChronicDuet_long_maybe.wav');
    Avideo = Svideo(:,1); % Use only one channel of the audio

rango = [-3.0248, 9.9904]; % Range from the w(11).tim that matches the video file

figure(1); clf; hold on;

    tim = 1/Fs:1/Fs:length(Avideo)/Fs;
    plot(tim, Avideo);
    
    tt = find(w(11).tim > rango(1) & w(11).tim <= rango(2));
    plot(w(11).tim(tt) - w(11).tim(tt(1)), w(11).duet(tt));

%% Setup the video read and the video output
v = VideoReader('~/Sync/Wren/cVideo/ChronicDuet_long_maybe.mov');    % 391 frames
% 1920 x 1080, fps = 29.82, H.264, AC3 48000 Hz (REAL width 1440)

nFrames = ceil(v.FrameRate*v.Duration); 
s(nFrames) = struct('cdata',[],'colormap',[]);


%% Loop

hFig = figure('MenuBar','none', 'Units','pixels', 'InnerPosition',[100 100 1920 1080]);

hAx = axes('Parent',hFig,'Units','pixels','NextPlot','add','Visible','off','XTick',[],'YTick',[],'Position',[10 40 1980 1100]);
    im = readFrame(v);
    image(flipud(im));
    hold on;
    plot([1430 8 8 1430 1430], [10 10 157 157 10], 'k-', 'LineWidth', 6)

hBx = axes('Parent',hFig,'Units','pixels','NextPlot','add','Visible','off','XTick',[],'YTick',[],'Position',[20 50 1930 160]);
    specgram(w(11).duet, 512, w(11).Fs, [], round(0.95*512));
    ug = flipud(gray); colormap(ug); caxis ([-20 33]);
    hold on;
    curtim = v.CurrentTime;
    
    plot([4 4], [0 5000], 'r-', 'Linewidth', 3);


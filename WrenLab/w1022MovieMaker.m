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

%% 
v = VideoReader('~/Sync/Wren/cVideo/ChronicDuet_long_maybe.mov');    % 391 frames
% 1920 x 1080, fps = 29.82, H.264, AC3 48000 Hz (REAL width 1440)

nFrames = ceil(v.FrameRate*v.Duration); 
s(nFrames) = struct('cdata',[],'colormap',[]);

hFig = figure('MenuBar','none', 'Units','pixels', 'InnerPosition',[100 100 1920 1080]);

% hAx = axes('Parent',hFig,'Units','pixels','NextPlot','add','Visible','off','XTick',[],'YTick',[],'Position',[0 0 1920 1080]);
hAx = axes('Parent',hFig,'Units','pixels','NextPlot','add','Visible','off','XTick',[],'YTick',[],'Position',[10 40 1980 1100]);

im = readFrame(v);
image(flipud(im));
hold on;
plot([10 1430], [155 155], 'k-', 'LineWidth', 8)
plot([10 1430], [10 10], 'k-', 'LineWidth', 6)
plot([10 10], [10 155], 'k-', 'LineWidth', 7)
plot([1430 1430], [10 155], 'k-', 'LineWidth', 6)

hBx = axes('Parent',hFig,'Units','pixels','NextPlot','add','Visible','off','XTick',[],'YTick',[],'Position',[20 50 1930 160]);
specgram(w(11).duet, 512, w(11).Fs, [], round(0.95*512));
ug = flipud(gray); colormap(ug); caxis ([-20 33]);


%% Make the video

f = readFrame(v);
fp = figure(1); clf; set(fp,'position',[150 150 1920 1080]);
imshow(f,'InitialMagnification','fit');
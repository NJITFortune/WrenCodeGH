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
% 1920 x 1080, fps = 29.82, H.264, AC3 48000 Hz

nFrames = ceil(v.FrameRate*v.Duration); 
s(nFrames) = struct('cdata',[],'colormap',[]);

hFig = figure('MenuBar','none', 'Units','pixels', 'Position',[100 100 1920 1080]);

hAx = axes('Parent',hFig,...
    'Units','pixels',...
    'Position',[0 0 1920 1080],...
    'NextPlot','add',...
    'Visible','off',...
    'XTick',[],...
    'YTick',[]);
hIm = image(uint8(zeros(1080,1920,3)),...
    'Parent',hAx);
hLine(1) = plot(hAx,[1 v.Width],[3 3],'-b','LineWidth',2);
hLine(2) = plot(hAx,[1 1],[3 3],'-b','LineWidth',4);
hLine(3) = plot(hAx,1,3,'ob','MarkerSize',10,'MarkerFaceColor','b');
im = readFrame(v);
hIm.CData = im; hIm.CData = flipud(im);


%% Make the video

f = readFrame(v);
fp = figure(1); clf; set(fp,'position',[150 150 1920 1080]);
imshow(f,'InitialMagnification','fit');
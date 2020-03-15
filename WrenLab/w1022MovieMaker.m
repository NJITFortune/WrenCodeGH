%% Load the data and plot an initial confirmation

load /Users/eric/Sync/Wren/ChronicCompleat2019f.mat

% [Svideo, Fs] = audioread('~/Sync/Wren/cVideo/ChronicDuet_long_maybe.wav');
%     Avideo = Svideo(:,1); % Use only one channel of the audio

rango = [-3.0248, 9.9904]; % Range from the w(11).tim that matches the video file
specpos = -rango(1);

% figure(1); clf; hold on;

%     tim = 1/Fs:1/Fs:length(Avideo)/Fs;
%     plot(tim, Avideo);
    
    tt = find(w(11).tim > rango(1) & w(11).tim <= rango(2));
%     plot(w(11).tim(tt) - w(11).tim(tt(1)), w(11).duet(tt));

%% Setup the video read and the video output
v = VideoReader('~/Sync/Wren/cVideo/ChronicDuet_long_maybe.mov');    % 391 frames
% 1920 x 1080, fps = 29.82, H.264, AC3 48000 Hz (REAL width 1440)

nFrames = ceil(v.FrameRate * v.Duration); 
s(nFrames) = struct('cdata',[],'colormap',[]);

hFig = figure('MenuBar','none', 'Units','pixels', 'InnerPosition', [100 100 1920 1080]);

   writerObj = VideoWriter('mymovie2.avi', 'Uncompressed AVI');
   writerObj.FrameRate = v.FrameRate;
   open(writerObj);


%% Loop

% for jjj = 1:5
while hasFrame(v) 

    im = readFrame(v);


hAx = axes('Parent',hFig,'Units','pixels','NextPlot','add','Visible','off','XTick',[],'YTick',[],'Position',[10 40 1980 1100]);
    image(flipud(im));
    hold on;
    plot([1430 8 8 1430 1430], [10 10 157 157 10], 'k-', 'LineWidth', 6)
    
% hBx = axes('Parent',hFig,'Units','pixels','NextPlot','add','Visible','off','XTick',[],'YTick',[],'Position',[20 50 1930 160]);
hBx = axes('Parent',hFig,'Units','pixels','NextPlot','add','Visible','off','XTick',[],'YTick',[],'Position',[20 50 2020 160]);

    specgram(w(11).duet(tt), 512, w(11).Fs, [], round(0.95*512));
    ug = flipud(gray); colormap(ug); caxis ([-20 33]);
    hold on;
    curtim = v.CurrentTime;
    plot([curtim curtim], [0, 5000], 'r-', 'LineWidth', 5); % Red progress line

    % Plot transparent color boxes
    
    for f = 1:length(w(11).syl)
        
        if w(11).syl(f).tim(1) < curtim - specpos
           if w(11).sylsex(f) == 1 % Male
           maxX = min([w(idx(11)).syl(f).tim(2)+specpos, curtim]);    
fill([specpos+w(11).syl(f).tim(1), maxX, maxX, specpos+w(11).syl(f).tim(1)], ...
    [750, 750, 4000, 4000], 'c', 'FaceAlpha', 0.1, 'LineStyle', 'none');
           end
           if w(11).sylsex(f) == 2 % Female
           maxX = min([w(idx(12)).syl(f).tim(2)+specpos, curtim]);    
fill([specpos+w(12).syl(f).tim(1), maxX, maxX, specpos+w(12).syl(f).tim(1)], ...
    [750, 750, 4000, 4000], 'm', 'FaceAlpha', 0.1, 'LineStyle', 'none');
           end
        
        end
    end
    
    % Plot spikes
    
        for j = 1:4
            
           malespkidx = find(w(11).Cspikes{j} < curtim - specpos & w(11).Cspikes{j} > -specpos);
           femalespkidx = find(w(12).Cspikes{j} < curtim - specpos & w(12).Cspikes{j} > -specpos);
           
           for k = 1:length(malespkidx)
               plot([specpos+w(11).Cspikes{j}(malespkidx(k)), specpos+w(11).Cspikes{j}(malespkidx(k))], [100+(j*200), 200+(j*200)+180], 'b-', 'LineWidth', 1);
           end
           for k = 1:length(femalespkidx)
               plot([specpos+w(12).Cspikes{j}(femalespkidx(k)), specpos+w(12).Cspikes{j}(femalespkidx(k))], [3900+(j*200), 3900+(j*200)+180], 'm-', 'LineWidth', 1);
           end
        end
        
    drawnow
    pause(0.1);
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
    pause(0.1);
 clf;

end
 
close(writerObj);

 
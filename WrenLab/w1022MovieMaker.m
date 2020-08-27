%% Load the data and plot an initial confirmation

load /Users/eric/Sync/Wren/ChronicCompleat2019f.mat


rango = [-3.0248, 9.9904]; % Range from the w(11).tim that matches the video file
specpos = -rango(1);
    
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
    plot([1430 8 8 1430 1430], [10 10 188 188 10], 'k-', 'LineWidth', 6);
    text(150, 1020, 'Female', 'Color', 'm', 'FontSize', 64);
    text(750, 1020, 'Male', 'Color', 'b', 'FontSize', 64);
    
% hBx = axes('Parent',hFig,'Units','pixels','NextPlot','add','Visible','off','XTick',[],'YTick',[],'Position',[20 50 1930 160]);
hBx = axes('Parent',hFig,'Units','pixels','NextPlot','add','Visible','off','XTick',[],'YTick',[],'Position',[20 20 2020 220]);

    specgram(w(11).duet(tt), 512, w(11).Fs, [], round(0.95*512));
    ug = flipud(gray); colormap(ug); caxis ([-20 33]);
    hold on;
    curtim = v.CurrentTime;
    
    fff = (w(11).syl(4).tim(1)-w(11).syl(3).tim(2))/2;
    plot([w(11).syl(3).tim(2)+specpos+fff w(11).syl(3).tim(2)+specpos+fff], [0, 5000], 'k--', 'LineWidth', 1);

    
    plot([curtim curtim], [0, 5000], 'r-', 'LineWidth', 5); % Red progress line

    % Plot transparent color boxes
    
    for f = 1:length(w(11).syl)
        
        if w(11).syl(f).tim(1) < curtim - specpos
           if w(11).sylsex(f) == 1 % Male
           maxX = min([w(11).syl(f).tim(2)+specpos, curtim]);    
fill([specpos+w(11).syl(f).tim(1), maxX, maxX, specpos+w(11).syl(f).tim(1)], ...
    [750, 750, 4000, 4000], 'c', 'FaceAlpha', 0.2, 'LineStyle', 'none');
           end
           if w(11).sylsex(f) == 2 % Female
           maxX = min([w(12).syl(f).tim(2)+specpos, curtim]);    
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
               plot([specpos+w(11).Cspikes{j}(malespkidx(k)), specpos+w(11).Cspikes{j}(malespkidx(k))], ...
                   [(j*200)-180, (j*200)+0], 'b-', 'LineWidth', 1);
           end
           for k = 1:length(femalespkidx)
               plot([specpos+w(12).Cspikes{j}(femalespkidx(k)), specpos+w(12).Cspikes{j}(femalespkidx(k))], ...
                   [3900+(j*200), 3900+(j*200)+180], 'm-', 'LineWidth', 1);
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


%% Make the audio track

rango = [-3.0248, 9.9904]; % Range from the w(11).tim that matches the video file
specpos = -rango(1);

% Read the audio file
[Svideo, Fs] = audioread('~/Sync/Wren/cVideo/ChronicDuet_long_maybe.wav');
     Avideo = Svideo(:,1); % Use only one channel of the audio

% Truncate 
    tim = 1/Fs:1/Fs:length(Avideo)/Fs;
    tim = tim - specpos;

%    tt = find(tim > rango(1) & tim <= rango(2));
%    Avideo = Avideo(tt);
%    tim = tim(tt);

% Make the Fake Spikes

    spiketim = 1/Fs:1/Fs:0.004; % 4 msec duration for our fake spikes
    len = length(spiketim);
    rmp = 1/floor(len/2):1/floor(len/2):1;

    for kk = 1:4
        fspike(:,kk) = sin(2*pi*(3000)*spiketim) * 0.5 ; % 1000 Hz for females
            fspike(1:length(rmp),kk) = fspike(1:length(rmp),kk)' .* rmp;
            fspike(end+1-length(rmp):end,kk) = fspike(end+1-length(rmp):end,kk)' .* rmp(end:-1:1);
        mspike(:,kk) = sin(2*pi*(1000)*spiketim) * 0.5 ; % 200 Hz for males
            mspike(1:length(rmp),kk) = mspike(1:length(rmp),kk)' .* rmp;
            mspike(end+1-length(rmp):end,kk) = mspike(end+1-length(rmp):end,kk)' .* rmp(end:-1:1);
    end
        
% Make the spike sound trains
    fem = zeros(1,length(tim));
    mal = zeros(1,length(tim));

for k = 1:4
    spkidx = find(w(11).Cspikes{k} > rango(1) & w(11).Cspikes{k} < rango(2));
    for j = 1:length(spkidx)   
        curidx = find(tim >= w(11).Cspikes{k}(spkidx(j)), 1, 'first');    
    if curidx+len-1 < length(tim)  % Need the if not to go over the end.  
        mal(curidx:curidx+len-1) = mal(curidx:curidx+len-1) + mspike(:,k)';
    end
    
    end
end

% idx = 2; % This is the female (even)

for k = 1:4
    spkidx = find(w(12).Cspikes{k} > rango(1) & w(12).Cspikes{k} < rango(2));
    for j = 1:length(spkidx)   
        curidx = find(tim >= w(12).Cspikes{k}(spkidx(j)), 1, 'first');    
    if curidx+len-1 < length(tim)  % Need the if not to go over the end.   
        fem(curidx:curidx+len-1) = fem(curidx:curidx+len-1) + fspike(:,k)';
    end
    
    end
end

spks = fem' + mal';

combo = spks + (Avideo*0.4); 
    combo = 0.95 * (combo / max(abs(combo)));
femonly = (Avideo*0.4) + fem';
    femonly = 0.95 * (femonly / max(abs(femonly)));
malonly = (Avideo*0.4) + mal';
    malonly = 0.95 * (malonly / max(abs(malonly)));
audiowrite('combo.wav', combo, Fs);
audiowrite('malonly.wav', malonly, Fs);
audiowrite('femonly.wav', femonly, Fs);
 
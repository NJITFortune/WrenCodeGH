%% Make the audio file
idx = [1, 2] ; % This is the Male (odd)
rango = [0.1, 5.8];
%idx = [5, 6] ; % This is the Male (odd)
%rango = [0, 9];

    specpos = 0;
    if rango(1) < 0; specpos = abs(rango(1)); end
    if rango(1) > 0; specpos = -rango(1); end

Fs = w(idx).Fs;

% Make the Fake Spikes

    spiketim = 1/Fs:1/Fs:0.005; % 4 msec duration for our fake spikes
    len = length(spiketim);
    rmp = 1/floor(len/2):1/floor(len/2):1;

    for kk = 1:4
        fspike(:,kk) = sin(2*pi*(1200+(50*kk))*spiketim) * 0.5 ; % 800 Hz for females
            fspike(1:length(rmp),kk) = fspike(1:length(rmp),kk)' .* rmp;
            fspike(end+1-length(rmp):end,kk) = fspike(end+1-length(rmp):end,kk)' .* rmp(end:-1:1);
        mspike(:,kk) = sin(2*pi*(800+(80*kk))*spiketim) * 0.25 ; % 300 Hz for males
            mspike(1:length(rmp),kk) = mspike(1:length(rmp),kk)' .* rmp;
            mspike(end+1-length(rmp):end,kk) = mspike(end+1-length(rmp):end,kk)' .* rmp(end:-1:1);
    end
    
    outim = w(idx(1)).tim(w(idx(1)).tim > rango(1) & w(idx(1)).tim < rango(2));
    outduet = w(idx(1)).duet(w(idx(1)).tim > rango(1) & w(idx(1)).tim < rango(2));
    
% Make the spike sound trains
    fem = zeros(1,length(outim));
    mal = zeros(1,length(outim));

for k = 1:4
    spkidx = find(w(idx(1)).Cspikes{k} > rango(1) & w(idx(1)).Cspikes{k} < rango(2));
    for j = 1:length(spkidx)   
        curidx = find(outim >= w(idx(1)).Cspikes{k}(spkidx(j)), 1, 'first');    
    if curidx+len-1 < length(outim)  % Need the if not to go over the end.   
        mal(curidx:curidx+len-1) = mal(curidx:curidx+len-1) + mspike(:,k)';
        % mal(curidx:curidx+len-1) = mspike;
    end
    
    end
end

% idx = 2; % This is the female (even)

for k = 1:4
    spkidx = find(w(idx(2)).Cspikes{k} > rango(1) & w(idx(2)).Cspikes{k} < rango(2));
    for j = 1:length(spkidx)   
        curidx = find(outim >= w(idx(2)).Cspikes{k}(spkidx(j)), 1, 'first');    
    if curidx+len-1 < length(outim)  % Need the if not to go over the end.   
        fem(curidx:curidx+len-1) = fspike;
    end
    
    end
end

spks = fem' + mal';

combo = spks + (outduet*0.8);
femonly = (outduet*0.8) + fem';
malonly = (outduet*0.8) + mal';
audiowrite('combo.wav', combo, 10000);
audiowrite('malonly.wav', malonly, 10000);
audiowrite('femonly.wav', femonly, 10000);

%% Quick plot in prep for the video

    vFs = 30; % frames per second

figure(1); clf; 
    set(gcf,'position',[100,100,640,480]);
    specgram(outduet, 512, Fs);
    colormap('HOT');
    % truesize([960 1280]); 
    % truesize([480 640]);
    % specgram(outduet, 512, Fs);
    caxis([-25 25])
    cmp = flipud(gray);
    %colormap(cmp);
    hold on;

    %% Generate the video
    
% Initialize the "object" that will be the final movie
    writerObj = VideoWriter('mymovie2.avi', 'Uncompressed AVI');
    writerObj.FrameRate = 30;
%writerObj.VideoFormat = 'RGB24';
    %writerObj.Quality = 100;
    open(writerObj);

for vtim = outim(1):1/30:outim(end)
    clf; 
    specgram(outduet, 512, Fs, [], round(512*0.95));
    colormap(cmp);
    caxis([-25 25])
    hold on;
    
    plot([specpos+vtim,specpos+vtim], [0, 5000], 'r-', 'LineWidth', 5); % Red progress line
    
    % Plot transparent color boxes
    
    for f = 1:length(w(idx(1)).syl)
%    for f = 1:13
        
        if w(idx(1)).syl(f).tim(1) >= rango(1) && w(idx(1)).syl(f).tim(1) < vtim
           maxX = min([w(idx(1)).syl(f).tim(2), vtim]);
           if w(idx(1)).sylsex(f) == 1 % Male
           fill([specpos+w(idx(1)).syl(f).tim(1), specpos+maxX, specpos+maxX, specpos+w(idx(1)).syl(f).tim(1)],...
               [750, 750, 4000, 4000], 'c', 'FaceAlpha', 0.2, 'LineStyle', 'none');
           end
           if w(idx(1)).sylsex(f) == 2 % Female
           fill([specpos+w(idx(1)).syl(f).tim(1), specpos+maxX, specpos+maxX, specpos+w(idx(1)).syl(f).tim(1)],...
               [750, 750, 4000, 4000], 'm', 'FaceAlpha', 0.1, 'LineStyle', 'none');
           end
        
        end
    end
    
    % Plot spikes
    
        for j = 1:4
           malespkidx = find(w(idx(1)).Cspikes{j} > rango(1) & w(idx(1)).Cspikes{j} < vtim);
           femalespkidx = find(w(idx(2)).Cspikes{j} > rango(1) & w(idx(2)).Cspikes{j} < vtim);
           
           for k = 1:length(malespkidx)
               plot([specpos+w(idx(1)).Cspikes{j}(malespkidx(k)), specpos+w(idx(1)).Cspikes{j}(malespkidx(k))], [200+(j*100), 200+(j*100)+90], 'b-', 'LineWidth', 1);
           end
           for k = 1:length(femalespkidx)
               plot([specpos+w(idx(2)).Cspikes{j}(femalespkidx(k)), specpos+w(idx(2)).Cspikes{j}(femalespkidx(k))], [4000+(j*100), 4000+(j*100)+90], 'm-', 'LineWidth', 1);
           end
        end
    pause(0.1);
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
end


 close(writerObj);
 
 % ffmpeg -i mymovie.avi -i combo.wav -c:v copy output.avi
 % ffmpeg -f concat -i mylist.txt -c copy combination.mp4
 % ffmpeg -i combination.mp4 -filter:v "setpts=2.0*PTS" -filter:a "atempo=0.5" slow.mp4
 % ffmpeg -i output2.avi -strict -2 -pix_fmt yuv420p m17spikes.mp4

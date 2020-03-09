%% Make the audio file
idx = [1, 2] ; % This is the Male (odd)
rango = [-2, 7];

    specpos = 0;
    if rango(1) < 0; specpos = abs(rango(1)); end

Fs = w(idx).Fs;

% Make the Fake Spikes

    spiketim = 1/Fs:1/Fs:0.005; % 5 msec duration for our fake spikes
    len = length(spiketim);

    fspike = sin(2*pi*4000*spiketim) * 0.2 ; % 4 kHz for females
    mspike = sin(2*pi*3000*spiketim) * 0.2 ; % 3 kHz for males
    
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
        mal(curidx:curidx+len-1) = mspike;
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
audiowrite('combo.wav', combo, 10000);

%% Make the video

    vFs = 30; % frames per second

figure(1); clf; 
    specgram(outduet, 512, Fs);
    colormap('HOT');
    % truesize([960 1280]); 
    truesize([480 640]);
    % specgram(outduet, 512, Fs);
    caxis([-25 25])
    cmp = flipud(gray);
    %colormap(cmp);
    hold on;
    
% % Initialize the "object" that will be the final movie
%     writerObj = VideoWriter('mymovie.avi');
%     writerObj.FrameRate = 30;
% %writerObj.VideoFormat = 'RGB24';
%     writerObj.Quality = 90;
%     open(writerObj);

for vtim = outim(1):1:outim(end)
    clf; 
    specgram(outduet, 512, Fs, [], round(512*0.95));
    colormap('HOT');
    caxis([-25 25])
    hold on;
    
    plot([specpos+vtim,specpos+vtim], [0, 5000], 'r-', 'LineWidth', 5); % Red progress line
    
        for j = 1:4
           malespkidx = find(w(idx(1)).Cspikes{j} > rango(1) & w(idx(1)).Cspikes{j} < vtim);
           femalespkidx = find(w(idx(2)).Cspikes{j} > rango(1) & w(idx(2)).Cspikes{j} < vtim);
           
           for k = 1:length(malespkidx)
               plot([specpos+w(idx(1)).Cspikes{j}(malespkidx(k)), specpos+w(idx(1)).Cspikes{j}(malespkidx(k))], [4000+(j*100), 4000+(j*100)+((j-1)*110)], 'b-', 'LineWidth', 1);
           end
           for k = 1:length(femalespkidx)
               plot([specpos+w(idx(2)).Cspikes{j}(femalespkidx(k)), specpos+w(idx(2)).Cspikes{j}(femalespkidx(k))], [100+(j*100), 100+(j*100)+((j-1)*110)], 'm-', 'LineWidth', 1);
           end
        end
    
%    frame = getframe(gcf);
%    writeVideo(writerObj, frame);
end


% close(writerObj);

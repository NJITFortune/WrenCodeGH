
idx = 1; % This is the Male (odd)
rango = [-2, 7];
Fs = w(idx).Fs;

% Make the Fake Spikes

    spiketim = 1/Fs:1/Fs:0.005; % 5 msec duration for our fake spikes
    len = length(spiketim);

    fspike = sin(2*pi*4000*spiketim) * 0.2 ; % 4 kHz for females
    mspike = sin(2*pi*3000*spiketim) * 0.2 ; % 3 kHz for males
    
    outim = w(idx).tim(w(idx).tim > rango(1) & w(idx).tim < rango(2));
    outduet = w(idx).duet(w(idx).tim > rango(1) & w(idx).tim < rango(2));
    
% Make the spike sound trains
    fem = zeros(1,length(outim));
    mal = zeros(1,length(outim));

for k = 1:4
    spkidx = find(w(idx).Cspikes{k} > rango(1) & w(idx).Cspikes{k} < rango(2));
    for j = 1:length(spkidx)   
        curidx = find(outim >= w(idx).Cspikes{k}(spkidx(j)), 1, 'first');    
    if curidx+len-1 < length(outim)  % Need the if not to go over the end.   
        mal(curidx:curidx+len-1) = mspike;
    end
    
    end
end

idx = 2; % This is the female (even)

for k = 1:4
    spkidx = find(w(idx).Cspikes{k} > rango(1) & w(idx).Cspikes{k} < rango(2));
    for j = 1:length(spkidx)   
        curidx = find(outim >= w(idx).Cspikes{k}(spkidx(j)), 1, 'first');    
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
    colormap(cmp);
    hold on;
    
% Initialize the "object" that will be the final movie
    writerObj = VideoWriter('mymovie.avi');
    writerObj.FrameRate = 30;
%writerObj.VideoFormat = 'RGB24';
    writerObj.Quality = 90;
    open(writerObj);

for i = outim(1):outim/30:outtim(end)
    plot(
    hold on; 
    plot(data.x(i), data.y(i), 'm*', data.x(i-90:i), data.y(i-90:i), 'm-'); 
    hold off;
%    frame = getframe(gcf);
%    writeVideo(writerObj, frame);
end


close(writerObj);

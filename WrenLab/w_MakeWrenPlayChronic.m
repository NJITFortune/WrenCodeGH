
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

figure(1); clf; hold on;
    specgram(outduet, 512, Fs);
    
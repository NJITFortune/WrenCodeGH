
idx = 1; % This is the Male (odd)

Fs = w(idx).Fs;

spiketim = 1/Fs:1/Fs:0.005; % 5 msec duration for our fake spikes

fspike = sin(2*pi*4000*spiketim) * 0.2 ; % 4 kHz for females
mspike = sin(2*pi*3000*spiketim) * 0.2 ; % 3 kHz for males


fem = zeros(1,length(w(idx).tim));
mal = zeros(1,length(w(idx).tim));

len = length(mspike);

for k = 1:4
for j = 1:length(w(idx).Cspikes{k})
   
    curidx = find(w(idx).tim >= w(idx).Cspikes{k}(j), 1, 'first');
    
    if curidx+len-1 < length(w(idx).tim)
         
        mal(curidx:curidx+len-1) = mspike;
        
    end
    
end
end

idx = 2;

for k = 1:4
for j = 1:length(w(idx).Cspikes{k})
   
    curidx = find(w(idx).tim >= w(idx).Cspikes{k}(j), 1, 'first');
    
    if curidx+len-1 < length(w(idx).tim)
        
        fem(curidx:curidx+len-1) = fspike;
        
    end
    
end
end

spks = fem' + mal';

combo = spks + (w(idx).duet*0.8);
audiowrite('combo.wav', combo, 10000);

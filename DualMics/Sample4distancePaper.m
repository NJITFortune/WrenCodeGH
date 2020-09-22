% load('/Users/eric/Sync/Wren/DistanceData-V10.mat')


figure(1); clf
subplot(211); 
    specgram(dd(12).femMic, 1024, dd(12).Fs, [], 1000);
    yarg = flipud(gray);
    colormap(yarg);
    ylim([0 5000]);

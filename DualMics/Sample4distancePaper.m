% load('/Users/eric/Sync/Wren/DistanceData-V10.mat')

specgram(dd(12).femMic, 512, dd(12).Fs);
subplot(211); specgram(dd(12).femMic, 512, dd(12).Fs);

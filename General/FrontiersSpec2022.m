% Spectrogram figure for Frontiers Review 2022

[s, Fs] = audioread('/Users/eric/SparkleShareOLD/WrenMS/Science2011/PlainTail_DuetSwitch.wav');

%%
tt = [34.9 42.82];
figure(1); specgram(s(tim > tt(1) & tim < tt(2)), 2048, Fs, [], 2000);
figure(1); colormap(flipud(gray)); caxis([-10 30]); ylim([500 4200])
set(gcf, 'Renderer', 'painters');

%% 

load /Users/eric/Sync/Wren/ChronicCompleat2021CodesG.mat

%%




% Spectrogram figure for Frontiers Review 2022

[s, Fs] = audioread('/Users/eric/SparkleShareOLD/WrenMS/Science2011/PlainTail_DuetSwitch.wav');

%%
tim = 1/Fs:1/Fs:length(s)/Fs;
tt = [34.9 42.82];
figure(1); specgram(s(tim > tt(1) & tim < tt(2)), 2048, Fs, [], 2000);
figure(1); colormap(flipud(gray)); caxis([-10 30]); ylim([500 4200])
set(gcf, 'Renderer', 'painters');

%% 

load /Users/eric/Sync/Wren/ChronicCompleat2021CodesG.mat

%%

spectim = 1/w(11).Fs:1/w(11).Fs:length(w(11).duet)/w(11).Fs;
tt = find(spectim > 3.6 & spectim < 10.6);

specgram(w(11).duet(tt), 512, w(11).Fs, [], 500); 
    ylim([500 4200]); colormap(flipud(gray)); 
    caxis([-13 22]);

set(gcf, 'Renderer', 'painters');

%%
tim = 1/Fs:1/Fs:length(s)/Fs;
tt = [0 10];
figure(1); specgram(s(tim > tt(1) & tim < tt(2)), 2048, Fs, [], 2000);
figure(1); colormap(flipud(gray)); caxis([-10 30]); ylim([500 4200])
set(gcf, 'Renderer', 'painters');

%%
[k, Fs] = audioread('~/Downloads/plaintail-1.mp3');
%%
tim = 1/Fs:1/Fs:length(k)/Fs;
tt = [0 10];
figure(1); specgram(k(tim > tt(1) & tim < tt(2)), 2048, Fs, [], 2000);
figure(1); colormap(flipud(gray)); caxis([-10 30]); ylim([500 4200])
set(gcf, 'Renderer', 'painters');



figure(2); specgram(newK(tim < 6), 2048, Fs, [], 2000);
figure(2); colormap(flipud(gray)); caxis([-10 30]); ylim([500 4200])
set(gcf, 'Renderer', 'painters');

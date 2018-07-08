

% [sng, Fs] = audioread('filename.wav');
% load filename.mat

% Plot the analysis on top of the spectrogram before we start just for kicks.
numcolors = 20;
oclrs = hsv(numcolors); % 20 colors for plotting up to 20 different fishes
clrs = zeros(numcolors, 3);
% figure(3); clf; hold on; for kk=1:20; plot(kk,kk, '*', 'MarkerEdgeColor', oclrs(kk,:)); end;
    shuff = [11 2 7 13 18 1 16 12 10 17 3 8 14 19 20 4 16 5 15 9];
    for i=1:20; clrs(i,:) = oclrs(shuff(i),:); end; 

figure(2); clf;

specgram(sng, 1024, Fs, [], 1000); ylim([200 5500]);
caxis([0 20]); colormap('HOT');

hold on;

for i = 1:length(asdf)
    
    plot(asdf(i).trace_tim + asdf(i).syltim(1), asdf(i).trace_freq, 'g', 'LineWidth', 2);

    
end
figure(1); clf; 
qwer = slicer(asdf);

figure(3); clf;

specgram(sng, 1024, Fs, [], 1000); ylim([200 5500]);
caxis([-10 20]); colormap('GRAY');

hold on;


for i = 1:length(qwer) % This is for each syllable type
    
    for j = 1:length(qwer(i).num)
        plot(asdf(qwer(i).num(j)).trace_tim + asdf(qwer(i).num(j)).syltim(1), asdf(qwer(i).num(j)).trace_freq, 'Color', clrs(i,:), 'LineWidth', 2);

    end
end

function out = sylsorter(clks, wav, Fs)
% Usage: out = sylsorter(clickdata, wavdata, samplerate)

numcolors = 50;
    clrs = hsv(numcolors); % 50 colors for plotting up to 50 syllables
    arandperm = [49 33 41 26 46 18 42 4 27 10 24 28 22 12 20 1 36 50 8 5 31 39 48 30 6 13 43 17 34 40 47 19 9 29 3 2 14 11 38 37 44 23 32 21 15 35 25 45 16 7];
    clrs = clrs(arandperm,:);
    
out = wrencleaner(clks);

%% Plot the data over the specgram with numbers

figure(1); clf; hold on;
specgram(wav, 1024, Fs, [], 1000); ylim([500 5500]); colormap('GRAY'); caxis([-10 20]); 

for k = 1:length(clks)
   
    plot(clks(k).trace_tim + clks(k).syltim(1), out(k).trace_freq, 'Color', clrs(k,:), 'LineWidth', 2);
    text(mean([clks(k).syltim]), 5100, int2str(k), 'FontSize', 14, 'Color', 'y');
    text(mean([clks(k).syltim]), 600, int2str(k), 'FontSize', 14, 'Color', 'y');
    
end

%% Run Slicer to get the different syllables and replot

    clusts = slicer(out);
    
    figure(1); clf; hold on;
    specgram(wav, 1024, Fs, [], 1000); ylim([500 5500]); colormap('GRAY'); caxis([-10 20]); 

for cc = 1:length(clusts)
   
    for ss = 1:length(clusts(cc).num)
        idx = clusts(cc).num(ss);
        plot(clks(idx).trace_tim + clks(idx).syltim(1), out(idx).trace_freq, 'Color', clrs(cc,:), 'LineWidth', 2);
        text(mean([clks(idx).syltim]), 5100, int2str(cc), 'FontSize', 14, 'Color', clrs(cc,:));
        text(mean([clks(idx).syltim]), 600, int2str(idx), 'FontSize', 14, 'Color', 'y');
    end
end

%% Check and resolve errors 

numsyls = length(clusts);

    for j=1:numsyls
        clustlen(j) = unique(diff(clusts(j).num));
    	if length(clustlen(j)) == 1
            fprintf('Syllable %i OK.\n', j); pause(0.1);
        end;
    end

    if length(unique(clustlen)) == 1 
        fprintf('\n Looks clean with %i syllables. \n', unique(clustlen));
    end
    
    
% load DistanceData-V10.mat

idx = 1;

%% Make the primary plot (4 windows)

figure(1); clf;

tim = 1/dd(idx).Fs:1/dd(1).Fs:length(dd(idx).maleMic)/dd(idx).Fs;

ax(1) = subplot(411); hold on;
specgram(dd(idx).femMic, 1024, dd(idx).Fs, [], 1000);
colormap('HOT'); caxis([-50 30]);
ylim([0 8000]);

ax(2) = subplot(412); hold on;
plot(tim, dd(idx).femMic, 'm-');

ax(3) = subplot(413); hold on;
specgram(dd(idx).maleMic, 1024, dd(idx).Fs, [], 1000);
colormap('HOT'); caxis([-50 30]);
ylim([0 8000]);

ax(4) = subplot(414); hold on;
plot(tim, dd(idx).maleMic, 'b-');
linkaxes(ax, 'x');

%% Cycle to plot each to show each syllable


% Female Microphone first

for j = 1:length(dd(idx).fsyl)
    
    ax(1) = subplot(411); 
        plot([dd(idx).fsyl(j).syltim(1), dd(idx).fsyl(j).syltim(1)], [0 8000], '-g', 'LineWidth', 2);
        plot([dd(idx).fsyl(j).syltim(2), dd(idx).fsyl(j).syltim(2)], [0 8000], '-r', 'LineWidth', 2);
        if dd(idx).fsyl(j).sexsyltype > 25; colr = 'Magenta'; else; colr = 'Blue'; end
        text(dd(idx).fsyl(j).syltim(1)+0.02, 7500, num2str(dd(idx).fsyl(j).sexsyltype), 'Color', 'White', 'FontWeight', 'bold', 'FontSize', 14);   
        text(dd(idx).fsyl(j).syltim(1)+0.02, 7500, num2str(dd(idx).fsyl(j).sexsyltype), 'Color', colr, 'FontSize', 14);   
    ax(2) = subplot(412);
        plot([dd(idx).fsyl(j).syltim(1), dd(idx).fsyl(j).syltim(1)], [-1 1], '-g', 'LineWidth', 2);
        plot([dd(idx).fsyl(j).syltim(2), dd(idx).fsyl(j).syltim(2)], [-1 1], '-r', 'LineWidth', 2);
        text(dd(idx).fsyl(j).syltim(1)+0.02, 0.80, num2str(dd(idx).fsyl(j).sexsyltype), 'Color', colr);   
    
end

% Now the male Microphone 

for j = 1:length(dd(idx).msyl)
    
    ax(3) = subplot(413); 
        plot([dd(idx).msyl(j).syltim(1), dd(idx).msyl(j).syltim(1)], [0 8000], '-g', 'LineWidth', 2);
        plot([dd(idx).msyl(j).syltim(2), dd(idx).msyl(j).syltim(2)], [0 8000], '-r', 'LineWidth', 2);
        if dd(idx).msyl(j).sexsyltype > 25; colr = 'Magenta'; else; colr = 'Blue'; end
        text(dd(idx).msyl(j).syltim(1)+0.02, 7500, num2str(dd(idx).msyl(j).sexsyltype), 'Color', 'White', 'FontWeight', 'bold', 'FontSize', 14);   
        text(dd(idx).msyl(j).syltim(1)+0.02, 7500, num2str(dd(idx).msyl(j).sexsyltype), 'Color', colr, 'FontSize', 14);   
    ax(4) = subplot(414);
        plot([dd(idx).msyl(j).syltim(1), dd(idx).msyl(j).syltim(1)], [-1 1], '-g', 'LineWidth', 2);
        plot([dd(idx).msyl(j).syltim(2), dd(idx).msyl(j).syltim(2)], [-1 1], '-r', 'LineWidth', 2);
        text(dd(idx).msyl(j).syltim(1)+0.02, 0.80, num2str(dd(idx).fsyl(j).sexsyltype), 'Color', colr);   
    
end




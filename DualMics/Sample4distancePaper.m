% load('/Users/eric/Sync/Wren/DistanceData-V10.mat')


figure(1); clf
subplot(211); 
    specgram(dd(12).femMic, 1024, dd(12).Fs, [], 1000);
    yarg = flipud(gray);
    colormap(yarg);
    ylim([0 5000]);
    xlim([0.6 5.6]);
    caxis([-30 35]);
    
hold on;
    
 for j=1:length(dd(12).fsyl)
    
     plot([dd(12).fsyl(j).syltim(1), dd(12).fsyl(j).syltim(2)], [100, 100], 'm-', 'LineWidth', 8);
     
 end
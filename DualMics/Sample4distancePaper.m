% load('/Users/eric/Sync/Wren/DistanceData-V10.mat')


figure(1); clf;
subplot(211); 
    specgram(dd(12).femMic, 1024, dd(12).Fs, [], 1000);
    yarg = flipud(gray);
    colormap(yarg);
    ylim([250 4500]);
    xlim([0.6 5.6]);
    caxis([-30 35]);
    
hold on;
    
fidx = find([dd(12).fsyl.sexsyltype] > 50);
midx = find([dd(12).fsyl.sexsyltype] < 50);

 for j=1:length(fidx)
    
     plot([dd(12).fsyl(fidx(j)).syltim(1), dd(12).fsyl(fidx(j)).syltim(2)], [100, 100], 'm-', 'LineWidth', 4);
     
 end
 
  for j=1:length(midx)
    
     plot([dd(12).fsyl(midx(j)).syltim(1), dd(12).fsyl(midx(j)).syltim(2)], [100, 100], 'b-', 'LineWidth', 4);
     
 end
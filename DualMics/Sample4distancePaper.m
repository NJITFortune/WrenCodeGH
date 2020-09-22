% load('/Users/eric/Sync/Wren/DistanceData-V10.mat')

[b,a] = butter(5, 500/(dd(12).Fs/2), 'high');
figure(1); clf;
subplot(211); 
    specgram(filtfilt(b,a,dd(12).femMic), 1024, dd(12).Fs, [], 1000);
    yarg = flipud(gray);
    colormap(yarg);
    ylim([0 4500]);
    xlim([0.6 5.6]);
    caxis([-30 35]);
    
hold on;
    
fidx = find([dd(12).fsyl.sexsyltype] > 50);
midx = find([dd(12).fsyl.sexsyltype] < 50);

 for j=1:length(fidx)
    
     plot([dd(12).fsyl(fidx(j)).syltim(1), dd(12).fsyl(fidx(j)).syltim(2)], [200, 200], 'm-', 'LineWidth', 4);
     plot([dd(12).fsyl(fidx(j)).syltim(1), dd(12).fsyl(fidx(j)).syltim(1)], [200, 4300], 'm-', 'LineWidth', 0.5);
     plot([dd(12).fsyl(fidx(j)).syltim(2), dd(12).fsyl(fidx(j)).syltim(2)], [200, 4300], 'm-', 'LineWidth', 0.5);
     
 end
 
  for j=1:length(midx)
    
     plot([dd(12).fsyl(midx(j)).syltim(1), dd(12).fsyl(midx(j)).syltim(2)], [200, 200], 'b-', 'LineWidth', 4);
     plot([dd(12).fsyl(midx(j)).syltim(1), dd(12).fsyl(midx(j)).syltim(1)], [200, 4300], 'b-', 'LineWidth', 0.5);
     plot([dd(12).fsyl(midx(j)).syltim(2), dd(12).fsyl(midx(j)).syltim(2)], [200, 4300], 'b-', 'LineWidth', 0.5);
     
 end
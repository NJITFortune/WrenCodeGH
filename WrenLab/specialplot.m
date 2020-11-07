
% 
subplot(3,1,1)
    specgram(w(13).duet, 1024, w(1).Fs); ylim([0 5000]);
    colormap(flipud(gray));
    caxis([-20 20]);
subplot(3,1,1)

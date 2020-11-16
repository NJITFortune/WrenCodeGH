% This plots the spectrogram and oscillogram in Figure 1 of the 
% Manuscript. It highlights the relations between solo and duet singing
% In plain-tailed wrens.

% Construct a time sequence for the signal
    tim = 1/w(13).Fs:1/w(13).Fs:length(w(13).duet)/w(13).Fs;

%% Row one - full specgram
subplot(3,1,1)
    specgram(w(13).duet, 1024, w(1).Fs); ylim([600 4300]);
    colormap(flipud(gray));
    caxis([-20 20]);
    
%% Rows two (spectrograms) and three (oscillograms) closeups (6-second duration clips)

% Pick first range
tt = find(tim > 0.6 & tim < 6.6);

    subplot(3,3,4);
    specgram(w(13).duet(tt), 1024, w(1).Fs); ylim([600 4300]);
    colormap(flipud(gray)); caxis([-20 20]);

    subplot(3,3,7);
    plot(tim(tt),w(13).duet(tt)); 
        xlim([tim(tt(1)), tim(tt(end))]); ylim([-1 1]);
    
% Pick second range
tt = find(tim > 26 & tim < 32);

    subplot(3,3,5); 
    specgram(w(13).duet(tt), 1024, w(1).Fs); ylim([600 4300]);
    colormap(flipud(gray)); caxis([-20 20]); set(gca, 'yticklabel', []);
    
    subplot(3,3,8);
    plot(tim(tt),w(13).duet(tt)); 
        xlim([tim(tt(1)), tim(tt(end))]); ylim([-1 1]);   
    
% Pick third range
tt = find(tim > 33 & tim < 39);

    subplot(3,3,6);
    specgram(w(13).duet(tt), 1024, w(1).Fs); ylim([600 4300]);
    colormap(flipud(gray)); caxis([-20 20]); set(gca, 'yticklabel', []);
    
    subplot(3,3,9);
    plot(tim(tt),w(13).duet(tt)); 
        xlim([tim(tt(1)), tim(tt(end))]); ylim([-1 1]); 
        
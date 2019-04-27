% load ChronicCompleat2019f.mat

sig = w(1).duet(w(1).tim > -9.5 & w(1).tim < -4.5)';
sig = [sig, w(1).duet(w(1).tim > 0.5 & w(1).tim < 6)'];
Fs = 10000;
timtim = 1/Fs:1/Fs:length(sig)/Fs;


figure(1); clf;
ax(1) = subplot(211);
%    specgram(w(1).duet, 512, 10000, [], 500); 
    specgram(sig, 512, 10000, [], 500);
    Yarg = flipud(gray);
    colormap(Yarg); caxis([-35 30]); ylim([0 4200]);
ax(2) = subplot(212); 
%    plot(w(1).tim(1:5:end) - w(1).tim(1), w(1).duet(1:5:end), 'k'); ylim([-1.1, 1.1]);
    plot(timtim, sig, 'k'); ylim([-1.1, 1.1]);
    hold on;
    
for j=1:3 %length(w(1).syl)

    if w(1).sylsex(j) == 2
        plot([w(1).syl(j).tim(1) - w(1).tim(1) - 0.5, w(1).syl(j).tim(2) - w(1).tim(1) - 0.5], [1, 1], 'm-', 'LineWidth', 2);
%        plot([w(1).syl(j).tim(1) - w(1).tim(1), w(1).syl(j).tim(2) - w(1).tim(1)], [1, 1], 'm-', 'LineWidth', 2);
    elseif w(1).sylsex(j) == 1
        plot([w(1).syl(j).tim(1) - w(1).tim(1) - 0.5, w(1).syl(j).tim(2) - w(1).tim(1) - 0.5], [1, 1], 'b-', 'LineWidth', 2);
%        plot([w(1).syl(j).tim(1) - w(1).tim(1), w(1).syl(j).tim(2) - w(1).tim(1)], [1, 1], 'b-', 'LineWidth', 2);
    end
    
end
for j=4:length(w(1).syl)

    if w(1).sylsex(j) == 2
        plot([w(1).syl(j).tim(1) - w(1).tim(1) - 5.5, w(1).syl(j).tim(2) - w(1).tim(1) - 5.5], [1, 1], 'm-', 'LineWidth', 2);
%        plot([w(1).syl(j).tim(1) - w(1).tim(1), w(1).syl(j).tim(2) - w(1).tim(1)], [1, 1], 'm-', 'LineWidth', 2);
    elseif w(1).sylsex(j) == 1
        plot([w(1).syl(j).tim(1) - w(1).tim(1) - 5.5, w(1).syl(j).tim(2) - w(1).tim(1) - 5.5], [1, 1], 'b-', 'LineWidth', 2);
%        plot([w(1).syl(j).tim(1) - w(1).tim(1), w(1).syl(j).tim(2) - w(1).tim(1)], [1, 1], 'b-', 'LineWidth', 2);
    end
    
end
    
linkaxes(ax, 'x');
xlim([0 10.5]);

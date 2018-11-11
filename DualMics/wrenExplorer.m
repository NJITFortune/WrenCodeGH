function wrenExplorer(duet, idx)
% Plot the spectrograms for each microphone with lines indicating starts 
% and ends for each syllable with other information.

figure(1); clf; 

Fs = duet(idx).Fs;
prpspc = flipud(gray);
plotlims = [-10 30];

ax(1) = subplot(211); hold on; 
    specgram(duet(idx).femMic, 1024, Fs, [], 1000); ylim([0 6000]); colormap(prpspc); caxis(plotlims);
    text(0.1, 5000, 'Female Microphone', 'Color', 'm');
    for j=1:length(duet(idx).fsyl) 
        if duet(idx).fsyl(j).sexsyltype < 49 % Male Syllable
           plot([duet(idx).fsyl(j).syltim(1) duet(idx).fsyl(j).syltim(2)], [5500, 5500], 'b', 'LineWidth', 2);
           text(mean([duet(idx).fsyl(j).syltim(1) duet(idx).fsyl(j).syltim(2)]), 5500, num2str(duet(idx).fsyl(j).sexsyltype), 'Color', 'k');
        else
           plot([duet(idx).fsyl(j).syltim(1) duet(idx).fsyl(j).syltim(2)], [500, 500], 'm', 'LineWidth', 2);
           text(mean([duet(idx).fsyl(j).syltim(1) duet(idx).fsyl(j).syltim(2)]), 500, num2str(duet(idx).fsyl(j).sexsyltype-50), 'Color', 'k');            
        end
        plot([duet(idx).fsyl(j).syltim(1) duet(idx).fsyl(j).syltim(1)], [500 5500], 'g', 'LineWidth', 1);
        plot([duet(idx).fsyl(j).syltim(2) duet(idx).fsyl(j).syltim(2)], [500 5500], 'r', 'LineWidth', 1);
    end

ax(2) = subplot(212); hold on;
    specgram(duet(idx).maleMic, 1024, Fs, [], 1000); ylim([0 6000]); colormap(prpspc); caxis(plotlims);
    text(0.1, 5000, 'Male Microphone', 'Color', 'b');
    for j=1:length(duet(idx).msyl) 
        if duet(idx).msyl(j).sexsyltype < 49 % Male Syllable
           plot([duet(idx).msyl(j).syltim(1) duet(idx).msyl(j).syltim(2)], [500, 500], 'b', 'LineWidth', 2);
           text(mean([duet(idx).msyl(j).syltim(1) duet(idx).msyl(j).syltim(2)]), 200, num2str(duet(idx).msyl(j).sexsyltype), 'Color', 'k', 'HorizontalAlignment', 'center');
        else
           plot([duet(idx).msyl(j).syltim(1) duet(idx).msyl(j).syltim(2)], [5500, 5500], 'm', 'LineWidth', 2);
           text(mean([duet(idx).msyl(j).syltim(1) duet(idx).msyl(j).syltim(2)]), 5200, num2str(duet(idx).msyl(j).sexsyltype-50), 'Color', 'k', 'HorizontalAlignment', 'center');            
        end
        plot([duet(idx).msyl(j).syltim(1) duet(idx).msyl(j).syltim(1)], [500 5500], 'g', 'LineWidth', 1);
        plot([duet(idx).msyl(j).syltim(2) duet(idx).msyl(j).syltim(2)], [500 5500], 'r', 'LineWidth', 1);
    end

linkaxes(ax, 'xy'); xlim([0 duet(idx).msyl(end).syltim(2)+0.5]);

%     %% Get lists of syllables
%     
%     malsyls = [];
%     for j=1:length(duet(1).anal.maltype)
%         malsyls = [malsyls duet(1).anal.maltype{j}];
%     end
%     malsyls = sort(malsyls);
%     
%     femsyls = [];
%     for j=1:length(duet(1).anal.femtype)
%         femsyls = [femsyls duet(1).anal.femtype{j}];
%     end
%     femsyls = sort(femsyls);
%     
%     speedofsound = 0.0033;
%     
%     distance = 3 * 2;
%     
%     sndelay = distance * speedofsound;
%     
%     
%         for j=1:length(malsyls)
%             subplot(211); 
%             %plot([in(1).m(malsyls(j)).syltim(1)+sndelay, in(1).m(malsyls(j)).syltim(1)+sndelay], [1000 4000], 'k', 'LineWidth', 1);
%             %plot([in(1).m(malsyls(j)).syltim(2)+sndelay, in(1).m(malsyls(j)).syltim(2)+sndelay], [1000 4000], 'k', 'LineWidth', 1);
%             plot(duet(1).f(malsyls(j)).trace_tim + duet(1).f(malsyls(j)).syltim(1), duet(1).f(malsyls(j)).trace_freq, 'k');
%             subplot(212);
%             plot(duet(1).m(malsyls(j)).trace_tim + duet(1).m(malsyls(j)).syltim(1), duet(1).m(malsyls(j)).trace_freq, 'k');
%         end;
%     
%         for j=1:length(femsyls)
%             subplot(212); 
%             %plot([in(1).f(femsyls(j)).syltim(1)+0, in(1).f(femsyls(j)).syltim(1)+0], [1000 4000], 'k', 'LineWidth', 1);
%             %plot([in(1).f(femsyls(j)).syltim(2)+0, in(1).f(femsyls(j)).syltim(2)+0], [1000 4000], 'k', 'LineWidth', 1);
%             subplot(211);
%             %plot(in(1).f(femsyls(j)).trace_tim + in(1).f(femsyls(j)).syltim(1), in(1).f(femsyls(j)).trace_freq, 'k');
%         end;
%         
%         linkaxes(ax, 'x');
%         
%         
%         
% %% Pre- Post- Analysis
% 
% % MALE SYLLABLES (use recordings aligned to female syllables)
% 
% for pp = 1:5
%     Mpre{pp} = []; Mpost{pp} = [];
% end
% 
% for j=1:length(duet) % For each duet
% 
%     numsyltypes = length(duet(j).anal.maltype); % Get number of syllables for the male in the duet
%         malsyls = [];
%         for jj=1:length(duet(1).anal.maltype)
%             malsyls = [malsyls duet(1).anal.maltype{jj}];
%         end
%     malsyls = sort(malsyls);
% 
%         femsyls = [];
%         for jj=1:length(duet(1).anal.femtype)
%             femsyls = [femsyls duet(1).anal.femtype{jj}];
%         end
%     femsyls = sort(femsyls);
%     
%     for k=1:numsyltypes % For each syllable type
%        
%         for i=1:length(duet(j).anal.maltype{k}) % For each syllable of that type
%         
%             currsyllableIDX = duet(j).anal.maltype{k}(i);
% 
%             if find(femsyls == currsyllableIDX-1)
%                 endofpresyllable = duet(j).f(currsyllableIDX-1).syltim(2);
%                 Mpre{k}(end+1) = duet(j).m(currsyllableIDX).syltim(1) - endofpresyllable;
%             end
% 
%             if find(femsyls == currsyllableIDX+1)
%                 startofpostsyllable = duet(j).f(currsyllableIDX+1).syltim(1);
%                 Mpost{k}(end+1) = startofpostsyllable - duet(j).m(currsyllableIDX).syltim(2);
%             end
% 
%         end
%         
%     end
% 
% end
% 
% out.Mpre = Mpre;
% out.Mpost = Mpost;

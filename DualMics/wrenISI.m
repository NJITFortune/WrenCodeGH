function out = wrenISI(in)
% Usage: out = dd_getISI(in)
% 'in' is a structure for two-microphone recordings of wrens
% 'in' is composed of 
% day month year        Date of the recording 
% Location              Name of the study site
% vision                1 Visual cues present, 0 Visual cues absent
% distance              Distance between the two cages in meters  
% femMic                Audio recording from microphone at female cage
% maleMic               Audio recording from microphone at male cage
% Fs                    Sample rate of audio recordings
% fsyl msyl             Structure with syllable data
%   fftFreqs            Frequencies for plotting with fftPower
%   fftPower            Power at frequencies in fftFreqs for the syllable
%   sex                 'F' for female syllable, 'M' for male
%   sexsyltype          Numeric identifiers for individual syllables, males 1-20, females 50-70   
%   


% out includes the following:
% out.Fmf = []; Female microphone, ISI male-to-female 
% out.Fmfd = []; Female microphone, distance for current ISI male-to-female
% out.Ffm = []; Female microphone, ISI female-to-male
% out.Ffmd = []; Female microphone, distance for current ISI female-to-male
% out.Mmf = []; Male microphone, ISI male-to-female
% out.Mmfd = []; Male microphone, distance for current ISI male-to-female
% out.Mfm = []; Male microphone, ISI female-to-male
% out.Mfmd = []; Male microphone, distance for current ISI female-to-male

out.Fmf = []; out.Fmfd = [];
out.Ffm = []; out.Ffmd = [];
out.Mmf = []; out.Mmfd = [];
out.Mfm = []; out.Mfmd = [];

% out.mm = []; out.ff = [];

spdosnd = 1/331.2; % Speed of sound is 331.2 meters per second


%% Cycle through every duet in the structure to retrieve ISIs
for d = 1:length(in)
    
    numsyls = length(in(d).fsyl);
    
    % Cycle through every syllable except the last (avoid last because it
    % has no inter-syllable-interval.
    
    for s=1:numsyls-1
        
        % Get the ISI for the current syllable on each microphone
        currFisi = in(d).fsyl(s+1).syltim(1) - in(d).fsyl(s).syltim(2);
        currMisi = in(d).msyl(s+1).syltim(1) - in(d).msyl(s).syltim(2);
            
        if in(d).fsyl(s).sexsyltype < 49 && in(d).fsyl(s+1).sexsyltype > 49 % Male to Female
                out.Fmf(end+1) = currFisi;
                out.Fmfd(end+1) = in(d).distance;
                out.Mmf(end+1) = currMisi;
                out.Mmfd(end+1) = in(d).distance;
        end
        if in(d).fsyl(s).sexsyltype > 49 && in(d).fsyl(s+1).sexsyltype < 49 % Female to Male
                out.Ffm(end+1) = currFisi;
                out.Ffmd(end+1) = in(d).distance;
                out.Mfm(end+1) = currMisi;
                out.Mfmd(end+1) = in(d).distance;
        end        
                        
    end
                    
end

%% Plot the raw data
figure(1); clf; subplot(211); hold on; subplot(212); hold on;

figure(1); subplot(211); plot([0 8], [0 0], 'k-'); ylim([-0.05 0.25]);
figure(1); subplot(212); plot([0 8], [0 0], 'k-'); ylim([-0.05 0.25]);
            subplot(211); text(3, 0.2, 'Female Microphone', 'Color', 'm');
            subplot(211); plot(out.Fmfd+0.1, out.Fmf, 'mo', 'MarkerSize', 10);
            subplot(211); plot(out.Ffmd-0.1, out.Ffm, 'bd', 'MarkerSize', 6);
            subplot(212); text(3, 0.2, 'Male Microphone', 'Color', 'b');
            subplot(212); plot(out.Mmfd-0.1, out.Mmf, 'md', 'MarkerSize', 6);
            subplot(212); plot(out.Mfmd+0.1, out.Mfm, 'bo', 'MarkerSize', 10);


%% Perform statistics

distances = sort(unique([in.distance]));

for jj = 1:length(distances)
    
    Ffm(jj) = mean(out.Ffm([out.Ffmd] == distances(jj)));
    Fmf(jj) = mean(out.Fmf([out.Fmfd] == distances(jj)));
    Ffmstd(jj) = std(out.Ffm([out.Ffmd] == distances(jj)));
    Fmfstd(jj) = std(out.Fmf([out.Fmfd] == distances(jj)));

    Mfm(jj) = mean(out.Mfm([out.Ffmd] == distances(jj)));
    Mmf(jj) = mean(out.Mmf([out.Fmfd] == distances(jj)));
    Mfmstd(jj) = std(out.Mfm([out.Ffmd] == distances(jj)));
    Mmfstd(jj) = std(out.Mmf([out.Fmfd] == distances(jj)));

end


%% Plot the analyzed data
figure(2); clf; subplot(211); hold on; subplot(212); hold on;
ax(1) = subplot(211); 
    plot(distances-0.1, Ffm(1)+(2*(distances-1)*spdosnd), 'k-');
    errorbar(distances-0.1, Ffm, Ffmstd, 'db');
    errorbar(distances+0.1, Fmf, Fmfstd, 'om', 'LineWidth', 2);
    text(3, 0.15, 'Female Microphone', 'Color', 'm');
    text(3, 0.0, 'Expected effect of distance', 'Color', 'k');

ax(2) = subplot(212); 
    plot(distances-0.1, Mmf(1)+(2*(distances-1)*spdosnd), 'k-');
    errorbar(distances+0.1, Mfm, Ffmstd, 'ob', 'LineWidth', 2);
    errorbar(distances-0.1, Mmf, Fmfstd, 'dm');
    text(3, 0.15, 'Male Microphone', 'Color', 'b');
    text(3, 0.0, 'Expected effect of distance', 'Color', 'k');

linkaxes(ax, 'xy'); xlim([0 8]); ylim([-0.02 0.18]);
    
%    subplot(211); errorbar(jj, mean(out.Ffm([out.Ffmd] == jj)) , std(out.Ffm([out.Ffmd] == jj)), 'om');
%    subplot(212); errorbar(jj, mean(out.Fmf([out.Fmfd] == jj)) , std(out.Fmf([out.Fmfd] == jj)), '*m');
%    subplot(211); plot(jj, mean(out.Ffm([out.Ffmd] == jj)),'mo');
%    subplot(212); plot(jj, mean(out.Fmf([out.Ffmd] == jj)),'m*');

%    subplot(211); errorbar(jj+0.2, mean(out.Mmf([out.Mmfd] == jj)) , std(out.Mmf([out.Mmfd] == jj)), 'ob');
%    subplot(212); errorbar(jj+0.2, mean(out.Mfm([out.Mfmd] == jj)) , std(out.Mfm([out.Mfmd] == jj)), '*b');
%    subplot(211); plot(jj+0.2, mean(out.Mmf([out.Mmfd] == jj)),'bo');
%    subplot(212); plot(jj+0.2, mean(out.Mfm([out.Mfmd] == jj)),'b*');

%% Autogenous ISI - Start F1 through M to Start F2 and Start M1 through F to Start M2
% 
% F2F = []; F2Fd = []; F2Fm = [];
% M2M = []; M2Md = []; M2Mf = [];
% 
% % Cycle through every duet in the structure
% for d = 1:length(in)
%     
%     numsyls = length(in(d).fsyl);
%     
%     % Cycle through every syllable except the last 2
%     for s=1:numsyls-2
%        
%         if in(d).fsyl(s).sexsyltype > 49 % We have a female syllable
%             if in(d).fsyl(s+1).sexsyltype < 49 % Next syllable is a male
% %                 nextfemsylidx = find(find(in(d).fsy(s).sexsyltype > 49) > s+1, 1); % Find next female syllable, if it exists
% %                 if nextfemsylidx > 1 % We found one
%                 if in(d).fsyl(s+2).sexsyltype > 49 % The subsequent syllable is also a female
%                     F2F(end+1) = in(d).fsyl(s+2).syltim(1) - in(d).fsyl(s).syltim(2); 
%                     F2Fm(end+1) = in(d).fsyl(s+2).syltim(1) - in(d).fsyl(s).syltim(2) - in(d).fsyl(s+1).sylen;
%                     F2Fd(end+1) = in(d).distance;
%                 end
%             end
%         end
% 
%         if in(d).msyl(s).sexsyltype < 49 % We have a male syllable
%             if in(d).msyl(s+1).sexsyltype > 49 % Next syllable is a female
%                 if in(d).msyl(s+2).sexsyltype < 49 % The subsequent syllable is also a male syllable
%                     M2M(end+1) = in(d).msyl(s+2).syltim(1) - in(d).msyl(s).syltim(2); 
%                     M2Mf(end+1) = in(d).msyl(s+2).syltim(1) - in(d).msyl(s).syltim(2) - in(d).msyl(s+1).sylen; 
%                     M2Md(end+1) = in(d).distance;
%                 end
%             end
%         end  
%         
%     end
% end
% 
% out.F2F = F2F; out.F2Fd = F2Fd; out.F2Fm = F2Fm;
% out.M2M = M2M; out.M2Md = M2Md; out.M2Mf = M2Mf;
% 
% 
% figure(3); clf; 
% subplot(211); hold on; 
%     for j=1:length(out.F2Fd); plot(out.F2Fd(j), out.F2F(j), '*m');end;
%     for j=1:length(out.M2Md); plot(out.M2Md(j)+0.1, out.M2M(j), '*b');end;   
% subplot(212); hold on;
%     for jj = 1:length(distances)
%         FFF(jj) = mean(out.F2F([out.F2Fd] == distances(jj)));
%         MMM(jj) = mean(out.M2M([out.M2Md] == distances(jj)));
%         FFFstd(jj) = std(out.F2F([out.F2Fd] == distances(jj)));
%         MMMstd(jj) = std(out.M2M([out.M2Md] == distances(jj)));
%     end
% 
%     errorbar(distances-0.1, FFF, FFFstd, 'om');
%     errorbar(distances+0.1, MMM, MMMstd, 'ob');
% 
% figure(4); clf; 
% subplot(211); hold on; 
%     for j=1:length(out.F2Fd); plot(out.F2Fd(j), out.F2Fm(j), '*m');end;
%     for j=1:length(out.M2Md); plot(out.M2Md(j)+0.1, out.M2Mf(j), '*b');end;   
% subplot(212); hold on;
%     for jj = 1:length(distances)
%         FFFm(jj) = mean(out.F2Fm([out.F2Fd] == distances(jj)));
%         MMMf(jj) = mean(out.M2Mf([out.M2Md] == distances(jj)));
%         FFFmstd(jj) = std(out.F2Fm([out.F2Fd] == distances(jj)));
%         MMMfstd(jj) = std(out.M2Mf([out.M2Md] == distances(jj)));
%     end
% 
%     errorbar(distances-0.1, FFFm, FFFmstd, 'om');
%     errorbar(distances+0.1, MMMf, MMMfstd, 'ob');
    
    

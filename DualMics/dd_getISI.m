function out = dd_getISI(in)
% out = dd_getISI(in)
% Where "in" is a structure for two-microphone recordings of wrens
% "in" is composed of 
% day
% distance
% femMic    
% Fs        
% fsyl      
% Location  
% maleMic   
% month     
% msyl      
% vision    
% year      

% This is not a general tool.
% out includes the following:
% out.Fmf = []; Female microphone, ISI male-to-female 
% out.Fmfd = []; Female microphone, distance for current ISI male-to-female
% out.Ffm = []; Female microphone, ISI female-to-male
% out.Ffmd = []; Female microphone, distance for current ISI female-to-male
% out.Mmf = []; Male microphone, ISI male-to-female
% out.Mmfd = []; Male microphone, distance for current ISI male-to-female
% out.Mfm = []; Male microphone, ISI female-to-male
% out.Mfmd = []; Male microphone, distance for current ISI female-to-male

% mmfcnt = 0; mfmcnt = 0; fmfcnt = 0; ffmcnt = 0;

out.Fmf = []; out.Fmfd = [];
out.Ffm = []; out.Ffmd = [];
out.Mmf = []; out.Mmfd = [];
out.Mfm = []; out.Mfmd = [];

% out.mm = []; out.ff = [];

spdosnd = 1/331.2; % Speed of sound is 331.2 meters per second

figure(1); clf; subplot(211); hold on; subplot(212); hold on;

%% Cycle through every duet in the structure
for d = 1:length(in)
    
    numsyls = length(in(d).fsyl);
    
    % Cycle through every syllable except the last (avoid last because it
    % has no inter-syllable-interval.
    
    for s=1:numsyls-1
        
        % Get the ISI for the current syllable on each microphone
        currFisi = in(d).fsyl(s+1).syltim(1) - in(d).fsyl(s).syltim(2);
        currMisi = in(d).msyl(s+1).syltim(1) - in(d).msyl(s).syltim(2);

        % These are very long ISIs
%         if currFisi > 0.22 
%             if in(d).fsyl(s).sexsyltype < 49 && in(d).fsyl(s+1).sexsyltype > 49
%             out.ff(end+1) = d;
%             end
%         end
%         if currMisi > 0.22
%             if in(d).fsyl(s).sexsyltype < 49 && in(d).fsyl(s+1).sexsyltype > 49
%             out.mm(end+1) = d;
%             end
%         end
            
        figure(1);

        if in(d).fsyl(s).sexsyltype < 49 && in(d).fsyl(s+1).sexsyltype > 49 % Male to Female
            subplot(211); plot(in(d).distance+0.1, currFisi, 'k*');
                out.Fmf(end+1) = currFisi;
                out.Fmfd(end+1) = in(d).distance;
            subplot(212); plot(in(d).distance, currMisi, 'bo');
                out.Mmf(end+1) = currMisi;
                out.Mmfd(end+1) = in(d).distance;
        end
        if in(d).fsyl(s).sexsyltype > 49 && in(d).fsyl(s+1).sexsyltype < 49; % Female to Male
            subplot(211); plot(in(d).distance, currFisi, 'mo');
                out.Ffm(end+1) = currFisi;
                out.Ffmd(end+1) = in(d).distance;
            subplot(212); plot(in(d).distance+0.1, currMisi, 'k*');
                out.Mfm(end+1) = currMisi;
                out.Mfmd(end+1) = in(d).distance;
        end        
                        
    end
            
        
end

figure(1); subplot(211); plot([0 8], [0 0], 'k-'); ylim([-0.05 0.25]);
figure(1); subplot(212); plot([0 8], [0 0], 'k-'); ylim([-0.05 0.25]);



figure(2); clf; subplot(211); hold on; subplot(212); hold on;

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

axxx(1) = subplot(211); 
    plot(distances-0.1, Ffm(1)+(2*(distances-1)*spdosnd), 'k-');
    errorbar(distances-0.1, Ffm, Ffmstd, 'ob');
    errorbar(distances+0.1, Fmf, Fmfstd, '*m', 'LineWidth', 2);
    text(10, 0.0, 'Female Microphone', 'Color', 'm');

axxx(2) = subplot(212); 
    plot(distances-0.1, Mmf(1)+(2*(distances-1)*spdosnd), 'k-');
    errorbar(distances+0.1, Mfm, Ffmstd, '*b', 'LineWidth', 2);
    errorbar(distances-0.1, Mmf, Fmfstd, 'om');
    text(10, 0.0, 'Male Microphone', 'Color', 'b');

linkaxes(axxx, 'xy'); xlim([-1 12]); ylim([-0.02 0.18]);

fprintf('Mean male AH 1 and 7 meters: %1.5f %1.5f \n', Mmf(1), Mmf(2));
fprintf('STD male AH 1 and 7 meters: %1.5f %1.5f \n', Mmfstd(1), Mmfstd(2));

fprintf('Mean male HA 1 and 7 meters: %1.5f %1.5f \n', Mfm(1), Mfm(2));
fprintf('STD male HA 1 and 7 meters: %1.5f %1.5f \n', Mfm(1), Mfm(2));

fprintf('Mean female AH 1 and 7 meters: %1.5f %1.5f \n', Ffm(1), Ffm(2));
fprintf('Mean female HA 1 and 7 meters: %1.5f %1.5f \n', Fmf(1), Fmf(2));

% Male HA 1 to 7 meters
[aa,bb,cc,dd] = ttest2(out.Mfm([out.Mfmd] == 1),  out.Mfm([out.Mfmd] == 7), 'vartype', 'unequal');
fprintf('Male HA 1 vs 7 meters, p = %1.5f, tstat = %1.5f, df = %1.5f. \n', bb, dd.tstat, dd.df);

% Male AH 1 to 7 meters
[aa,bb,cc,dd] = ttest2(out.Mmf([out.Mmfd] == 1),  out.Mmf([out.Mmfd] == 7), 'vartype', 'unequal');
fprintf('Male AH 1 vs 7 meters, p = %1.5f, tstat = %1.5f, df = %1.5f. \n', bb, dd.tstat, dd.df);

% FeMale AH 1 to 7 meters
[aa,bb,cc,dd] = ttest2(out.Ffm([out.Ffmd] == 1),  out.Ffm([out.Ffmd] == 7), 'vartype', 'unequal');
fprintf('Female AH 1 vs 7 meters, p = %1.5f, tstat = %1.5f, df = %1.5f. \n', bb, dd.tstat, dd.df);

% FeMale HA 1 to 7 meters
[aa,bb,cc,dd] = ttest2(out.Fmf([out.Fmfd] == 1),  out.Fmf([out.Fmfd] == 7), 'vartype', 'unequal');
fprintf('Female HA 1 vs 7 meters, p = %1.5f, tstat = %1.5f, df = %1.5f. \n', bb, dd.tstat, dd.df);

    
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
    
    

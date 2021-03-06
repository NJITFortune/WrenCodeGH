function out = dd_ISI_chron(in)
% out = dd_ISI_chron(in)
% Where "in" is a structure for two-microphone recordings of wrens.

spdosnd = 1/331.2; % Speed of sound is 331.2 meters per second

out.F.ha = []; out.F.had = []; % Female microphone, Het-Aut ISI
out.F.ah = []; out.F.ahd = []; % Female microphone, Aut-Het ISI
out.F.aha = []; out.F.ahad = []; % Female microphone, Aut-Het-Aut ISI
out.M.ah = []; out.M.ahd = []; % Male microphone, Aut-Het ISI
out.M.ha = []; out.M.had = []; % Male microphone, Het-Aut ISI
out.M.aha = []; out.M.ahad = []; % Male microphone, Aut-Het-Aut ISI

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
            
        figure(1);

        if in(d).fsyl(s).sexsyltype < 49 && in(d).fsyl(s+1).sexsyltype > 49 % Male to Female
            subplot(211); plot(in(d).distance+0.1, currFisi, 'k*');
                out.F.ha(end+1) = currFisi;
                out.F.had(end+1) = in(d).distance;
            subplot(212); plot(in(d).distance, currMisi, 'bo');
                out.M.ah(end+1) = currMisi;
                out.M.ahd(end+1) = in(d).distance;
                
                if s+2 <= numsyls
                if in(d).fsyl(s+2).sexsyltype < 49 % Male to Female to Male (AHA interval for male)
                    out.M.aha(end+1) = in(d).msyl(s+2).syltim(1) - in(d).msyl(s).syltim(2);
                    out.M.ahad(end+1) = in(d).distance;
                end
                end
        end
        
        if in(d).fsyl(s).sexsyltype > 49 && in(d).fsyl(s+1).sexsyltype < 49 % Female to Male
            subplot(211); plot(in(d).distance, currFisi, 'mo');
                out.F.ah(end+1) = currFisi;
                out.F.ahd(end+1) = in(d).distance;
            subplot(212); plot(in(d).distance+0.1, currMisi, 'k*');
                out.M.ha(end+1) = currMisi;
                out.M.had(end+1) = in(d).distance;
                
                if s+2 <= numsyls
                if in(d).fsyl(s+2).sexsyltype > 49 % Female to Male to Female (AHA interval for female)
                    out.F.aha(end+1) = in(d).fsyl(s+2).syltim(1) - in(d).fsyl(s).syltim(2);
                    out.F.ahad(end+1) = in(d).distance;
                end
                end
        end        
                        
    end
            
        
end

figure(1); subplot(211); plot([0 8], [0 0], 'k-'); ylim([-0.05 0.25]);
figure(1); subplot(212); plot([0 8], [0 0], 'k-'); ylim([-0.05 0.25]);



figure(2); clf; subplot(211); hold on; subplot(212); hold on;

distances = sort(unique([in.distance]));

for jj = length(distances):-1:1
    
    Fah(jj) = mean(out.F.ah([out.F.ahd] == distances(jj)));
    Fha(jj) = mean(out.F.ha([out.F.had] == distances(jj)));
    Fahstd(jj) = std(out.F.ah([out.F.ahd] == distances(jj)));
    Fhastd(jj) = std(out.F.ha([out.F.had] == distances(jj)));

    Mha(jj) = mean(out.M.ha([out.M.had] == distances(jj)));
    Mah(jj) = mean(out.M.ah([out.M.ahd] == distances(jj)));
    Mhastd(jj) = std(out.M.ha([out.M.had] == distances(jj)));
    Mahstd(jj) = std(out.M.ah([out.M.ahd] == distances(jj)));
    
    Faha(jj) = mean(out.F.aha([out.F.ahad] == distances(jj)));
    Fahastd(jj) = std(out.F.aha([out.F.ahad] == distances(jj)));
    Maha(jj) = mean(out.M.aha([out.M.ahad] == distances(jj)));
    Mahastd(jj) = std(out.M.aha([out.M.ahad] == distances(jj)));

end

axxx(1) = subplot(311); 
    plot(distances-0.1, Fah(1)+(2*(distances-1)*spdosnd), 'k-');
    errorbar(distances-0.1, Fah, Fahstd, 'ob');
    errorbar(distances+0.1, Fha, Fhastd, '*m', 'LineWidth', 2);
    text(10, 0.0, 'Female Microphone', 'Color', 'm');

axxx(2) = subplot(312); 
    plot(distances-0.1, Mah(1)+(2*(distances-1)*spdosnd), 'k-');
    errorbar(distances+0.1, Mha, Mhastd, '*b', 'LineWidth', 2);
    errorbar(distances-0.1, Mah, Mahstd, 'om');
    text(10, 0.0, 'Male Microphone', 'Color', 'b');

axxx(3) = subplot(313); 
    plot(distances-0.1, Mah(1)+(2*(distances-1)*spdosnd), 'k-');
    errorbar(distances+0.1, Mha, Mhastd, '*b', 'LineWidth', 2);
    errorbar(distances-0.1, Mah, Mahstd, 'om');
    text(10, 0.0, 'Male Microphone', 'Color', 'b');

linkaxes(axxx, 'xy'); xlim([-1 12]); ylim([-0.02 0.18]);
    
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
    
    

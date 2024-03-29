function out = dd_getISidx(in)
% out = dd_getISI(in)
% Where in is a structure for two-microphone recordings of wrens
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

idxpairs = sort(unique([in.pairnum]));

    blues = winter(length(idxpairs)*2); % Colors for male birds
    reds = copper(length(idxpairs)*2); % Colors for female birds

% mmfcnt = 0; mfmcnt = 0; fmfcnt = 0; ffmcnt = 0;

out(length(idxpairs)).Fmf = []; out(length(idxpairs)).Fmfd = []; out(length(idxpairs)).Fmfpredur = [];
out(length(idxpairs)).Ffm = []; out(length(idxpairs)).Ffmd = []; out(length(idxpairs)).Ffmpredur = [];
out(length(idxpairs)).Mmf = []; out(length(idxpairs)).Mmfd = []; out(length(idxpairs)).Mmfpredur = [];
out(length(idxpairs)).Mfm = []; out(length(idxpairs)).Mfmd = []; out(length(idxpairs)).Mfmpredur = [];

mm = []; ff = [];

%% Cycle through every duet in the structure

for jj = length(idxpairs):-1:1 % For each pair

    idx = find([in.pairnum] == idxpairs(jj));

for d = length(idx):-1:1 % Go through each duet

    numsyls = length(in(idx(d)).fsyl);
    
    % Cycle through every syllable except the last
    for s=1:numsyls-1
        
        % Get the ISI for the current syllable on each microphone
        currFisi = in(idx(d)).fsyl(s+1).syltim(1) - in(idx(d)).fsyl(s).syltim(2);
        currMisi = in(idx(d)).msyl(s+1).syltim(1) - in(idx(d)).msyl(s).syltim(2);

        if currFisi > 0.22 
            if in(idx(d)).fsyl(s).sexsyltype < 49 && in(idx(d)).fsyl(s+1).sexsyltype > 49 % MF
                ff(end+1) = d;
            end
        end
        if currMisi > 0.22 
            if in(idx(d)).fsyl(s).sexsyltype < 49 && in(idx(d)).fsyl(s+1).sexsyltype > 49 % MF
                mm(end+1) = d;
            end
        end
            
        % Figure out male and female syllable sequences

        if in(idx(d)).fsyl(s).sexsyltype < 49 && in(idx(d)).fsyl(s+1).sexsyltype > 49 % Male to Female
            %subplot(211); plot(in(idx(d)).distance+0.1, currFisi, 'k*');
                out(jj).Fmf(end+1) = currFisi;
                out(jj).Fmfd(end+1) = in(idx(d)).distance;
                out(jj).Fmfpredur(end+1) = in(idx(d)).msyl(s).sylen;
            %subplot(212); plot(in(idx(d)).distance, currMisi, 'bo');
                out(jj).Mmf(end+1) = currMisi;
                out(jj).Mmfd(end+1) = in(idx(d)).distance;
                out(jj).Mmfpredur(end+1) = in(idx(d)).msyl(s).sylen;
        end
        if in(idx(d)).fsyl(s).sexsyltype > 49 && in(idx(d)).fsyl(s+1).sexsyltype < 49 % Female to Male
            %subplot(211); plot(in(idx(d)).distance, currFisi, 'mo');
                out(jj).Ffm(end+1) = currFisi;
                out(jj).Ffmd(end+1) = in(idx(d)).distance;
                out(jj).Ffmpredur(end+1) = in(idx(d)).fsyl(s).sylen;
            %subplot(212); plot(in(idx(d)).distance+0.1, currMisi, 'k*');
                out(jj).Mfm(end+1) = currMisi;
                out(jj).Mfmd(end+1) = in(idx(d)).distance;
                out(jj).Mfmpredur(end+1) = in(idx(d)).fsyl(s).sylen;
        end        
                        
    end
        
end

end

%figure(1); subplot(211); plot([0 11], [0 0], 'k-');
%figure(1); subplot(212); plot([0 11], [0 0], 'k-');
%out.ff = ff;
%out.mm = mm;

for kk = 1:length(idxpairs)

figure(kk*2); clf; 
    subplot(211); hold on; plot([0, 10], [0, 0.033] + 0.05, 'k-', 'LineWidth', 2);
        text(0, 0.2, in(find([in.pairnum] == idxpairs(kk), 1, 'first')).pairname);
    subplot(212); hold on; plot([0, 10], [0, 0.033] + 0.05, 'k-', 'LineWidth', 2);

figure((kk*2) + 1); clf; 
    subplot(211); hold on; 
        text(0, 0.2, in(find([in.pairnum] == idxpairs(kk), 1, 'first')).pairname);
        plot([0, 10], [0, 0.033] + 0.05, 'k-', 'LineWidth', 2);
        plot([0, 10], [0, 0] + 0.05, 'r-');

    subplot(212); hold on; 
        plot([0, 10], [0, 0.033] + 0.05, 'k-', 'LineWidth', 2);
        plot([0, 10], [0, 0] + 0.05, 'r-');
%    subplot(411); hold on; plot([0, 10], [0, 0.033] + 0.05, 'k-');
%    subplot(412); hold on; plot([0, 10], [0, 0.033] + 0.05, 'k-');
%    subplot(413); hold on; plot([0, 10], [0, 0.033] + 0.05, 'k-');
%    subplot(414); hold on; plot([0, 10], [0, 0.033] + 0.05, 'k-');

    idx = find([in.pairnum] == idxpairs(kk));

    distances = sort(unique([in(idx).distance]));

    for dn = 1:length(distances)
        
        currFMdIDX = find([out(kk).Ffmd] == distances(dn));
        currMFdIDX = find([out(kk).Fmfd] == distances(dn));

        dad(kk).Ffm(dn) = mean(out(kk).Ffm(currFMdIDX));
        dad(kk).Fmf(dn) = mean(out(kk).Fmf(currMFdIDX));
        dad(kk).Ffmstd(dn) = std(out(kk).Ffm(currFMdIDX));
        dad(kk).Fmfstd(dn) = std(out(kk).Fmf(currMFdIDX));
   
        dad(kk).Mfm(dn) = mean(out(kk).Mfm(currFMdIDX));
        dad(kk).Mmf(dn) = mean(out(kk).Mmf(currMFdIDX));
        dad(kk).Mfmstd(dn) = std(out(kk).Mfm(currFMdIDX));
        dad(kk).Mmfstd(dn) = std(out(kk).Mmf(currMFdIDX));
    

figure(kk*2); % Raw dots

subplot(211); hold on;
    jt = rand(1,length(currFMdIDX));
    plot((distances(dn) * ones(1,length(currFMdIDX))) + (jt/5), out(kk).Ffm(currFMdIDX), '.', 'Color', blues(kk*2,:), 'MarkerSize', 16);
    ylim([0 0.2]); xlim([-1 11]); ylim([-0.05 0.22]);
subplot(212); hold on;
    jt = rand(1,length(currMFdIDX));
    plot((distances(dn) * ones(1,length(currMFdIDX))) + (jt/5), out(kk).Mmf(currMFdIDX), '.', 'Color', reds(kk*2,:), 'MarkerSize', 16);
    ylim([0 0.2]); xlim([-1 11]); ylim([-0.05 0.22]);
    end

figure((kk*2) + 1); % Error bars

subplot(211); hold on;
    errorbar(distances, dad(kk).Ffm, dad(kk).Ffmstd, 'o', 'Color', blues(kk*2,:), 'LineWidth', 2);
    %text(1, 0.0, 'Female Microphone', 'Color', 'm');
    %xlim([-1 11]); ylim([-0.1 0.22]);

    errorbar(distances + 0.2, dad(kk).Fmf, dad(kk).Fmfstd, '*', 'Color', reds(kk*2,:));
    text(1, 0.0, 'Female Microphone', 'Color', 'm');
    xlim([-1 11]); ylim([-0.05 0.22]);

subplot(212); hold on;
    errorbar(distances, dad(kk).Mmf, dad(kk).Fmfstd, 'o', 'Color', reds(kk*2,:), 'LineWidth', 2);
%    text(1, 0.0, 'Male Microphone', 'Color', 'b');
%    xlim([-1 11]); ylim([-0.1 0.22]);

%subplot(414); hold on;
    errorbar(distances + 0.2, dad(kk).Mfm, dad(kk).Ffmstd, '*', 'Color', blues(kk*2,:));
    text(1, 0.0, 'Male Microphone', 'Color', 'b');
    xlim([-1 11]); ylim([-0.05 0.22]);

end
    
%    subplot(211); errorbar(jj, mean(out.Ffm([out.Ffmd] == jj)) , std(out.Ffm([out.Ffmd] == jj)), 'om');
%    subplot(212); errorbar(jj, mean(out.Fmf([out.Fmfd] == jj)) , std(out.Fmf([out.Fmfd] == jj)), '*m');
%    subplot(211); plot(jj, mean(out.Ffm([out.Ffmd] == jj)),'mo');
%    subplot(212); plot(jj, mean(out.Fmf([out.Ffmd] == jj)),'m*');

%    subplot(211); errorbar(jj+0.2, mean(out.Mmf([out.Mmfd] == jj)) , std(out.Mmf([out.Mmfd] == jj)), 'ob');
%    subplot(212); errorbar(jj+0.2, mean(out.Mfm([out.Mfmd] == jj)) , std(out.Mfm([out.Mfmd] == jj)), '*b');
%    subplot(211); plot(jj+0.2, mean(out.Mmf([out.Mmfd] == jj)),'bo');
%    subplot(212); plot(jj+0.2, mean(out.Mfm([out.Mfmd] == jj)),'b*');

%% Autogenous ISI - Start F1 through M to Start F2 and Start M1 through F to Start M2

F2F = []; F2Fd = []; F2Fm = [];
M2M = []; M2Md = []; M2Mf = [];

% Cycle through every duet in the structure
for d = 1:length(idx)
    
    numsyls = length(in(idx(d)).fsyl);
    
    % Cycle through every syllable except the last 2
    for s=1:numsyls-2
       
        if in(idx(d)).fsyl(s).sexsyltype > 49 % We have a female syllable
            if in(idx(d)).fsyl(s+1).sexsyltype < 49 % Next syllable is a male
%                 nextfemsylidx = find(find(in(idx(d)).fsy(s).sexsyltype > 49) > s+1, 1); % Find next female syllable, if it exists
%                 if nextfemsylidx > 1 % We found one
                if in(idx(d)).fsyl(s+2).sexsyltype > 49 % The subsequent syllable is also a female
                    F2F(end+1) = in(idx(d)).fsyl(s+2).syltim(1) - in(idx(d)).fsyl(s).syltim(2); 
                    F2Fm(end+1) = in(idx(d)).fsyl(s+2).syltim(1) - in(idx(d)).fsyl(s).syltim(2) - in(idx(d)).fsyl(s+1).sylen;
                    F2Fd(end+1) = in(idx(d)).distance;
                end
            end
        end

        if in(idx(d)).msyl(s).sexsyltype < 49 % We have a male syllable
            if in(idx(d)).msyl(s+1).sexsyltype > 49 % Next syllable is a female
                if in(idx(d)).msyl(s+2).sexsyltype < 49 % The subsequent syllable is also a male syllable
                    M2M(end+1) = in(idx(d)).msyl(s+2).syltim(1) - in(idx(d)).msyl(s).syltim(2); 
                    M2Mf(end+1) = in(idx(d)).msyl(s+2).syltim(1) - in(idx(d)).msyl(s).syltim(2) - in(idx(d)).msyl(s+1).sylen; 
                    M2Md(end+1) = in(idx(d)).distance;
                end
            end
        end  
        
    end
end
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
    
    

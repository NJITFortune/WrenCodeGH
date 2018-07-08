function out = duetsylsanity(in)


Fs = 1 / (in(1).fc.ftim(2) - in(1).fc.ftim(1));


%% Plot the spectrograms with the lines.

figure(1); clf; 

ax(1) = subplot(211); specgram(in(1).fc.fem, 1024, Fs, [], 1000); ylim([0 6000]);
subplot(211); hold on; 

    for j=1:length(in(1).f) 
        plot([in(1).f(j).syltim(1) in(1).f(j).syltim(1)], [500 5500], 'g', 'LineWidth', 2);
        plot([in(1).f(j).syltim(2) in(1).f(j).syltim(2)], [500 5500], 'r', 'LineWidth', 2);
    end;

ax(2) = subplot(212); specgram(in(1).fc.mal, 1024, Fs, [], 1000); ylim([0 6000]);
subplot(212); hold on; 

    for j=1:length(in(1).m) 
        plot([in(1).m(j).syltim(1) in(1).m(j).syltim(1)], [500 5500], 'g', 'LineWidth', 2);
        plot([in(1).m(j).syltim(2) in(1).m(j).syltim(2)], [500 5500], 'r', 'LineWidth', 2);
    end;

    %% Get lists of syllables
    
    malsyls = [];
    for j=1:length(in(1).anal.maltype)
        malsyls = [malsyls in(1).anal.maltype{j}];
    end
    malsyls = sort(malsyls);
    
    femsyls = [];
    for j=1:length(in(1).anal.femtype)
        femsyls = [femsyls in(1).anal.femtype{j}];
    end
    femsyls = sort(femsyls);
    
    speedofsound = 0.0033;
    
    distance = 3 * 2;
    
    sndelay = distance * speedofsound;
    
    
        for j=1:length(malsyls)
            subplot(211); 
            %plot([in(1).m(malsyls(j)).syltim(1)+sndelay, in(1).m(malsyls(j)).syltim(1)+sndelay], [1000 4000], 'k', 'LineWidth', 1);
            %plot([in(1).m(malsyls(j)).syltim(2)+sndelay, in(1).m(malsyls(j)).syltim(2)+sndelay], [1000 4000], 'k', 'LineWidth', 1);
            plot(in(1).f(malsyls(j)).trace_tim + in(1).f(malsyls(j)).syltim(1), in(1).f(malsyls(j)).trace_freq, 'k');
            subplot(212);
            plot(in(1).m(malsyls(j)).trace_tim + in(1).m(malsyls(j)).syltim(1), in(1).m(malsyls(j)).trace_freq, 'k');
        end;
    
        for j=1:length(femsyls)
            subplot(212); 
            %plot([in(1).f(femsyls(j)).syltim(1)+0, in(1).f(femsyls(j)).syltim(1)+0], [1000 4000], 'k', 'LineWidth', 1);
            %plot([in(1).f(femsyls(j)).syltim(2)+0, in(1).f(femsyls(j)).syltim(2)+0], [1000 4000], 'k', 'LineWidth', 1);
            subplot(211);
            %plot(in(1).f(femsyls(j)).trace_tim + in(1).f(femsyls(j)).syltim(1), in(1).f(femsyls(j)).trace_freq, 'k');
        end;
        
        linkaxes(ax, 'x');
        
        
        
%% Pre- Post- Analysis

% MALE SYLLABLES (use recordings aligned to female syllables)

for pp = 1:5
    Mpre{pp} = []; Mpost{pp} = [];
end

for j=1:length(in) % For each duet

    numsyltypes = length(in(j).anal.maltype); % Get number of syllables for the male in the duet
        malsyls = [];
        for jj=1:length(in(1).anal.maltype)
            malsyls = [malsyls in(1).anal.maltype{jj}];
        end
    malsyls = sort(malsyls);

        femsyls = [];
        for jj=1:length(in(1).anal.femtype)
            femsyls = [femsyls in(1).anal.femtype{jj}];
        end
    femsyls = sort(femsyls);
    
    for k=1:numsyltypes % For each syllable type
       
        for i=1:length(in(j).anal.maltype{k}) % For each syllable of that type
        
            currsyllableIDX = in(j).anal.maltype{k}(i);

            if find(femsyls == currsyllableIDX-1)
                endofpresyllable = in(j).f(currsyllableIDX-1).syltim(2);
                Mpre{k}(end+1) = in(j).m(currsyllableIDX).syltim(1) - endofpresyllable;
            end

            if find(femsyls == currsyllableIDX+1)
                startofpostsyllable = in(j).f(currsyllableIDX+1).syltim(1);
                Mpost{k}(end+1) = startofpostsyllable - in(j).m(currsyllableIDX).syltim(2);
            end

        end
        
    end

end

out.Mpre = Mpre;
out.Mpost = Mpost;

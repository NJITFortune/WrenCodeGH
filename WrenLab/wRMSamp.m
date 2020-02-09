function [m, f] = wRMSamp(in)
% Compares solo and duet syllable amplitudes.
% Use only recordings in which both solo and duet syllables are present
% as a form of internal control.
% Requires ChronicCompleat2019f.mat, wData.m

msoloAmp = []; mduetAmp = []; msoloDur = []; mduetDur = [];
fsoloAmp = []; fduetAmp = []; fsoloDur = []; fduetDur = [];

% Get our list of data
[msolosyls, mduetsyls, fsolosyls, fduetsyls, ~] = wData;

%% Male syllables
for j=1:length(msolosyls)
    
  if ~isempty(msolosyls{j}) && ~isempty(mduetsyls{j})
      
    % Get amplitudes and durations of each solo syllable
    for k=1:length(msolosyls{j})
        msoloAmp(end+1) = rms(in(j*2).duet(in(j*2).tim > in(j*2).syl(msolosyls{j}(k)).tim(1) & in(j*2).tim < in(j*2).syl(msolosyls{j}(k)).tim(2)));
        msoloDur(end+1) = abs(in(j*2).syl(msolosyls{j}(k)).tim(2) - in(j*2).syl(msolosyls{j}(k)).tim(1));
    end
    
    % Get amplitudes and duration of each duet syllable
    for k=1:length(mduetsyls{j})
        mduetAmp(end+1) = rms(in(j*2).duet(in(j*2).tim > in(j*2).syl(mduetsyls{j}(k)).tim(1) & in(j*2).tim < in(j*2).syl(mduetsyls{j}(k)).tim(2)));
        mduetDur(end+1) = abs(in(j*2).syl(mduetsyls{j}(k)).tim(2) - in(j*2).syl(mduetsyls{j}(k)).tim(1));
    end
            
  end
    
end

%% Female syllables
for j=1:length(fsolosyls)
    
  if ~isempty(fsolosyls{j}) && ~isempty(fduetsyls{j})
      
    % Get amplitudes and durations of each solo syllable
    for k=1:length(fsolosyls{j})
        fsoloAmp(end+1) = rms(in(j*2).duet(in(j*2).tim > in(j*2).syl(fsolosyls{j}(k)).tim(1) & in(j*2).tim < in(j*2).syl(fsolosyls{j}(k)).tim(2)));
        fsoloDur(end+1) = abs(in(j*2).syl(fsolosyls{j}(k)).tim(2) - in(j*2).syl(fsolosyls{j}(k)).tim(1));
    end
    
    % Get amplitudes and duration of each duet syllable
    for k=1:length(fduetsyls{j})
        fduetAmp(end+1) = rms(in(j*2).duet(in(j*2).tim > in(j*2).syl(fduetsyls{j}(k)).tim(1) & in(j*2).tim < in(j*2).syl(fduetsyls{j}(k)).tim(2)));
        fduetDur(end+1) = abs(in(j*2).syl(fduetsyls{j}(k)).tim(2) - in(j*2).syl(fduetsyls{j}(k)).tim(1));
    end
            
  end
    
end

%% Plot and such

figure(1); clf; 
subplot(211); hold on; 
    yyaxis left; plot(msoloAmp, '-*b'); ylim([0 0.3]); ylabel('Solo Amp');
    yyaxis right; plot(msoloDur); ylabel('Solo Duration');
subplot(212); hold on; 
    yyaxis left; plot(mduetAmp, '-*b'); ylim([0 0.3]); ylabel('Duet Amp');
    yyaxis right; plot(mduetDur); ylabel('Duet Duration');

m.amp = 20*log(mean(mduetAmp)/mean(msoloAmp));
    fprintf('Mean amplitude increase in Males: %f2.4 \n', m.amp);
[m.H,m.P,m.CI,m.Stats] = ttest2(msoloAmp, mduetAmp)

figure(2); clf; 
subplot(211); hold on; 
    yyaxis left; plot(fsoloAmp, '-*m'); ylim([0 0.3]); ylabel('Solo Amp');
    yyaxis right; plot(fsoloDur); ylabel('Solo Duration');
subplot(212); hold on; 
    yyaxis left; plot(fduetAmp, '-*m'); ylim([0 0.3]); ylabel('Duet Amp');
    yyaxis right; plot(fduetDur); ylabel('Duet Duration');

f.amp = 20*log(mean(fduetAmp(1:end-2))/mean(fsoloAmp(1:end-1)))
[f.H,f.P,f.CI,f.Stats] = ttest2(fsoloAmp(1:end-1), fduetAmp(1:end-2))

end
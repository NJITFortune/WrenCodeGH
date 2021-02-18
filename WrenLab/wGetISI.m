function out = wGetISI(in)
% Usage out = wGetISI(in)
% load ChronicCompleat2020c.mat
% out = wGetISI(w); % Is the command you should use
% This provides behavioral measurements for intersyllable intervals, and
% syllable durations.


% Get the list of syllables
[msolosyls, mduetsyls, fsolosyls, fduetsyls, ~, ~] = wData;



% Initialize variables 
    out.FM = []; out.MF = [];
    out.Fdur = []; out.Mdur = [];

%% Cycle for each duet syllable

% MALE DUET SYLLABLES
for j=1:length(mduetsyls)    % For each duet
    for k=1:length(mduetsyls{j}) % For each male duet syllable in each duet
        
        % Durations of male syllables
        out.Mdur(end+1) = in(j*2).syl(mduetsyls{j}(k)).tim(2) - in(j*2).syl(mduetsyls{j}(k)).tim(1);
        
        % Intersyllable interval between end of female syllable and start
        % of subsequent male syllable
        if ~isempty(find(fduetsyls{j} == mduetsyls{j}(k)-1, 1)) % If there was a prior female syllable
            out.FM(end+1) = in(j*2).syl(mduetsyls{j}(k)).tim(1) - in(j*2).syl(mduetsyls{j}(k)-1).tim(2);
        end
    end
end

% FEMALE DUET SYLLABLES
for j=1:length(fduetsyls)    
    for k=1:length(fduetsyls{j})        
        
        % Durations of female syllables
        out.Fdur(end+1) = in(j*2).syl(fduetsyls{j}(k)).tim(2) - in(j*2).syl(fduetsyls{j}(k)).tim(1);
        
        % Intersyllable interval between end of male syllable and start
        % of subsequent female syllable
        if ~isempty(find(mduetsyls{j} == fduetsyls{j}(k)-1, 1))
            out.MF(end+1) = in(j*2).syl(fduetsyls{j}(k)).tim(1) - in(j*2).syl(fduetsyls{j}(k)-1).tim(2);
        end
    end
end

%% Cycle for each solo syllable

% Initialize variables 
    out.FS = []; out.MS = [];
    out.FSdur = []; out.MSdur = [];

% MALE SOLO SYLLABLES
for j=1:length(msolosyls)    % For each duet
    for k=1:length(msolosyls{j}) % For each male solo syllable 
        
        % Durations of male syllables
        out.MSdur(end+1) = in(j*2).syl(msolosyls{j}(k)).tim(2) - in(j*2).syl(msolosyls{j}(k)).tim(1);
        
        % Intersyllable interval between end of male solo syllable and start
        % of subsequent syllable
            out.MS(end+1) = in(j*2).syl(msolosyls{j}(k)+1).tim(1) - in(j*2).syl(msolosyls{j}(k)).tim(2);
    end
end

% FEMALE SOLO SYLLABLES
for j=1:length(fsolosyls)    
    for k=1:length(fsolosyls{j})        
        
        % Durations of female syllables
        out.FSdur(end+1) = in(j*2).syl(fsolosyls{j}(k)).tim(2) - in(j*2).syl(fsolosyls{j}(k)).tim(1);
        
        % Intersyllable interval between end of male syllable and start
        % of subsequent female syllable
            out.FS(end+1) = in(j*2).syl(fsolosyls{j}(k)+1).tim(1) - in(j*2).syl(fsolosyls{j}(k)).tim(2);
    end
end



%% Output the measurements to the screen
fprintf('The mean and std for M2F ISI is  %1.5f %1.5f. n= %i \n', mean(out.MF), std(out.MF), length(out.MF));
fprintf('The mean and std for F2M ISI is  %1.5f %1.5f. n= %i \n', mean(out.FM), std(out.FM), length(out.FM));
[~, pVal, ~, ~] =  ttest2(out.FM, out.MF, 'vartype', 'unequal');
fprintf('Are these different or not p = %1.5f \n', pVal);

fprintf('************\n');

fprintf('Male Duet Syllable Duration: mean %1.3f, std %1.3f, n %i. \n', mean(out.Mdur), std(out.Mdur), length(out.Mdur));
fprintf('Female Duet Syllable Duration: mean %1.3f, std %1.3f, n %i. \n', mean(out.Fdur), std(out.Fdur), length(out.Fdur));
[~, pVal, ~, ~] =  ttest2(out.Mdur, out.Fdur, 'vartype', 'unequal');
fprintf('Are these different or not p = %1.5f \n', pVal);

fprintf('************\n');

fprintf('The mean and std for Male after ISI is  %1.5f %1.5f. n= %i \n', mean(out.MS), std(out.MS), length(out.MS));
fprintf('The mean and std for Female after ISI is  %1.5f %1.5f. n= %i \n', mean(out.FS), std(out.FS), length(out.FS));
[~, pVal, ~, ~] =  ttest2(out.FS, out.MS, 'vartype', 'unequal');
fprintf('Are these different or not p = %1.5f (which does not matter at all). \n', pVal);

fprintf('************\n');

fprintf('Male Solo Syllable Duration: mean %1.3f, std %1.3f, n %i. \n', mean(out.MSdur), std(out.MSdur), length(out.MSdur));
fprintf('Female Solo Syllable Duration: mean %1.3f, std %1.3f, n %i. \n', mean(out.FSdur), std(out.FSdur), length(out.FSdur));
[~, pVal, ~, ~] =  ttest2(out.MSdur, out.FSdur, 'vartype', 'unequal');
fprintf('Are these different or not p = %1.5f \n', pVal);



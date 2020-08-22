function out = wGetISI(in)
% Usage out = wGetISI(in)

% Get the list of duet syllables
[~, mduetsyls, ~, fduetsyls, ~, ~] = wData;

out.FM = []; out.MF = [];
out.Fdur = []; out.Mdur = [];
for j=1:length(mduetsyls)    % For each duet
    for k=1:length(mduetsyls{j}) % For each male duet syllable in each duet
        
        out.Mdur(end+1) = in(j*2).syl(mduetsyls{j}(k)).tim(2) - in(j*2).syl(mduetsyls{j}(k)).tim(1);
        
        if ~isempty(find(fduetsyls{j} == mduetsyls{j}(k)-1, 1)) % If there was a prior female syllable
            out.FM(end+1) = in(j*2).syl(mduetsyls{j}(k)).tim(1) - in(j*2).syl(mduetsyls{j}(k)-1).tim(2);
        end
    end
end

for j=1:length(fduetsyls)    
    for k=1:length(fduetsyls{j})        
        
        out.Fdur(end+1) = in(j*2).syl(fduetsyls{j}(k)).tim(2) - in(j*2).syl(fduetsyls{j}(k)).tim(1);
        
        if ~isempty(find(mduetsyls{j} == fduetsyls{j}(k)-1, 1))
            out.MF(end+1) = in(j*2).syl(fduetsyls{j}(k)).tim(1) - in(j*2).syl(fduetsyls{j}(k)-1).tim(2);
        end
    end
end


fprintf('The mean and std for M2F ISI is  %1.5f %1.5f \n', mean(out.MF), std(out.MF));
fprintf('The mean and std for F2M ISI is  %1.5f %1.5f \n', mean(out.FM), std(out.FM));
[~, pVal, ~, ~] =  ttest2(out.FM, out.MF, 'vartype', 'unequal');
fprintf('Are these different or not p = %1.5f \n', pVal);

fprintf('************\n');

fprintf('Male Syllable Duration: mean %1.3f, std %1.3f, n %i. \n', mean(out.Mdur), std(out.Mdur), length(out.Mdur));
fprintf('Female Syllable Duration: mean %1.3f, std %1.3f, n %i. \n', mean(out.Fdur), std(out.Fdur), length(out.Fdur));
[~, pVal, ~, ~] =  ttest2(out.Mdur, out.Fdur, 'vartype', 'unequal');
fprintf('Are these different or not p = %1.5f \n', pVal);



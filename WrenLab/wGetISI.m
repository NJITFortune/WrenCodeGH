function out = wGetISI(in)
% Usage out = wGetISI(in)

% Get the list of duet syllables
[~, mduetsyls, ~, fduetsyls, ~, ~] = wData;

out.FM = []; out.MF = [];


for j=1:length(mduetsyls)    % For each duet
    for k=1:length(mduetsyls{j}) % For each male duet syllable in each duet
        if ~isempty(find(fduetsyls{j} == mduetsyls{j}(k)-1, 1)) % If there was a prior female syllable
            out.FM(end+1) = in(j*2).syl(mduetsyls{j}(k)).tim(1) - in(j*2).syl(mduetsyls{j}(k)-1).tim(2);
        end
    end
end

for j=1:length(fduetsyls)    
    for k=1:length(fduetsyls{j})        
        if ~isempty(find(mduetsyls{j} == fduetsyls{j}(k)-1, 1))
            out.MF(end+1) = in(j*2).syl(fduetsyls{j}(k)).tim(1) - in(j*2).syl(fduetsyls{j}(k)-1).tim(2);
        end
    end
end



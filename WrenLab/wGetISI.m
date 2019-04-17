function out = wGetISI(in)
% Usage out = wGetISI(in)

% Get the list of duet syllables
[~, mduetsyls, ~, fduetsyls, ~, ~] = wData;

out.m = []; out.f = [];


for j=1:length(mduetsyls)    
    j
    for k=1:length(mduetsyls{j})
        k
        if ~isempty(find(fduetsyls{j} == mduetsyls{j}(k)+1, 1))
            out.m(end+1) = in(j*2).syl(mduetsyls{j}(k)+1).tim(1) - in(j*2).syl(mduetsyls{j}(k)).tim(2);
        end
    end
end

for j=1:length(fduetsyls)    
    for k=1:length(fduetsyls{j})        
        if ~isempty(find(mduetsyls{j} == fduetsyls{j}(k)+1, 1))
            out.f(end+1) = in(j*2).syl(fduetsyls{j}(k)+1).tim(1) - in(j*2).syl(fduetsyls{j}(k)).tim(2);
        end
    end
end



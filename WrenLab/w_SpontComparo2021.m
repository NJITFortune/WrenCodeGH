%% Get spontaneous windows

[~, ~, ~, ~, Cspon, ~, ~] = wData;

%% Male

malecnt = zeros(1,8);

for j=1:2:15
    
    spontims = Cspon(:,ceil(j/2));

    for k = length(w(j).Cspikes)
        malecnt(ceil(j/2)) = malecnt(ceil(j/2)) + (length(find(w(j).Cspikes{k} > spontims(1) & w(j).Cspikes{k} < spontims(2))));
    end
    malecnt(ceil(j/2)) =  malecnt(ceil(j/2)) / (spontims(2) - spontims(1));
end



%% Female

femcnt = zeros(1,8);

for j=2:2:16
    
    spontims = Cspon(:,j/2);

    for k = length(w(j).Cspikes)
        femcnt(j/2) = femcnt(j/2) + length(find(w(j).Cspikes{k} > spontims(1) & w(j).Cspikes{k} < spontims(2)));
    end
    femcnt(ceil(j/2)) =  femcnt(ceil(j/2)) / (spontims(2) - spontims(1));
    
end


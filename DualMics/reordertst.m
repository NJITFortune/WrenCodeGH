sylnum = length(out);

% We need to identify the first occurance of each syllable type
for i = 1:sylnum;
        pp = sort([out(i).num]);
        id(i) = pp(1);
end;

idi = sort(id);

for j=1:sylnum;

    sylord(j) = find(id == idi(j));
    gout(j).num = sort(out(sylord(j)).num);

end;

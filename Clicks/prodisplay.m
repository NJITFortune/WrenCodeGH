function prodisplay(joy, ind);

boo = length(ind);

rows = ceil(boo / 5);

for i = 1:length(ind);
subplot(rows,5,i);
specgram(joy.song(joy.ind(ind(i),1):joy.ind(ind(i),2)),1024,joy.Fs,[],1000);
ylim([500 5000]);
end;


[sng, Fs] = audioread('/Users/eric/Documents/29_Dec_2012_012_04_sample_Bellavista.wav');

sng = sng(:,1); % One channel of a stereo recording
Fsylstart = []; Fsylend = []; Msylstart = []; Msylend = [];
figure(1); clf; hold on;

specgram(sng, 1024, Fs, [], 1000);
ylim([500 7500]);
colormap(flipud(gray));
caxis([-20 20]);

% Syllable 1, Male 
Msylstart(end+1) = 0.0081; Msylend(end+1) = 0.2113;
% Syllable 2, Female 
Fsylstart(end+1) = 0.2754; Fsylend(end+1) = 0.6254;
% Syllable 3, Male 
Msylstart(end+1) = 0.6463; Msylend(end+1) = 1.0887;
% Syllable 4, Female 
Fsylstart(end+1) = 1.1048; Fsylend(end+1) = 1.6565;
% Syllable 5, Male 
Msylstart(end+1) = 1.6823; Msylend(end+1) = 1.8855;
% Syllable 6, Female 
Fsylstart(end+1) = 1.9629; Fsylend(end+1) = 2.3145;
% Syllable 7, Male 
Msylstart(end+1) = 2.3371; Msylend(end+1) = 2.7726;

xlim([0.5 2.5]);

% Motif start and end
mstart = 0.6344; mend = 2.3247;
plot([mstart, mstart], [500, 7500], '-g', 'LineWidth', 0.5);
plot([mend, mend], [500, 7500], '-g', 'LineWidth', 0.5);


for j=1:length(Fsylstart)
    plot([Fsylstart(j), Fsylstart(j)], [500 7500], '-m', 'LineWidth', 0.2)
    plot([Fsylend(j), Fsylend(j)], [500 7500], '-m', 'LineWidth', 0.2)
end

for j=1:length(Msylstart)
    plot([Msylstart(j), Msylstart(j)], [500 7500], '-b', 'LineWidth', 0.2)
    plot([Msylend(j), Msylend(j)], [500 7500], '-b', 'LineWidth', 0.2)
end

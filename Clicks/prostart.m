function prostart(joy, ind, started);

ind(1)

colrs{5} = [1 1 0]; % yellow
colrs{4} = [1 0 1]; % magenta 
colrs{3} = [0 0 1]; % blue
colrs{6} = [0 1 1]; % cyan
colrs{2} = [1 0 0]; % red
colrs{1} = [0 0 0]; % black
colrs{3} = [0 1 0]; % green

colrs{7} = [0.3 0.3 0.3]; % grey (light black)
colrs{8} = [0.5 0 0]; % pink (light red)
colrs{9} = [0 0 0.5]; % light blue

colrs{10} = [0 0 1]; % blue [only for background]

% cr{1}='k*'; cr{2}='r*'; cr{3}='g*';
% cr{4}='m*'; cr{5}='y*'; cr{6}='c*';
% cr{7}='k+'; cr{8}='r+'; cr{9}='g+';

figure(2)

tim = joy.ind(1,1)/joy.Fs:1/joy.Fs:joy.ind(10,2)/joy.Fs;
%ns = ones(length(tim),1);
ns = ones(length(tim),1)*500;

if (started == 1);
%	subplot(2,1,1);
	specgram(joy.song(1:joy.ind(10,2)), 1024, joy.Fs, [], 1000); 
	ylim([0 5000]); caxis([-100 50]);

hold on;

%	subplot(2,1,2);
	area(tim, ns, 'FaceColor', colrs{10});

rea = joy.ind(ind(1),1)/joy.Fs:1/joy.Fs:joy.ind(ind(1),2)/joy.Fs;
%reao = ones(length(rea),1);
reao = ones(length(rea),1)*500;
area(rea, reao, 'FaceColor', colrs{started});
end

if (started > 1);
%	subplot(2,1,2);
	hold on;
rea = joy.ind(ind(1),1)/joy.Fs:1/joy.Fs:joy.ind(ind(1),2)/joy.Fs;
%reao = ones(length(rea),1);
reao = ones(length(rea),1)*500;
area(rea, reao, 'FaceColor', colrs{started});
end


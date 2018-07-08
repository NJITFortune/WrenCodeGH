function out = humanpicker(in);

close all;
adv = 1;
Fs = in.Fs;

% use a 5 second window

wind = 5;

numwin = ceil( max(in.time) / wind);

for i = 1:numwin;

	tt = find(in.time > wind*(i-1) & in.time <= wind*i);
	figure(i); set(gcf,'Position',[200+(i*10) 200+(i*10) 1200 400]);
	oscson(in.song(tt), Fs); xlim([0 wind]);

	sylind{i} = find(in.syl(:,2) > wind*(i-1) & in.syl(:,2) <= wind*i);

	figure(i); hold on;
	for j = 1:length(sylind{i})
		plot([ (in.syl(sylind{i}(j),1)-wind*(i-1)) (in.syl(sylind{i}(j),2)-wind*(i-1)) ], [4900 4900], 'w', 'LineWidth',2);
	end; hold off;

end

j = 1; 
while j < 99;

out.xs{j} = 0;

	p = 1;
	while p <= numwin;

		figure(p); [x y] = ginput;

		        figure(p); hold on;
			for k = 1:length(x);
				idx(k) = find(in.syl(:,1) < x(k)+(wind*(p-1)) & in.syl(:,2) > x(k)+(wind*(p-1)))
		                plot([ (in.syl(idx(k),1)-wind*(p-1)) (in.syl(idx(k),2)-wind*(p-1)) ], [4900 4900], 'g', 'LineWidth',2);
        		end; hold off;

		yy = input('If OK, type 0 or hit return, if not OK then any other number: ');
		if isempty(yy) == 1; yy = 0; end;

		if yy == 0; 
			out.xs{j} = [out.xs{j} idx]; 
                        figure(p); hold on;
			for k = 1:length(x);
                                idx(k) = find(in.syl(:,1) < x(k)+(wind*(p-1)) & in.syl(:,2) > x(k)+(wind*(p-1)))
                                plot([ (in.syl(idx(k),1)-wind*(p-1)) (in.syl(idx(k),2)-wind*(p-1)) ], [4900 4900], 'r', 'LineWidth',2);
                        end; hold off;
			p = p + 1; 
		end;
		
		if yy ~= 0;                         
			figure(p); hold on;
			for k = 1:length(x);
                                idx(k) = find(in.syl(:,1) < x(k)+(wind*(p-1)) & in.syl(:,2) > x(k)+(wind*(p-1)));
                                plot([ (in.syl(idx(k),1)-wind*(p-1)) (in.syl(idx(k),2)-wind*(p-1)) ], [4900 4900], 'w', 'LineWidth',2);
                        end; hold off;
		end;

		clear idx x
	end;

	endme = input('If you have labeled all of the syllables, type 99: ');
		if endme == 99; j = 99; end;

	j=j+1;
end

for f = 1:length(out.xs);
	out.xs{f} = out.xs{f}(2:end);
end

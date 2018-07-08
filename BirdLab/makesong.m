function song = makesong(syls, nops, seq);

% Example: seq = ['f1';'n1';'r2';'n2';'f11';'n3';'s5'];
% syntax for seq is "type""num"
% f for forward syllable
% r is for reverse syllable
% s is for a silence equal to the normal length of that syllable
% n for the silence following syllable num

song = [0 0];

for i = 1:length(seq)

%% Reverse case
	if seq(i,1) == 'r'
		ind = int8(seq(i,2:end)) - 48;
                if ind(2) < 0
                        ind = ind(1);
                else
                        ind = ind(2) + ind(1)*10;
                end
		syl = syls{ind(1)}';
		song = [song syl(end:-1:1)];
		clear syl;
	end

%% Forward case
	if seq(i,1) == 'f'
		ind = int8(seq(i,2:end)) - 48;
                if ind(2) < 0
                        ind = ind(1);
                else
                        ind = ind(2) + ind(1)*10;
                end
		syl = syls{ind(1)}';
		song = [song syl];
		clear syl;
	end

%% NOP case
	if seq(i,1) == 'n'
		ind = int8(seq(i,2:end)) - 48;
		if ind(2) < 0
			ind = ind(1);
		else
			ind = ind(2) + ind(1)*10;
		end
		syl = nops{ind(1)}';
		song = [song syl];
		clear syl;
	end

%% Silent Syllable case
	if seq(i,1) == 's'
                ind = int8(seq(i,2:end)) - 48;
                if ind(2) < 0
                        ind = ind(1);
                else
                        ind = ind(2) + ind(1)*10;
                end
		syl = 0*syls{ind(1)}';
		song = [song syl];
		clear syl;
	end

end
subplot(2,1,1); specgram(song,512,[],512,500);
subplot(2,1,2); plot(song);


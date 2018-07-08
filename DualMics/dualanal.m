function [ anal ] = dualanal( md, fd, Fs )
%[ anal ] = dualanal( md, fd, Fs )
%   The Kitchen Sink of analyses

% In general, we want to use data on the female microphone for analyzing
% data for females, and vice-versa for males.

% Get the syllable indices for male and female syllables
midx = find([md.sex] == 1); fidx = find([fd.sex] == 0);

% Get the indices for the syllable-to-syllable transitions
mTf = find(diff(midx) == 1); 
fTm = find(diff(fidx) == -1);
mTm = find(diff(midx) == 0 & midx(1:end-1) == 1);
fTf = find(diff(fidx) == 0 & fidx(1:end-1) == 0);

%% FIGURE 1
figure(1);
% Plot syllable durations when the animals sang them
subplot(221);
for i=length(midx):-1:1; mtmp(i) = md(midx(i)).syltim(1); end;
for i=length(fidx):-1:1; ftmp(i) = fd(fidx(i)).syltim(1); end;
hold off; plot(mtmp, [md(midx).sylen],'b-*'); 
hold on; plot(ftmp, [fd(fidx).sylen],'m-*'); 
% Plot syllable GCs when the animals sang them
subplot(222);
hold off; plot(mtmp, [md(midx).ofer_GC],'b-*'); 
hold on; plot(ftmp, [fd(fidx).ofer_GC],'m-*'); 
% Plot simple ISIs for both birds
for i=1:length(md); m2tmp(i) = md(i).syltim(2); end;
for i=1:length(fd); f2tmp(i) = fd(i).syltim(2); end;
subplot(223);
hold off; plot(m2tmp(1:end-1),[md.ISI], 'b-x');
hold on; plot(f2tmp(1:end-1),[fd.ISI], 'm-x');
retr = 0; if md(end).sex == 1; retr = 1; end;
for i=1:length(midx(1:end-retr)); 
    plot(md(midx(i)).syltim(2), md(midx(i)).ISI, 'bo');
end;
retr = 0; if fd(end).sex == 0; retr = 1; end;
for i=1:length(fidx(1:end-retr)); 
    plot(fd(fidx(i)).syltim(2), fd(fidx(i)).ISI, 'mo');
end;
% Plot syllable slope STD when the animals sang them
subplot(224);
hold off;
hold on;
for var=1:length(midx);
    entr = mean(md(midx(var)).ofer_mEntropy);
    plot(mtmp(var), entr,'b-*');
end;
for var=1:length(fidx);
    entr = mean(fd(fidx(var)).ofer_mEntropy);
    plot(ftmp(var), entr,'m-*');
end;
hold off;

%% Save it for later
anal.midx = midx; anal.fidx = fidx;

end


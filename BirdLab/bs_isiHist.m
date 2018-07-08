function isih = bs_isiHist( spiketimes, binum )
% isih = bs_isiHist( spiketimes, bin-num )
% Default bin-num is 30
% spiketimes are cell array from the bs toolset

% Check to see if user specified the number of bins
if nargin < 2;  binum = 30; end;

reps=length(spiketimes);

% binlims = logspace(-3, 0.176, binum); % From 1 msec to 1.5 seconds in log steps
binlims = logspace(-3, 0.398, binum); % From 1 msec to 2.5 seconds in log steps
xc = zeros(length(binlims),reps);

for j=1:reps;    
    isis = diff( spiketimes{j} ); % ISIs for this rep
    xc(:,j) = histc(isis, binlims);
end;

for i = 1:length(xc)-1;
    isih.hist(i) = sum(xc(i,:));
end;
    isih.tims = binlims(1:end-1);

    % plot(isih.tims,isih.hist, '*-');
    area(isih.tims,isih.hist);
        

%% Old version that has linear bins
% for j=1:reps
%     
%     isis{j} = diff( spiketimes{j} ); % ISIs
%     
%     for i = 1:binum;
%         t_strt = ((i*5)+1);
%         t_stp = ((i*5)+6);
%         repbinisi(i) = sum ( isis{j}*1000 > t_strt & isis{j}*1000 <= t_stp );
%         a(i) = a(i) + repbinisi(i);
%         tim(i) = t_strt;
%     end
%     
% end

%isih.hist=(a)';
%isih.tims=(tim);
%plot(isih.tims,isih.hist);
%semilogx(isih.tims, isih.hist,'b*');
%hist(isih.hist, binum);
xlabel('t (ms)');
ylabel('Count');

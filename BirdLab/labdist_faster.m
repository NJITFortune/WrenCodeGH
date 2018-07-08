function lf = labdist_faster(SpikesA,la,SpikesB,lb,q,k)
%LABDIST_FASTER(SA,LA,SB,LB,Q,K)
% Calculates the multi-unit metric distance between two spike trains
% Uses a fast version of the algorithm 
% SA, SB - spike times on the two spike trains
% LA, LB - spike labels (positive integers) 
% Q - timing precision parameter 
% K - label reassigning parameter 
%
% Dmitriy Aronov, 6/20/01; modified by Thomas Kreuz, 10/19/08 

%Assign labels in the form 1,2,...,L and count spikes of each label 
lbs = unique([la lb]); 
L = size(lbs,2);
for c = 1:L 
    j = find(la==lbs(c));
    la(j) = c; 
    numa(c) = size(j,2); 
    j = find(lb==lbs(c)); 
    lb(j) = c; 
    numb(c) = size(j,2); 
end
%Choose the spike train to separate to subtrains 
if prod(numb+1)*(sum(numa)+1) > prod(numa+1)*(sum(numb)+1)
    t = la; 
    la = lb; 
    lb = t; 
    t = SpikesA; 
    SpikesA = SpikesB; 
    SpikesB = t; 
    t = numa; 
    % numa=numb;
    numb = t; 
end
tb=zeros(L,max(numb));
for c = 1:L 
    tb(c,1:numb(c))=SpikesB(logical(lb==c));
end
%Set up an indexing system 
ind = []; 
for c = 1:L 
    j = repmat(0:numb(c),prod(numb(c+1:end)+1),1);
    j = repmat(reshape(j,numel(j),1),prod(numb(1:c-1)+1),1);
    ind = [ind j];
end
ind = sortrows([sum(ind,2) ind]); 
ind = ind(:,2:end); 
%Initialize the array 
m = zeros(size(ind,1),size(SpikesA,2)+1); 
m(1,:) = 0:size(SpikesA,2); 
m(:,1) = sum(ind,2); 
%Perform the calculation
for v = 2:size(m,1)
    fa2=find(m(:,1)==m(v,1)-1);
    fa=fa2(logical(sum(ind(fa2,:)-repmat(ind(v,:),length(fa2),1)==0,2)==L-1));
    fth=find(ind(v,:)>0)';
    bsv=diag(tb(fth,ind(v,fth)));
    for w = 2:size(m,2) 
        m(v,w)=min([m(v,w-1)+1; m(fa,w)+1; m(fa,w-1)+q*abs(SpikesA(w-1)-bsv)+k*not(la(w-1)==fth)]);
    end
end 
lf = m(end,end); 

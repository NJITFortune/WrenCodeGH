function H = JointEntropy(x, y, binSize)
 
%I suggest binSize=[range(x)/std(x)*5,range(y)/std(y)*5]
[N,C] = hist3([x y],binSize);
 
hx = C{1,1};
hy = C{1,2};
 
%xyPos=meshgrid(hx,hy);
 
% Normalize the area of the histogram to make it a pdf
N = N ./ sum(sum(N));
b=hx(2)-hx(1);
l=hy(2)-hy(1);
 
% Calculate the entropy
indices = N ~= 0;
H = -b*l*sum(N(indices).*log2(N(indices)));
 
end
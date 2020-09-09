%% Setup

distances = [0, 1, 2, 3, 5, 7, 9, 10];

vdists = [1, 3, 5, 7];

%% Get the indices


for j=1:length(vdists)
   
    distIDX = find([dd.distance] == vdists(j));
    blindidx = distIDX([dd(distIDX).vision] == 0);
    visionidx = distIDX([dd(distIDX).vision] == 1);
    
    v(j) = dd_getISidx(dd, visionidx);
    b(j) = dd_getISidx(dd, blindidx);
    
end

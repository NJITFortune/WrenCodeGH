function [colorMat,colorsToChoose]=colorChoose(numColors)

% colorRGB=readtable('/Users/daynf/dayta/Clustered/RColorRGB.csv','ReadVariableNames',0);
% colorRGB=table2cell(colorRGB);
load colorRGB.mat
colorsToChoose={'red','blue','blueviolet','darkgoldenrod1','chartreuse4','deeppink','cyan2','tan4','yellowgreen','thistle4','gray60','chocolate4','skyblue4','aquamarine4','forestgreen','mediumpurple1','azure3','rosybrown','salmon3','seashell2','cornflowerblue','maroon','peru','hotpink','gold4','darkorange'};

colorMat=zeros(numColors,3);

if numColors>length(colorsToChoose);
  d=numColors-length(colorsToChoose);
  r=randperm(length(f),d);
  colorsToChoose=[colorsToChoose colorRGB(f(r))];
end

for i=1:numColors
    
     clusterColorIdx = find(ismember(colorRGB(:,1), colorsToChoose{i}));
     colorMat(i,:)=[colorRGB{clusterColorIdx,2:4}]/255; %rgb
        
end
colorsToChoose=colorsToChoose(1:numColors);
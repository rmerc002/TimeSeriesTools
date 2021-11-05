function [canvas] = addContrast(panMP)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
numPercentiles = 10;
ps = linspace(0,100,numPercentiles+1)';
percentileValues = zeros(1,size(ps,1));

for index=2:size(ps,1)
    percentileValues(index) = prctile(panMP,ps(index),'all');
end

canvas = panMP;
            
%adjust contrast
tempCanvas = canvas;
for index=2:numPercentiles+1
    selection = percentileValues(index-1) < canvas & canvas <= percentileValues(index);
    tempPMP = canvas(selection);
    tempPMP = normalize(tempPMP,'range');
    tempCanvas(selection) = tempPMP*(1/(numPercentiles-1)) + (index-2)/(numPercentiles-1);
end
canvas = tempCanvas;

end


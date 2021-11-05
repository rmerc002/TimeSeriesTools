function plotDatasetFeatures(dataMatrix, outputPath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%assume number of features is less than the number of datapoints per sample
% if size(dataMatrix,1) < size(dataMatrix,2)
%    dataMatrix = dataMatrix';
% end

numFeatures = size(dataMatrix,2);
numSamples = size(dataMatrix,1);

fig = figure;
hold on;
for plotIndex = 1:numFeatures
   tempTS = dataMatrix(:,plotIndex);
   tempMin = nanmin(tempTS);
   tempMax = nanmax(tempTS);
   tempRange = nanmax(1e-8, tempMax - tempMin);
   tempS = -plotIndex + 0.9*(tempTS - tempMin)/tempRange;
   
   plot(tempS);
end
hold off;
set(gca, 'TickDir','out');
    box off;

if nargin == 2
    fileName = "platos";
    filePath = fullfile(outputPath, fileName + ".fig");
    savefig(fig, filePath);
    
    filePath = fullfile(outputPath, fileName + ".png");
    saveas(gcf, filePath);

    filePath = fullfile(outputPath, fileName + ".emf");
    print(filePath,'-dmeta'); 
end

end


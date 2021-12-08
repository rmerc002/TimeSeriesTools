function plotDatasetFeatures(dataMatrix, outputPath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%assume number of features is less than the number of datapoints per sample
% if size(dataMatrix,1) < size(dataMatrix,2)
%    dataMatrix = dataMatrix';
% end

numFeatures = size(dataMatrix,2);
numSamples = size(dataMatrix,1);

inset = 0.9;
fig = figure;
hold on;
for plotIndex = 1:numFeatures
   tempTS = dataMatrix(:,plotIndex);
   tempMin = min(tempTS, [], 'omitnan');
   tempMax = max(tempTS, [], 'omitnan');
   tempRange = tempMax - tempMin;
   tempS = -plotIndex + inset*(tempTS - tempMin)/tempRange;
   
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


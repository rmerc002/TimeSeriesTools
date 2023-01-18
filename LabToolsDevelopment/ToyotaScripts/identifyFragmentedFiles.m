inputPath = "C:\Users\Ryan\Dropbox\Education\UCR\Research\0_Projects\ActiveResearch\Toyota\Data\Car-following Data_2022-05-09\Trajectory";

inputFiles = dir(inputPath);

validFiles = false(length(inputFiles),1);
for ii = 1:length(inputFiles)
    if inputFiles(ii).name(1) ~= '.' && strcmp(inputFiles(ii).name(end-3:end), '.txt') == 1
        validFiles(ii) = true;
    end
end

inputFiles(~validFiles)= [];

startIndex = 1;
endIndex = 1;
for ii = 1:length(inputFiles)
    filePath = fullfile(inputPath, inputFiles(ii).name);
    dataTemp = importdata(filePath);
    
    [maxGap, maxIndex] = max(dataTemp(2:end,1) - dataTemp(1:end-1,1));
    if maxGap > 0.5
        fprintf("Fragmented timestamp: %s, maxGap: %.2f seconds, index:%d\n", inputFiles(ii).name, maxGap, maxIndex + 1);
    end
end
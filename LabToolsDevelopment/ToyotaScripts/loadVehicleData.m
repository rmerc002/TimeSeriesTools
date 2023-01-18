function dataConcat = loadVehicleData(inputPath)
inputFiles = dir(inputPath);

validFiles = false(length(inputFiles),1);
for ii = 1:length(inputFiles)
    if inputFiles(ii).name(1) ~= '.' && strcmp(inputFiles(ii).name(end-3:end), '.txt') == 1
        validFiles(ii) = true;
    end
end

inputFiles(~validFiles)= [];

fileBytes = zeros(length(inputFiles),1);
for ii = 1:length(fileBytes)
    fileBytes(ii) = inputFiles(ii).bytes;
end
% 
% [sortValues,sortIndices] = sort(fileBytes,'descend');
% 
% inputFiles = inputFiles(sortIndices);

filePath1 = fullfile(inputPath, inputFiles(1).name);
file1 = importdata(filePath1);
fileBytes1 = inputFiles(1).bytes;
fileRows1 = size(file1,1);

rowsPerBytes = fileRows1/fileBytes1;

totalBytes = sum(fileBytes);
totalRows = totalBytes*rowsPerBytes;
totalRowsBuffer = ceil(totalRows*1.05 + length(inputFiles));

dataConcat = nan(totalRowsBuffer, size(file1,2));

startIndex = 1;
endIndex = 1;
for ii = 1:length(inputFiles)
    if true || strcmp(inputFiles(ii).name(1:4),'base') == 1
        filePath = fullfile(inputPath, inputFiles(ii).name);
        dataTemp = importdata(filePath);
        endIndex = startIndex + size(dataTemp,1) - 1;
        dataConcat(startIndex:endIndex,:) = dataTemp;
        startIndex = endIndex + 2;
    end
end
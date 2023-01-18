% inputPath = "C:\Users\Ryan\Dropbox\Education\UCR\Research\0_Projects\ActiveResearch\Toyota\Data\Car-following Data_2022-05-09\Trajectory";

inputPath = "C:\Users\Ryan\Dropbox\Education\UCR\Research\0_Projects\ActiveResearch\Toyota\Data\Car-following Data_2022-05-17\Data_anom";
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

[sortValues,sortIndices] = sort(fileBytes,'descend');

inputFiles = inputFiles(sortIndices);

filePath1 = fullfile(inputPath, inputFiles(1).name);
file1 = importdata(filePath1);
fileBytes1 = inputFiles(1).bytes;
fileRows1 = size(file1,1);

rowsPerBytes = fileRows1/fileBytes1;

totalBytes = sum(sortValues);
totalRows = totalBytes*rowsPerBytes;
totalRowsBuffer = ceil(totalRows*1.05 + length(inputFiles));

dataConcatBase = nan(totalRowsBuffer, size(file1,2));

startIndex = 1;
endIndex = 1;
for ii = 1:length(inputFiles)
    if strcmp(inputFiles(ii).name(1:4),'base') == 1
        filePath = fullfile(inputPath, inputFiles(ii).name);
        dataTemp = importdata(filePath);
        endIndex = startIndex + size(dataTemp,1) - 1;
        dataConcatBase(startIndex:endIndex,:) = dataTemp;
        startIndex = endIndex + 2;
    end
end

dataConcatAnom = nan(totalRowsBuffer, size(file1,2));

startIndex = 1;
endIndex = 1;
for ii = 1:length(inputFiles)
    if strcmp(inputFiles(ii).name(1:4),'anom') == 1
        filePath = fullfile(inputPath, inputFiles(ii).name);
        dataTemp = importdata(filePath);
        endIndex = startIndex + size(dataTemp,1) - 1;
        dataConcatAnom(startIndex:endIndex,:) = dataTemp;
        startIndex = endIndex + 2;
    end
end

columnNumber = 5;

% upsampleRate = 1;%5;
% posnan = dataConcatAnom(:,columnNumber);
% % posnan = resample(posnan,upsampleRate,1);
% posnan = posnan + randn(length(posnan),1)*1e-5;
% 
% negnan = dataConcatBase(:,columnNumber);
% % negnan = resample(negnan,upsampleRate,1);
% negnan = negnan + randn(length(negnan),1)*1e-5;
% 
% pos = posnan;
% pos(isnan(posnan)) = mean(posnan,'omitnan');
% % pos = pos + randn(length(pos),1)*1e-5;
% 
% neg = negnan;
% neg(isnan(negnan)) = mean(negnan,'omitnan');
% % neg = neg + randn(length(neg),1)*1e-5;
% 
% mm = 10*upsampleRate;
% threshold = 0.25;
% ep = OnlineEmergenceProfile(posnan, negnan, mm, threshold);
% ep.plot();
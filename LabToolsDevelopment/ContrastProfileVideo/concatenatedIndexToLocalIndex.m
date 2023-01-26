function [index, tsLocal] = concatenatedIndexToLocalIndex(index, ts, fileLengths)
    if nargin < 3
        fileLengths = [inf];
    end
    %%% assume files are separated by a nan
    fileStartIndex = 1;
    fileEndIndex = fileStartIndex + fileLengths(1) -1;
    for ii = 1:length(fileLengths)
        if index <= fileLengths(ii)
            break;
        else
            index = index - fileLengths(ii) - 1;%-1 for nan
            fileStartIndex = fileStartIndex + fileLengths(ii) + 1;
            fileEndIndex = fileStartIndex + fileLengths(ii) -1;
        end
    end
    tsLocal = ts(fileStartIndex:fileEndIndex);
end
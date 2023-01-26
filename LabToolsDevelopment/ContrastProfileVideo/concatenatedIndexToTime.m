function [timeStr, fileIndex] = concatenatedIndexToTime(behaviorIndex, sampleRate, fileLengths)
    %%% assume files are separated by a nan
    if nargin < 3
        fileLengths = [inf];
    end
    
    for ii = 1:length(fileLengths)
        if behaviorIndex <= fileLengths(ii)
            fileIndex = ii;
            break;
        else
            behaviorIndex = behaviorIndex - fileLengths(ii) - 1;%-1 for nan
        end
    end

    totalSeconds = seconds(behaviorIndex/sampleRate);
    totalSeconds.Format = "hh:mm:ss.SS";
    timeStr = char(totalSeconds);
end
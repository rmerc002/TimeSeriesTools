function [segLengths, csld] = continuousSegmentDistribution(data, maxSegLength)
%%%Finds distribution of segments
%%%expects 1D data
segLengths = zeros(1,sum(isnan(data)));
csld = zeros(1,maxSegLength);
startIndex = 1;
while isnan(data(startIndex))
    startIndex = startIndex + 1;
end
segmentIndex = 1;
for endIndex = startIndex + 1:length(data)
    if startIndex > length(data)
        break;
    end
    if endIndex > startIndex && isnan(data(endIndex))
        tempLength = endIndex -startIndex;

        indexLength = min(tempLength, maxSegLength);
        csld(tempLength) = csld(indexLength) + 1;

        segLengths(segmentIndex) = tempLength;
        segmentIndex = segmentIndex+1;
        
        startIndex = endIndex + 1;
        while isnan(data(startIndex)) && startIndex < length(data)
            startIndex = startIndex + 1;
        end
    end
    
end

end
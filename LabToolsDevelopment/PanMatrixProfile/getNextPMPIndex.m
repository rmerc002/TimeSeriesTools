function nextIndex = getNextPMPIndex(PMP, prevIndices)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    prevIndices = sort(prevIndices);
    maxDist = 0;
    nextIndex = 0;
    for metaIdx = 2:length(prevIndices)
        idx1 = prevIndices(metaIdx-1);
        idx2 = prevIndices(metaIdx);
        if idx1 + 1 < idx2 %can't search between integers
            thisDist = PMP(idx2,:)-PMP(idx1,:);
            thisDist = thisDist(~isnan(thisDist));
            thisDist = norm(thisDist)/sqrt(2*length(thisDist));
            if thisDist > maxDist
                maxDist = thisDist;
                nextIndex = ceil((idx2+idx1)/2);
            end
        end
    end
end


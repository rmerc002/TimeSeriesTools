%%
% Returns a binary tree traversal of the integers from 1 to n
function [IDX] = MaxSplit(Range, PMP)
IDX = [1,length(Range)];
prevIndices = [1,length(Range)];
for i = 3:length(Range)
    metaIdx = getNextIndex(PMP, Range, prevIndices);
    IDX = [IDX,metaIdx];
    prevIndices = sort([prevIndices,metaIdx]);
end
end

function nextIndex = getNextIndex(PMP, range, prevIndices)
    prevIndices = sort(prevIndices);
    maxDist = 0;
    nextIndex = 0;
%     nextIndex = find((prevIndices(2:end)-prevIndices(1:end-1))-1>0);
%     if length(nextIndex) == 0
%         return
%     end
%     nextIndex = nextIndex(1)+1;
    for metaIdx = 2:length(prevIndices)
        idx1 = prevIndices(metaIdx-1);
        idx2 = prevIndices(metaIdx);
        if idx1 + 1 < idx2 %can't search between integers
            thisDist = PMP(range(idx2),:)-PMP(range(idx1),:);
            thisDist = thisDist(~isnan(thisDist));
            thisDist = norm(thisDist)/sqrt(2*length(thisDist));
            if thisDist > maxDist
                maxDist = thisDist;
                nextIndex = ceil((idx2+idx1)/2);
            end
        end
    end
end

function [sortedLowestIndices] = allLowestDistanceIndices_V06(ts, subLength, exclusionLength)
%%%by Ryan Mercer
%%%For Contrast Profile Experiments
%%%   Sara Alaee, PLATO Verificaiton Tool
%%%Specific Alterations:
%%%   Finds all possible subsequences with low value
%%%   Returns indices rather than subsequences
%%%   Exclusion zone is subLength
%%%   Oriented to select the lowest value first.

if isempty(exclusionLength)
    exclusionLength = ceil(subLength/2);
end

N = 2*floor(length(ts)/(subLength));

A = zeros(2,length(ts));
A(1,:) = 1:length(ts);
A(2,:) = (-ts)';
A = A';
B = flip(sortrows(A,2),1);

exclusionZone = zeros(1,length(ts));

sortedLowestIndices = nan(N,1);
index = 1;
iterNum = 1; %for debugging
while index <= N && iterNum < size(B,1)
    trialIndex = B(iterNum,1);
    if exclusionZone(trialIndex) == 0
        sortedLowestIndices(index) = trialIndex;
        index = index + 1;
        leftOffset = exclusionLength;
        rightOffset = exclusionLength;
        exclusionZone(max(1,trialIndex - leftOffset):min(length(ts), trialIndex + rightOffset)) = 1;
    end
    iterNum = iterNum + 1;
end

sortedLowestIndices = sortedLowestIndices(~isnan(sortedLowestIndices));
end


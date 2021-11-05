function [bestIndices] = KLowestDistanceIndices_V05(ts, subLength, K, exclusionLength)
%%%by Ryan Mercer
%%%For Contrast Profile Experiments
%%%   Sara Alaee, PLATO Verificaiton Tool
%%%Specific Alterations:
%%%   Uses top K
%%%   Returns indices rather than subsequences
%%%   Exclusion zone is subLength
%%%   Oriented to select the lowest value first.

if isempty(exclusionLength)
    exclusionLength = ceil(subLength);
end

A = zeros(2,length(ts));
A(1,:) = 1:length(ts);
A(2,:) = (-ts)';
A = A';
B = flip(sortrows(A,2),1);

exclusionZone = zeros(1,length(ts));

bestIndices = nan(K,1);
index = 1;
iterNum = 1; %for debugging
while index <= K && iterNum < size(B,1)
    trialIndex = B(iterNum,1);
    if exclusionZone(trialIndex) == 0
        bestIndices(index) = trialIndex;
        index = index + 1;
        leftOffset = exclusionLength;
        rightOffset = exclusionLength;
        exclusionZone(max(1,trialIndex - leftOffset):min(length(ts), trialIndex + rightOffset)) = 1;
    end
    iterNum = iterNum + 1;
end

bestIndices = bestIndices(~isnan(bestIndices));
end


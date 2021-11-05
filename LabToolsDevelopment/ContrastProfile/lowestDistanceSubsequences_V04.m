function [bestSubsequences] = lowestDistanceSubsequences_V04(ts, subLength, K)
%%%by Ryan Mercer
%%%For Contrast Profile Experiments
%%%   Alireza Chicken data
%%%Specific Alterations:
%%%   Uses top K
%%%   Returns subsequences rather than indices
%%%   Exclusion zone is subLength/4
%%%   Oriented to select the lowest value first.


exclusionLength = ceil(subLength/4);

A = zeros(2,length(ts));
A(1,:) = 1:length(ts);
A(2,:) = (-ts)';
A = A';
B = flip(sortrows(A,2),1);

exclusionZone = zeros(1,length(ts));

bestSubsequences = nan(K,subLength);
index = 1;
iterNum = 1; %for debugging
while index <= K && iterNum < size(B,1)
    trialIndex = B(iterNum,1);
    if exclusionZone(trialIndex) == 0
        startIndex = trialIndex;
        endIndex = min(length(ts), trialIndex+subLength-1);
        sampleLength = endIndex - startIndex + 1;
        bestSubsequences(index,1:sampleLength) = zscore(ts(startIndex:endIndex));
        index = index + 1;
        leftOffset = exclusionLength;
        rightOffset = exclusionLength;
        exclusionZone(max(1,trialIndex - leftOffset):min(length(ts), trialIndex + rightOffset)) = 1;
    end
    iterNum = iterNum + 1;
end
end


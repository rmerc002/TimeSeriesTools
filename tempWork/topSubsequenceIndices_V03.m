function [topIndices] = topSubsequenceIndices_V03(ts, subsequenceLength)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

exclusionLength = ceil(subsequenceLength);
N = 2*floor(length(ts)/(subsequenceLength));

A = zeros(2,length(ts));
A(1,:) = 1:length(ts);
A(2,:) = ts';
A = A';
B = flip(sortrows(A,2),1);

exclusionZone = zeros(1,length(ts));

topIndices = nan(1,N);
index = 1;
iterNum = 1; %for debugging
while index <= N && iterNum < size(B,1)
    trialIndex = B(iterNum,1);
    if exclusionZone(trialIndex) == 0
        topIndices(index) = trialIndex;
        index = index + 1;
        leftOffset = exclusionLength;
        rightOffset = exclusionLength;
        exclusionZone(max(1,trialIndex - leftOffset):min(length(ts), trialIndex + rightOffset)) = 1;
    end
    iterNum = iterNum + 1;
end

topIndices = topIndices(~isnan(topIndices));
end


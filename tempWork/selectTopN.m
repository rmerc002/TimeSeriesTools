function [topNIndices] = selectTopN(N, distances, exclusionLength)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%Choose Top 11
A = zeros(2,length(distances));
A(1,:) = 1:length(distances);
A(2,:) = distances';
A = A';
B = sortrows(A,2);

exclusionZone = zeros(1,length(distances));

topNIndices = nan(1,N);
index = 1;
iterNum = 1;
while index <= N
    trialIndex = B(iterNum,1);
    if exclusionZone(trialIndex) == 0
        topNIndices(index) = trialIndex;
        index = index + 1;
        exclusionZone(max(1,ceil(trialIndex-exclusionLength)):ceil(trialIndex+exclusionLength)) = 1;
    end
    iterNum = iterNum + 1;
end

end


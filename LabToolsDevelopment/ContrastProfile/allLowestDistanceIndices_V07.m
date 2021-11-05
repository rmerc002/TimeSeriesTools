function [sortedLowestIndices, correspondingDistances] = allLowestDistanceIndices_V07(distanceProfile, subLength, exclusionLength, K)
%%%by Ryan Mercer
%%%For Contrast Profile Experiments
%%%   Sara Alaee, PLATO Verificaiton Tool
%%%Specific Alterations:
%%%   Finds all possible subsequences with low value
%%%   Returns indices rather than subsequences
%%%   Exclusion zone is subLength
%%%   Oriented to select the lowest value first.

if nargin == 2
    exclusionLength = ceil(subLength/2);
end
if nargin <= 3
    K = 2*floor(length(distanceProfile)/(subLength));
elseif nargin < 2 || nargin > 4
    error('incorrect number of input arguments');
end



A = zeros(2,length(distanceProfile));
A(1,:) = 1:length(distanceProfile);
A(2,:) = (-distanceProfile)';
A = A';
B = flip(sortrows(A,2),1);

exclusionZone = zeros(1,length(distanceProfile));

sortedLowestIndices = nan(K,1);
correspondingDistances = nan(K,1);
index = 1;
iterKum = 1; %for debugging
while index <= K && iterKum < size(B,1)
    trialIndex = B(iterKum,1);
    if exclusionZone(trialIndex) == 0
        sortedLowestIndices(index) = trialIndex;
        correspondingDistances(index) = distanceProfile(trialIndex);
        index = index + 1;
        leftOffset = exclusionLength;
        rightOffset = exclusionLength;
        exclusionZone(max(1,trialIndex - leftOffset):min(length(distanceProfile), trialIndex + rightOffset)) = 1;
    end
    iterKum = iterKum + 1;
end

sortedLowestIndices = sortedLowestIndices(~isnan(sortedLowestIndices));
correspondingDistances = correspondingDistances(~isnan(sortedLowestIndices));
end


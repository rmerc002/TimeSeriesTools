function [sortedLowestIndices, correspondingDistances] = NearestNeighborSelection(distanceProfile, subLength, K)
%%% Returns the indices of the lowest values
%%%   sorted by lowest values
%%%   using excluson zone around matches
%%%   K (optional): return Top-K min indices. Default return all possible indices, 
%%%   WARNING: may not
exclusionLength = ceil(subLength/2);

KMax = 2*floor(length(distanceProfile)/(subLength));

if nargin == 2
    K = KMax;
elseif nargin < 2 || nargin > 4
    error('incorrect number of input arguments');
end
if K > KMax
   fprintf("Warning: K is %d, which is larger than the max possible %d. Setting K=%d\n",K,KMax,KMax);
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
    if exclusionZone(trialIndex) == 0 && ~isnan(distanceProfile(trialIndex))
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


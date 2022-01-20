function [peakIndices, peakValues] = exclusionZonePeaks(profile, m, exclusionLength, K)
%%% Returns the indices of the lowest values
%%%   sorted by lowest values
%%%   using excluson zone around matches
%%%   K (optional): return Top-K min indices. Default return all possible indices, 
%%%   WARNING: may not
if nargin < 3
    exclusionLength = ceil(m/2);
end

numSubsequences = length(profile) - m + 1;
maxNumPeaks = min(K, ceil(2*numSubsequences/exclusionLength));


A = zeros(2,length(profile));
A(1,:) = 1:length(profile);
A(2,:) = (profile)';
A = A';
B = flip(sortrows(A,2),1);

exclusionZone = zeros(1,length(profile));

peakIndices = nan(maxNumPeaks,1);
peakValues = nan(maxNumPeaks,1);
index = 1;
Kindex = 1; %for debugging
while index <= maxNumPeaks && Kindex < size(B,1)
    trialIndex = B(Kindex,1);
    if exclusionZone(trialIndex) == 0 && ~isnan(profile(trialIndex))
        peakIndices(index) = trialIndex;
        peakValues(index) = profile(trialIndex);
        index = index + 1;
        leftOffset = exclusionLength;
        rightOffset = exclusionLength;
        exclusionZone(max(1,trialIndex - leftOffset):min(length(profile), trialIndex + rightOffset)) = 1;
    end
    Kindex = Kindex + 1;
end

peakIndices = peakIndices(~isnan(peakIndices));
peakValues = peakValues(~isnan(peakIndices));

[peakIndices, sortIndices] = sort(peakIndices);
peakValues = peakValues(sortIndices);
end


function [PMPKNN, PMPiKNN] = mpxABKNN(tsA, tsB, subLength, K)
exclusionLength = subLength;

dataOrientation = 0; %0 for row, 1 for column
%Change to row vector for internal use
if size(tsA,1) == 1 
    %save the data orientation for matching output format to input
    %give orientation priority to positiveTS
    dataOrientation = 1; 
%     tsA = tsA';
else
    tsA = tsA';
end
if size(tsB,1) ~= 1 
    %give orientation priority to positiveTS, do not save negative
    %orientation if different than positive
    tsB = tsB';
end
    
    
tsANonNan = tsA;
tsANonNan(isnan(tsA)) = nanmean(tsA);

tsBNonNan = tsB;
tsBNonNan(isnan(tsB)) = nanmean(tsB);

%%%This is just for time refence. The output is not used
% tic;
% [matrixProfile_AB, matrixProfileIdx_AB] = mpx_ABBA_v2(tsANonNan, tsBNonNan, subLength);
% endTime = toc;
% fprintf("Matrix Profile completed in %.2f seconds\n",endTime);

% tic;
% K = 100;
PMPKNN = nan(K, length(tsANonNan));
PMPiKNN = nan(K, length(tsANonNan));
for indexTime = 1:length(tsANonNan) - subLength + 1
   subsequence = tsANonNan(indexTime:indexTime + subLength - 1);
   if sum(isnan(subsequence)) > 0
       continue;
   end
   distProfile = real(MASS_V2(tsBNonNan, subsequence));
   
   %%%TODO: get candidates so exclusion zones are used
   [sortedIndices, sortedDist] = allLowestDistanceIndices(distProfile, subLength, exclusionLength);
   KMin = min(K, length(sortedIndices));
   PMPKNN(1:KMin, indexTime) = sortedDist(1:KMin);
   PMPiKNN(1:KMin, indexTime) = sortedIndices(1:KMin);
end

% endTime = toc;
% fprintf("MP Frequency completed in %.2f seconds\n",endTime);


% visualizeMMPAB(tsANonNan, PMPKNN, PMPiKNN, [], [], [], 100, 1:100, "KNN","");

%%% Now output the best candidate for each number of nearest neighbors
% [bestCandidateDistances, minIndices] = min(PMPKNN, [],2);
% plotSortedBinnedSamples(minIndices, tsANonNan, [], subLength, 10,10,"");




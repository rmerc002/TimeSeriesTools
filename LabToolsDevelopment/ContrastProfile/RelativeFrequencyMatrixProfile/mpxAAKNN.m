function [PMPKNN, PMPiKNN] = mpxAAKNN(ts, subLength, K)
exclusionLength = subLength;

% ts = tsMale(1:100000);
% ts = randn(10000,1);
% sinIndices = 1:1000:10000-subLength;
% for i = 1:length(sinIndices)
%     si = sinIndices(i);
%    ts(si:si+subLength-1) = sin(linspace(0,2*pi,subLength)); 
% end

dataOrientation = 0; %0 for row, 1 for column
%Change to row vector for internal use
if size(ts,1) == 1 
    %save the data orientation for matching output format to input
    %give orientation priority to positiveTS
    dataOrientation = 1; 
%     ts = ts';
else
    ts = ts';
end


tsNonNan = ts;
tsNonNan(isnan(ts)) = nanmean(ts);

%%%This is just for time refence. The output is not used
% tic;
% [matrixProfile_AA, matrixProfileIdx_AB, isvalidwindow, motifsIdx, discordsIdx] = mpx_v3(tsNonNan, subLength, subLength, false);
% endTime = toc;
% fprintf("Matrix Profile completed in %.2f seconds\n",endTime);

% tic;
% K = 100;
PMPKNN = nan(K, length(tsNonNan));
PMPiKNN = nan(K, length(tsNonNan));
for indexTime = 1:length(tsNonNan) - subLength + 1
    subsequence = tsNonNan(indexTime:indexTime + subLength - 1);
    if sum(isnan(subsequence)) > 0
       continue;
   end
   distProfile = MASS_V2(tsNonNan, subsequence);
   
   %%%Exlusion zone around the current subsequence
   startIndex = max(1, indexTime-subLength);
   endIndex = min(length(tsNonNan), indexTime + subLength - 1);
   distProfile(startIndex:endIndex) = inf;
   %%%TODO: get candidates so exclusion zones are used
   [sortedIndices, sortedDist] = allLowestDistanceIndices(distProfile, subLength, exclusionLength);
   KMin = min(K, length(sortedIndices));
   PMPKNN(1:KMin, indexTime) = sortedDist(1:KMin);
   PMPiKNN(1:KMin, indexTime) = sortedIndices(1:KMin);
end
% PMPMean = nan(size(PMPKNN));
% for indexK = 1:K
%    PMPMean(indexK, :) = mean(PMPKNN(1:indexK,:), 1); 
% end
% endTime = toc;
% fprintf("MP Frequency completed in %.2f seconds\n",endTime);


% visualizeMMPAB(tsNonNan, PMPKNN, PMPiKNN, [], [], [], 100, 1:100, "KNN","");
% plotDatasetFeatures(motifCandidateScores(1:100,:)')

%%% Now output the best candidate for each number of nearest neighbors
% [bestCandidateDistances, minIndices] = min(PMPKNN, [],2);
% plotSortedBinnedSamples(minIndices, tsNonNan, [], subLength, 10,10,"");




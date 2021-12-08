function [RFMP, RFMPi] = RelativeFrequencyMatrixProfile(tsA, tsB, m, maxFreq)

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


    
selfFlag = false;
if isempty(tsB)
   selfFlag = true;
   tsB = tsA;
end

if ~(size(tsB,1) == 1 )
%give orientation priority to positiveTS, do not save negative
%orientation if different than positive
    tsB = tsB';
end

if  length(tsA) == length(tsB) && mean(tsA == tsB)>0.9
    selfFlag = true;
end
    
tsANonNan = tsA;
tsANonNan(isnan(tsA)) = nanmean(tsA);

tsBNonNan = tsB;
tsBNonNan(isnan(tsB)) = nanmean(tsB);



RFMP = nan(maxFreq, length(tsANonNan));
RFMPi = nan(maxFreq, length(tsANonNan));
for indexTime = 1:length(tsANonNan) - m + 1
   subsequence = tsA(indexTime:indexTime + m - 1);
   if sum(isnan(subsequence)) > 0
       continue;
   end
   DP = real(MASS_V2(tsBNonNan, subsequence));
   for nanIndex = 1:length(tsB)-m+1
       if isnan(tsB(nanIndex))
        startIndex = max(1,nanIndex-m+1);
        DP(startIndex:nanIndex) = nan;
       end
   end
   
   %%%TODO: get candidates so exclusion zones are used
   [sortedIndices, sortedDist] = NearestNeighborSelection(DP, m, maxFreq+selfFlag);
   KMin = min(maxFreq+selfFlag, length(sortedIndices));
   RFMP(1:KMin-selfFlag, indexTime) = sortedDist(1+selfFlag:KMin);
   RFMPi(1:KMin-selfFlag, indexTime) = sortedIndices(1+selfFlag:KMin);
end

% endTime = toc;
% fprintf("MP Frequency completed in %.2f seconds\n",endTime);


% visualizeMMPAB(tsANonNan, PMPKNN, PMPiKNN, [], [], [], 100, 1:100, "KNN","");

%%% Now output the best candidate for each number of nearest neighbors
% [bestCandidateDistances, minIndices] = min(PMPKNN, [],2);
% plotSortedBinnedSamples(minIndices, tsANonNan, [], subLength, 10,10,"");




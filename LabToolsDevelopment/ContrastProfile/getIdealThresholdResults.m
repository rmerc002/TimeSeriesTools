function [tpIndices, fpIndices, fnIndices, idealThreshold] = getIdealThresholdResults(distanceProfile, groundTruthIndices, subLength, figureSavePath, K, queryOffset)
%%% I'm going to temporarily run this independently of the test code
%%% I will save results in the test level directory
%%% 
%%% Then this will be the output of running the code natively in the test
%%% code

%%% I can then generate another path script that fetches all the results
%%%   at the test directory level and calculates statistics (mean, median,
%%%   std)

%%% How will I write this function
%%% 1) get candidates
%%% 2) keep candidates(1:numSignals)
%%% 3) calculate the tp fp fn results

    if nargin < 4
        figureSavePath = "";
    end
    if nargin < 5
        K = length(groundTruthIndices);
    end
    if nargin < 6
        queryOffset = 0;
    end
    if nargin < 3 || nargin > 6
        error('incorrect number of input arguments');
    end
    
    %%%Bug warning: There may be an issue if fewer than K indices present
%     candidateIndices = topSubsequenceIndices_V03(-profileDistances, subLength);
    exclusionLength = ceil(subLength/2);
%     KLowestDistanceIndices_V05(profileDistances, subLength, K, exclusionLength);
    candidateIndices = allLowestDistanceIndices_V06(distanceProfile, subLength, exclusionLength);
    topCandidateIndices = candidateIndices(1:K);
    topCandidateIndices = max(1,min(length(distanceProfile),topCandidateIndices));

    idealThreshold = (distanceProfile(candidateIndices(K)) + distanceProfile(candidateIndices(K+1)))/2;

    tpIndices = [];
    fpIndices = [];
    fnIndices = [];
    

    groundTruthFound = zeros(length(groundTruthIndices),1);

    for i1 = 1:length(topCandidateIndices)
        ri = topCandidateIndices(i1);
        isResultMatched = false;
       for i2 = 1:length(groundTruthIndices)
           gti = groundTruthIndices(i2);
           %%% If the position of the pulse's max amplitude is in the middle
            %%%   of plato, then the offset should be set to this index.
            %%% The queryOffset is unchanging for an identified plato.
            %%% One alternative is to add queryOffset to topCandidateIndices,
            %%%   but I thought that may be misleading for those used to the
            %%%   convention of reading distanceProfiles.
           absDistanceBetweenCandidateAndGroundTruth = abs(ri + queryOffset - gti);
           if groundTruthFound(i2) == 0 && absDistanceBetweenCandidateAndGroundTruth < subLength
                   tpIndices = [tpIndices,ri];
                   groundTruthFound(i2) = 1;
                   isResultMatched = true;
                   break;
           end
%            if groundTruthFound(i2) == 1 && abs(ri - gti) < ceil(subLength/2)
%                    fpIndices = [fpIndices,ri];
%                    isResultMatched = true;
%                    break;
%            end
           %%%if groundTruthFound(i2) == 0 && abs(ri - gti) > subLength
           %%% then I don't care

           %%% if groundTruthFound(i2) == 1 && abs(ri - gti) > subLength
           %%% then I don't care
       end
       if isResultMatched == false
           fpIndices = [fpIndices, ri];
       end
    end

    fnIndices = groundTruthIndices(groundTruthFound == 0);
    [~,sortedIndices] = sort(distanceProfile(fnIndices));
    fnIndices = fnIndices(sortedIndices);
    
    
    
    if ~isempty(figureSavePath) || figureSavePath == ""
        dataSetFileName = "internVars_idealThresholdVals";
        dataSetFilePath = fullfile(figureSavePath, dataSetFileName + ".mat");
        save(dataSetFilePath, 'tpIndices','fpIndices','fnIndices','idealThreshold');
    end
end
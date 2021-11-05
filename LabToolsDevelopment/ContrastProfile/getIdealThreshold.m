function [idealThreshold] = getIdealThreshold(queryShapelet,  testTS, subLength, numSignals, figureSavePath)
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

    profileDistances = real(mpx_AB_shapeletDiscovery(testTS, queryShapelet, subLength));
    
    candidateIndices = topSubsequenceIndices_V03(-profileDistances, subLength);

    idealThreshold = (profileDistances(candidateIndices(numSignals)) + profileDistances(candidateIndices(numSignals+1)))/2;

    %%% Save Variables
    
    dataSetFileName = "internVars_idealThresholdVals";
    dataSetFilePath = fullfile(figureSavePath, dataSetFileName + ".mat");
    save(dataSetFilePath, 'idealThreshold');
    
end
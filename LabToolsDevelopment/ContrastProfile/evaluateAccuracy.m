function [tpIndices, fpIndices, fnIndices] = evaluateAccuracy(positiveTS, subLength, groundTruthIndices, queryShapelet, queryThreshold, thresholdConfidence, figureSavePath)

    if size(positiveTS,1) == 1
        positiveTS = positiveTS';
    end
    if size(queryShapelet,1) == 1
       queryShapelet = queryShapelet'; 
    end
    
    numSignals = length(groundTruthIndices);
    
    activeLength = length(positiveTS) - subLength + 1;
    defaultRate = (activeLength - numSignals * 2*subLength)/activeLength;
    fprintf("Defualt Rate: single hit: %.0f%%, single miss: %.0f%%\n", (1-defaultRate)*100, defaultRate*100);

    MP_AQ = real(mpx_AB_shapeletDiscovery(positiveTS, queryShapelet, length(queryShapelet)));
    [tpCount,fpCount,fnCount, tpIndices, fpIndices, fnIndices] = getAccuracyMetrics(positiveTS, queryShapelet, MP_AQ, queryThreshold, subLength, groundTruthIndices, figureSavePath);
    fprintf("True Positives: %d, False Positives: %d, False Negatives: %d\n", tpCount, fpCount, fnCount);
    
    %%%This is redudant for test 18, should be ignored afterwards
    %%%TODO: remove after test 18
    fileName = "internVars_findingIdealThreshold";
    filePath = fullfile(figureSavePath, fileName + ".mat");
    save(filePath, 'groundTruthIndices', 'MP_AQ');
    
    
    
    %%% assure var path exists
    varsDirPath = fullfile(figureSavePath, 'Vars');
    if ~isfolder(varsDirPath)
    	mkdir(varsDirPath);
    end

    %%% I did not think any internal vars were worthwhile
    %%% internal vars
    dataSetFileName = "evaluateAccuracy_internal";
    internalVarsFilePath = fullfile(varsDirPath,dataSetFileName + ".mat");
    save(internalVarsFilePath, 'MP_AQ','groundTruthIndices');

    %%% output vars
    dataSetFileName = "evaluateAccuracy_output";
    outputVarsFilePath = fullfile(varsDirPath,dataSetFileName + ".mat");
    save(outputVarsFilePath, 'tpIndices', 'fpIndices', 'fnIndices');
end
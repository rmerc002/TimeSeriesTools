function [queryShapelet, queryThreshold, thresholdConfidence] = shapeletDiscovery(positiveTS, negativeTS, subLength, groundTruthIndices, figureSavePath)
showPlot = false;
    rng(1);
    %  minLength = 3;
    %  maxLength = ceil(positiveTS/3);
     %%% Temporary parameters that should be taken care of automatically:
     %%% 1) subLength
     %%% 2) numSignals

     numSignals = length(groundTruthIndices);
     if size(positiveTS,1) == 1
         positiveTS = positiveTS';
     end
     if size(negativeTS,1) == 1
        negativeTS = negativeTS'; 
     end
      %%% For debugging and visualizaiton
    %  groundTruthIndices = [14943       39910       45067       82151      101176];
     if showPlot == true
        figure;
        subplot(2,1,1);
        plot(negativeTS);

        subplot(2,1,2);
        plot(positiveTS);
        hold on;
        for i = groundTruthIndices
        plot(i:i+subLength-1,positiveTS(i:i+subLength-1),'green');
        end
        hold off;
     end

%     fprintf("\t#$#$FigureSavePath: %s\nS",figureSavePath);
     %%% Main Work
     shapeletCandidateIndices = getFilteredShapeletCandidates(positiveTS, negativeTS, subLength, groundTruthIndices, figureSavePath); %two steps in one function
     queryShapelet = getCentralSubsequenceOfTopCluster(positiveTS, shapeletCandidateIndices, subLength, numSignals, figureSavePath); %two steps in one function
     [queryThreshold, thresholdConfidence] = getQueryThresholdAndConfidence(positiveTS, queryShapelet, numSignals, figureSavePath); %two steps in one function

     fprintf("Threshold: %.2f, Confidence: %.2f\n",queryThreshold, thresholdConfidence);
     
     dataSetFileName = "z_outputVars_shapeletDiscovery";
     dataSetFilePath = fullfile(figureSavePath, 'Vars');
     if ~isfolder(dataSetFilePath)
         mkdir(dataSetFilePath);
     end
     dataSetFilePath = fullfile(dataSetFilePath,dataSetFileName + ".mat");
     save(dataSetFilePath, 'queryShapelet', 'queryThreshold', 'thresholdConfidence');
end

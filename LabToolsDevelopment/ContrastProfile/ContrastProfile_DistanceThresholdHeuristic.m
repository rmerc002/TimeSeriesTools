function ContrastProfile_DistanceThresholdHeuristic(trainPos, trainNeg, subLength, validateTS, validateLabels, testTS, testLabels, binSize, queryOffset)
%     close all;
    rng(1);
    
    colorGreen = [0.1, 0.8, 0.3];
    colorRed = [0.8, 0.3, 0.1];
    colorBlue = [0.1, 0.3, 0.8];
    colorPurple = [0.5, 0.1, 0.5];

    exclusionLength = floor(subLength*1.1);
    
    
    validateLabels = [validateLabels(queryOffset+1:end); zeros(queryOffset,1)];
    testLabels = [testLabels(queryOffset+1:end); zeros(queryOffset,1)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                            Train/Test
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Train Plato
    [plato, plato_twin, CP] = ContrastProfile(trainPos, trainNeg, subLength, false);

    figure;
    title("Plato");
    plot(plato);
    xlim([0, length(plato)]);
    xticks([0, length(plato)]);
    yticks([]);
    %%% Train Distance Threshold
    distanceProfile = MASS_V2(validateTS, plato);

    [sortedIndices, sortedDistances] = allLowestDistanceIndices(distanceProfile, subLength, exclusionLength);

    classGroundTruth = validateLabels(sortedIndices);
    classResults = -1*ones(size(classGroundTruth));
    distanceThreshold = inf;

    for sampleIndex = sortedIndices
        sample = validateTS(sampleIndex:sampleIndex + subLength - 1);
        nonNanIndices = ~isnan(sample);
        sample = sample(nonNanIndices);
        sample = sample + 1e-5*(rand(length(sample),1)*2-1);

        if length(plato) > length(sample)
           classGroundTruth(sampleIndex) = -1;
        end
    end
    %%% Results based on the number of ground truth class positive
    numClassPos = sum(classGroundTruth==1);
    classResults(1:numClassPos) = 1;
    classResults(numClassPos+1:end) = 0;

    %%% Using FScore as heuristic
    growingPrecision = [];
    growingRecall = [];
    classResultSums = [0,0]; %FP and TP counts
    for sortedOrderIndex = 1:length(sortedIndices)
       orderIndex = sortedIndices(sortedOrderIndex);
       classGT = classGroundTruth(sortedOrderIndex);

       result = (classGT == 1);
       classResultSums(result+1) = classResultSums(result+1) + 1;

       tempPrecision = classResultSums(2) / sum(classResultSums);
       growingPrecision(end+1) = tempPrecision;

       tempRecall = classResultSums(2) / sum(classGroundTruth == 1);
       growingRecall(end+1) = tempRecall;
    end
    beta2 = 0.25.^2;
    growingFScore = (1+beta2) * growingPrecision .* growingRecall./((beta2 .* growingPrecision) + growingRecall);
    [maxFScore,sortedMaxFScoreIndex] = max(growingFScore);
    maxFScoreIndex = sortedIndices(sortedMaxFScoreIndex);
    distanceThreshold = sortedDistances(sortedMaxFScoreIndex);

    fprintf("Distance Threshold: %.4f, max F-Score Index: %d\n",distanceThreshold, maxFScoreIndex);

    figure;
    hold on;
    plot(sortedDistances, growingPrecision,'r');
    plot(sortedDistances, growingRecall,'b');
    plot(sortedDistances, growingFScore,'Color',[0.8, 0.3, 0.8]);
    hold off;
    title("Precision, Recall, FScore");

%     plotName = "Validate";
%     defaultRateComparisonPlot(classGroundTruth, sortedDistances, binSize, plotName);
%     
%     plotSortedBinnedSamples(sortedIndices, validateTS, classGroundTruth, subLength, binSize, plotName);
%     
%     sortedCandidateSubsequences = zeros(length(sortedIndices), subLength);
%     for sortedOrderIndex = 1:length(sortedIndices)
%         baseIndex = sortedIndices(sortedOrderIndex);
%         sortedCandidateSubsequences(sortedOrderIndex,:) = zscore(validateTS(baseIndex: baseIndex + subLength -1 ));
%     end
%     mapcaplot(sortedCandidateSubsequences);
%     sortedCandidateSubsequencesValidate = sortedCandidateSubsequences;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Test distance threshold on test data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    distanceProfile = MASS_V2(testTS, plato);
    [sortedIndices, sortedDistances] = allLowestDistanceIndices(distanceProfile, subLength, exclusionLength);
%     plotDistanceProfileCandidates(testTS, distanceProfile, sortedIndices, testLabels);
     
    classGroundTruth = testLabels(sortedIndices);
    classResults = sortedDistances <= distanceThreshold;

    truePositiveClass1 = sum(classGroundTruth(classResults) == 1);
    falsePositiveClass1 = sum(classGroundTruth(classResults) == 0);
    numClass1 = sum(classGroundTruth == 1);
    precision = truePositiveClass1/(truePositiveClass1 + falsePositiveClass1);
    recall = truePositiveClass1/numClass1;
    fScore = 1.25*precision*recall/(0.25*precision + recall);
    
    %%% Ideal FScore
    growingPrecision = [];
    growingRecall = [];
    classResultSums = [0,0]; %FP and TP counts
    for sortedOrderIndex = 1:length(sortedIndices)
       orderIndex = sortedIndices(sortedOrderIndex);
       classGT = classGroundTruth(sortedOrderIndex);

       result = (classGT == 1);
       classResultSums(result+1) = classResultSums(result+1) + 1;

       tempPrecision = classResultSums(2) / sum(classResultSums);
       growingPrecision(end+1) = tempPrecision;

       tempRecall = classResultSums(2) / sum(classGroundTruth == 1);
       growingRecall(end+1) = tempRecall;
    end
    beta2 = 1.^2;
    growingFScore = (1+beta2) * growingPrecision .* growingRecall./((beta2 .* growingPrecision) + growingRecall);
    [maxFScore,sortedMaxFScoreIndex] = max(growingFScore);
    maxFScoreIndex = sortedIndices(sortedMaxFScoreIndex);
    oracleDistanceThreshold = sortedDistances(sortedMaxFScoreIndex);
    
    
    %%% Random guess
    defaultPrecision = sum(classGroundTruth == 1)/length(classGroundTruth);
    
    fprintf("Test Results:\n");
    fprintf("\t%.2f: Default Precision\n", defaultPrecision);
    fprintf("\t-------------------------------\n");
    fprintf("\t%.2f: Precision\n", precision);
    fprintf("\t%.2f: Recall\n", recall);
    fprintf("\t%.2f: FScore\n", fScore);
    fprintf("\t%.4f: distanceThreshold (heuristic)\n", distanceThreshold);
    fprintf("\t-------------------------------\n");
    fprintf("\t%.2f: FScore (oracle)\n", maxFScore);
    fprintf("\t%.4f: distanceThreshold (oracle)\n", oracleDistanceThreshold);

%     groupSize = 10;
    plotName = "Test";
    defaultRateComparisonPlot(classGroundTruth, sortedDistances, binSize, plotName);
    numExemplars = 10;
    plotSortedBinnedSamples(sortedIndices, testTS, classGroundTruth, subLength, binSize, numExemplars, plotName);
    
    sortedCandidateSubsequences = zeros(length(sortedIndices), subLength);
    for sortedOrderIndex = 1:13000%length(sortedIndices)
        baseIndex = sortedIndices(sortedOrderIndex);
        sortedCandidateSubsequences(sortedOrderIndex,:) = zscore(testTS(baseIndex: baseIndex + subLength -1 ));
    end
%     mapcaplot(sortedCandidateSubsequences);
    sortedCandidateSubsequencesTest = sortedCandidateSubsequences;
    
end
  

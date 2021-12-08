function [platoBinPrecision, randomBinPrecision] = plotRandomSubsequenceComparison(testTS,trainTS, platos, labels, m, binSize, KNN, numBins, plotName, outputPath)
    platoBinPrecision = getPlatoBinPrecision(platos, testTS, labels, binSize);
    
    numIterations = 5;
    numEffectiveBins = length(platoBinPrecision);
    randomBinPrecision = zeros(1, numEffectiveBins);
    
    for ii = 1:numIterations
        randomSubsequences = getRandomSubsequences(trainTS, m, KNN); 
        binPrecision = getPlatoBinPrecision(randomSubsequences, testTS, labels, binSize);
        numEffectiveBins = min(length(randomBinPrecision), length(binPrecision));
        randomBinPrecision(1:numEffectiveBins) = randomBinPrecision(1:numEffectiveBins) + binPrecision(1:numEffectiveBins);
    end
    randomBinPrecision = randomBinPrecision/numIterations;
    
    plotBinImprovement(platoBinPrecision(1:numBins), randomBinPrecision(1:numBins), binSize, plotName,outputPath);
end




function binPrecision = getBinPrecision(DP, labels, m, binSize)
    maxCandidates = length(DP)/m;
    NN = NearestNeighborSelection(DP, m, maxCandidates);

    sortedLabels = labels(NN);

    numCandidates = length(NN);
        
    numSamples = length(NN);
    
    numBins = ceil(numSamples/binSize);
    binPrecision = zeros(1,numBins);
    for binIndex = 1:numBins
       startIndex = (binIndex-1)*binSize + 1;
       endIndex = min(numSamples, startIndex + binSize-1);
       
       binPrecision(binIndex) = nanmean(sortedLabels(startIndex:endIndex));
    end
end
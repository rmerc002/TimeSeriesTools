function alignmentDist = getAlignmentDist(subsequenceWithContext, query)
    %%%assume equal context lengths padding a learned shape with the same
    %%%length as the query
    m = length(query);
    contextLength = floor((length(subsequenceWithContext)-m)/2);

    startIndex = contextLength + 1;
    endIndex = startIndex + m - 1;
    subsequence = subsequenceWithContext(startIndex:endIndex);
    centralDist = norm(zscore(subsequence) - zscore(query));
    DP = real(MASS_V2(subsequenceWithContext, query));

    %%% ramp
%     tempWeights = linspace(0,1,m+1);
    %%% sigmoid weights
    tempWeights = linspace(-5,5,contextLength+1);
    tempWeights = 1./(1+exp(-tempWeights));

    tempWeights = tempWeights(1:end-1);
    weights = [tempWeights, 1, tempWeights(end:-1:1)]';
    DPGain = DP - centralDist;
    weightedDPGain = weights.*DPGain;

    alignmentDist = min(weightedDPGain) + centralDist;
end
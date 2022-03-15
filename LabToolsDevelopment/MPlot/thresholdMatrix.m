function simmatNorm = thresholdMatrix(simmat, distPercentile, minOrMax)    
    %%% Suppress distance above some threshold
    simmatThresh = simmat;
    simmatFlat = reshape(simmat,1,[]);
    if lower(minOrMax) == "max"
        distanceThreshold = prctile(simmatFlat, 100-distPercentile);
        simmatThresh(simmatThresh < distanceThreshold) = distanceThreshold;
    else
        distanceThreshold = prctile(simmatFlat, distPercentile);
        simmatThresh(simmatThresh > distanceThreshold) = distanceThreshold;
    end
    
    %%% Maximize the colormap
    simmatNorm = simmatThresh - min(simmatThresh,[],'all');
    simmatNorm = simmatNorm./max(simmatNorm,[],'all');
end
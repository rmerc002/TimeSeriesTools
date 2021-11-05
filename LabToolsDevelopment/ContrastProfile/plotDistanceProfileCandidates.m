function plotDistanceProfileCandidates(timeseries, labels, query)
    grayColor = [0.7, 0.7, 0.7];
    colorGreen = [0.1, 0.8, 0.3];
    colorRed = [0.8, 0.3, 0.1];
    
    if size(query,1) == 1
        query = query';
    end
    if size(timeseries,1) == 1
        timeseries = timeseries';
    end
    
    m = length(query);
    exclusionLength = m;
    
    DP = MASS_V2(timeseries, query);
    
    timeseriesNoNan = timeseries;
    timeseriesNoNan(isnan(timeseries)) = 0;
    DP = MASS_V2(timeseriesNoNan, query);

    maxVal = sqrt(4*m);
    for nanIndex = 1:length(timeseriesNoNan)-m+1
        if sum(isnan(timeseries(nanIndex)))
            startIndex = max(1, nanIndex-m+1);
            DP(startIndex:nanIndex) = maxVal;
        end
    end
    
    [candidateIndices, correspondingDistances] = allLowestDistanceIndices(DP, m, exclusionLength);
    
    figure;
    subplot(2,1,1);
    hold on;
    plot(timeseries);
    xlim([0, length(timeseries)]);
    
    for orderIndex = 1:length(candidateIndices)
        candidateIndex = candidateIndices(orderIndex);
        colorLabel = colorGreen;
        if labels(candidateIndex) == 0
            colorLabel = colorRed;
        end
        candidate = timeseries(candidateIndex:candidateIndex+m-1);
        plot(candidateIndex:candidateIndex+m-1, candidate,'Color',colorLabel,'LineWidth',3)
    end
    hold off;
    
    subplot(2,1,2);
    hold on;
    plot(DP);
    hold off;
    xlim([0, length(timeseries)]);
end
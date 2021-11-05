function plotDefaultRateComparison(DP, labels, m, binSize, plotName)
    grayColor = [0.3, 0.3, 0.3];
    colorGreen = [0.1, 0.8, 0.3];
    colorRed = [0.8, 0.3, 0.1];
    boxAlpha = 0.3;
    
    maxCandidates = length(DP)/m;
    NN = NearestNeighborSelection(DP, m, maxCandidates);

    sortedClassGroundTruth = labels(NN);
    sortedDistances = DP(NN);
    
    numSamples = length(sortedClassGroundTruth);
    defaultRate = sum(sortedClassGroundTruth)/numSamples;
    
    figure('Name',plotName);
    hold on;
    titleFormat = sprintf("Binned (%d) default rate comparison. Sorted samples by Euclidean distance.", binSize);
    title(titleFormat);
    xlim([0,length(sortedDistances)]);
    ylim([-0.1, 1.1]);
    
    xStart = 0;
    xEnd = length(sortedDistances);
    
    plot([xStart, xEnd], [defaultRate, defaultRate],'--', 'Color', grayColor);
    plot([xStart, xStart], [0, 1], 'Color', grayColor, 'LineWidth', 3);
    xlabel("Sorted test sample indices (by ED)");
    ylabel("Precision per bin");
    
    for groupIndex = 1:binSize:numSamples
       startIndex = groupIndex;
       endIndex = min(numSamples, groupIndex + binSize-1);
%        xStart = sortedDistances(startIndex);
%        xEnd = sortedDistances(endIndex);
       xStart = startIndex;
       xEnd = endIndex;
       xLength = xEnd - xStart;
       
       groupPrecision = sum(sortedClassGroundTruth(startIndex:endIndex))/binSize;
       groupColor = colorGreen;
       if groupPrecision < defaultRate
           groupColor = colorRed;
       end
       
       yStart = defaultRate;
       yEnd = groupPrecision;
       yLength = yEnd - yStart;
       if yLength == 0
           continue;
       elseif yLength < 0
           yStart = yEnd;
           yLength = abs(yLength);
       end
       rectangle('Position',[xStart, yStart, xLength, yLength],'FaceColor', [groupColor,boxAlpha], 'EdgeColor',[groupColor,0.8]);
    end
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  X-Axis is Distance %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Name',plotName);
    hold on;
    titleFormat = sprintf("Binned (%d) default rate comparison. Sorted samples by Euclidean distance.", binSize);
    title(titleFormat);
    xlim([0,sortedDistances(end)]);
    ylim([-0.1, 1.1]);
    
    xStart = sortedDistances(1);
    xEnd = sortedDistances(end);
    
    plot([xStart, xEnd], [defaultRate, defaultRate],'--', 'Color', grayColor);
    plot([xStart, xStart], [0, 1], 'Color', grayColor, 'LineWidth', 3);
    xlabel("Euclidean Distance of sorted test samples (by ED)");
    ylabel("Precision per bin");
    
    for groupIndex = 1:binSize:numSamples
       startIndex = groupIndex;
       endIndex = min(numSamples, groupIndex + binSize-1);
       xStart = sortedDistances(startIndex);
       xEnd = sortedDistances(endIndex);
%        xStart = startIndex;
%        xEnd = endIndex;
       xLength = xEnd - xStart;
       
       groupPrecision = sum(sortedClassGroundTruth(startIndex:endIndex))/binSize;
       groupColor = colorGreen;
       if groupPrecision < defaultRate
           groupColor = colorRed;
       end
       
       yStart = defaultRate;
       yEnd = groupPrecision;
       yLength = yEnd - yStart;
       if yLength == 0
           continue;
       elseif yLength < 0
           yStart = yEnd;
           yLength = abs(yLength);
       end
       rectangle('Position',[xStart, yStart, xLength, yLength],'FaceColor', [groupColor,boxAlpha], 'EdgeColor',[groupColor,0.8]);
    end
end
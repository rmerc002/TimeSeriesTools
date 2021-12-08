function plotDefaultRateComparison(DP, labels, m, binSize, plotName, numBinsPlot, outputPath)
    grayColor = [0.3, 0.3, 0.3];
    colorGreen = [0.1, 0.8, 0.3];
    colorRed = [0.8, 0.3, 0.1];
    boxAlpha = 0.3;
    
    maxCandidates = length(DP)/m;
    NN = NearestNeighborSelection(DP, m, maxCandidates);

    sortedLabels = labels(NN);
    sortedDistances = DP(NN);
    
    numSamples = length(sortedLabels);
    defaultRate = nanmean(sortedLabels);
    
    fig = figure('Name',plotName);
    hold on;
    titleFormat = sprintf("Binned (%d) default rate comparison.\nSorted samples by Euclidean distance.", binSize);
    title(titleFormat);
    xlim([0,binSize*numBinsPlot]);
    ylim([-0.1, 1.1]);
    
    xStart = 0;
    xEnd = length(sortedDistances);
    
    plot([xStart, xEnd], [defaultRate, defaultRate],'--', 'Color', grayColor);
%     plot([xStart, xStart], [0, 1], 'Color', grayColor, 'LineWidth', 3);
    xlabel("Sorted test sample indices (by ED)");
    ylabel("Precision per bin");
    set(gca, 'TickDir','out');
    box off;
    
    binCount = 0;
    for binIndex = 1:binSize:numSamples
        binCount = binCount + 1;
        if binCount > numBinsPlot
           break; 
        end
       startIndex = binIndex;
       endIndex = min(numSamples, binIndex + binSize-1);
%        xStart = sortedDistances(startIndex);
%        xEnd = sortedDistances(endIndex);
       xStart = startIndex;
       xEnd = endIndex;
       xLength = xEnd - xStart;
       
       groupPrecision = nanmean(sortedLabels(startIndex:endIndex));
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
    
    if ~isempty(outputPath)
        fileName =  "improvedOverDefault_xSample";
        filePath = fullfile(outputPath, fileName + ".fig");
        savefig(fig, filePath);
        
        filePath = fullfile(outputPath, fileName + ".png");
        saveas(gcf, filePath);

        filePath = fullfile(outputPath, fileName + ".emf");
        print(filePath,'-dmeta'); 
    end
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  X-Axis is Distance %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure('Name',plotName);
    hold on;
    titleFormat = sprintf("Binned (%d) default rate comparison.\nSorted samples by Euclidean distance.", binSize);
    title(titleFormat);
    xlim([0,sortedDistances(end)]);
    ylim([-0.1, 1.1]);
    
    xStart = sortedDistances(1);
    xEnd = sortedDistances(end);
    
    plot([xStart, xEnd], [defaultRate, defaultRate],'--', 'Color', grayColor);
    plot([xStart, xStart], [0, 1], 'Color', grayColor, 'LineWidth', 3);
    xlabel("Euclidean Distance of sorted test samples (by ED)");
    ylabel("Precision per bin");
    set(gca, 'TickDir','out');
    box off;
    
    binCount = 0;
    for binIndex = 1:binSize:numSamples
        
        binCount = binCount + 1;
        if binCount > numBinsPlot
           break; 
        end
       
       startIndex = binIndex;
       endIndex = min(numSamples, binIndex + binSize-1);
       xStart = sortedDistances(startIndex);
       xEnd = sortedDistances(endIndex);
%        xStart = startIndex;
%        xEnd = endIndex;
       xLength = xEnd - xStart;
       
       groupPrecision = nanmean(sortedLabels(startIndex:endIndex));
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
    
    if ~isempty(outputPath)
        fileName =  "improvedOverDefault_xED";
        filePath = fullfile(outputPath, fileName + ".fig");
        savefig(fig, filePath);
        
        filePath = fullfile(outputPath, fileName + ".png");
        saveas(gcf, filePath);

        filePath = fullfile(outputPath, fileName + ".emf");
        print(filePath,'-dmeta'); 
    end
end
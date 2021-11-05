function plotQueryNN(ts, query, KNN)
    m = length(query);
    %%%Generating unique colors for classes when there can be many classes
    %%% I will assume no more than 1000 classes
    numColors = 1000;
    colors = lines(numColors);
    
    %plot with context
    contextLength = m;
    
    
    DP = getPlatoDistanceProfile(ts, query);
    
    NN = NearestNeighborSelection(DP, m, KNN);

    sampleLength = contextLength + m + contextLength;
    subsequences = nan(1+KNN,sampleLength);
    
    %%%Store the query
    startIndex = contextLength+1;
    endIndex = startIndex + m-1;
    subsequences(1,startIndex:endIndex) = query;
    
    %%%Store the nearest neighbors
    ts = [nan(1,contextLength),ts,nan(1,contextLength)]; %easy way to handles edge cases
    for ssIndex = 1:KNN
        nnIndex = NN(ssIndex) + contextLength;
        startIndex = nnIndex - contextLength;
        endIndex = startIndex + sampleLength - 1;
        subsequences(ssIndex+1,:) = ts(startIndex:endIndex);
    end
    
    inset = 0.9;
    figure;
    hold on;
    xVals = 1-contextLength:m+contextLength;
    %Plot individually
    for plotIndex = 1:size(subsequences,1)
        ssIndex = plotIndex;
        %%%Min/Max normalize just NN, not context
        startIndex = contextLength+1;
        endIndex = startIndex + m-1;
        tempTS = subsequences(plotIndex, startIndex:endIndex);
        tempMin = nanmin(tempTS);
        tempMax = nanmax(tempTS);
        tempRange = tempMax - tempMin;
        %%%Now apply normalization to context as well
        tempTS = subsequences(ssIndex, :);
        tempPlot = -plotIndex + inset*(tempTS - tempMin)/tempRange;
        
        plot(xVals, tempPlot,'Color',[0.5,0.5,0.5]);
        color = colors(ssIndex,:);
        startIndex = contextLength+1;
        endIndex = startIndex + m-1;
        plot(xVals(startIndex:endIndex),tempPlot(startIndex:endIndex),'Color',color, 'LineWidth',2);
    end
    plot([xVals(1),xVals(end)], [-1.1, -1.1], '--', 'Color', [0.5, 0.5, 0.5]);
    plot([xVals(1),xVals(end)], [-plotIndex-0.1, -plotIndex-0.1], '--', 'Color', [0.5, 0.5, 0.5]);
    
    plotIndex = plotIndex + 2;
    %Now all overlapping
    %%%z-normalize just NN, not context
    for ssIndex = 1:size(subsequences,1)
        startIndex = contextLength+1;
        endIndex = startIndex + m-1;
        tempTS = subsequences(ssIndex, startIndex:endIndex);
        tempMean = nanmean(tempTS);
        tempSTD = nanstd(tempTS);
        %%%Now apply normalization to context as well
        tempTS = subsequences(ssIndex, :);
        tempPlot = -plotIndex + inset*(tempTS - tempMean)/(3*tempSTD);
        
        plot(xVals, tempPlot,'Color',[0.5,0.5,0.5]);
        color = colors(ssIndex,:);
        startIndex = contextLength+1;
        endIndex = startIndex + m-1;
        if ssIndex == 1
            lw = 1;
        else
            lw = 0.5;
        end
        plot(xVals(startIndex:endIndex),tempPlot(startIndex:endIndex),'Color',color,'LineWidth',lw);
    end
    hold off;
    
    titleFormat = sprintf("Plato with %d NN", KNN);
    title(titleFormat);
    xlabel("Time Index");
    set(gca,'xtick',[1-contextLength,1,m,m+contextLength],'ytick',[], 'TickDir','out');
    box off;
end
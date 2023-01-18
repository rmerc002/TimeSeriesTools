function [subsequences, ssIndices] = plotQueryNN(ts, query, KNN, contextLength)
    
    query = reshape(query,length(query),1);
    ts = reshape(ts,length(ts),1);
    m = length(query);
    epsilon = 1e-6;
    %%%Generating unique colors for classes when there can be many classes
    %%% I will assume no more than 1000 classes
    numColors = 1000;
    colors = lines(numColors);
    
    %plot with context
    if nargin < 4
        contextLength = m;
    end
    
    
    DP = getPlatoDistanceProfile(ts, query);
    
    [NN, NNDist] = NearestNeighborSelection(DP, m, KNN+1);
    
    %%%Ignore the first if match is exact, it is likely finding query from
    %%%the time series it came from
    if NNDist(1) < epsilon
       NN(1) = [];
       NNDist(1) = [];
    end
    
    sampleLength = contextLength + m + contextLength;
    subsequences = nan(1+KNN,sampleLength);
    
    %%%Store the query
    startIndex = contextLength+1;
    endIndex = startIndex + m-1;
    subsequences(1,startIndex:endIndex) = query;
    ssIndices = zeros(KNN,1);
    %%%Store the nearest neighbors
    ts = [nan(contextLength,1);ts;nan(contextLength,1)]; %easy way to handles edge cases
    for ssIndex = 1:KNN
        nnIndex = NN(ssIndex) + contextLength;
        startIndex = nnIndex - contextLength;
        endIndex = startIndex + sampleLength - 1;

        ssIndices(ssIndex) = startIndex;

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
        plot(xVals(startIndex:endIndex),tempPlot(startIndex:endIndex),'Color',color, 'LineWidth',1);
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
            lw = 2;
        else
            lw = 0.5;
        end
        plot(xVals(startIndex:endIndex),tempPlot(startIndex:endIndex),'Color',color,'LineWidth',lw);
    end
    hold off;
    
    titleFormat = sprintf("Query with %d NN", KNN);
    title(titleFormat);
    xlabel("Time Index");
    set(gca,'xtick',[1-contextLength,1,m,m+contextLength],'ytick',[], 'TickDir','out');
    box off;

    %%% Distance Profile Locations

    figure;
    hold on;
    plot(DP);
    NNDist = DP(NN);
    scatter(NN,NNDist,'filled');
    hold off;
end
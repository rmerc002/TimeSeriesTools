function plotSortedBinnedSamples(sortedIndices, timeSeries, classLabels, subLength, binSize, numExemplars, plotName)
    grayColor = [0.7, 0.7, 0.7];
    colorGreen = [0.1, 0.8, 0.3];
    colorRed = [0.8, 0.3, 0.1];
    boxAlpha = 0.3;
    
    if isempty(classLabels)
       classLabels = zeros(length(timeSeries),1); 
    end
    
    lengthContext = floor(subLength);
    lengthInstance = lengthContext + subLength + lengthContext;
    
    inset = 0.9;
    
    figure('Name',plotName);
    xticks([]);
    yticks([]);
    
    ylim([-numExemplars - 1, 2]);
    
    hold on;
    for indexSortedOrder = 1:length(sortedIndices)
        indexSample = sortedIndices(indexSortedOrder);
        
        
        xOffset = 1.1*lengthInstance*floor((indexSortedOrder-1)/binSize);
        yOffset = -mod(indexSortedOrder-1, binSize);
        if -yOffset > numExemplars
            continue;
        end
            
        
        %%% Annotations
        if yOffset == 0
            binNumber = floor((indexSortedOrder-1)/binSize) + 1;
            binAnnotation = sprintf("Bin %d",binNumber);
            text(xOffset, 2, binAnnotation);
        end
        
        binAnnotation = sprintf("Idx:%d",indexSample);
        text(xOffset, yOffset+0.9, binAnnotation,'FontSize', 8);

        startIndexContext = indexSample - lengthContext;
        endIndexContext = indexSample + subLength + lengthContext - 1;
        
        nanPrefix = NaN(max(0, -1*startIndexContext + 1),1);
        nanSuffix = NaN(max(0, endIndexContext - length(timeSeries)),1);
        
        startIndexContext = max(1, startIndexContext);
        endIndexContext = min(length(timeSeries), endIndexContext);
        
        startIndexCandidate = indexSample;
        endIndexCandidate = indexSample + subLength - 1;
        
        tempTSCandidate = timeSeries(startIndexCandidate:endIndexCandidate);
        
        tempTSContext = timeSeries(startIndexContext:endIndexContext);
        tempMin = min(tempTSContext);
        tempMax = max(tempTSContext);
        tempRange = max(1e-5, tempMax-tempMin);
        
        tempClass = classLabels(indexSortedOrder);
        classColor = colorGreen;
        if tempClass == 0
           classColor = colorRed; 
        end

        

        xSeriesWithContext = (xOffset:xOffset + lengthContext + subLength + lengthContext - 1)';
        ySeriesWithContext = yOffset + inset * (tempTSContext - tempMin) / tempRange;
        ySeriesWithContext = [nanPrefix; ySeriesWithContext(:); nanSuffix];
        
        xSeriesOnlyCandidate = (xOffset + lengthContext:xOffset + lengthContext + subLength -1)';
        ySeriesOnlyCandidate = yOffset + inset * (tempTSCandidate - tempMin) / tempRange;
        
        

        plot(xSeriesWithContext, ySeriesWithContext, 'Color', grayColor);
        plot(xSeriesOnlyCandidate, ySeriesOnlyCandidate, 'Color', classColor, 'LineWidth', 3);
    
    end
    hold off;
    
end
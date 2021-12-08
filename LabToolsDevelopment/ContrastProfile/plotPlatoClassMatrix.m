function plotPlatoClassMatrix(classPlatos)
    numColors = 1000;
    colors = lines(numColors);  

    numClasses = size(classPlatos,1);
    figure;
    hold on;
    xOffset = 0;
    for negClassIndex = 1:numClasses
        maxLength = 0;
        for posClassIndex = 1:numClasses
            if posClassIndex == negClassIndex
                continue;
            end
            tempTS = classPlatos{posClassIndex, negClassIndex};
            tempMin = nanmin(tempTS);
            tempMax = nanmax(tempTS);
            tempRange = tempMax - tempMin;
    
            maxLength = max(maxLength, length(tempTS));
    
            tempPlot = -posClassIndex + 0.9*(tempTS - tempMin)/tempRange;
            xAxis = 1:length(tempTS);
            xAxis = xAxis + xOffset;
            plot(xAxis, tempPlot, 'Color', colors(posClassIndex,:));
        end
        xOffset = xOffset + maxLength + 20;
    end
    hold off;
end
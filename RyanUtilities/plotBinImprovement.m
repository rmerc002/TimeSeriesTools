function plotBinImprovement(platoBinPrecision, randomBinPrecision, binSize, plotName, outputPath)
    grayColor = [0.3, 0.3, 0.3];
    colorGreen = [0.1, 0.8, 0.3];
    colorRed = [0.8, 0.3, 0.1];
    boxAlpha = 0.3;
    
    fig = figure('Name',plotName);
    hold on;
    titleFormat = sprintf("Binned (%d) random subsequence comparison.\nSorted samples by Euclidean distance.", binSize);
    title(titleFormat);
    xlim([0,binSize*length(platoBinPrecision)]);
    ylim([-1.1, 1.1]);
   
    xlabel("Sorted test sample indices (by ED)");
    ylabel("Precision gain per bin");
    
    numBins = length(platoBinPrecision);
    
    xStart = 0;
    xEnd = binSize*numBins;
    plot([xStart, xEnd], [0,0],'--', 'Color', grayColor);
%     plot([xStart, xStart], [0, 1], 'Color', grayColor, 'LineWidth', 3);
    
    
    for binIndex = 1:numBins
       startIndex = (binIndex-1)*binSize + 1;
       endIndex = startIndex + binSize - 1;

       xStart = startIndex;
       xEnd = endIndex;
       xLength = xEnd - xStart;
       
       binImprovement = platoBinPrecision(binIndex) - randomBinPrecision(binIndex);
       groupColor = colorGreen;
       if binImprovement <= 0
           groupColor = colorRed;
       end
       
       yStart = 0;
       yEnd = binImprovement;
       yLength = yEnd - yStart;
       if yLength == 0
           continue;
       elseif yLength < 0
           yStart = yEnd;
           yLength = abs(yLength);
       end
       rectangle('Position',[xStart, yStart, xLength, yLength],'FaceColor', [groupColor,boxAlpha], 'EdgeColor',[groupColor,0.8]);
    end

    set(gca, 'TickDir','out');
    box off;
    
    if ~isempty(outputPath)
        fileName = "improvedOverRandSubSeq_xSample";
        filePath = fullfile(outputPath, fileName + ".fig");
        savefig(fig, filePath);
        
        filePath = fullfile(outputPath, fileName + ".png");
        saveas(gcf, filePath);

        filePath = fullfile(outputPath, fileName + ".emf");
        print(filePath,'-dmeta'); 
    end
end
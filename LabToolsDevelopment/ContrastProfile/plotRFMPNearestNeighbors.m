function plotRFMPNearestNeighbors(targetIndex, RFMP_AA_Indices)
    maxFreq = size(RFMP_AA_Indices, 1); %Number of nearest neighbors
    m = 200; %subsequence length
    figure;
    hold on;
    inset = 0.9;
    for plotIndex = 1:maxFreq
        nnIndex = RFMP_AA_Indices(plotIndex,targetIndex);
        tempTS = pos(nnIndex:nnIndex+m-1);
        tempMin = nanmin(tempTS);
        tempMax = nanmax(tempTS);
        tempRange = tempMax - tempMin;
        tempPlot = -plotIndex + inset*(tempTS - tempMin)/tempRange;

        plot(tempPlot);
    end
    hold off;
    formattedTitle = sprintf("Plato's Nearest Neighbors");
    title(formattedTitle);
    set(gca,'xtick',[1,m],'ytick',[], 'TickDir','out');
    xlim([0,m]);
    ylim([-plotIndex-0.5,1]);
    box off;
end
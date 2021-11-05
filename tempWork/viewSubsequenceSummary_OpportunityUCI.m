mask = S2Drill{1:10:end,250} == 407521;

data = S2Drill{1:10:end,64};
figure;
plotCount = 0;
startIndex = 1;
for i = 2:length(data)
    if mask(i) > mask(i-1)
        startIndex = i;
    elseif mask(i) < mask(i-1)
        plotOffset = floor(plotCount/2)+1+(1-mod(plotCount,2))*20;
        plot(zscore(data(startIndex:i))/2+plotOffset);
        hold on;
        plotCount = plotCount + 1;
    end
end
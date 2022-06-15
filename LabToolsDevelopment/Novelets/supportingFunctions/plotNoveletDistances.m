redColor = [0.73,0.05,0];
%             greenColor = [0,0.73,0.41]; 
blueColor = [0,0.29,0.73];
grayColor = [0.75,0.75,0.75];
lightGrayColor = [0.9, 0.9, 0.9];
%             lightBlueColor = [0.01, 0.83,0.99];
noveletColor = [143/255, 226/255, 227/255];
noveletLabelColor = [56/255, 226/255, 235/255];
noveletNNColor = [227/255, 143/255, 219/255];
noveletNNLabelColor = [235/255, 56/255, 217/255];

splitPos = [0.2 0.8];
emptyPosition = [0 splitPos(2) splitPos(1) 1-splitPos(2)];
posPosition = [splitPos(1) splitPos(2) 1-splitPos(1) 1-splitPos(2)];
noveletPosition = [0 0 splitPos(1) splitPos(2)];
distProfPosition = [splitPos(1) 0 1-splitPos(1) splitPos(2)];

figure;
ax1 = subplot('Position', posPosition);
plot(np.positiveTS);
set(gca,'xtick',[],'ytick',[], 'TickDir','out');
xlim([1,length(np.positiveTS)]);
box off;



nargin = 1;
if nargin <= 1
    %%% Options to sort by. First is at the top.
    %%% "descendingScore": plot novelets with highest novety
    %%%     score at the top
    %%% "ascendingChronology": plot novelets in order
    %%%     discovered
    sortMode = "ascendingChronology";
end

if sortMode == "descendingScore"
    [~,sortedOrder] = sort(np.noveletScores,"descend");
    sortedIndices = np.noveletIndices(sortedOrder);
    sortedNNIndices = np.noveletNNIndices(sortedOrder);
    sortTitle = "Descending Novelty Score";
else %%% "ascendingChronology"
    sortedOrder = 1:length(np.noveletIndices);
    sortedIndices = np.noveletIndices;
    sortedNNIndices = np.noveletNNIndices;
    sortTitle = "Order Discovered";
end


        
%%%Generating unique colors for classes when there can be many classes
%%% I will assume no more than 1000 results
numColors = 1000;
colors = lines(numColors);

ax2 = subplot('Position', noveletPosition);
hold on;
numNovelets = size(np.novelets,1);
startIndex = np.contextLength+1;
endIndex = startIndex + np.mm - 1;
for ki = 1:numNovelets
    %%%Nearest neighbor plotted under
    startIndexNN = sortedNNIndices(ki);
    endIndexNN = startIndexNN + np.mm - 1;
    tempTS = np.positiveTS(startIndexNN: endIndexNN);
    tempMin = min(tempTS);
    tempMax = max(tempTS);
    tempRange = tempMax-tempMin;
    plot(1:np.mm,-ki + 0.9*(tempTS - tempMin)/tempRange,'Color', 1-0.3*(1-colors(ki,:)));

    %%% prepare the Novelet without context
    ni = sortedOrder(ki);
    tempTS = np.novelets(ni,startIndex:endIndex);
    tempMin = min(tempTS);
    tempMax = max(tempTS);
    tempRange = tempMax-tempMin;

    %%% plot the Novelet on top
    plot(1:np.mm,-ki + 0.9*(tempTS - tempMin)/tempRange,'Color',colors(ki,:));
end
hold off;
% formattedTitle = sprintf("%d Novelets > %.2f, Sorted by %s", numNovelets, np.noveltyThreshold, sortTitle);
% title(formattedTitle);
set(gca,'xtick',[1,np.mm],'ytick',[], 'TickDir','out');
xlim([1,np.mm]);
box off;

ax3 = subplot('Position', distProfPosition);
numNovelets = size(np.novelets,1);
hold on;
for ki = 1:numNovelets
    %%% prepare the Novelet without context
    ni = sortedOrder(ki);
    startIndex = sortedIndices(ki);
    endIndex = startIndex + np.mm - 1;
    novelet = np.positiveTS(startIndex:endIndex);
    tempTS = real(MASS_V2(np.positiveTS, novelet));
    tempMin = min(tempTS);
    tempMax = max(tempTS);
    tempRange = tempMax-tempMin;

    %%% plot the distance threshold used to classify the behavior
    noveltyThreshold = np.noveletDistanceThresholds(ni);
    offsetThreshold = -ki + 0.9*(noveltyThreshold-tempMin)/tempRange;
    plot([0,length(np.positiveTS)],[offsetThreshold, offsetThreshold],'--','Color',[0.3, 0.3, 0.3]);
    %%% plot the Novelet on top
    plot(-ki + 0.9*(tempTS - tempMin)/tempRange,'Color',colors(ki,:));

    [peakIndices, peakValues] = exclusionZonePeaks(-tempTS, np.mm, np.mm, length(tempTS), -noveltyThreshold);
    for pi = 1:length(peakIndices)
        offsetMarker = -ki + 0.9*1;
        scatter(peakIndices, offsetMarker*ones(length(peakIndices),1),'v','filled','MarkerEdgeColor', redColor,'MarkerFaceColor', redColor);
    end

end
hold off;
% formattedTitle = sprintf("%d Novelets > %.2f, Sorted by %s", numNovelets, np.noveltyThreshold, sortTitle);
% title(formattedTitle);
set(gca,'xtick',[],'ytick',[], 'TickDir','out');
xlim([1,length(np.positiveTS)]);
box off;

linkaxes([ax1, ax3],'x');
linkaxes([ax2 ax3],'y');
function  [magicMP, profileIndices] = magicMatrixProfileVisualizer(data, magicMP, profileIndices, subLenSeries)
% close all;

%improves contrast
ps = linspace(0,100,11)';
percentileValues = zeros(1,size(ps,1));

for index=2:size(ps,1)
    percentileValues(index) = prctile(magicMP,ps(index),'all');
end
tempCanvas = magicMP;
numPercentiles = size(percentileValues,2);
for index=2:numPercentiles
    selection = percentileValues(index-1) < magicMP & magicMP <= percentileValues(index);
    tempMMP = magicMP(selection);
    tempMMP = normalize(tempMMP,'range');
    tempCanvas(selection) = tempMMP*(1/(numPercentiles-1)) + (index-2)/(numPercentiles-1);
end
magicMP = tempCanvas;


data = data/max(max(data));
f1 = figure;

ax1 = subplot(10,1,1);

ax2 = subplot(10,1,2);
plot(data,'Color',[.5,.5,.5]);
xlim([1,length(data)]);

%give an initial guess 
ax3 = subplot(10,1,3:10);
threshold = 0.1;%max(max(magicMP))*(1/3);
magicMP_contrast = centerOnThreshold(magicMP, threshold);
plotSurface(magicMP_contrast,subLenSeries);

linkaxes([ax2,ax3],'x');
adjustContrast = true;

while true
    [timeseriesIndex, subLength] = ginput(1);
    
    if checkContrastButtonClicked(timeseriesIndex, subLength, magicMP, subLenSeries)
        adjustContrast = ~adjustContrast;
        continue
    end
    
    timeseriesIndex = uint32(timeseriesIndex);
    [~,subLengthIndex] = max(subLenSeries((subLenSeries-subLength)<= 0));
    subLength = subLenSeries(subLengthIndex);
    
    subplot(10,1,3:10);
    if adjustContrast == true
        threshold = magicMP(subLengthIndex, timeseriesIndex);
        magicMP_contrast = centerOnThreshold(magicMP, threshold);
        plotSurface(magicMP_contrast,subLenSeries);
    end
    
    %%% Plot nearest neighbor subsequence
    startIndex1 = timeseriesIndex;
    endIndex1 = min(length(data), startIndex1 + subLength);
    startIndex2 = profileIndices(subLengthIndex,startIndex1);
    endIndex2 = min(length(data), startIndex2 + subLength);
    plotNearestNeighborSubsequences(data, profileIndices, startIndex1, startIndex2, subLength);
    
    %%% Plot original data with highlighted subsequences
    subplot(10,1,2);
    plot(data,'Color',[.5,.5,.5]);
    xlim([1,length(data)]);
    hold on;
    plot(startIndex1:startIndex1 + subLength,data(startIndex1:startIndex1 + subLength),'g');
    plot(startIndex2:startIndex2 + subLength,data(startIndex2:startIndex2 + subLength),'r');
    hold off;
    
    %give some descriptive data
    subplot(10,1,3:10);
    xlabel(sprintf('index1 = %d, index2 = %d, subLength = %d, dist = %f',startIndex1, startIndex2, subLenSeries(subLengthIndex), magicMP(subLengthIndex,timeseriesIndex)));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   HELPER FUNCTIONS   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSurface(magicMP, subLenSeries)
    magicMP = drawContrastButton(magicMP, subLenSeries);
    
    
    startIndex = 1;%2000;
    endIndex = size(magicMP,2);%4000;
%     endIndex = size(magicMP,2)
    step = 1;
    if size(magicMP,2)>1024
        step = floor((endIndex-startIndex)/1024);
    end
    plotIndices = startIndex:step:endIndex;
    [X,Y] = meshgrid(plotIndices, subLenSeries);
    h = surf(X,Y,magicMP(:,plotIndices));
%     shading('interp');
    view(2);
    set(h,'LineStyle','none');
%     c1 = hot(150);
%     c1 = c1(1:100,:);
%     c2 = flipud(winter(100));
%     customColor = [c1;ones(5,3);c2];
%     colormap(customColor);
%     cm = load('MMPColormap.mat');
%     cm = load('MyColormap9.mat');
%         cm = cell2mat(struct2cell(cm));
    cm = gray();
    colormap(gca, cm);
    
%     colormap(flipud(jet));
    % colormap(jet);
    xlim([startIndex,endIndex]);
    ylim([1,subLenSeries(end)]);
%     colorbar;
    xlabel('Timeseries indices');
    ylabel('SubLength');
    
end

function magicMP = drawContrastButton(magicMP, subLenSeries)
    [x,y] = getButtonBounds(size(magicMP,2), subLenSeries(end));
    series = 1:length(subLenSeries);
    index = min(series(subLenSeries >= y));
    magicMP(index:end, x:end) = 0;
end

function buttonClicked = checkContrastButtonClicked(xInput, yInput, magicMP, subLenSeries)
    buttonClicked = false;
    [x,y] = getButtonBounds(size(magicMP,2), subLenSeries(end));
    if xInput >= x && yInput >= y
        buttonClicked = true;
    end
end

function [x, y] = getButtonBounds(xMax,yMax)
    x = ceil(xMax - xMax/25);
    y = ceil(yMax - yMax/10);
    
end

function [magicMP_contrast] = centerOnThreshold(magicMP, threshold)

    magicMP_contrast = magicMP;% - threshold;
    %matching indices
    indices = magicMP <= threshold;
    magicMP_contrast(indices) = magicMP(indices) - min(min(magicMP(indices)));
    magicMP_contrast(indices) = magicMP_contrast(indices)./max(max(magicMP_contrast(indices)));
    magicMP_contrast(indices) = magicMP_contrast(indices)./2;

    %differing indices
    indices = ~indices;
    magicMP_contrast(indices) = magicMP(indices) - min(min(magicMP(indices)));
    magicMP_contrast(indices) = magicMP_contrast(indices)./max(max(magicMP_contrast(indices)));
    magicMP_contrast(indices) = magicMP_contrast(indices)./2 + 0.5;
    
end

function plotNearestNeighborSubsequences(data, profileIndices, index1, index2, subLength)
    subplot(10,1,1);
    
    startIndex = index1;
    endIndex = startIndex+subLength;

    startIndex2 = index2;
    endIndex2 = min(length(data), startIndex2 + subLength);
    
    plot(zscore(data(startIndex:endIndex)),'g');
    hold on;
    plot(zscore(data(startIndex2:endIndex2)),'r');
    hold off;
    xlim([1,subLength]);
end
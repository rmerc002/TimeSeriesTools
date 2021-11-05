function profileIndexVisualizer(data, magicMP, profileIndices, subLenSeries)
% close all;

%add contrast using percentiles
tempMMP0 = magicMP;
ps = linspace(0,100,10)';
pV = zeros(1,size(ps,1));
for index=2:size(ps,1)
    pV(index) = prctile(magicMP,ps(index),'all');
end
for index=2:size(ps,1)
    selection = pV(index-1) < magicMP & magicMP <= pV(index);
    tempMMP = magicMP(selection);
    tempMMP = normalize(tempMMP,'range');
    tempMMP0(selection) = tempMMP*(1/(size(ps,1)-1)) + (index-2)/(size(ps,1)-1);
end
magicMP = tempMMP0;
max(magicMP,[],'all')
min(magicMP,[],'all')

data = data/max(max(data));
f1 = figure;

ax1 = subplot(10,1,1);

ax2 = subplot(10,1,2);
plot(data,'Color',[.5,.5,.5]);
xlim([1,length(data)]);

%give an initial guess 
ax3 = subplot(10,1,3:10);
selectedIndex = ceil(size(data,2)/2);

piContrast = centerOnThreshold(profileIndices, selectedIndex);
plotSurface(magicMP, piContrast,subLenSeries);

linkaxes([ax2,ax3],'x');
adjustContrast = true;

while true
    [timeseriesIndex, subLength] = ginput(1);
    
    if checkContrastButtonClicked(timeseriesIndex, subLength, profileIndices, subLenSeries)
        adjustContrast = ~adjustContrast;
        continue
    end
    
    timeseriesIndex = uint32(timeseriesIndex);
    [~,subLengthIndex] = max(subLenSeries((subLenSeries-subLength)<= 0));
    subLength = subLenSeries(subLengthIndex);
    
    subplot(10,1,3:10);
    if adjustContrast == true
        piContrast = centerOnThreshold(profileIndices, timeseriesIndex);
        plotSurface(magicMP, piContrast,subLenSeries);
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
%     xlabel(sprintf('index1 = %d, index2 = %d, subLength = %d, dist = %f',startIndex1, startIndex2, subLenSeries(subLengthIndex), profileIndices(subLengthIndex,timeseriesIndex)));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   HELPER FUNCTIONS   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSurface(magicMP, profileIndices, subLenSeries)
    canvas = drawContrastButton(magicMP, subLenSeries);
    indices = profileIndices ~= 0;
    piNorm = profileIndices./size(profileIndices,2);
    canvas(indices) = piNorm(indices)+1;
    
    

    
    startIndex = 1;%2000;
    endIndex = size(canvas,2);%4000;
%     endIndex = size(profileIndices,2)
    step = floor((endIndex-startIndex)/1024)
    plotIndices = startIndex:step:endIndex;
    
    %To keep colorbar consistent
    canvas = canvas(:,plotIndices);
    
    canvas = [linspace(1,2,size(canvas,2)); canvas];
    subLenSeries = [0,subLenSeries];
    [X,Y] = meshgrid(plotIndices, subLenSeries);
    h = surf(X,Y,canvas);
    view(2);
    set(h,'LineStyle','none');
%     c1 = hot(150);
%     c1 = c1(1:100,:);
%     c2 = flipud(winter(100));
%     customColor = [c1;ones(5,3);c2];
%     colormap(customColor);
%     cm = load('MMPColormap.mat');
%     cm = load('MyColormap9.mat');
%     cm = cell2mat(struct2cell(cm));

    cm = [gray(100);hsv(100)];
    colormap(gca, cm);
%     colorbar('Location','south');
    
%     colormap(flipud(jet));
    % colormap(jet);
    xlim([startIndex,endIndex]);
    ylim([1,subLenSeries(end)]);
%     colorbar;
    xlabel('Timeseries indices');
    ylabel('SubLength');
    
end

function profileIndices = drawContrastButton(profileIndices, subLenSeries)
    [x,y] = getButtonBounds(size(profileIndices,2), subLenSeries(end));
    series = 1:length(subLenSeries);
    index = min(series(subLenSeries >= y));
    profileIndices(index:end, x:end) = 0;
end

function buttonClicked = checkContrastButtonClicked(xInput, yInput, profileIndices, subLenSeries)
    buttonClicked = false;
    [x,y] = getButtonBounds(size(profileIndices,2), subLenSeries(end));
    if xInput >= x && yInput >= y
        buttonClicked = true;
    end
end

function [x, y] = getButtonBounds(xMax,yMax)
    x = ceil(xMax - xMax/25);
    y = ceil(yMax - yMax/10);
    
end

function [piContrast] = centerOnThreshold(profileIndices, selectedIndex)
%     piContrast = profileIndices;
%     return;
    
    piMap = zeros(size(profileIndices));
    piMap(:,selectedIndex) = 1;
    
    sequencialIndices = abs(profileIndices(:,1:end-1) - profileIndices(:,2:end))<=5;
    sequencialIndices = [sequencialIndices(:,1:selectedIndex-1), zeros(size(profileIndices,1),1), sequencialIndices(:,selectedIndex:end)];
    
    %select sequential indices, center to right
    cont = ones(size(profileIndices,1),1);
    for i=selectedIndex+1:size(profileIndices,2)
        cont = cont & sequencialIndices(:,i);
        piMap(:,i) = cont;
    end
    %select sequential indices, center to left
    cont = ones(size(profileIndices,1),1);
    for i=selectedIndex-1:-1:1
        cont = cont & sequencialIndices(:,i);
        piMap(:,i) = cont;
    end
    
    piContrast = piMap .* profileIndices;
    
%     piNan = profileIndices;
%     piNan(piNan==0) = nan;
% %     nonzeroIndices = piMap == 1;
%     piMap(nonzeroIndices .* piNan(nonzeroIndices,i)) = 1;
    
%     for i=1:size(profileIndices,2)
%         nonzeroIndices = piMap(:,i) == 1 & profileIndices(:,i) ~= 0;
%         piMap(nonzeroIndices,profileIndices(nonzeroIndices,i)') = 1;
%     end
    
    for i=1:size(profileIndices,1)
        for j=1:size(profileIndices,2)
            if piContrast(i,j) ~= 0
                piMap(i,profileIndices(i,j)) = 1;
            end
        end
    end
   
    piContrast = piMap .* profileIndices;
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
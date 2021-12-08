numColors = 1000;
colors = lines(numColors);  

shapeletData = ChlorineConcentrationTEST{:,:};
labels = shapeletData(:,1);
samples = shapeletData(:,2:end);


[labels, sortedLabelIndices] = sort(labels);
samples = samples(sortedLabelIndices,:);

classes = sort(unique(labels));
numClasses = length(classes);

%%%Plot the samples which are sorted by class
%%%Each class is plotted with a unique color
fig = figure('Name','Shapelet Classes','NumberTitle','off');
fig.Position = [0 0 300 1000];
inset = 0.9;
hold on;
prevClassIndex = -inf;
textLength = 0.3*size(samples,2);
for plotIndex = 1:size(samples,1)
    tempTS = samples(plotIndex,:);
    tempMin = min(tempTS);
    tempMax = max(tempTS);
    tempRange = max(1e-5, tempMax - tempMin);
    tempPlot = -plotIndex + inset*(tempTS - tempMin)/tempRange;

    classPosIndex = find(labels(plotIndex) == classes,1);
    plot(tempPlot,'Color',colors(classPosIndex,:));

    if classPosIndex > prevClassIndex
        text(-textLength, -plotIndex+1,sprintf("Class %d",labels(plotIndex))); 
        prevClassIndex = classPosIndex;  
    end

end
hold off;
formattedTitle = sprintf("Training samples sorted by class");
title(formattedTitle);
set(gca,'xtick',[1,size(samples,2)],'ytick',[], 'TickDir','out');
xlim([-textLength,size(samples,2)]);
ylim([-plotIndex-0.5,1]);
box off;
hold off;

%%% Isolated Classes
for ii = 1:numClasses
    className = classes(ii);

    numInstances = sum(labels == className);
    
    classShapes = nan(numInstances, size(samples,2));
    count = 0;
    
    for li = 1:length(labels)
        if labels(li) ~= className
            continue;
        end
        count = count + 1;
        classShapes(count,:) = zscore(samples(li,:));
    end
    
    meanShape = nanmean(classShapes,1);
    stdShape = nanstd(classShapes,1);
    
    figure;
    plot(classShapes');
    
    figure; 
    hold on;
    plot(meanShape);
    plot(meanShape + stdShape,'--');
    plot(meanShape - stdShape,'--');
    hold off;
    
    plotDatasetFeatures(classShapes');
end
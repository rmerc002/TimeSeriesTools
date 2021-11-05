

for classIndex = 1:length(classes)
    plotDataTransforms(classTS{classIndex}(1:10000));
end

%%%Plot Plato with nearestest neighbors in context
numNN = 5;
plotQueryNN(trainClassTS{classIndex}, classPlatos{classIndex}(2,:), numNN);
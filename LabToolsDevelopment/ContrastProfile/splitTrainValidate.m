function [dataFoldIndices] = splitTrainValidate(labels, nFolds, maxTrainSamples)
    
    if length(nargin) == 2
        maxTrainSamples = inf;
    end
    classes = unique(labels);
    numClasses = length(classes);

    classStruct = {};
    numSamples = length(labels);
    xaxis= 1:numSamples;

    for classIndex = 1:numClasses
       class = classes(classIndex);
       classStruct{end+1} = xaxis(labels==class);
    end

    dataFoldIndices = table({},{},'VariableNames',{'Train','Validate'});

    for foldIndex = 1:nFolds
        tempSampleStructTrain = {};
        tempSampleStructValidate = {};
        for classIndex = 1:numClasses
            class = classes(classIndex);
            tempIndices = classStruct{classIndex}(randperm(length(classStruct{classIndex})));
            numTrain = floor(length(tempIndices)/2);
            numTrain = min(numTrain, maxTrainSamples);
            tempSampleStructTrain{classIndex} = tempIndices(1:numTrain);
            tempSampleStructValidate{classIndex} = tempIndices(numTrain+1:end);
        end

        tempTable = table(tempSampleStructTrain, tempSampleStructValidate,'VariableNames',{'Train','Validate'});
     
        dataFoldIndices = [dataFoldIndices;tempTable];
    end
end
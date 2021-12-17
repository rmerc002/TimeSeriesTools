function [platoModel, error] = validatePlatosOnMulticlassSegments(samples, labels, classPlatos, outputPath)

    if nargin == 3
        outputPath = "";
    end

    %%%validation errors for each length, 
    %%%
    platoModel = struct;

    classes = sort(unique(labels));
    numClasses = length(classes);
        

    lengthSortScores = inf(numClasses,numClasses);
    bestLengthIndices = zeros(numClasses,numClasses);

    classSortScores = cell(numClasses,numClasses);

    for posClassIndex = 1:numClasses
        posClass = classes(posClassIndex);
        posSamples = samples(labels==posClass,:);
        for negClassIndex = 1:numClasses
            if posClassIndex == negClassIndex
                continue;
            end
            negClass = classes(negClassIndex);
            negSamples = samples(labels==negClass,:);
            numLengths = length(classPlatos{posClassIndex, negClassIndex});
            classSortScores{posClassIndex, negClassIndex} = cell(numLengths,1);
            for lengthIndex = 1:numLengths
                platosStruct = classPlatos{posClassIndex, negClassIndex}{lengthIndex};
                platosPos = platosStruct.platos;
                [tempSortScore, bestK] = sortScoreOnTwoClassSegmented_OR(posSamples, negSamples, platosPos);
                validationHeuristicStruct = struct;
                validationHeuristicStruct.sortScore = tempSortScore;
                validationHeuristicStruct.K = bestK;
                classSortScores{posClassIndex, negClassIndex}{lengthIndex} = validationHeuristicStruct;
            end
        end
    end

    platoModel = constructPlatoModel(classPlatos, classSortScores);
    [ypred] = classifyWithPlatoModel_OR(samples, platoModel);
    error = 1-mean(ypred == labels);
    fprintf("validation error: %.2f\n", error);
end
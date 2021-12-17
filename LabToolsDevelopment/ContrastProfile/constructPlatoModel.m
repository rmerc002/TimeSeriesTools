function platoModel = constructPlatoModel(classPlatos, classSortScores)
    numClasses = size(classPlatos,1);
    platoModel = cell(numClasses,numClasses);
    for posClassIndex = 1:numClasses
        for negClassIndex = 1:numClasses
            if posClassIndex == negClassIndex
                continue;
            end
            
            numLengths = length(classSortScores{posClassIndex, negClassIndex});
            bestLengthIndex = 0;
            bestScore = 0;
            bestK = 1;
            for lengthIndex = 1:numLengths
                validationHeuristicStruct = classSortScores{posClassIndex, negClassIndex}{lengthIndex};
                sortScore = validationHeuristicStruct.sortScore;
                if sortScore > bestScore
                    bestScore = sortScore;
                    bestLengthIndex = lengthIndex;
                    bestK = validationHeuristicStruct.K;
                end
            end
            model = classPlatos{posClassIndex, negClassIndex}{bestLengthIndex};
            model.operator = "OR";
            model.bestK = bestK;
            platoModel{posClassIndex, negClassIndex} = model;
        end
    end
end
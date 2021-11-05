function [positiveTS,negativeTS,subLength,class,groundTruthIndices,dataChoiceIndices,localPath] = extractVariables(dataFileName)
    data = load(dataFileName);
    positiveTS = data.positiveTS;
    negativeTS = data.negativeTS;
    subLength = data.subLength;
    class = data.class;
    groundTruthIndices = data.groundTruthIndices;
    dataChoiceIndices = data.dataChoiceIndices;
    localPath = "";%data.localPath; %TODO: fix
end
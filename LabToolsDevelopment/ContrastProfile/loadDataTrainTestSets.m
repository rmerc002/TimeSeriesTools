function [dataTrainTestFileNames, numDataSet] = loadDataTrainTestSets(versionedCuratedDataPath)
    %%% Output:
    %%%   The dataset should be saved in a table
    %%%      First entry is the name of the data set
    %%%      Next is the [ {'dataset_name',{{train,test},...}, ...}
    numDataSet = 0;
    dataTrainTestFileNames = {};
    %%% Get all dataset directory names
    dirinfo1 = dir(versionedCuratedDataPath);
    dirinfo1(~[dirinfo1.isdir]) = [];  %remove non-directories

    for K = 1 : length(dirinfo1)
        trainTestFiles = {};
        thisDataSetDir = dirinfo1(K).name;
        if thisDataSetDir(1) == '.'
            continue
        end
        if thisDataSetDir(2) == '.'
            continue
        end


        %%% Iterate through all the classes
        dataSetPath = fullfile(versionedCuratedDataPath,thisDataSetDir);
        dirinfo2 = dir(dataSetPath);
        dirinfo2(~[dirinfo2.isdir]) = [];  %remove non-directories

        for K = 1 : length(dirinfo2)
            thisClassDir = dirinfo2(K).name;
            if thisClassDir(1) == '.' || thisClassDir(1) == '_'
                continue
            end
            if thisClassDir(2) == '.'
                continue
            end
        
            %%% Iterate through all the data gen iterations
            classPath = fullfile(dataSetPath, thisClassDir);
            dirinfo3 = dir(classPath);
            dirinfo3(~[dirinfo3.isdir]) = [];  %remove non-directories

            for K = 1 : length(dirinfo3)
                thisIterationDir = dirinfo3(K).name;
                if thisIterationDir(1) == '.' 
                    continue
                end
                if thisIterationDir(2) == '.'
                    continue
                end
            
                trainFileName = fullfile(versionedCuratedDataPath, thisDataSetDir, thisClassDir, thisIterationDir, 'train','dataVariables.mat');
                testFileName = fullfile(versionedCuratedDataPath, thisDataSetDir, thisClassDir, thisIterationDir, 'test','dataVariables.mat');
                
                if isfile(trainFileName) && isfile(testFileName)
                    trainTestFiles{end+1} = {trainFileName, testFileName};
                    numDataSet = numDataSet + 1;
                end
            end
        end
        dataTrainTestFileNames{end+1} = {thisDataSetDir, trainTestFiles};
    end
end
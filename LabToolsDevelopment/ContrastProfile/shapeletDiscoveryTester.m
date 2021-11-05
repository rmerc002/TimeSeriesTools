function tpfpfnIndices = shapeletDiscoveryTester(testVersion, dataVersionName, continueAnytime)
    rng(1);
    %%% This should eventually check previously created datasets
    %%% In the meantime, run this manually
    curatedDataPath = "C:\Users\Ryan\Dropbox\Education\UCR\Research\0_Projects\ActiveResearch\shapeletDiscovery\curatedData";

    versionedCuratedDataPath = fullfile(curatedDataPath, dataVersionName);

    testVersionString = sprintf('TestVersion_%04d',testVersion);
    


    %%% I want to show results at the dataset level. 
    %%% shuffle generated datasets at the data set level
    %%% Final plot showing number of false negatives, false positives 
    %%% at data set level
    [dataTrainTestFileNames, numGenDataSets] = loadDataTrainTestSets(versionedCuratedDataPath);

    for dataSetIndex = 1:length(dataTrainTestFileNames)
        dataSetName = dataTrainTestFileNames{dataSetIndex}{1};
        trainTestSets = dataTrainTestFileNames{dataSetIndex}{2};
        
        dataSetResultsPath = fullfile(versionedCuratedDataPath, dataSetName, '_Results');
        if ~isfolder(dataSetResultsPath)
            mkdir(dataSetResultsPath);
        end

        numTrainTestSets = length(trainTestSets);
        shuffledDataSetIndices = randperm(numTrainTestSets);
        
        fpfnCounts = ones(numTrainTestSets,2)*-1;
        tpfpfnIndices = cell(numTrainTestSets,1);
        orderedFileList = cell(numTrainTestSets,1);
        
        for shuffledIndex = 1:numTrainTestSets
            %%%stop shuffling
            dataIndex = shuffledIndex;%shuffledDataSetIndices(shuffledIndex);
            
            %%% Anytime-testing: check if already completed
            %%% save the completion file in the test versioned directory
            testFileName = trainTestSets{dataIndex}{2};
            [figureSavePath,~,~] = fileparts(testFileName);
            figureSavePath = fullfile(figureSavePath, testVersionString);
            completionFlagFileName = fullfile(figureSavePath,'completedTest.txt');
            if continueAnytime
                if isfile(completionFlagFileName)
                    %%% Load the stored fpfnCounts and tpfpfnIndices
                    %%% variables
                    previousOutputData = load(fullfile(figureSavePath, 'outputData.mat'));
                    fpfnCounts(dataIndex,:) = previousOutputData.fpfnCounts_single;
                    tpfpfnIndices{dataIndex} = previousOutputData.tpfpfnIndices_single;
                    orderedFileList{dataIndex} = figureSavePath;
                   continue; 
                end
            end
            
            %%% Train
            trainFileName = trainTestSets{dataIndex}{1};
            [figureSavePath,~,~] = fileparts(trainFileName);
            figureSavePath = fullfile(figureSavePath, testVersionString);
            
            if ~isfolder(figureSavePath)
                fprintf('figureSavePath: %s\n', figureSavePath);
                mkdir(figureSavePath);
            end
            
            
            
            fprintf("%d of %d: train file name: %s\n",shuffledIndex, numTrainTestSets, trainFileName);
            [positiveTSTrain,negativeTSTrain,subLength,class,groundTruthIndicesTrain,dataChoiceIndices,~] = extractVariables(trainFileName);
            [queryShapelet, queryThreshold, thresholdConfidence] = shapeletDiscovery(positiveTSTrain, negativeTSTrain, subLength, groundTruthIndicesTrain, figureSavePath);

            %%% Test
            testFileName = trainTestSets{dataIndex}{2};
            [figureSavePath,~,~] = fileparts(testFileName);
            figureSavePath = fullfile(figureSavePath, testVersionString);
            if ~isfolder(figureSavePath)
                mkdir(figureSavePath);
            end
            
            
            
            fprintf("%d of %d: test  file name: %s\n",shuffledIndex, numTrainTestSets, testFileName);
            [testTS, negativeTSTest, subLength, class, groundTruthIndicesTest, dataChoiceIndices, ~] = extractVariables(testFileName);
            [tpIndices, fpIndices,fnIndices] = evaluateAccuracy(testTS, subLength, groundTruthIndicesTest, queryShapelet, queryThreshold, thresholdConfidence, figureSavePath);

            %%%TODO: need to save the results in the parent directory containing train and test
            
            fpfnCounts_single = [length(fpIndices),length(fnIndices)]; %TODO: store the indices as well
            tpfpfnIndices_single = {tpIndices, fpIndices,fnIndices}; %TODO: fix ordering. needs to be in order of file list, not order of execution
            
            save(fullfile(figureSavePath, 'outputData.mat'), 'fpfnCounts_single','tpfpfnIndices_single');
            orderedFileList{dataIndex} = figureSavePath;
            fpfnCounts(dataIndex,:) = fpfnCounts_single;
            tpfpfnIndices{dataIndex} = tpfpfnIndices_single;
            
            
            %%% Save global results after every iteration for anytime
            %%%    testing purposes
            plotFPvsFN(fpfnCounts, dataSetResultsPath);

            fileName = "performanceVariables";
            filePath = fullfile(dataSetResultsPath, fileName + ".mat");
            save(filePath, 'tpfpfnIndices', 'fpfnCounts','orderedFileList');
        
            %%% Anytime-testing: check if already completed
            fclose(fopen(completionFlagFileName, 'w'));
        end

    end
    %%% Save the last completed test version
    %%% This help with anytime testing
    versionFileName = "lastVersionCompleted.mat";
    versionFilePath = fullfile(curatedDataPath, versionFileName);
    fid = fopen(versionFilePath, 'w');
    fprintf(fid, "%04d", testVersion);
    fclose(fid);
end
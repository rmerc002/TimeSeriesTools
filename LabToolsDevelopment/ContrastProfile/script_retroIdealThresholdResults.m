%%% I will save results in the test level directory

%%% get list of files
fileList = {}; %TODO

versionedCuratedDataPath = "C:\Users\Ryan\Dropbox\Education\UCR\Research\0_Projects\ActiveResearch\shapeletDiscovery\curatedData\datasetCreation_V08";
[~,dataVersion,~] = fileparts(versionedCuratedDataPath);
testVersion = 'TestVersion_0026';
%%% save results at 3 levels
%%%    1) overall, test version
%%%    2) data set
%%%    3) class 

%maybe this should be redundant, where the overall level can see the
%categorization of the previous levels, for ease of retreiving results


 

dirinfo1 = dir(versionedCuratedDataPath);
dirinfo1(~[dirinfo1.isdir]) = [];  %remove non-directories

testVersionResults = {versionedCuratedDataPath};
numDataSets = 0;
dataSetAggregate = [0,0,0];
for K = 1 : length(dirinfo1)
    thisDataSetDir = dirinfo1(K).name;
    
    if thisDataSetDir(1) == '.' || thisDataSetDir(1) == '_'
        continue
    end
    if thisDataSetDir(2) == '.'
        continue
    end
    fprintf('### %s\n',thisDataSetDir);
    %%% Iterate through all the classes
    dataSetPath = fullfile(versionedCuratedDataPath,thisDataSetDir);
    dirinfo2 = dir(dataSetPath);
    dirinfo2(~[dirinfo2.isdir]) = [];  %remove non-directories

    dataSetResults = {dataSetPath};
    numClasses = 0;
    classAggregate = [0,0,0];
    for K = 1 : length(dirinfo2)
        thisClassDir = dirinfo2(K).name;
        
        if thisClassDir(1) == '.' || thisClassDir(1) == '_'
            continue
        end
        if thisClassDir(2) == '.'
            continue
        end
        fprintf('###\t %s\n', thisClassDir);
        numClasses = numClasses + 1;
        %%% Iterate through all the data gen iterations
        classPath = fullfile(dataSetPath, thisClassDir);
        dirinfo3 = dir(classPath);
        dirinfo3(~[dirinfo3.isdir]) = [];  %remove non-directories

        classResults = {classPath};
        numIterations = 0;
        iterationAggregate = [0,0,0];
        for K = 1 : length(dirinfo3)
            thisIterationDir = dirinfo3(K).name;
            
            if thisIterationDir(1) == '.' || thisIterationDir(1) == '_'
                continue
            end
            if thisIterationDir(2) == '.'
                continue
            end
            fprintf('###\t\t %s\n',thisIterationDir);
            numIterations = numIterations +1;
            
            inputDataPath = fullfile(versionedCuratedDataPath, thisDataSetDir, thisClassDir, thisIterationDir, 'test','dataVariables.mat');
            threshDataPath = fullfile(versionedCuratedDataPath, thisDataSetDir, thisClassDir, thisIterationDir, 'test',testVersion,'internVars_findingIdealThreshold.mat');
            
            if isfile(inputDataPath) && isfile(threshDataPath)
                inputVars = load(inputDataPath);
                subLength = inputVars.subLength;
                
                threshVars = load(threshDataPath);
                profileDistances = threshVars.MP_AQ;
                groundTruthIndices = threshVars.groundTruthIndices;
                numSignals = length(groundTruthIndices);

                [dataSavePath,~,~] = fileparts(threshDataPath);
                [tpIndices, fpIndices, fnIndices, idealThreshold] = getIdealThresholdResults(profileDistances, groundTruthIndices, subLength, dataSavePath);
                classResults{end+1} = {tpIndices, fpIndices, fnIndices,idealThreshold};
                iterationAggregate = iterationAggregate + [length(tpIndices),length(fpIndices),length(fnIndices)];
            end
        end
        %%% class level
        iterationMean = iterationAggregate/numIterations;
        
        tempRemoveResults = fullfile(versionedCuratedDataPath, thisDataSetDir, thisClassDir,'Results');
        if isfolder(tempRemoveResults)
            rmdir(tempRemoveResults,'s');
        end
        
        classPathDirs = [versionedCuratedDataPath, thisDataSetDir, thisClassDir,'_Results'];
        classPath = mkdirRecursive(classPathDirs);
        dataSetFileName = testVersion;
        dataSetFilePath = fullfile(classPath, dataSetFileName + ".mat");
        save(dataSetFilePath, 'classResults','iterationAggregate');
        
        dataSetResults{end+1} = classResults;
        classAggregate = classAggregate + iterationAggregate;
    end
    %%% data set level
      classMean = classAggregate/numClasses;
      
      tempRemoveResults = fullfile(versionedCuratedDataPath, thisDataSetDir,'Results');
      if isfolder(tempRemoveResults)
        rmdir(tempRemoveResults,'s');
      end
        
      dataSetPathDirs = [versionedCuratedDataPath, thisDataSetDir, '_Results'];
      dataSetPath = mkdirRecursive(dataSetPathDirs);
      dataSetFileName = testVersion;
      dataSetFilePath = fullfile(dataSetPath, dataSetFileName + ".mat");
      save(dataSetFilePath, 'dataSetResults','classAggregate');

      testVersionResults{end+1} = dataSetResults;
      dataSetAggregate = dataSetAggregate + classAggregate;
end
%%% test version level
dataSetMean = dataSetAggregate/numDataSets;

tempRemoveResults = fullfile(versionedCuratedDataPath, 'Results');
if isfolder(tempRemoveResults)
    rmdir(tempRemoveResults, 's');
end
      
testVersionPathDirs = [versionedCuratedDataPath, '_Results'];
testVersionPath = mkdirRecursive(testVersionPathDirs);
testVersionFileName = testVersion;
testVersionFilePath = fullfile(testVersionPath, testVersionFileName + ".mat");
save(testVersionFilePath, 'testVersionResults','dataSetAggregate');

function proposedPath = mkdirRecursive(proposedDirs)
    proposedPath = ""; 
    for dir = proposedDirs
       proposedPath = fullfile(proposedPath,dir);
       if exist(proposedPath,'dir')~=7
           mkdir(proposedPath);
       end
    end
end

            
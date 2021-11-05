function datasetCreation_V04(datasets, curatedDataPath)
%%% Purpose:
%%%    Make the time series stationary


    numDataSetsPerClass = 2; %%%TODO: make this 5
    numDataSamples = 5;
    %%% temp for testing
%     curatedDataPath = "C:\Users\Ryan\Dropbox\Education\UCR\Research\0_Projects\ActiveResearch\shapeletDiscovery\curatedData";
%     datasets = {};
%     datasets{length(datasets)+1} = {"D:\Datasets\UCRArchive_2018\UCRArchive_2018\Mallat\Mallat_TRAIN.tsv","Mallat","train","segmented",false,0,1,2,1025};
%     datasets{length(datasets)+1} = {"D:\Datasets\UCRArchive_2018\UCRArchive_2018\Mallat\Mallat_TEST.tsv","Mallat","test","segmented",false,0,1,2,1025};
    %%% Needs to save all generated datasets in curatedDataPath
    
    [~,thisVersionName,~] = fileparts(mfilename('fillpath'));
    proposedDirsLevel2 = [curatedDataPath,thisVersionName]; %start at level2
    
    %handle each entry in datasets independently
    for di = 1:length(datasets)
        path = char(datasets{di}{1});
        dataSetName = datasets{di}{2};
        trainTest = datasets{di}{3};
        type = datasets{di}{4};
        needsTranspose = datasets{di}{5};
        headerlines = datasets{di}{6};
        classIndex = datasets{di}{7};
        sampleStartIndex = datasets{di}{8};
        sampleEndIndex = datasets{di}{9};
        
        
        
        %%% {"D:\Datasets\UCRArchive_2018\UCRArchive_2018\Mallat\Mallat_TRAIN.tsv","Mallat","train","segmented",false,1}; 
        pathLength = strlength(path);
        proposedDirsLevel3 = [proposedDirsLevel2,dataSetName];
        
        %%% Specific load based on extension
        data = [];
        if strcmp(path(pathLength-3:pathLength), ".csv")
            opts = detectImportOptions(path,'FileType','text');
            data = readtable(path,opts);%,'FileType','text','HeaderLines',headerlines);
        elseif strcmp(path(pathLength-3:pathLength),".tsv")
            opts = detectImportOptions(path,'FileType','text');
            data = readtable(path,opts);%,'FileType','text','HeaderLines',headerlines);
        else
            fprintf("Only .csv and .tsv types understood. Deal with %s to make it so.\n", path);
            continue;
        end
        
        %%% Set to expected shape
        data = data{:,:};
        if needsTranspose
           data = data'; 
        end
    
        
        classes = data(:,classIndex);
        if sampleEndIndex == -1
           sampleEndIndex = size(data,2); 
        end
        samples = data(:,sampleStartIndex:sampleEndIndex);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Generate Time Series %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        classLabels = unique(classes);
        numClasses = length(classLabels);
        for classIndex = 1:numClasses %%%TODO: comment this back in
            
            class = classLabels(classIndex);
            proposedDirsLevel4 = [proposedDirsLevel3,sprintf("Class_%d",class)];
            numClassInstances = sum(classes == class);
            %use this to look up the index of the i-th instance of a
            %of the current class
            intraClassIndexToDataIndex = 1:length(data);
            intraClassIndexToDataIndex = intraClassIndexToDataIndex(classes == class);
            %%% if there are not enough samples, do not bother with a
            %%% dataset
            if length(intraClassIndexToDataIndex) < numDataSamples
                continue;
            end
            for classIteration = 1:numDataSetsPerClass
                proposedDirsLevel5 = [proposedDirsLevel4, sprintf("Iteration_%d",classIteration)];
                proposedDirsLevel6 = [proposedDirsLevel5, trainTest];
                localPath = mkdirRecursive(proposedDirsLevel6);
                %%% plot the samples chosen
                
                numSamples = size(data,1);

                

                %sample without replacement
                intraClassChoiceIndices = randsample(numClassInstances,numDataSamples);
                dataChoiceIndices = intraClassIndexToDataIndex(intraClassChoiceIndices);
                fprintf('The following indices are class %d, data type = %s\n',class, trainTest);
                disp(dataChoiceIndices);
                disp(classes(dataChoiceIndices));
                
                figRandomSignals = figure;
                plot(0,0);
                hold on;
                %%% plot the chosen samples
                for plotIndex = 1:numDataSamples
                   classSampleIndex = intraClassChoiceIndices(plotIndex);
                   dataIndex = intraClassIndexToDataIndex(classSampleIndex);
                   tempTS = samples(dataIndex,:);
                   tempMin = min(tempTS);
                   tempMax = max(tempTS);
                   tempRange = max(1e-5, tempMax-tempMin);
                   tempTSNorm = (tempTS-tempMin)/tempRange;

                   plot(-plotIndex+0.9*tempTSNorm);
                end
                hold off;
                
                set(gcf, 'Position', [0,100,200,400]);
                set(gca,'xtick',[],'ytick',[]);
                formattedTitle = sprintf("%d Randomly Chosen Indices for Positive Class %d: ",numDataSamples, class);
                for dci = dataChoiceIndices
                   formattedTitle = formattedTitle + sprintf(", %d", dci); 
                end
                xlim([1,size(samples,2)]);
                title(formattedTitle);
                
                figureFileName = "randomSignalSelection";
                saveas(figRandomSignals, fullfile(localPath, figureFileName + ".fig"));
                saveas(figRandomSignals, fullfile(localPath, figureFileName + ".png"));
                close(figRandomSignals);
%tested above (quickly)

                positiveSamples = samples(dataChoiceIndices,:);
                signalLength = size(positiveSamples,2);
                desiredDefaultRate = 0.90; %for negative class
                sampleSectionRange = ceil(signalLength/(1-desiredDefaultRate));
                fillerLength = sampleSectionRange - signalLength;

                pregeneratedRandomWalk = getRandWalk(ceil(fillerLength*numDataSamples*1.1));%extra just so I don't have to be precise in the calculation below
                
                subLength = signalLength;
                randOffset = randi(fillerLength, 1,numDataSamples,1);
                randOffset(1) = max(signalLength,randOffset(1));%I use the last <signalLength> portion to make the time series stationary

                groundTruthIndices = zeros(numDataSamples,1);
                fprintf("randomWalkLength: %d, preGenLength: %d\n",randOffset(1), length(pregeneratedRandomWalk));

                positiveTS = pregeneratedRandomWalk(1:randOffset(1));
                pregeneratedRandomWalk(1:randOffset(1)) = [];
                for i = 1:numDataSamples
%                     part1 = zscore(getRandWalk(ro));
%                     part1 = part1 - part1(1) + positiveTS(end);
                    groundTruthIndices(i) = length(positiveTS)+1;
                    signalPart = zscore(positiveSamples(i,:));
                    prevRandSTD = std(positiveTS(end-signalLength+1:end));
%                     figure; plot(positiveTS(end-signalLength+1:signalLength));
                    fprintf("prevRandSTD: %f\n", prevRandSTD);
                    signalPart = signalPart * (prevRandSTD);
                    signalPart = signalPart - signalPart(1) + positiveTS(end);
                    
                    signalPart2 = signalPart;
                    signalPart2 = signalPart2 - signalPart(2) + signalPart(end);
                    
                    currentRandomOffset = randOffset(i);
                    if i < numDataSamples
                        nextRandomOffset = randOffset(i+1);
                        randomWalkLength = fillerLength - currentRandomOffset + nextRandomOffset;
                        fprintf("randomWalkLength: %d, preGenLength: %d\n",randomWalkLength, length(pregeneratedRandomWalk));
                        randomPart = pregeneratedRandomWalk(1:randomWalkLength);
                        pregeneratedRandomWalk(1:randomWalkLength) = [];
                        randomPart = randomPart - randomPart(1) + signalPart(end);
                    else
                        randomWalkLength = fillerLength - currentRandomOffset;
                        randomPart = pregeneratedRandomWalk(1:randomWalkLength);
                        pregeneratedRandomWalk(1:randomWalkLength) = [];
                        randomPart = randomPart - randomPart(1) + signalPart(end);
                    end
                    positiveTS = [positiveTS, signalPart, signalPart2, randomPart];
                    fprintf("length positiveTS: %d\n",length(positiveTS));
                end
                
                if trainTest == "train"
                    negativeTS = getRandWalk(length(positiveTS));
                else
                   netativeTS = getRandWalk(signalLength);%just some dummy data, should not be used. 
                end
                
                %%%%%%%%%%%%%%%%%%%%
                %%% Save Results %%%
                %%%%%%%%%%%%%%%%%%%%
                
                %%%[positiveTSTrain, 
                %%% negativeTSTrain, 
                %%% subLength, 
                %%% groundTruthIndicesTrain, 
                %%% pathToOriginalData,
                %%% indices of samples used in order,
                %%% class type of each sample]
                
                %%% save the generated data
                dataSetFileName = "dataVariables";
                dataSetFilePath = fullfile(localPath, dataSetFileName + ".mat");
                save(dataSetFilePath, 'positiveTS','negativeTS','subLength','class','groundTruthIndices','dataChoiceIndices','localPath');
                
                %%% append the newly created filename to a list for
                %%% convenience
                datasetVersionPath = mkdirRecursive(proposedDirsLevel2);
                if strcmp(trainTest, 'train')
                    trainTestFilePath = fullfile(datasetVersionPath, "trainList.txt");
                else
                    trainTestFilePath = fullfile(datasetVersionPath, "testList.txt");
                end
                fid = fopen(trainTestFilePath, 'a+');
                fprintf(fid, "%s\n",dataSetFilePath);
                fclose(fid);
            
                %%% plot the generated time series
                figTimeSeries = figure;
                set(gcf, 'Position', [0,100,2050,200]);
                subplot(2,1,1);
                plot(negativeTS);
                formattedTitle = sprintf("Negative Time Series, Length=%d",length(negativeTS));
                title(formattedTitle);
                set(gca,'xtick',[],'ytick',[]);
                xlim([1,max(length(negativeTS), length(positiveTS))]);

                subplot(2,1,2);
                plot(positiveTS);
                hold on;
                for i = 1:numDataSamples
                    index = groundTruthIndices(i);
                   plot(index:index+subLength-1, positiveTS(index:index+subLength-1),'green');
                end
                hold off;
                formattedTitle = sprintf("Positive Time Series, Length=%d",length(positiveTS));
                title(formattedTitle);
                set(gca,'xtick',[],'ytick',[]);
                xlim([1,max(length(negativeTS), length(positiveTS))]);
                
                figureFileName = "negativeAndPositiveTimeSeries";
                saveas(figTimeSeries, fullfile(localPath, figureFileName + ".fig"));
                saveas(figTimeSeries, fullfile(localPath, figureFileName +".png"));
                close(figTimeSeries);
                
                
                %%% Plot showing stationary property of positive time
                %%% series
                figTimeSeries = figure;
                subplot(2,1,1);
                plot(positiveTS);
                hold on;
                for i = 1:numDataSamples
                    index = groundTruthIndices(i);
                   plot(index:index+subLength-1, positiveTS(index:index+subLength-1),'green');
                end
                hold off;
                formattedTitle = sprintf("Positive Time Series, Length=%d",length(positiveTS));
                title(formattedTitle);
                set(gca,'xtick',[],'ytick',[]);
                xlim([1,length(positiveTS)]);
                
                subplot(2,1,2);
                slidingSTD = movstd(positiveTS,signalLength);
                plot(slidingSTD);
                formattedTitle = sprintf("Sliding STD, Window Length=%d",signalLength);
                title(formattedTitle);
                set(gca,'xtick',[],'ytick',[]);
                xlim([1,length(positiveTS)]);
                
                figureFileName = "stationaryPositiveTimeSeries";
                saveas(figTimeSeries, fullfile(localPath, figureFileName + ".fig"));
                saveas(figTimeSeries, fullfile(localPath, figureFileName +".png"));
                close(figTimeSeries);
            end
        end
    end 
end


function proposedPath = mkdirRecursive(proposedDirs)
    proposedPath = ""; 
    for dir = proposedDirs
       proposedPath = fullfile(proposedPath,dir);
       if exist(proposedPath,'dir')~=7
           mkdir(proposedPath);
       end
    end
end
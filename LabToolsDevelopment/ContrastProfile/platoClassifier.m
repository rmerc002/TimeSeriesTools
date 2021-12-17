function model = platoClassifier(pos, neg)
    rng(1);
%     archivePath = "D:\Datasets\UCRArchive_2018\UCRArchive_2018";
%     outputPath = "C:\Users\Ryan\Dropbox\Education\UCR\Research\0_Projects\ActiveResearch\ContrastProfile\Experiments\2021-11-17_UCRArchive_Attempt2";
    experimentName = "DatasetResults_RFCP2_NegIndividual";
    baseOutputPath = fullfile(outputPath, experimentName);
    if ~exist(baseOutputPath, 'dir')
           mkdir(baseOutputPath);
    end

    forcePlot = false;
    subLengthCoefs = [1/4];%, 1/4, 1/2, 3/4];

    %%%Generating unique colors for classes when there can be many classes
    %%% I will assume no more than 1000 classes
    numColors = 1000;
    colors = lines(numColors);

    % Get a list of all files in the folder with the desired file name pattern.
    filePattern = fullfile(archivePath); % Change to whatever pattern you need.
    theFiles = dir(filePattern);
    %%%The shapelets papers did not run on every dataset due to time
    %%%contraints. For this reason, we limited our comparison to only
    %%%dataset which we could compare to
%     shapeletDatasets = [4,6,9,11,12,13,18,22,24,28,29,30,27,33,40,42,43,44,46,50,54,53,66,67,70,71,82,74,76,75,86,87];%2015 archive indices
    shapeletDatasets = [4, 10, 13, 16, 17, 18, 24, 31, 33, 39, 40, 41, 43, 54, 67, 69, 70, 71, 73, 81, 84, 85, 109, 110, 113, 114, 115, 118, 120, 119, 126, 131];
    for fileIndex = 16%shapeletDatasets%1:size(theFiles,1)
        close all;
        datasetName = theFiles(fileIndex).name;
        if datasetName(1) == '.'
            continue;
        end
        baseFileParts = strsplit(theFiles(fileIndex).folder,'\');
        subjectName = baseFileParts(end);

        fullTrainFileName = fullfile(theFiles(fileIndex).folder, datasetName,sprintf("%s_TRAIN.tsv",datasetName));
        fullTestFileName = fullfile(theFiles(fileIndex).folder, datasetName,sprintf("%s_TEST.tsv",datasetName));
        fprintf(1, 'Now reading %s\n', datasetName);
        % Now do whatever you want with this file name,
        % such as reading it in as an image array with imread()
        if ~exist(fullTrainFileName, 'file') || ~exist(fullTestFileName, 'file')
            fprintf("\n%%%%%%%%%% ERROR: skipping: %s\n\n", datasetName);
           continue; 
        end

        dataTrain = importdata(fullTrainFileName);
        dataTest = importdata(fullTestFileName);
        datasetOutputPath = fullfile(baseOutputPath, datasetName);
        if ~exist(datasetOutputPath, 'dir')
           mkdir(datasetOutputPath);
        end

%         fileName = "completed.txt";
%         filePath = fullfile(datasetOutputPath, fileName);
%         if exist(filePath)
%            continue; 
%         end

        %%% Training
        labelsTrain = dataTrain(:,1);
        classes = sort(unique(labelsTrain));
        numClasses = length(classes);

        samples = dataTrain(:,2:end);
        
%         if size(samples,1) > 500
%            continue; 
%         end

        [labelsTrain, sortedLabelIndices] = sort(labelsTrain);
        samples = samples(sortedLabelIndices,:);

        %%%Plot the samples which are sorted by class
        %%%Each class is plotted with a unique color
        fig = figure('Name','UCR Contrast Profile: Top-K Platos','NumberTitle','off');
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

            classPosIndex = find(labelsTrain(plotIndex) == classes,1);
            plot(tempPlot,'Color',colors(classPosIndex,:));

            if classPosIndex > prevClassIndex
                text(-textLength, -plotIndex+1,sprintf("Class %d",labelsTrain(plotIndex))); 
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

        fileName = "Samples";
        filePath = fullfile(datasetOutputPath, fileName + ".fig");
        savefig(fig, filePath);

        filePath = fullfile(datasetOutputPath, fileName + ".emf");
        print(filePath,'-dmeta');

        ts = {};
        for classPosIndex =1:numClasses
           ts{classPosIndex} = []; 
        end

        for sampleIndex = 1:size(samples,1)
            classPosIndex = find(labelsTrain(sampleIndex) == classes,1);
            ts{classPosIndex} = [ts{classPosIndex}, NaN, samples(sampleIndex,:)];
        end

        %%%Search set of subLengthCoefs.
        %%%This is a guided Pan Contrast Profile
        classPlatos = cell(numClasses, numClasses);
        classContrast = nan(numClasses, numClasses);
    %     classPlatoDistances = {};
        for classPosIndex = 1:numClasses
            positiveTS = ts{classPosIndex};
            
            for classNegIndex = 1:numClasses
                if classPosIndex == classNegIndex
                    continue;
                end
                classPlatos{classPosIndex, classNegIndex} = {};

                negativeTS = ts{classNegIndex};
                
            
                bestContrast = -1;
                bestCoefIndex = 1;
                for coefIndex = 1:length(subLengthCoefs)
                    subLengthCoef = subLengthCoefs(coefIndex);
                    subLength = round(subLengthCoef*size(samples,2));
%                   [platos,platoIndices,CP] = ContrastProfile_TopK(positiveTS, negativeTS, subLength, maxK, forcePlot);
%                     [plato, ~, CP] = ContrastProfile(positiveTS, negativeTS, subLength, forcePlot);
%                     [maxContrast, maxCPIndex] = max(CP);
                    KPos = ceil(sum(labelsTrain == classes(classPosIndex))); %%number of samples for the target class
                    KNeg = ceil(sum(labelsTrain == classes(classNegIndex))); %%number of samples for the target class
                    K = min(KPos, KNeg);
%                     K = min(5, K);
                    [plato, RFCP, RFMP_AA_Indices] = RelativeFrequencyContrastProfile(positiveTS, negativeTS, subLength, K, forcePlot);
                    CPmean = zeros(1,size(RFCP,2));
                    for ti = 1:length(positiveTS)-subLength+1
                        subsequence = positiveTS(ti:ti+subLength-1);
                        if sum(isnan(subsequence)) > 0
                            continue
                        end
                       CPmean(ti) = norm(RFCP(:,ti))/sqrt(K); %% root mean squared
                    end
                    [maxContrast, maxCPIndex] = max(CPmean);
                    startIndex = maxCPIndex;
                    endIndex = startIndex + subLength - 1;
                    plato = positiveTS(startIndex:endIndex);
                    fprintf("classPosIndex: %d, classNegIndex: %d, \tcoef: %f, plato contrast: %.4f\n",classPosIndex, classNegIndex, subLengthCoefs(coefIndex), maxContrast);
                    if maxContrast > bestContrast
                        bestContrast = maxContrast;
                        bestCoefIndex = coefIndex;
                        classPlatos{classPosIndex, classNegIndex} = plato;
                        classContrast(classPosIndex, classNegIndex) = maxContrast;
                    end
                end
                fprintf("classPosIndex: %d, classNegIndex: %d, bestCoef: %.3f, bestContrast: %.4f\n",classPosIndex, classNegIndex, subLengthCoefs(bestCoefIndex),bestContrast);
            end
        end

        %%% Plot heatmap of class contrast
        figure;
        heatmap(classContrast);

        fileName = "Platos.mat";
        filePath = fullfile(datasetOutputPath, fileName);
        save(filePath, 'classPlatos', 'classContrast', 'classes');
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Training Best SubLength for Given K %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%indices into length coefficients
%         %%%each selection of K is assigned a learned best coef
% 
%         model = {};
%         for classIndex = 1:numClasses
%            model{classIndex} = classPlatos{classIndex}; 
%         end
% 
%         [ypred] = classify(dataTrain, model);
% 
%         ypredClassLabels = nan(size(ypred));
%         for classIndex = 1:numClasses
%            ypredClassLabels(ypred==classIndex) = classes(classIndex);
%         end
%         %%%
%         results = labelsTrain == ypredClassLabels;
%         error = 1-mean(results);
%                 errorsKCoef(KIndex, coefIndex) = error;
% 
%         
%         [minKErrors,bestSubLengthCoefIndices] = min(errorsKCoef,[],2);
%         trainedKsubLengthCoef = subLengthCoefs(bestSubLengthCoefIndices);
%         [minError, bestKIndex] = min(minKErrors);
%         bestK = Ks(bestKIndex);
%         bestCoefIndex = bestSubLengthCoefIndices(bestKIndex);
%         bestCoef = subLengthCoefs(bestCoefIndex);
%         fileName = "trainedKsubLengthCoef.mat";
%         filePath = fullfile(datasetOutputPath, fileName);
%         save(filePath,'trainedKsubLengthCoef','Ks', 'subLengthCoefs', 'errorsKCoef','bestK','bestCoef','bestCoefIndex','bestKIndex','bestSubLengthCoefIndices');
%         
%         fileName = "trainedKsubLengthCoef.csv";
%         filePath = fullfile(datasetOutputPath, fileName);
%         writematrix(trainedKsubLengthCoef,filePath);
%         
%         fileName = "trainingErrors.csv";
%         filePath = fullfile(datasetOutputPath, fileName);
%         writematrix(errorsKCoef,filePath);
        
%%%%%%%%%%%%%%%
%%% Testing %%%
%%%%%%%%%%%%%%%
        labels = dataTest(:,1);
        samples = dataTest(:,2:end);
        
        %%%prepare the model based on desired K and subLength
%         model = {};
%         for classPosIndex = 1:numClasses
%            model{classPosIndex} = classPlatos{classPosIndex,:}; 
%         end
        model = classPlatos;
        %%%ADD classifier here
        [ypred] = classify(samples, model);
        ypredClassLabels = nan(size(ypred));
        for classPosIndex = 1:numClasses
           ypredClassLabels(ypred==classPosIndex) = classes(classPosIndex);
        end
        %%%
        results = labels == ypredClassLabels;

        fileName = "individualClassifications.mat";
        filePath = fullfile(datasetOutputPath, fileName);
        save(filePath,'labels','ypredClassLabels');

        %%% Confusion matrix
        fig = figure('Name','Plato As Classifier: Confusion Matrix','NumberTitle','off');
        confusionchart(labels, ypredClassLabels);

        fileName = sprintf("ConfusionMatrix");
        filePath = fullfile(datasetOutputPath, fileName + ".fig");
        savefig(fig, filePath);

        filePath = fullfile(datasetOutputPath, fileName + ".emf");
        print(filePath,'-dmeta');

        %%%
        error = 1-mean(results);

        
        fprintf("error: %.4f,",error);
        
        fprintf("\n");
        fileName = "error.csv";
        filePath = fullfile(datasetOutputPath, fileName);
        writematrix(error, filePath);

        %%plot the learned platos
        for classPosIndex = 1:numClasses
            %%%Plot Platos with class color
            fig = figure('Name','UCR Contrast Profile: RFCP','NumberTitle','off');
            fig.Position = [300 300 200 200];
            inset = 0.9;
            hold on;
            maxPlatoLength = 0;
            for plotIndex = 1:size(classPlatos,2)
                if classPosIndex == plotIndex
                    continue;
                end
                tempTS = classPlatos{classPosIndex, plotIndex};
                tempMin = min(tempTS);
                tempMax = max(tempTS);
                tempRange = max(1e-5, tempMax - tempMin);
                tempPlot = -plotIndex + inset*(tempTS - tempMin)/tempRange;

                plot(tempPlot,'Color',colors(classPosIndex,:));

                platoLength = length(tempTS);
                maxPlatoLength = max(maxPlatoLength, length(tempTS));
            end
            formattedTitle = sprintf("Plato, Class %d",classes(classPosIndex));
            title(formattedTitle);
            set(gca,'xtick',[1,maxPlatoLength],'ytick',[], 'TickDir','out');
            xlim([1,platoLength]);
            box off;
            hold off;

            fileName = sprintf("LearnedPlatos_ClassIndex%03d",classPosIndex);
            filePath = fullfile(datasetOutputPath, fileName + ".fig");
            savefig(fig, filePath);

            filePath = fullfile(datasetOutputPath, fileName + ".emf");
            print(filePath,'-dmeta');
        end
        
        

        
        fileName = "completed.txt";
        filePath = fullfile(datasetOutputPath, fileName);
        fclose(fopen(filePath, 'w'));
    end
end

function [ypred] = classify(X, model)
    ypred = zeros(size(X,1),1);
    for sampleIndex = 1:size(X,1)
        sample = X(sampleIndex,:);
        closestClassIndex = NaN;
        minDist = inf;
        %%%test each plato of each class
        for classIndex = 1:length(model)
           numPlatos = size(model,2);
           classRMSE = 0;
           for platoIndex = 1:numPlatos
               if classIndex == platoIndex
                   continue;
               end
               plato = model{classIndex, platoIndex};
               distanceProfile = real(MASS_V2(sample,plato));
               dist = min(distanceProfile);
               classRMSE = classRMSE + dist^2;
           end
           classRMSE = sqrt(classRMSE/numPlatos);
           if classRMSE < minDist
               minDist = classRMSE;
               closestClassIndex = classIndex;
           end
        end
        ypred(sampleIndex) = closestClassIndex; 
    end
end
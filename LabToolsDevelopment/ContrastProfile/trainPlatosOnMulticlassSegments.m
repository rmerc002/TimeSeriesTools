function classPlatos = trainPlatosOnMulticlassSegments(samples, labels, subLengths, K, maxFreq, outputPath)

    if nargin == 3
        K = 1;
    elseif nargin == 4
        maxFreq = 1;
    elseif nargin == 5
        outputPath = "";
    end

    classes = unique(labels);
    numClasses = length(classes);

    classPlatos = cell(numClasses, numClasses);
        
    ts = {};
    for classPosIndex =1:numClasses
       ts{classPosIndex} = []; 
    end

    for sampleIndex = 1:size(samples,1)
        classPosIndex = find(labels(sampleIndex) == classes,1);
        ts{classPosIndex} = [ts{classPosIndex}, NaN, samples(sampleIndex,:)];
    end

    for classPosIndex = 1:numClasses
        positiveTS = ts{classPosIndex};
        
        for classNegIndex = 1:numClasses
            if classPosIndex == classNegIndex
                continue;
            end
            fprintf("Learning Plato for classPosIndex: %d, classNegIndex: %d\n", classPosIndex, classNegIndex);

            classPlatos{classPosIndex, classNegIndex} = {};

            negativeTS = ts{classNegIndex};
            
            for lengthIndex = 1:length(subLengths)
                subLength = subLengths(lengthIndex);

                numPos = ceil(sum(labels == classes(classPosIndex))); %%number of samples for the target class
                numNeg = ceil(sum(labels == classes(classNegIndex))); %%number of samples for the target class

                if maxFreq == 1
                    if K == 1
                        [platos, ~, CP] = ContrastProfile(positiveTS, negativeTS, subLength);
                        [platoContrasts, ~] = max(CP);
                    else
                        [platos, ~, platoContrasts] = ContrastProfile_TopK(positiveTS, negativeTS, subLength, K, false);
                    end
                else
                    maxFreq = min([numPos, numNeg, maxFreq]);
                    if K == 1
                        [platos, ~, platoContrasts] = RelativeFrequencyContrastProfile(positiveTS, negativeTS, subLength, maxFreq, false);
                    else
                        [platos, ~, platoContrasts] = RelativeFrequencyContrastProfile_TopK(positiveTS, negativeTS, subLength, maxFreq, K, false);
                    end
%                     CPmean = zeros(1,size(RFCP,2));
%                     for ti = 1:length(positiveTS)-subLength+1
%                         subsequence = positiveTS(ti:ti+subLength-1);
%                         if sum(isnan(subsequence)) > 0
%                             continue
%                         end
%                        CPmean(ti) = norm(RFCP(:,ti))/sqrt(K); %% root mean squared
%                     end
%                     [maxContrast, maxCPIndex] = max(CPmean);
%                     startIndex = maxCPIndex;
%                     endIndex = startIndex + subLength - 1;
%                     platos = positiveTS(startIndex:endIndex);
                end
                
                tempStruct = struct;
                tempStruct.subLength = subLength;
                tempStruct.platos = platos;
                tempStruct.contrasts = platoContrasts;
                classPlatos{classPosIndex, classNegIndex}{end+1} = tempStruct; 
            end
        end
    end

%     %%% Plot heatmap of class contrast
%     figure;
%     heatmap(classContrast);
% 
%     fileName = "Platos.mat";
%     %         filePath = fullfile(datasetOutputPath, fileName);
%     %         save(filePath, 'classPlatos', 'classContrast', 'classes');
end
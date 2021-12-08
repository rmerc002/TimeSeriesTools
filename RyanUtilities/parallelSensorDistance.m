function parallelSensorDistance(data)
    %%%dimsension: sensors x samples
    featureWindow = 20*60;
    for rowIndex = 1:size(data,1)
        rowIndex
       data(rowIndex,:) = movcomplexity(data(rowIndex,:), featureWindow);
    end
    plotDatasetFeatures(data(:,featureWindow:end)');
    
    output = nan(size(data));
    m = 10*20;
    for colIndex = 1:size(data,2)-m+1
        startIndex = colIndex
        endIndex = startIndex + m - 1;
        for rowIndex1 = 4%:size(data,1)
            subSeq1 = data(rowIndex1, startIndex:endIndex);
            subSeq1Norm = normalize(subSeq1,'range');
            for rowIndex2 = 1:size(data,1)
                subSeq2 = data(rowIndex2, startIndex:endIndex);
                subSeq2Norm = normalize(subSeq2,'range');
            
                output(rowIndex2, startIndex) = norm(subSeq1Norm - subSeq2Norm);
            end
        end
    end
    
    outputNorm = nanmin(output, sqrt(2*m));
    outputNorm = outputNorm./sqrt(2*m);
    outputNorm = 1-outputNorm;
    outputNorm(isnan(output)) = nan;
    plotDatasetFeatures(outputNorm');
end
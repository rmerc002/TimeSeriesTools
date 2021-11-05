function averageError = getSmoothingError(smoothLength, subLength)
    data = movmean(randn(1,subLength*100), smoothLength);
    averageError = 0;
%     numSums = 100000;
%     maxIndex = length(data)-subLength+1;
%     for index = 1:numSums
%         index1 = randi([1 maxIndex]);
%         index2 = randi([1 maxIndex]);
%         tempError = norm(data(index1:index1+subLength-1) - data(index2:index2+subLength-1)) / sqrt(2*subLength);
%         averageError = averageError + tempError;
%     end
%     averageError = averageError/numSums;
    numSums = 50;
    for index = 1:numSums
       tempError = mean(mpx(movmean(randn(subLength*10,1),smoothLength),ceil(subLength/2), subLength)); 
       averageError = averageError + tempError;
    end
    averageError = averageError/numSums;
end
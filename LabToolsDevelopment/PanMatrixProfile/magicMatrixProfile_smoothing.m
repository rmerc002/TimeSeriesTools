function  [dataSeries, magicMP, profileIndices,smoothSeries] = magicMatrixProfile_smoothing(data, subLength, startLength, endLength, numLengths)
if size(data,1) > size(data,2)
    data = data';
end

smoothSeries = getSubLenSeries(startLength, endLength, numLengths);

dataSeries = zeros(length(smoothSeries), length(data));
magicMP = zeros(length(smoothSeries),length(data)); %preallocate
profileIndices = zeros(length(smoothSeries), length(dataSeries)); %preallocate

for smoothIndex = 1:length(smoothSeries)
    smoothSeries(smoothIndex)
    smoothData = movmean(data, smoothSeries(smoothIndex));
    dataSeries(smoothIndex,:) = smoothData;
    
    
    [matrixProfile,profileIndex] = mpx(smoothData',ceil(subLength/2),subLength);

    
    
    matrixProfile(matrixProfile < 0) = 0;
%     averageNoiseError = getSmoothingError(smoothSeries(smoothIndex), subLength);
    magicMP(smoothIndex,1:length(matrixProfile)) = matrixProfile;%./averageNoiseError;
    magicMP(smoothIndex, length(smoothData)-subLength:end) = nan;

    profileIndices(smoothIndex,1:length(profileIndex)) = profileIndex; 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   HELPER FUNCTIONS   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function subLenSeries = getSubLenSeries(startLength, endLength, numLengths)
%%% Purpose: Reduce space. matrix profile distances are pretty similar from
%%% from one subLen to the next. By trial and error, adding the square root
%%% of the current subLen seems to produce a good distribution.
%     subLenSeries = [];
%     index = startLength;
%     
%     while index < endLength
%         subLenSeries = [subLenSeries,index];
%         index = index + max(1,ceil(sqrt(index)));
%     end
%     
%     subLenSeries = [subLenSeries,endLength];

powerMin = log10(startLength);
powerMax = log10(endLength);
powerStep = (powerMax-powerMin)/numLengths;
powers = powerMin:powerStep:powerMax;
subLenSeries = unique(ceil(power(10,powers)));

end

function averageError = getSmoothingError(smoothLength, subLength)
%     data = movmean(randn(1,subLength*100), smoothLength);
%     averageError = 0;
%     numSums = 100000;
%     maxIndex = length(data)-subLength+1;
%     for index = 1:numSums
%         index1 = randi([1 maxIndex]);
%         index2 = randi([1 maxIndex]);
%         tempError = norm(data(index1:index1+subLength-1) - data(index2:index2+subLength-1)) / sqrt(2*subLength);
%         averageError = averageError + tempError;
%     end
%     averageError = averageError/numSums;
    data = movmean(randn(1,max(subLength,smoothLength)*100), smoothLength);
    averageError = 0;
    numSums = 25;
    for index = 1:numSums
       tempError = max(mpx(movmean(randn(subLength*10,1),smoothLength),ceil(subLength/2), subLength),[],'all'); 
       averageError = averageError + tempError;
    end
    averageError = averageError/numSums;
end
